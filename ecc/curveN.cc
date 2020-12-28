// +build curve,!sm2p

/*
 * Copyright (c) 2013, 2014 Kenneth MacKay. All rights reserved.
 * Copyright (c) 2019 Vitaly Chikunov <vt@altlinux.org>
 * Copyright(c) 2020 Jesse Kuang  <jkuang@21cn.com>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <errno.h>
#ifdef	WITH_SYS_RANDOM
#include <sys/random.h>
#endif
#include "ecc.h"
#include "curve_defs.hpp"
#include "mont.hpp"

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#pragma GCC pop_options

#ifndef	ARRAY_SIZE
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))
#endif


static forceinline const ecc_curve *ecc_get_curve(uint curve_id)
{
	switch (curve_id) {
	/* In FIPS mode only allow P256 and higher */
	case ECC_CURVE_SECP256K1:
		return &secp256k1;
	case ECC_CURVE_NIST_P256:
		return &nist_p256;
	case ECC_CURVE_SM2:
		return &sm2_p256;
	default:
		return nullptr;
	}
}


CURVE_HND	get_curve(uint curve_id)
{
	return (CURVE_HND)ecc_get_curve(curve_id);
}


/**
 * get_curve_params	--	get curve params
 * p, n, b, gx, gy	--	bn_t 256 Bits
 */
void	get_curve_params(u64 *p, u64 *n, u64 *b, u64 *gx, u64 *gy,
				CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	ecc_curve	*curve=(ecc_curve *)curveH;
	if (curve->name[0] == 0 || curve->ndigits != 4) return;
	vli_set<4>(p, curve->p);
	vli_set<4>(n, curve->n);
	vli_set<4>(b, curve->b);
	vli_set<4>(gx, curve->gx);
	vli_set<4>(gy, curve->gy);
}


/* ------ Point operations ------ */

/* Returns true if p_point is the point at infinity, false otherwise. */
static forceinline
bool ecc_point_is_zero(const Point *p)
{
	if (/*p->isZero ||*/ vli_is_zero<4>(p->z)) return true;
	return (vli_is_zero<4>(p->x) && vli_is_zero<4>(p->y));
}

template<uint ndigits> forceinline
static void
vli_mod_mult_f(u64 *result, const u64 *left, const u64 *right,
				const ecc_curve *curve)
{
#ifdef	WITH_BARRETT
	if ( ! curve->use_barrett ) return;	// ERROR NOOP
	u64 product[2 * ECC_MAX_DIGITS];

	vli_mult<ndigits * 2>(product, left, right);
	vli_mmod_barrett<ndigits>(result, product, curve->p);
#else
	if ( curve->rr_p == nullptr ) return;	// ERROR NOOP
	u64		xp[ECC_MAX_DIGITS];
	u64		yp[ECC_MAX_DIGITS];
	u64		r[ECC_MAX_DIGITS];

	mont_mult<ndigits>(xp, left, curve->rr_p, curve->p, curve->k0_p);
	mont_mult<ndigits>(yp, right, curve->rr_p, curve->p, curve->k0_p);
	mont_mult<ndigits>(r, xp, yp, curve->p, curve->k0_p);
	mont_reduction<ndigits>(result, r, curve->p, curve->k0_p);
#endif
}

/* Computes result = left^2 % curve_prime. */
template<uint ndigits> forceinline
static void
vli_mod_square_f(u64 *result, const u64 *left, const ecc_curve *curve) noexcept
{
#ifdef	WITH_BARRETT
	if ( ! curve->use_barrett ) return;	// ERROR NOOP
	u64 product[2 * ECC_MAX_DIGITS];

	vli_square<ndigits>(product, left);
	vli_mmod_barrett<ndigits>(result, product, curve->p);
#else
	if ( curve->rr_p == nullptr ) return;	// ERROR NOOP
	u64		xp[ECC_MAX_DIGITS];
	u64		r[ECC_MAX_DIGITS];
	mont_mult<ndigits>(xp, left, curve->rr_p, curve->p, curve->k0_p);
	mont_mult<ndigits>(r, xp, xp, curve->p, curve->k0_p);
	mont_reduction<ndigits>(result, r, curve->p, curve->k0_p);
#endif
}

/* Point multiplication algorithm using Montgomery's ladder with co-Z
 * coordinates. From http://eprint.iacr.org/2011/338.pdf
 */

/*  RESULT = 2 * POINT  (Weierstrass version). */
/* Double in place */
template<uint ndigits> forceinline
static void ecc_point_double_jacobian(u64 *x3, u64 *y3, u64 *z3,
					u64 *x1, u64 *y1, u64 *z1, const ecc_curve *curve)
{
#define x3 (result->x)
#define y3 (result->y)
#define z3 (result->z)
#define t1 (ctx->t.scratch[0])
#define t2 (ctx->t.scratch[1])
#define t3 (ctx->t.scratch[2])
#define l1 (ctx->t.scratch[3])
#define l2 (ctx->t.scratch[4])
#define l3 (ctx->t.scratch[5])

	if (!mpi_cmp_ui(point->y, 0) || !mpi_cmp_ui(point->z, 0)) {
		/* P_y == 0 || P_z == 0 => [1:1:0] */
		mpi_set_ui(x3, 1);
		mpi_set_ui(y3, 1);
		mpi_set_ui(z3, 0);
	} else {
		if (ec_get_a_is_pminus3(ctx)) {
			/* Use the faster case.  */
			/* L1 = 3(X - Z^2)(X + Z^2) */
			/*                          T1: used for Z^2. */
			/*                          T2: used for the right term. */
			ec_pow2(t1, point->z, ctx);
			ec_subm(l1, point->x, t1, ctx);
			ec_mulm(l1, l1, mpi_const(MPI_C_THREE), ctx);
			ec_addm(t2, point->x, t1, ctx);
			ec_mulm(l1, l1, t2, ctx);
		} else {
			/* Standard case. */
			/* L1 = 3X^2 + aZ^4 */
			/*                          T1: used for aZ^4. */
			ec_pow2(l1, point->x, ctx);
			ec_mulm(l1, l1, mpi_const(MPI_C_THREE), ctx);
			ec_powm(t1, point->z, mpi_const(MPI_C_FOUR), ctx);
			ec_mulm(t1, t1, ctx->a, ctx);
			ec_addm(l1, l1, t1, ctx);
		}
		/* Z3 = 2YZ */
		ec_mulm(z3, point->y, point->z, ctx);
		ec_mul2(z3, z3, ctx);

		/* L2 = 4XY^2 */
		/*                              T2: used for Y2; required later. */
		ec_pow2(t2, point->y, ctx);
		ec_mulm(l2, t2, point->x, ctx);
		ec_mulm(l2, l2, mpi_const(MPI_C_FOUR), ctx);

		/* X3 = L1^2 - 2L2 */
		/*                              T1: used for L2^2. */
		ec_pow2(x3, l1, ctx);
		ec_mul2(t1, l2, ctx);
		ec_subm(x3, x3, t1, ctx);

		/* L3 = 8Y^4 */
		/*                              T2: taken from above. */
		ec_pow2(t2, t2, ctx);
		ec_mulm(l3, t2, mpi_const(MPI_C_EIGHT), ctx);

		/* Y3 = L1(L2 - X3) - L3 */
		ec_subm(y3, l2, x3, ctx);
		ec_mulm(y3, y3, l1, ctx);
		ec_subm(y3, y3, l3, ctx);
	}

#undef x3
#undef y3
#undef z3
#undef t1
#undef t2
#undef t3
#undef l1
#undef l2
#undef l3
}

/* Modify (x1, y1) => (x1 * z^2, y1 * z^3) */
template<uint ndigits> forceinline
static void apply_z(u64 *x1, u64 *y1, u64 *z, const ecc_curve *curve) noexcept
{
	u64 t1[ECC_MAX_DIGITS];

	vli_mod_square_f<ndigits>(t1, z, curve);    /* z^2 */
	vli_mod_mult_f<ndigits>(x1, x1, t1, curve); /* x1 * z^2 */
	vli_mod_mult_f<ndigits>(t1, t1, z, curve);  /* z^3 */
	vli_mod_mult_f<ndigits>(y1, y1, t1, curve); /* y1 * z^3 */
}


/* RESULT = P1 + P2  (Weierstrass version).*/
template<uint ndigits> forceinline
static void
point_add_jacobian(POINT result,
		POINT p1, POINT p2,
		const ecc_curve *curve) noexcept
{
#define x1 (p1->x)
#define y1 (p1->y)
#define z1 (p1->z)
#define x2 (p2->x)
#define y2 (p2->y)
#define z2 (p2->z)
#define x3 (result->x)
#define y3 (result->y)
#define z3 (result->z)
#define l1 (ctx->t.scratch[0])
#define l2 (ctx->t.scratch[1])
#define l3 (ctx->t.scratch[2])
#define l4 (ctx->t.scratch[3])
#define l5 (ctx->t.scratch[4])
#define l6 (ctx->t.scratch[5])
#define l7 (ctx->t.scratch[6])
#define l8 (ctx->t.scratch[7])
#define l9 (ctx->t.scratch[8])
#define t1 (ctx->t.scratch[9])
#define t2 (ctx->t.scratch[10])

	if ((!mpi_cmp(x1, x2)) && (!mpi_cmp(y1, y2)) && (!mpi_cmp(z1, z2))) {
		/* Same point; need to call the duplicate function.  */
		mpi_ec_dup_point(result, p1, ctx);
	} else if (!mpi_cmp_ui(z1, 0)) {
		/* P1 is at infinity.  */
		mpi_set(x3, p2->x);
		mpi_set(y3, p2->y);
		mpi_set(z3, p2->z);
	} else if (!mpi_cmp_ui(z2, 0)) {
		/* P2 is at infinity.  */
		mpi_set(x3, p1->x);
		mpi_set(y3, p1->y);
		mpi_set(z3, p1->z);
	} else {
		int z1_is_one = !mpi_cmp_ui(z1, 1);
		int z2_is_one = !mpi_cmp_ui(z2, 1);

		/* l1 = x1 z2^2  */
		/* l2 = x2 z1^2  */
		if (z2_is_one)
			mpi_set(l1, x1);
		else {
			ec_pow2(l1, z2, ctx);
			ec_mulm(l1, l1, x1, ctx);
		}
		if (z1_is_one)
			mpi_set(l2, x2);
		else {
			ec_pow2(l2, z1, ctx);
			ec_mulm(l2, l2, x2, ctx);
		}
		/* l3 = l1 - l2 */
		ec_subm(l3, l1, l2, ctx);
		/* l4 = y1 z2^3  */
		ec_powm(l4, z2, mpi_const(MPI_C_THREE), ctx);
		ec_mulm(l4, l4, y1, ctx);
		/* l5 = y2 z1^3  */
		ec_powm(l5, z1, mpi_const(MPI_C_THREE), ctx);
		ec_mulm(l5, l5, y2, ctx);
		/* l6 = l4 - l5  */
		ec_subm(l6, l4, l5, ctx);

		if (!mpi_cmp_ui(l3, 0)) {
			if (!mpi_cmp_ui(l6, 0)) {
				/* P1 and P2 are the same - use duplicate function. */
				mpi_ec_dup_point(result, p1, ctx);
			} else {
				/* P1 is the inverse of P2.  */
				mpi_set_ui(x3, 1);
				mpi_set_ui(y3, 1);
				mpi_set_ui(z3, 0);
			}
		} else {
			/* l7 = l1 + l2  */
			ec_addm(l7, l1, l2, ctx);
			/* l8 = l4 + l5  */
			ec_addm(l8, l4, l5, ctx);
			/* z3 = z1 z2 l3  */
			ec_mulm(z3, z1, z2, ctx);
			ec_mulm(z3, z3, l3, ctx);
			/* x3 = l6^2 - l7 l3^2  */
			ec_pow2(t1, l6, ctx);
			ec_pow2(t2, l3, ctx);
			ec_mulm(t2, t2, l7, ctx);
			ec_subm(x3, t1, t2, ctx);
			/* l9 = l7 l3^2 - 2 x3  */
			ec_mul2(t1, x3, ctx);
			ec_subm(l9, t2, t1, ctx);
			/* y3 = (l9 l6 - l8 l3^3)/2  */
			ec_mulm(l9, l9, l6, ctx);
			ec_powm(t1, l3, mpi_const(MPI_C_THREE), ctx); /* fixme: Use saved value*/
			ec_mulm(t1, t1, l8, ctx);
			ec_subm(y3, l9, t1, ctx);
			ec_mulm(y3, y3, ec_get_two_inv_p(ctx), ctx);
		}
	}

#undef x1
#undef y1
#undef z1
#undef x2
#undef y2
#undef z2
#undef x3
#undef y3
#undef z3
#undef l1
#undef l2
#undef l3
#undef l4
#undef l5
#undef l6
#undef l7
#undef l8
#undef l9
#undef t1
#undef t2
}

#ifdef	ommit
template<uint ndigits> forceinline
static void ecc_point_mult(u64 *result_x, u64 *result_y,
			   const u64 *point_x, const u64 *point_y, const u64 *scalar,
			   u64 *initial_z, const struct ecc_curve *curve) noexcept
{
	/* R0 and R1 */
	u64 rx[2][ECC_MAX_DIGITS];
	u64 ry[2][ECC_MAX_DIGITS];
	u64 z[ECC_MAX_DIGITS];
	u64 sk[2][ECC_MAX_DIGITS];
	const u64 *curve_prime = curve->p;
	int i, nb;
	int num_bits;
	int carry;

	carry = vli_add<ndigits>(sk[0], scalar, curve->n);
	vli_add<ndigits>(sk[1], sk[0], curve->n);
	scalar = sk[!carry];
	num_bits = sizeof(u64) * ndigits * 8 + 1;

	vli_set<ndigits>(rx[1], point_x);
	vli_set<ndigits>(ry[1], point_y);

	xycz_initial_double<ndigits>(rx[1], ry[1], rx[0], ry[0], initial_z,
					curve);

	for (i = num_bits - 2; i > 0; i--) {
		nb = !vli_test_bit(scalar, i);
		xycz_add_c<ndigits>(rx[1 - nb], ry[1 - nb], rx[nb], ry[nb], curve);
		xycz_add<ndigits>(rx[nb], ry[nb], rx[1 - nb], ry[1 - nb], curve);
	}

	nb = !vli_test_bit(scalar, 0);
	xycz_add_c<ndigits>(rx[1 - nb], ry[1 - nb], rx[nb], ry[nb], curve);

	/* Find final 1/Z value. */
	/* X1 - X0 */
	vli_mod_sub<ndigits>(z, rx[1], rx[0], curve_prime);
	/* Yb * (X1 - X0) */
	vli_mod_mult_f<ndigits>(z, z, ry[1 - nb], curve);
	/* xP * Yb * (X1 - X0) */
	vli_mod_mult_f<ndigits>(z, z, point_x, curve);

	/* 1 / (xP * Yb * (X1 - X0)) */
	vli_mod_inv<ndigits>(z, z, curve_prime);

	/* yP / (xP * Yb * (X1 - X0)) */
	vli_mod_mult_f<ndigits>(z, z, point_y, curve);
	/* Xb * yP / (xP * Yb * (X1 - X0)) */
	vli_mod_mult_f<ndigits>(z, z, rx[1 - nb], curve);
	/* End 1/Z calculation */

	xycz_add<ndigits>(rx[nb], ry[nb], rx[1 - nb], ry[1 - nb], curve);

	apply_z<ndigits>(rx[0], ry[0], z, curve);

	vli_set<ndigits>(result_x, rx[0]);
	vli_set<ndigits>(result_y, ry[0]);
}
#endif

/* Computes R = P + Q mod p */
template<uint ndigits> forceinline
static void ecc_point_add(Point *result, const Point *p, const Point *q,
		   const struct ecc_curve *curve)
{
	u64 z[ECC_MAX_DIGITS];
	u64 px[ECC_MAX_DIGITS];
	u64 py[ECC_MAX_DIGITS];

	if (ndigits != curve->ndigits) return;	// NOOOO
	vli_set<ndigits>((u64 *)result->x, (u64 *)q->x);
	vli_set<ndigits>((u64 *)result->y, (u64 *)q->y);
	vli_mod_sub<ndigits>(z, result->x, p->x, curve->p);
	vli_set<ndigits>(px, p->x);
	vli_set<ndigits>(py, p->y);
	xycz_add<ndigits>(px, py, result->x, result->y, curve);
	vli_mod_inv<ndigits>(z, z, curve->p);
	apply_z<ndigits>(result->x, result->y, z, curve);
}

void    point_double_jacobian(Point *pt_r, const Point *pt, CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	ecc_curve	*curve=(ecc_curve *)curveH;
	if (curve->name[0] == 0 || curve->ndigits != 4) return;
	vli_set<4>(pt_r->x, pt->x);
	vli_set<4>(pt_r->y, pt->y);
	vli_set<4>(pt_r->z, pt->z);
	ecc_point_double_jacobian<4>(pt_r->x, pt_r->y, pt_r->z, curve);
}

void    point_double(Point *pt, const Point *p, CURVE_HND curveH)
{
	u64 z[ECC_MAX_DIGITS];
	if (curveH == nullptr) return;
	ecc_curve	*curve=(ecc_curve *)curveH;
	if (curve->name[0] == 0 || curve->ndigits != 4) return;
	vli_set<4>(pt->x, p->x);
	vli_set<4>(pt->y, p->y);
	vli_set<4>(pt->z, p->z);
	ecc_point_double_jacobian<4>(pt->x, pt->y, pt->z, curve);
	vli_mod_inv<4>(z, pt->z, curve->p);
	apply_z<4>(pt->x, pt->y, z, curve);
}

void    point_add_jacobian(Point *pt_r, const Point *pt1, const Point *pt2,
                CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	ecc_curve	*curve=(ecc_curve *)curveH;
	if (curve->name[0] == 0 || curve->ndigits != 4) return;
}

void    point_add(Point *pt, const Point *p, const Point *q,
                CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	ecc_curve	*curve=(ecc_curve *)curveH;
	if (curve->name[0] == 0 || curve->ndigits != 4) return;
	ecc_point_add<4>(pt, p, q, curve);
}

void    affine_from_jacobian(u64 *x, u64 *y, const Point *pt, CURVE_HND curveH)
{
	u64 z[ECC_MAX_DIGITS];
	if (curveH == nullptr) return;
	ecc_curve	*curve=(ecc_curve *)curveH;
	if (curve->name[0] == 0 || curve->ndigits != 4) return;
	vli_mod_inv<4>(z, pt->z, curve->p);
	apply_z<4>(x, y, z, curve);
}

#ifdef	ommit
void	point_mult(Point *pt_r, const Point *pt, const u64 *scalar,
				CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	ecc_curve	*curve=(ecc_curve *)curveH;
	if (curve->name[0] == 0 || curve->ndigits != 4) return;
	ecc_point_mult<4>(pt_r->x, pt_r->y, pt->x, pt->y, scalar, nullptr, curve);
}
#endif
