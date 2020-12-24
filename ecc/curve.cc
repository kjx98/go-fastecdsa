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
#include "vli.hpp"
#include "ecc.h"
#include "ecc_impl.hpp"
#include "ecc_curve_defs.h"

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


void mont_sm2_mod_mult_p(u64 *result, const u64 *x, const u64 *y)
{
	u64	xp[ECC_MAX_DIGITS];
	u64	yp[ECC_MAX_DIGITS];
	u64	r[ECC_MAX_DIGITS];
	auto rr = sm2_p256.rr_p;
	auto prime = sm2_p256.p;
	mont_multP(xp, x, rr, prime);
	mont_multP(yp, y, rr, prime);
	mont_multP(r, xp, yp, prime);
	mont_reductionP(result, r, prime);
}

void mont_sm2_mod_mult_n(u64 *result, const u64 *x, const u64 *y)
{
	u64	xp[ECC_MAX_DIGITS];
	u64	yp[ECC_MAX_DIGITS];
	u64	r[ECC_MAX_DIGITS];
	auto rr = sm2_p256.rr_n;
	auto prime = sm2_p256.n;
	auto k0 = sm2_p256.k0_n;
	mont_mult<4>(xp, x, rr, prime, k0);
	mont_mult<4>(yp, y, rr, prime, k0);
	mont_mult<4>(r, xp, yp, prime, k0);
	mont_reduction<4>(result, r, prime, k0);
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

/* Double in place */
template<uint ndigits> forceinline
static void ecc_point_double_jacobian(u64 *x1, u64 *y1, u64 *z1,
				      const ecc_curve *curve)
{
	/* t1 = x, t2 = y, t3 = z */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
	u64 t4[ECC_MAX_DIGITS];
	u64 t5[ECC_MAX_DIGITS];
	const u64 *curve_prime = curve->p;

	if (vli_is_zero<ndigits>(z1))
		return;

	/* t4 = y1^2 */
	vli_mod_square_f<ndigits>(t4, y1, curve);
	/* t5 = x1*y1^2 = A */
	vli_mod_mult_f<ndigits>(t5, x1, t4, curve);
	/* t4 = y1^4 */
	vli_mod_square_f<ndigits>(t4, t4, curve);
	/* t2 = y1*z1 = z3 */
	vli_mod_mult_f<ndigits>(y1, y1, z1, curve);
	/* t3 = z1^2 */
	vli_mod_square_f<ndigits>(z1, z1, curve);

	/* t1 = x1 + z1^2 */
	vli_mod_add<ndigits>(x1, x1, z1, curve_prime);
	/* t3 = 2*z1^2 */
	vli_mod_add<ndigits>(z1, z1, z1, curve_prime);
	/* t3 = x1 - z1^2 */
	vli_mod_sub<ndigits>(z1, x1, z1, curve_prime);
	/* t1 = x1^2 - z1^4 */
	vli_mod_mult_f<ndigits>(x1, x1, z1, curve);

	/* t3 = 2*(x1^2 - z1^4) */
	vli_mod_add<ndigits>(z1, x1, x1, curve_prime);
	/* t1 = 3*(x1^2 - z1^4) */
	vli_mod_add<ndigits>(x1, x1, z1, curve_prime);
	if (vli_test_bit(x1, 0)) {
		u64 carry = vli_add_to<ndigits>(x1, curve_prime);

		vli_rshift1<ndigits>(x1);
		x1[ndigits - 1] |= carry << 63;
	} else {
		vli_rshift1<ndigits>(x1);
	}
	/* t1 = 3/2*(x1^2 - z1^4) = B */

	/* t3 = B^2 */
	vli_mod_square_f<ndigits>(z1, x1, curve);
	/* t3 = B^2 - A */
	vli_mod_sub<ndigits>(z1, z1, t5, curve_prime);
	/* t3 = B^2 - 2A = x3 */
	vli_mod_sub<ndigits>(z1, z1, t5, curve_prime);
	/* t5 = A - x3 */
	vli_mod_sub<ndigits>(t5, t5, z1, curve_prime);
	/* t1 = B * (A - x3) */
	vli_mod_mult_f<ndigits>(x1, x1, t5, curve);
	/* t4 = B * (A - x3) - y1^4 = y3 */
	vli_mod_sub<ndigits>(t4, x1, t4, curve_prime);

	vli_set<ndigits>(x1, z1);
	vli_set<ndigits>(z1, y1);
	vli_set<ndigits>(y1, t4);
#pragma GCC diagnostic push
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

/* P = (x1, y1) => 2P, (x2, y2) => P' */
template<uint ndigits> forceinline
static void xycz_initial_double(u64 *x1, u64 *y1, u64 *x2, u64 *y2,
				u64 *p_initial_z, const ecc_curve *curve) noexcept
{
	u64 z[ECC_MAX_DIGITS];

	vli_set<ndigits>(x2, x1);
	vli_set<ndigits>(y2, y1);

	vli_clear<ndigits>(z);
	z[0] = 1;

	if (p_initial_z) {
		vli_set<ndigits>(z, p_initial_z);
		apply_z<ndigits>(x1, y1, z, curve);
	}

	ecc_point_double_jacobian<ndigits>(x1, y1, z, curve);

	apply_z<ndigits>(x2, y2, z, curve);
}

/* Input P = (x1, y1, Z), Q = (x2, y2, Z)
 * Output P' = (x1', y1', Z3), P + Q = (x3, y3, Z3)
 * or P => P', Q => P + Q
 */
template<uint ndigits> forceinline
static void
xycz_add(u64 *x1, u64 *y1, u64 *x2, u64 *y2, const ecc_curve *curve) noexcept
{
	/* t1 = X1, t2 = Y1, t3 = X2, t4 = Y2 */
	u64 t5[ECC_MAX_DIGITS];
	const u64 *curve_prime = curve->p;

	/* t5 = x2 - x1 */
	vli_mod_sub<ndigits>(t5, x2, x1, curve_prime);
	/* t5 = (x2 - x1)^2 = A */
	vli_mod_square_f<ndigits>(t5, t5, curve);
	/* t1 = x1*A = B */
	vli_mod_mult_f<ndigits>(x1, x1, t5, curve);
	/* t3 = x2*A = C */
	vli_mod_mult_f<ndigits>(x2, x2, t5, curve);
	/* t4 = y2 - y1 */
	vli_mod_sub<ndigits>(y2, y2, y1, curve_prime);
	/* t5 = (y2 - y1)^2 = D */
	vli_mod_square_f<ndigits>(t5, y2, curve);

	/* t5 = D - B */
	vli_mod_sub<ndigits>(t5, t5, x1, curve_prime);
	/* t5 = D - B - C = x3 */
	vli_mod_sub<ndigits>(t5, t5, x2, curve_prime);
	/* t3 = C - B */
	vli_mod_sub<ndigits>(x2, x2, x1, curve_prime);
	/* t2 = y1*(C - B) */
	vli_mod_mult_f<ndigits>(y1, y1, x2, curve);
	/* t3 = B - x3 */
	vli_mod_sub<ndigits>(x2, x1, t5, curve_prime);
	/* t4 = (y2 - y1)*(B - x3) */
	vli_mod_mult_f<ndigits>(y2, y2, x2, curve);
	/* t4 = y3 */
	vli_mod_sub<ndigits>(y2, y2, y1, curve_prime);

	vli_set<ndigits>(x2, t5);
}

/* Input P = (x1, y1, Z), Q = (x2, y2, Z)
 * Output P + Q = (x3, y3, Z3), P - Q = (x3', y3', Z3)
 * or P => P - Q, Q => P + Q
 */
template<uint ndigits> forceinline
static void xycz_add_c(u64 *x1, u64 *y1, u64 *x2, u64 *y2,
				const ecc_curve *curve) noexcept
{
	/* t1 = X1, t2 = Y1, t3 = X2, t4 = Y2 */
	u64 t5[ECC_MAX_DIGITS];
	u64 t6[ECC_MAX_DIGITS];
	u64 t7[ECC_MAX_DIGITS];
	const u64 *curve_prime = curve->p;

	/* t5 = x2 - x1 */
	vli_mod_sub<ndigits>(t5, x2, x1, curve_prime);
	/* t5 = (x2 - x1)^2 = A */
	vli_mod_square_f<ndigits>(t5, t5, curve);
	/* t1 = x1*A = B */
	vli_mod_mult_f<ndigits>(x1, x1, t5, curve);
	/* t3 = x2*A = C */
	vli_mod_mult_f<ndigits>(x2, x2, t5, curve);
	/* t4 = y2 + y1 */
	vli_mod_add<ndigits>(t5, y2, y1, curve_prime);
	/* t4 = y2 - y1 */
	vli_mod_sub<ndigits>(y2, y2, y1, curve_prime);

	/* t6 = C - B */
	vli_mod_sub<ndigits>(t6, x2, x1, curve_prime);
	/* t2 = y1 * (C - B) */
	vli_mod_mult_f<ndigits>(y1, y1, t6, curve);
	/* t6 = B + C */
	vli_mod_add<ndigits>(t6, x1, x2, curve_prime);
	/* t3 = (y2 - y1)^2 */
	vli_mod_square_f<ndigits>(x2, y2, curve);
	/* t3 = x3 */
	vli_mod_sub<ndigits>(x2, x2, t6, curve_prime);

	/* t7 = B - x3 */
	vli_mod_sub<ndigits>(t7, x1, x2, curve_prime);
	/* t4 = (y2 - y1)*(B - x3) */
	vli_mod_mult_f<ndigits>(y2, y2, t7, curve);
	/* t4 = y3 */
	vli_mod_sub<ndigits>(y2, y2, y1, curve_prime);

	/* t7 = (y2 + y1)^2 = F */
	vli_mod_square_f<ndigits>(t7, t5, curve);
	/* t7 = x3' */
	vli_mod_sub<ndigits>(t7, t7, t6, curve_prime);
	/* t6 = x3' - B */
	vli_mod_sub<ndigits>(t6, t7, x1, curve_prime);
	/* t6 = (y2 + y1)*(x3' - B) */
	vli_mod_mult_f<ndigits>(t6, t6, t5, curve);
	/* t2 = y3' */
	vli_mod_sub<ndigits>(y1, t6, y1, curve_prime);

	vli_set<ndigits>(x1, t7);
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
