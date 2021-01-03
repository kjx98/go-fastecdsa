// +build ncurve,!sm2p

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

using namespace vli;

static forceinline const ecc_curve<4> *ecc_get_curve(uint curve_id) noexcept
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
	auto	*curve=(ecc_curve<4> *)curveH;
	if ( !curve || curve->ndigits() != 4) return;
	curve->getP(p);
	curve->getN(n);
	curve->getB(b);
	curve->getGx(gx);
	curve->getGy(gy);
}


/* ------ Point operations ------ */

/* Point multiplication algorithm using Montgomery's ladder with co-Z
 * coordinates. From http://eprint.iacr.org/2011/338.pdf
 */

/*  RESULT = 2 * POINT  (Weierstrass version). */
/* Double in place */
template<const uint N> forceinline
static void
ecc_point_double_jacobian(u64 *x3, u64 *y3, u64 *z3, const u64 *x1,
			const u64 *y1, const u64 *z1, const ecc_curve<N> &curve) noexcept
{
	if (vli_is_zero<N>(y1) || vli_is_zero<N>(z1)) {
		/* P_y == 0 || P_z == 0 => [1:1:0] */
		vli_clear<N>(x3);
		vli_clear<N>(y3);
		vli_clear<N>(z3);
		x3[0] = 1;
		y3[0] = 1;
		return;
	}
	bool	z_is_one = vli_is_one<N>(z1);
	bignum<N>	t1, t2, l1, l2, l3;
	bignum<N>	xp, yp, zp;
	//bignum<N>	x3p, y3p, z3p;
	bignum<N>	*x3p = reinterpret_cast<bignum<N> *>(x3);
	bignum<N>	*y3p = reinterpret_cast<bignum<N> *>(y3);
	bignum<N>	*z3p = reinterpret_cast<bignum<N> *>(z3);
	curve.to_montgomery(xp, x1);
	curve.to_montgomery(yp, y1);
	if (z_is_one) {
		zp = curve.mont_one();
	} else {
		curve.to_montgomery(zp, z1);
	}
	if (curve.a_is_pminus3()) {
		/* Use the faster case.  */
		/* L1 = 3(X - Z^2)(X + Z^2) */
		/*                          T1: used for Z^2. */
		/*                          T2: used for the right term. */
		if (z_is_one) {
			// l1 = X - Z^2
			curve.mod_sub(l1, xp, curve.mont_one());
			// 3(X - Z^2) = 2(X - Z^2) + (X - Z^2)
			// t1 = 2 * l1
			curve.mont_mult2(t1, l1);
			// l1 = 2 * l1 + l1 = 3(X - Z^2)
			curve.mod_add_to(l1, t1);
			// t1 = X + Z^2
			curve.mod_add(t1, xp, curve.mont_one());
			// l1 = 3(X - Z^2)(X + Z^2)
			curve.mont_mult(l1, l1, t1);
		} else {
			// t2 = Z^2
			curve.mont_sqr(t2, zp);
			// l1 = X - Z^2
			curve.mod_sub(l1, xp, t2);
			// 3(X - Z^2) = 2(X - Z^2) + (X - Z^2)
			// t1 = 2 * l1
			curve.mont_mult2(t1, l1);
			// l1 = 2 * l1 + l1 = 3(X - Z^2)
			curve.mod_add_to(l1, t1);
			// t1 = X + Z^2
			curve.mod_add(t1, xp, t2);
			// l1 = 3(X - Z^2)(X + Z^2)
			curve.mont_mult(l1, l1, t1);
		}
	} else {
		/* Standard case. */
		/* L1 = 3X^2 + aZ^4 */
		/*                          T1: used for aZ^4. */
		// l1 = X^2
		curve.mont_sqr(l1, xp);
		curve.mont_mult2(t1, l1);
		curve.mod_add_to(l1, t1);	// l1 = 3X^2
		if (z_is_one) {
			curve.mod_add_to(l1, curve.paramA());
		} else {
			// t1 = Z^4
			curve.mont_sqr(t1, zp, 2);
			curve.mont_mult(t1, t1, curve.paramA());
			curve.mod_add_to(l1, t1);
		}
	}

	/* Z3 = 2YZ */
	if (z_is_one) {
		curve.mont_mult2(*z3p, yp);
	} else {
		// Z3 = YZ
		curve.mont_mult(*z3p, yp, zp);
		// Z3 *= 2
		curve.mont_mult2(*z3p, *z3p);
	}

	/* L2 = 4XY^2 */
	/*                              T2: used for Y2; required later. */
	curve.mont_sqr(t2, yp);	// t2 = Y^2
	// l2 = XY^2
	curve.mont_mult(l2, t2, xp);
	curve.mont_mult2(l2, l2);
	curve.mont_mult2(l2, l2);

	/* X3 = L1^2 - 2L2 */
	/*                              T1: used for 2L2. */
	curve.mont_sqr(*x3p, l1);
	curve.mont_mult2(t1, l2);
	curve.mod_sub_from(*x3p, t1);

	/* L3 = 8Y^4 */
	/*                              T2: taken from above. */
	curve.mont_sqr(t2, t2);
	curve.mont_mult2(l3, t2);
	curve.mont_mult2(l3, l3);
	curve.mont_mult2(l3, l3);

	/* Y3 = L1(L2 - X3) - L3 */
	curve.mod_sub(*y3p, l2, *x3p);
	curve.mont_mult(*y3p, *y3p, l1);
	curve.mod_sub_from(*y3p, l3);

	// montgomery reduction
	curve.from_montgomery(x3, *x3p);
	curve.from_montgomery(y3, *y3p);
	curve.from_montgomery(z3, *z3p);
}

/* Modify (x1, y1) => (x1 * z^2, y1 * z^3) */
template<const uint N> forceinline
static void
apply_z(u64 *x1, u64 *y1, u64 *z, const ecc_curve<N> &curve) noexcept
{
	bignum<N>	t1;

	if (vli_is_one<N>(z) || vli_is_zero<N>(z)) return; 
	bignum<N>	*xp = reinterpret_cast<bignum<N> *>(x1);
	bignum<N>	*yp = reinterpret_cast<bignum<N> *>(y1);
	bignum<N>	*zp = reinterpret_cast<bignum<N> *>(z);
	curve.to_montgomery(*zp, z);
	curve.mont_sqr(t1, *zp);	// t1 = z^2
	curve.to_montgomery(*xp, x1);
	curve.to_montgomery(*yp, y1);
	curve.mont_mult(*xp, *xp, t1);	// x1 * z^2
	curve.mont_mult(t1, t1, *zp);	// t1 = z^3
	curve.mont_mult(*yp, *yp, t1);	// y1 * z^3

	// montgomery reduction
	curve.from_montgomery(x1, *xp);
	curve.from_montgomery(y1, *yp);
}


/* RESULT = P1 + P2  (Weierstrass version).*/
template<const uint N> forceinline
static void
ecc_point_add_jacobian(u64 *x3, u64 *y3, u64 *z3, const u64 *x1,
			const u64 *y1, const u64 *z1, const u64 *x2,
			const u64 *y2, const u64 *z2, const ecc_curve<N> &curve) noexcept
{

	if (vli_cmp<N>(x1,x2) == 0 && vli_cmp<N>(y1, y2) == 0 &&
			vli_cmp<N>(z1, z2) == 0) {
		ecc_point_double_jacobian(x3, y3, z3, x1, y1, z1, curve);
		return;
	}
	if (vli_is_zero<N>(z1)) {
		vli_set<N>(x3, x2);
		vli_set<N>(y3, y2);
		vli_set<N>(z3, z2);
		return;
	} else if (vli_is_zero<N>(z2)) {
		vli_set<N>(x3, x1);
		vli_set<N>(y3, y1);
		vli_set<N>(z3, z1);
		return;
	}
	bool	z1_is_one = vli_is_one<N>(z1);
	bool	z2_is_one = vli_is_one<N>(z2);
	bignum<N>	l1, l2, l3, l4, l5, l6, l7, l8;
	bignum<N>	z1z1, z1p;
	bignum<N>	z2z2, z2p;
	bignum<N>	*x3p = reinterpret_cast<bignum<N> *>(x3);
	bignum<N>	*y3p = reinterpret_cast<bignum<N> *>(y3);
	bignum<N>	*z3p = reinterpret_cast<bignum<N> *>(z3);
	//curve.to_montgomery(x1p, x1);
	//curve.to_montgomery(y1p, y1);

	/* l1 = x1 z2^2  */
	/* l2 = x2 z1^2  */
	if (z2_is_one) {
		curve.to_montgomery(l1, x1);
	} else {
		curve.to_montgomery(z2p, z2);
		// z2z2 = z2^2
		curve.mont_sqr(z2z2, z2p);
		curve.to_montgomery(l1, x1);
		// l1 = x1 z2^2
		curve.mont_mult(l1, l1, z2z2);
	}
	if (z1_is_one) {
		curve.to_montgomery(l2, x2);
	} else {
		curve.to_montgomery(z1p, z1);
		// z1z1 = z1^2
		curve.mont_sqr(z1z1, z1p);
		curve.to_montgomery(l2, x2);
		// l2 = x2 z1^2
		curve.mont_mult(l2, l2, z1z1);
	}

	/* l3 = l1 - l2 */
	curve.mod_sub(l3, l1, l2);
	/* l4 = y1 z2^3  */
	// z2z2 = z2^3
	curve.mont_mult(z2z2, z2z2, z2p);
	// l4 = y1
	curve.to_montgomery(l4, y1);
	// l4 = y1 z2^3
	curve.mont_mult(l4, l4, z2z2);

	/* l5 = y2 z1^3  */
	// z1z1 = z1^3
	curve.mont_mult(z1z1, z1z1, z1p);
	// l5 = y2
	curve.to_montgomery(l5, y2);
	// l5 = y2 z1^3
	curve.mont_mult(l5, l5, z1z1);
	/* l6 = l4 - l5  */
	curve.mod_sub(l6, l4, l5);

	if (l3.is_zero()) {
		if (l6.is_zero()) {
			/* P1 and P2 are the same - use duplicate function. */
			ecc_point_double_jacobian(x3, y3, z3, x1, y1, z1, curve);
			return;
		}
		/* P1 is the inverse of P2.  */
		vli_clear<N>(x3);
		vli_clear<N>(y3);
		vli_clear<N>(z3);
		x3[0] = 1;
		y3[0] = 1;
		return;
	}
	/* l7 = l1 + l2  */
	curve.mod_add(l7, l1, l2);	// no use l1, l2 following
	/* l8 = l4 + l5  */
	curve.mod_add(l8, l4, l5);
	/* z3 = z1 z2 l3  */
	curve.mont_mult(*z3p, z1p, z2p);
	curve.mont_mult(*z3p, *z3p, l3);
	/* x3 = l6^2 - l7 l3^2  */
	// t1 = l3^2, reuse l1 for t1
	curve.mont_sqr(l1, l3);
	// x3 = l6^2
	curve.mont_sqr(*x3p, l6);
	// t2 = l7 l3^2, reuse l2 for t2
	curve.mont_mult(l2, l7, l1);
	// x3 = l6^2 - l7 l3^2
	curve.mod_sub_from(*x3p, l2);
	/* l9 = l7 l3^2 - 2 x3  */
	// t2 = 2 x3, reuse l2 for t2
	curve.mont_mult2(l2, *x3p);
	// l9 = l7 l3^2 = l7 t1, reuse l4 for l9
	curve.mont_mult(l4, l7, l1);
	//  l9 = l7 l3^2 - 2 x3 = l9 - t2
	curve.mod_sub_from(l4, l2);
	/* y3 = (l9 l6 - l8 l3^3)/2  */
	curve.mont_mult(l4, l4, l6);	// l9 = l9 l6
	curve.mont_mult(l1, l1, l3);	// t1 = l3^3
	curve.mont_mult(l2, l8, l1);	// t2 = l8 l3^3
	curve.mod_sub_from(l4, l2);	// l9 = l9 l6 - l8 l3^3
	curve.mont_mult(*y3p, l4, curve.mont_inv2());
	// montgomery reduction
	curve.from_montgomery(x3, *x3p);
	curve.from_montgomery(y3, *y3p);
	curve.from_montgomery(z3, *z3p);
}


void    point_double_jacobian(Point *pt, const Point *p, CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(ecc_curve<4> *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	ecc_point_double_jacobian<4>(pt->x, pt->y, pt->z, p->x, p->y, p->z, *curve);
}

// p->z MUST be one
void    point_double(Point *pt, const Point *p, CURVE_HND curveH)
{
	u64 z[4];
	if (curveH == nullptr) return;
	auto	*curve=(ecc_curve<4> *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	ecc_point_double_jacobian<4>(pt->x, pt->y, pt->z, p->x, p->y, p->z, *curve);
	vli_mod_inv<4>(z, pt->z, curve->paramP().data());
	apply_z<4>(pt->x, pt->y, z, *curve);
}

void    point_add_jacobian(Point *pt, const Point *p, const Point *q,
                CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(ecc_curve<4> *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	ecc_point_add_jacobian<4>(pt->x, pt->y, pt->z, p->x, p->y, p->z,
				q->x, q->y, q->z, *curve);
}

// p->z and q->z MUST be one
void    point_add(Point *pt, const Point *p, const Point *q,
                CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(ecc_curve<4> *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	ecc_point_add_jacobian<4>(pt->x, pt->y, pt->z, p->x, p->y, p->z,
				q->x, q->y, q->z, *curve);
	u64	z[4];
	vli_mod_inv<4>(z, pt->z, curve->paramP().data());
	apply_z<4>(pt->x, pt->y, z, *curve);
}

void    affine_from_jacobian(u64 *x, u64 *y, const Point *pt, CURVE_HND curveH)
{
	u64 z[4];
	if (curveH == nullptr) return;
	auto	*curve=(ecc_curve<4> *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	vli_mod_inv<4>(z, pt->z, curve->paramP().data());
	vli_set<4>(x, pt->x);
	vli_set<4>(y, pt->y);
	apply_z<4>(x, y, z, *curve);
}

#ifdef	ommit
void	point_mult(Point *pt, const Point *p, const u64 *scalar,
				CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	ecc_curve	*curve=(ecc_curve *)curveH;
	if (curve->name[0] == 0 || curve->ndigits != 4) return;
	ecc_point_mult<4>(pt->x, pt->y, p->x, p->y, scalar, nullptr, *curve);
}
#endif
