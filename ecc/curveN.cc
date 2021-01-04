// +build curve sm2p

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
		if (sm2_p256.init()) return &sm2_p256;
		return nullptr;
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
	if ( !(*curve) || curve->ndigits() != 4) return;
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
#ifndef	ommit
template<const uint N> forceinline
static void
ecc_point_double_jacobian(u64 *x3, u64 *y3, u64 *z3, const u64 *x1,
			const u64 *y1, const u64 *z1, const ecc_curve<N> &curve) noexcept
{
/* dbl-1998-cmo-2 algorithm
 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html
 * M ... l1
 * S ... l2
 * T ... x3
 */
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
	bignum<N>	t1, t2, l1, l2;
	bignum<N>	xp, yp, zp;
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
			// t1 = Z^2
			curve.mont_sqr(t1, zp);
			// l1 = X - Z^2
			curve.mod_sub(l1, xp, t1);
			// 3(X - Z^2) = 2(X - Z^2) + (X - Z^2)
			// t2 = 2 * l1
			curve.mont_mult2(t2, l1);
			// l1 = l1 + 2 * l1 = 3(X - Z^2)
			curve.mod_add_to(l1, t2);
			// t2 = X + Z^2
			curve.mod_add(t2, xp, t1);
			// l1 = 3(X - Z^2)(X + Z^2)
			curve.mont_mult(l1, l1, t2);
		}
	} else {
		/* Standard case. */
		/* L1 = 3X^2 + aZ^4 */
		/*                          T1: used for aZ^4. */
		// l1 = X^2
		curve.mont_sqr(l1, xp);
		curve.mont_mult2(t1, l1);
		// l1 = 3X^2
		curve.mod_add_to(l1, t1);
		if (curve.a_is_zero()) {
			/* Use the faster case.  */
			/* L1 = 3X^2 */
			// do nothing
		} else if (z_is_one) {
			// should be mont_paramA
			curve.mod_add_to(l1, curve.montParamA());
		} else {
			// t1 = Z^4
			curve.mont_sqr(t1, zp, 2);
			// t1 = a * Z^4
			curve.mont_mult(t1, t1, curve.montParamA());
			// l1 = 3 X^2 + a Z^4
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
	/* t2 = Y1^2 */
	curve.mont_sqr(t2, yp);
	// t2 = 2 Y^2
	curve.mont_mult2(t2, t2);
	// l2 =  2 XY^2
	curve.mont_mult(l2, t2, xp);
	// l2 = 4 X Y^2
	curve.mont_mult2(l2, l2);

	/* X3 = L1^2 - 2L2 */
	/*                              T1: used for 2L2. */
	curve.mont_sqr(*x3p, l1);
	curve.mont_mult2(t1, l2);
	curve.mod_sub_from(*x3p, t1);

	/* L3 = 8Y^4 */
	/*   L3 reuse t2, t2: taken from above. */
	curve.mont_sqr(t2, t2);		// t2 = t2^2 = 4Y^4
	curve.mont_mult2(t2, t2);	// t2 *= 2, t2 = 8Y^4

	/* Y3 = L1(L2 - X3) - L3 */
	curve.mod_sub(*y3p, l2, *x3p);
	curve.mont_mult(*y3p, l1, *y3p);
	curve.mod_sub_from(*y3p, t2);

	// montgomery reduction
	curve.from_montgomery(x3, *x3p);
	curve.from_montgomery(y3, *y3p);
	curve.from_montgomery(z3, *z3p);
}
#else
template<const uint N> forceinline
static void
ecc_point_double_jacobian(u64 *x3, u64 *y3, u64 *z3, const u64 *x1,
			const u64 *y1, const u64 *z1, const ecc_curve<N> &curve) noexcept
{
/* dbl-2001-b algorithm
 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html
 * only for a = -3
 */
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
	bignum<N>	alpha, beta, delta, gamma;
	bignum<N>	xp, yp, zp, t1;
	bignum<N>	*x3p = reinterpret_cast<bignum<N> *>(x3);
	bignum<N>	*y3p = reinterpret_cast<bignum<N> *>(y3);
	bignum<N>	*z3p = reinterpret_cast<bignum<N> *>(z3);
	// delta = z1^2
	if (z_is_one) {
		zp = curve.mont_one();
		delta = curve.mont_one();
	} else {
		curve.to_montgomery(zp, z1);
		curve.mont_sqr(delta, zp);
	}
	// gamma = y1^2
	curve.to_montgomery(yp, y1);
	curve.mont_sqr(gamma, yp);
   	curve.to_montgomery(xp, x1);
	// alpha = 3 (x1 - delta) (x1 + delta)
	// alpha = x1 - delta
	// t1 = x1 + delta
	curve.mod_sub(alpha, xp, delta);
	curve.mod_add(t1, xp, delta);
	// t1 = (x1 - delta) (x1 + delta)
	curve.mont_mult(t1, alpha, t1);
	alpha = t1;
	// t1 *= 2
	curve.mont_mult2(t1, t1);
	curve.mod_add_to(alpha, t1);

	// beta = x gamma
	curve.mont_mult(beta, xp, gamma);
	curve.mont_mult2(beta, beta);
	// beta = 4 beta = 4 x gamma
	curve.mont_mult2(beta, beta);
	// x3 = alpha^2 - 8 * beta
	curve.mont_sqr(*x3p, alpha);
	// t1 = 8 beta
	curve.mont_mult2(t1, beta);

	// x3 = alpha^2 - 8 beta
	curve.mod_sub_from(*x3p, t1);

	// z3 = (y1 + z1)^2 - gamma - delta
	curve.mod_add(t1, yp, zp);
	curve.mont_sqr(*z3p, t1);
	curve.mod_sub_from(*z3p, gamma);
	curve.mod_sub_from(*z3p, delta);

	// y3 = alpha ( 4 * beta - x3)
	// beta = 4 beta - x3
	curve.mod_sub_from(beta, *x3p);
	curve.mont_mult(*y3p, alpha, beta);

	// gamma = 8 gamma^2
	// t1 = 2 * gamma
	curve.mont_mult2(t1, gamma);
	// gamma = t1 ^ 2 = 4 gamma^2
	curve.mont_sqr(gamma, t1);
	curve.mont_mult2(gamma, gamma);
	// y3 = alpha ( 4 beta - x3) - 8 gamma^2
	curve.mod_sub_from(*y3p, gamma);


	// montgomery reduction
	curve.from_montgomery(x3, *x3p);
	curve.from_montgomery(y3, *y3p);
	curve.from_montgomery(z3, *z3p);
}
#endif


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

	bignum<N>	*x3p = reinterpret_cast<bignum<N> *>(x3);
	bignum<N>	*y3p = reinterpret_cast<bignum<N> *>(y3);
	bignum<N>	*z3p = reinterpret_cast<bignum<N> *>(z3);
#ifdef	WITH_Z1_Z2_EQ
	if ( vli_cmp<N>(z1, z2) == 0 ) {
		/* zadd-2007-m algorithm
		 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html
		 */
		if (vli_cmp<N>(x1,x2) == 0) {
			if (vli_cmp<N>(y1, y2) == 0) {
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
		bignum<N>		x1p, x2p, y1p, aa, bb,cc, dd;
		curve.to_montgomery(x1p, x1);
		curve.to_montgomery(x2p, x2);
		curve.to_montgomery(y1p, y1);
		bool	z1_is_one = vli_is_one<N>(z1);
		bignum<N>		t1;
		// t1 = x2 - x1
		curve.mod_sub(t1, x2p, x1p);
		// z3 = z1 * (x2 -x1)
		if (z1_is_one) {
			*z3p = t1;
		} else {
			curve.to_montgomery(*z3p, z1);
			curve.mont_mult(*z3p, *z3p, t1);
		}
		// aa = (x2 - x1)^2
		curve.mont_sqr(aa, t1);
		// bb = x1 * aa
		curve.mont_mult(bb, x1p, aa);
		// cc = x2 * aa
		curve.mont_mult(cc, x2p, aa);
		curve.to_montgomery(t1, y2);
		// t1 = y2 -y1
		curve.mod_sub_from(t1, y1p);
		// dd = (y2 -y1)^2
		curve.mont_sqr(dd, t1);
		// x3 = dd - bb - cc
		curve.mod_sub(*x3p, dd, bb);
		curve.mod_sub_from(*x3p, cc);
		// y3 = bb - x3
		curve.mod_sub(*y3p, bb, *x3p);
		// y3 = t1 * (bb - x3) = (y2 - y1) * (bb - x3)
		curve.mont_mult(*y3p, t1, *y3p);
		// t1 = c - b
		curve.mod_sub(t1, cc, bb);
		// t1 = y1 * (c - b)
		curve.mont_mult(t1, y1p, t1);
		// y3 = y3 - t1 = (y2 - y1)(bb - x3) - y1(c-b)
		curve.mod_sub_from(*y3p, t1);
	} else
#endif
	{
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
		// add-2007-bl
		/* add-2007-bl algorithm
		 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html
		 */
		// U1 ... l1
		// U2 ... l2
		// S1 ... l4
		// S2 ... l5
		bool	z1_is_one = vli_is_one<N>(z1);
		bool	z2_is_one = vli_is_one<N>(z2);
		bignum<N>	u1, u2, s1, s2, h, i, j, r, v;
		bignum<N>	t1;
		bignum<N>	z1z1, z1p;
		bignum<N>	z2z2, z2p;

		/* u1 = x1 z2^2  */
		/* u2 = x2 z1^2  */
		if (z2_is_one) {
			curve.to_montgomery(u1, x1);
		} else {
			curve.to_montgomery(z2p, z2);
			// z2z2 = z2^2
			curve.mont_sqr(z2z2, z2p);
			curve.to_montgomery(u1, x1);
			// u1 = x1 z2^2
			curve.mont_mult(u1, u1, z2z2);
		}
		if (z1_is_one) {
			curve.to_montgomery(u2, x2);
		} else {
			curve.to_montgomery(z1p, z1);
			// z1z1 = z1^2
			curve.mont_sqr(z1z1, z1p);
			curve.to_montgomery(u2, x2);
			// u2 = x2 z1^2
			curve.mont_mult(u2, u2, z1z1);
		}

		/* h = u2 - u1 */
		curve.mod_sub(h, u2, u1);
		/* s1 = y1 z2^3  */
		// s1 = y1
		curve.to_montgomery(s1, y1);
		if ( ! z2_is_one ) {
			// s1 = y1 z2^2
			curve.mont_mult(s1, s1, z2z2);
			curve.mont_mult(s1, s1, z2p);
		}

		/* s2 = y2 z1^3  */
		// s2 = y2
		curve.to_montgomery(s2, y2);
		if ( !z1_is_one ) {
			// s2 = y2 z1^3
			curve.mont_mult(s2, s2, z1z1);
			curve.mont_mult(s2, s2, z1p);
		}
		/* r = s2 - s1  */
		curve.mod_sub(r, s2, s1);

		if (h.is_zero()) {
			if (r.is_zero()) {
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
		// r = 2 * (s2 -s1)
		curve.mont_mult2(r, r);
		// i = (2*h)^2
		curve.mont_mult2(i, h);
		curve.mont_sqr(i, i);

		// j = h * i
		curve.mont_mult(j, h, i);
		// v = u1 * i
		curve.mont_mult(v, u1, i);

		// x3 = r^2 - j - 2*v
		curve.mont_sqr(*x3p, r);
		curve.mod_sub_from(*x3p, j);
		curve.mod_sub_from(*x3p, v);
		curve.mod_sub_from(*x3p, v);

		// y3 = v - x3
		curve.mod_sub(*y3p, v, *x3p);
		// y3 = r * (v - x3)
		curve.mont_mult(*y3p, r, *y3p);
		// t1 = 2 * s1 * j
		curve.mont_mult(t1, s1, j);
		curve.mont_mult2(t1, t1);
		// y3 = r * (v - x3) - 2 * s1 *j
		curve.mod_sub_from(*y3p, t1);

		if (z1_is_one && z2_is_one) {
			curve.mont_mult2(*z3p, h);
		} else {
			if (z1_is_one) {
				z1p = curve.mont_one();
				z1z1 = curve.mont_one();
			}
			if (z2_is_one) {
				z2p = curve.mont_one();
				z2z2 = curve.mont_one();
			}
			curve.mod_add(t1, z1p, z2p);
			curve.mont_sqr(t1, t1);
			curve.mod_sub_from(t1, z1z1);
			curve.mod_sub_from(t1, z2z2);
			curve.mont_mult(*z3p, t1, h);
		}
	}
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
	if (curveH == nullptr) return;
	auto	*curve=(ecc_curve<4> *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	ecc_point_double_jacobian<4>(pt->x, pt->y, pt->z, p->x, p->y, p->z, *curve);
	u64 z[4];
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
