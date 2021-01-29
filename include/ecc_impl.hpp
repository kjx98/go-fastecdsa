/*
 * Copyright (c) 2013, Kenneth MacKay
 * All rights reserved.
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
#pragma once
#ifndef	__ECC_IMPL_H__
#define __ECC_IMPL_H__

#include "cdefs.h"
#include <string>
#include <functional>
#include "vli.hpp"
#include "vli_bn.hpp"


namespace ecc {

using namespace vli;

// PublicKey struct also spoint_t
template<const uint N>
struct spoint_t {
	bignum<N>	x;
	bignum<N>	y;
	spoint_t() = default;
	spoint_t(const spoint_t &) = default;
	spoint_t(const bignum<N>& xx, const bignum<N>& yy) : x(xx), y(yy) {}
	explicit spoint_t(const u64 *xx, const u64 *yy) : x(xx), y(yy) { }
	void clear() {
		x.clear();
		y.clear();
	}
	bool operator==(const spoint_t& q) const noexcept {
		return (x ==  q.x && y == q.y);
	}
};


template<const uint N>
struct point_t {
	bignum<N>	x;
	bignum<N>	y;
	bignum<N>	z;
#ifdef	ommit
	u64 *xd() { return const_cast<u64 *>(x.data()); }
	u64 *yd() { return const_cast<u64 *>(y.data()); }
	u64 *zd() { return const_cast<u64 *>(z.data()); }
#endif
	point_t() = default;
	point_t(const point_t &) = default;
	//explicit point_t(const spoint_t<N>& pt) : x(pt.x), y(pt.y), z(1) {}
	explicit point_t(const bignum<N>& xx, const bignum<N>& yy,
		const bignum<N>& zz=bignum<N>(1)) : x(xx), y(yy), z(zz) {}
	explicit point_t(const u64 *xx, const u64 *yy, const u64 *zz) :
		x(xx), y(yy), z(zz) { }
	explicit point_t(const u64 *xx, const u64 *yy) : x(xx), y(yy), z(1) { }
	void clear() { x.clear(); y.clear(); z.clear(); }
	bool operator==(const point_t& q) const noexcept {
		if (x ==  q.x && y == q.y && z ==  q.z) return true;
		return false;
	}
	int is_zero() const noexcept {
		return z.is_zero(); // | y.is_zero();
	}
};


template<const uint N=4> forceinline static
bool operator==(const spoint_t<N>& p, const point_t<N>& q) noexcept {
	if (p.x ==  q.x && p.y == q.y) return true;
	return false;
}

// point add & double template
/* dbl-1998-cmo-2 algorithm
 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html
 * M ... l1
 * S ... l2
 * T ... x3
 */
// x,y,z in montgomery form
template<const bool A_is_n3=true, typename bnT, typename curveT>
forceinline static
void point_double_jacob(const curveT& curve, bnT& x3, bnT& y3, bnT& z3,
		const bnT &x1, const bnT &y1, const bnT &z1) noexcept
{
#ifdef	ommit
	if ( unlikely(z1.is_zero()) ) // | y1.is_zero())
	{
		/* P_y == 0 || P_z == 0 => [1:1:0] */
		x3 = bnT(1);
		y3 = bnT(1);
		z3.clear();
		return;
	}
#endif
	bool	z_is_one = (curve.mont_one() == z1);
	bnT	t1, t2, l1, l2;
	//if ( likely(curve.a_is_pminus3()) )
#if	__cplusplus >= 201703L
	if constexpr(A_is_n3)
#else
	if ( likely(A_is_n3) )
#endif
	{
		/* Use the faster case.  */
		/* L1 = 3(X - Z^2)(X + Z^2) */
		/*						T1: used for Z^2. */
		/*						T2: used for the right term. */
		if ( unlikely(z_is_one) ) {
			// l1 = X - Z^2
			curve.mod_sub(l1, x1, curve.mont_one());
			// 3(X - Z^2) = 2(X - Z^2) + (X - Z^2)
			// t1 = 2 * l1
			curve.mont_mult2(t1, l1);
			// l1 = 2 * l1 + l1 = 3(X - Z^2)
			curve.mod_add_to(l1, t1);
			// t1 = X + Z^2
			curve.mod_add(t1, x1, curve.mont_one());
			// l1 = 3(X - Z^2)(X + Z^2)
			curve.mont_mmult(l1, l1, t1);
		} else {
			// t1 = Z^2
			curve.mont_msqr(t1, z1);
			// l1 = X - Z^2
			curve.mod_sub(l1, x1, t1);
			// 3(X - Z^2) = 2(X - Z^2) + (X - Z^2)
			// t2 = 2 * l1
			curve.mont_mult2(t2, l1);
			// l1 = l1 + 2 * l1 = 3(X - Z^2)
			curve.mod_add_to(l1, t2);
			// t2 = X + Z^2
			curve.mod_add(t2, x1, t1);
			// l1 = 3(X - Z^2)(X + Z^2)
			curve.mont_mmult(l1, l1, t2);
		}
	} else {
		/* Standard case. */
		/* L1 = 3X^2 + aZ^4 */
		/*					T1: used for aZ^4. */
		// l1 = X^2
		curve.mont_msqr(l1, x1);
		curve.mont_mult2(t1, l1);
		// l1 = 3X^2
		curve.mod_add_to(l1, t1);
		if ( likely(curve.a_is_zero()) ) {
			/* Use the faster case.  */
			/* L1 = 3X^2 */
			// do nothing
		} else if ( unlikely(z_is_one) ) {
			// should be mont_paramA
			curve.mod_add_to(l1, curve.montParamA());
		} else {
			// t1 = Z^4
			curve.mont_msqr(t1, z1, 2);
			// t1 = a * Z^4
			curve.mont_mmult(t1, t1, curve.montParamA());
			// l1 = 3 X^2 + a Z^4
			curve.mod_add_to(l1, t1);
		}
	}

	/* Z3 = 2YZ */
	if ( unlikely(z_is_one) ) {
		curve.mont_mult2(z3, y1);
	} else {
		// Z3 = YZ
		curve.mont_mmult(z3, y1, z1);
		// Z3 *= 2
		curve.mont_mult2(z3, z3);
	}

	/* L2 = 4XY^2 */
	/* t2 = Y^2 */
	curve.mont_msqr(t2, y1);
	// t2 = 2 Y^2
	curve.mont_mult2(t2, t2);
	// l2 =  2 XY^2
	curve.mont_mmult(l2, t2, x1);
	// l2 = 4 X Y^2
	curve.mont_mult2(l2, l2);

	/* X3 = L1^2 - 2L2 */
	/*						T1: used for 2L2. */
	curve.mont_msqr(x3, l1);
	curve.mont_mult2(t1, l2);
	curve.mod_sub_from(x3, t1);

	/* L3 = 8Y^4 */
	/*   L3 reuse t2, t2: taken from above. */
	curve.mont_msqr(t2, t2);	// t2 = t2^2 = 4Y^4
	curve.mont_mult2(t2, t2);	// t2 *= 2, t2 = 8Y^4

	/* Y3 = L1(L2 - X3) - L3 */
	curve.mod_sub(y3, l2, x3);
	curve.mont_mmult(y3, l1, y3);
	curve.mod_sub_from(y3, t2);
}

// x,y in montgomery form
template<const bool A_is_n3=true, typename bnT, typename curveT>
forceinline static
void point_doublez_jacob(const curveT& curve, bnT& x3, bnT& y3, bnT& z3,
		const bnT &x1, const bnT &y1) noexcept
{
	bnT	t1, t2, l1, l2;
	/* Use the faster case.  */
	/* L1 = 3(X - Z^2)(X + Z^2) */
#if	__cplusplus >= 201703L
	if constexpr(A_is_n3)
#else
	if ( likely(A_is_n3) )
#endif
	{
		/* Use the faster case.  */
		/* L1 = 3(X - Z^2)(X + Z^2) */
		/*						T1: used for Z^2. */
		/*						T2: used for the right term. */
		// l1 = X - Z^2
		curve.mod_sub(l1, x1, curve.mont_one());
		// 3(X - Z^2) = 2(X - Z^2) + (X - Z^2)
		// t1 = 2 * l1
		curve.mont_mult2(t1, l1);
		// l1 = 2 * l1 + l1 = 3(X - Z^2)
		curve.mod_add_to(l1, t1);
		// t1 = X + Z^2
		curve.mod_add(t1, x1, curve.mont_one());
		// l1 = 3(X - Z^2)(X + Z^2)
		curve.mont_mmult(l1, l1, t1);
	} else {
		/* Standard case. */
		/* L1 = 3X^2 + aZ^4 */
		/*					T1: used for aZ^4. */
		// l1 = X^2
		curve.mont_msqr(l1, x1);
		curve.mont_mult2(t1, l1);
		// l1 = 3X^2
		curve.mod_add_to(l1, t1);
		if ( likely(curve.a_is_zero()) ) {
			/* Use the faster case.  */
			/* L1 = 3X^2 */
			// do nothing
		} else {
			// should be mont_paramA
			curve.mod_add_to(l1, curve.montParamA());
		}
	}

	/* Z3 = 2YZ */
	curve.mont_mult2(z3, y1);

	/* L2 = 4XY^2 */
	/* t2 = Y^2 */
	curve.mont_msqr(t2, y1);
	// t2 = 2 Y^2
	curve.mont_mult2(t2, t2);
	// l2 =  2 XY^2
	curve.mont_mmult(l2, t2, x1);
	// l2 = 4 X Y^2
	curve.mont_mult2(l2, l2);

	/* X3 = L1^2 - 2L2 */
	/*					T1: used for 2L2. */
	curve.mont_msqr(x3, l1);
	curve.mont_mult2(t1, l2);
	curve.mod_sub_from(x3, t1);

	/* L3 = 8Y^4 */
	/*   L3 reuse t2, t2: taken from above. */
	curve.mont_msqr(t2, t2);	// t2 = t2^2 = 4Y^4
	curve.mont_mult2(t2, t2);	// t2 *= 2, t2 = 8Y^4

	/* Y3 = L1(L2 - X3) - L3 */
	curve.mod_sub(y3, l2, x3);
	curve.mont_mmult(y3, l1, y3);
	curve.mod_sub_from(y3, t2);
}


// point add & double template
/* dbl-2004-hmv algorithm
 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html
 */
// x,y,z in montgomery form
template<typename bnT, typename curveT>
forceinline static
void point_double3n_jacob(const curveT& curve, bnT& x3, bnT& y3, bnT& z3,
		const bnT &x1, const bnT &y1, const bnT &z1) noexcept
{
#ifdef	ommit
	if ( unlikely(z1.is_zero()) ) // | y1.is_zero())
	{
		/* P_y == 0 || P_z == 0 => [1:1:0] */
		//x3 = bnT(1);
		x3.clear();
		y3 = bnT(1);
		z3.clear();
		return;
	}
#endif
	bool	z_is_one = (curve.mont_one() == z1);
	bnT	t1, t2, t3;
	// t2 = 3(x1 + z1^2)(x1 - z1^2)
	if ( unlikely(z_is_one) ) {
		// t2 = X - Z^2
		curve.mod_sub(t2, x1, curve.mont_one());
		curve.mod_add(t1, x1, curve.mont_one());
		curve.mont_mmult(t2, t2, t1);
		curve.mont_mult2(t1, t2);
		curve.mod_add_to(t2, t1);
	} else {
		// t1 = Z^2
		curve.mont_msqr(t1, z1);
		// t2 = X - Z^2
		curve.mod_sub(t2, x1, t1);
		// t1 = X + Z^2
		curve.mod_add_to(t1, x1);
		curve.mont_mmult(t2, t2, t1);
		// t2 = (x1 + z1^2) ( x1 - z1^2)
		// t1 = 2 * t2
		curve.mont_mult2(t1, t2);
		// t2 = t2 + 2 * t2 = 3(X - Z^2)
		curve.mod_add_to(t2, t1);
	}

	// y3 = 2 y1
	curve.mont_mult2(y3, y1);
	/* Z3 = 2YZ */
	if ( unlikely(z_is_one) ) {
		z3 = y3;
	} else {
		// Z3 = YZ
		curve.mont_mmult(z3, y3, z1);
	}

	// y3 = y3^2 = 4 y1^2
	curve.mont_msqr(y3, y3);
	/* t3 = 4XY^2 */
	curve.mont_mmult(t3, y3, x1);
	curve.mont_msqr(y3, y3);
	// y3 = y3^2 / 2 = 8 y1^4
	curve.mont_div2(y3);
	/* x3 = t2^2 */
	curve.mont_msqr(x3, t2);
	// t1 = 2 t3 = 8 x y^2
	curve.mont_mult2(t1, t3);
	// x3 = x3 - t1
	curve.mod_sub_from(x3, t1);
	// t1 = t3 - x3
	curve.mod_sub(t1, t3, x3);
	// t1 = t1 * t2
	curve.mont_mmult(t1, t1, t2);
	// y3 = t1 - y3
	curve.mod_sub(y3, t1, y3);
}

/* add-1998-cmo-2 algorithm
 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html
 */
// x,y,z in montgomery field
// point 1, point 2 not zero
template<const bool A_is_n3=true, typename bnT, typename curveT>
forceinline static
void point_add_jacob(const curveT& curve, bnT& x3, bnT& y3, bnT& z3,
			const bnT& x1, const bnT& y1, const bnT& z1, const bnT& x2,
			const bnT& y2, const bnT& z2) noexcept
{
#ifdef	ommit
	if ( unlikely(z1.is_zero()) )
	{
		x3 = x2;
		y3 = y2;
		z3 = z2;
		return;
	}
	else if (z2.is_zero()) {
		x3 = x1;
		y3 = y1;
		z3 = z1;
		return;
	}
#endif
	// add-2007-bl
	/* add-2007-bl algorithm
	 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html
	 */
	// U1 ... l1
	// U2 ... l2
	// S1 ... l4
	// S2 ... l5
	bool	z1_is_one = (z1 == curve.mont_one());
	bool	z2_is_one = (z2 == curve.mont_one());
#ifdef	WITH_ADD_2007bl
	bnT	u1, u2, s1, s2, h, i, j, r, v;
	bnT	t1;
#else
	bnT	u1, u2, s1, s2, h, hh, hhh, r, v;
#define	t1		z1z1
#endif
	bnT	z1z1;
	bnT	z2z2;

	/* u1 = x1 z2^2  */
	/* u2 = x2 z1^2  */
	u1 = x1;
	if ( likely(!z2_is_one) ) {
		// z2z2 = z2^2
		curve.mont_msqr(z2z2, z2);
		// u1 = x1 z2^2
		curve.mont_mmult(u1, u1, z2z2);
	}
	u2 = x2;
	if ( likely(!z1_is_one) ) {
		// z1z1 = z1^2
		curve.mont_msqr(z1z1, z1);
		// u2 = x2 z1^2
		curve.mont_mmult(u2, u2, z1z1);
	}

	/* h = u2 - u1 */
	curve.mod_sub(h, u2, u1);
	/* s1 = y1 z2^3  */
	// s1 = y1
	s1 = y1;
	if ( likely(!z2_is_one) ) {
		// s1 = y1 z2^3
		curve.mont_mmult(s1, s1, z2z2);
		curve.mont_mmult(s1, s1, z2);
	}

	/* s2 = y2 z1^3  */
	// s2 = y2
	s2 = y2;
	if ( likely(!z1_is_one) ) {
		// s2 = y2 z1^3
		curve.mont_mmult(s2, s2, z1z1);
		curve.mont_mmult(s2, s2, z1);
	}
	/* r = s2 - s1  */
	curve.mod_sub(r, s2, s1);

	if ( unlikely(h.is_zero()) ) {
		if (r.is_zero()) {
			/* P1 and P2 are the same - use duplicate function. */
			if ( unlikely(z1_is_one) )
				point_doublez_jacob<A_is_n3>(curve, x3, y3, z3, x1, y1);
			else
				point_double_jacob<A_is_n3>(curve, x3, y3, z3, x1, y1, z1);
			return;
		}
		/* P1 is the inverse of P2.  */
		x3 = bnT(1);
		y3 = bnT(1);
		z3.clear();
		return;
	}
#ifdef	WITH_ADD_2007bl
	// r = 2 * (s2 -s1)
	curve.mont_mult2(r, r);
	// i = (2*h)^2
	curve.mont_mult2(i, h);
	curve.mont_msqr(i, i);

	// j = h * i
	curve.mont_mmult(j, h, i);
	// v = u1 * i
	curve.mont_mmult(v, u1, i);

	// x3 = r^2 - j - 2*v
	curve.mont_msqr(x3, r);
	curve.mod_sub_from(x3, j);
	curve.mod_sub_from(x3, v);
	curve.mod_sub_from(x3, v);

	// y3 = v - x3
	curve.mod_sub(y3, v, x3);
	// y3 = r * (v - x3)
	curve.mont_mmult(y3, r, y3);
	// t1 = 2 * s1 * j
	curve.mont_mmult(t1, s1, j);
	curve.mont_mult2(t1, t1);
	// y3 = r * (v - x3) - 2 * s1 *j
	curve.mod_sub_from(y3, t1);

	if (z2_is_one) {
		// z3 = 2 z1 h
		if (z1_is_one) {
			curve.mont_mult2(z3, h);
		} else {
			curve.mont_mult2(t1, z1);
			curve.mont_mmult(z3, t1, h);
		}
	} else {
		// z3 = ((z1 + z2)^2 -z1z1 - z2z2) * h
		if (z1_is_one) {
			// t1 = 2 * z2
			// z3 = t1 * h = 2 * z2 * h
			curve.mont_mult2(t1, z2);
		} else {
			// t1 = (z1 + z2)^2 -z1z1 - z2z2
			curve.mod_add(t1, z1, z2);
			curve.mont_msqr(t1, t1);
			curve.mod_sub_from(t1, z1z1);
			curve.mod_sub_from(t1, z2z2);
		}
		curve.mont_mmult(z3, t1, h);
	}
#else
	// hh = h^2
	curve.mont_msqr(hh, h);
	// hhh = h * hh
	curve.mont_mmult(hhh, hh, h);
	// v = u1 * hh
	curve.mont_mmult(v, u1, hh);

	// x3 = r^2 - hhh - 2*v
	curve.mont_msqr(x3, r);
	curve.mod_sub_from(x3, hhh);
	curve.mod_sub_from(x3, v);
	curve.mod_sub_from(x3, v);

	// y3 = v - x3
	curve.mod_sub(y3, v, x3);
	// y3 = r * (v - x3)
	curve.mont_mmult(y3, r, y3);
	// t1 = s1 * hhh
	curve.mont_mmult(t1, s1, hhh);
	// y3 = r * (v - x3) - s1 *hhh
	curve.mod_sub_from(y3, t1);

	// z3 = z1 * z2 * h
	if ( unlikely(z2_is_one) ) {
		// z3 = z1 h
		if (z1_is_one) {
			z3 = h;
		} else {
			curve.mont_mmult(z3, z1, h);
		}
	} else {
		// z3 = z1 * z2 * h
		if (z1_is_one) {
			// t1 = z2
			t1 =  z2;
		} else {
			// t1 = z1 * z2
			curve.mont_mmult(t1, z1, z2);
		}
		curve.mont_mmult(z3, t1, h);
	}
#undef	t1
#endif
}

template<const bool A_is_n3=true, typename bnT, typename curveT>
forceinline static
void point_addz_jacob(const curveT& curve, bnT& x3, bnT& y3, bnT& z3,
			const bnT& x1, const bnT& y1, const bnT& z1, const bnT& x2,
			const bnT& y2) noexcept
{
#ifdef	ommit
	if ( unlikely(z1.is_zero()) ) {
		x3 = x2;
		y3 = y2;
		z3 = curve.mont_one();
		return;
	}
#endif
	bool	z1_is_one = (z1 == curve.mont_one());
	bnT	u1, u2, s1, s2, h, hh, hhh, r, v;
#define	t1		z1z1
	bnT	z1z1;

	/* u1 = x1 z2^2  */
	/* u2 = x2 z1^2  */
	u1 = x1;
	u2 = x2;
	if ( likely(!z1_is_one) ) {
		// z1z1 = z1^2
		curve.mont_msqr(z1z1, z1);
		// u2 = x2 z1^2
		curve.mont_mmult(u2, u2, z1z1);
	}

	/* h = u2 - u1 */
	curve.mod_sub(h, u2, u1);
	/* s1 = y1 z2^3  */
	// s1 = y1
	s1 = y1;

	/* s2 = y2 z1^3  */
	// s2 = y2
	s2 = y2;
	if ( likely(!z1_is_one) ) {
		// s2 = y2 z1^3
		curve.mont_mmult(s2, s2, z1z1);
		curve.mont_mmult(s2, s2, z1);
	}
	/* r = s2 - s1  */
	curve.mod_sub(r, s2, s1);

	if ( unlikely(h.is_zero()) ) {
		if (r.is_zero()) {
			/* P1 and P2 are the same - use duplicate function. */
			point_doublez_jacob<A_is_n3>(curve, x3, y3, z3, x2, y2);
			return;
		}
		/* P1 is the inverse of P2.  */
		x3 = bnT(1);
		y3 = bnT(1);
		z3.clear();
		return;
	}
	// hh = h^2
	curve.mont_msqr(hh, h);
	// hhh = h * hh
	curve.mont_mmult(hhh, hh, h);
	// v = u1 * hh
	curve.mont_mmult(v, u1, hh);

	// x3 = r^2 - hhh - 2*v
	curve.mont_msqr(x3, r);
	curve.mod_sub_from(x3, hhh);
	curve.mod_sub_from(x3, v);
	curve.mod_sub_from(x3, v);

	// y3 = v - x3
	curve.mod_sub(y3, v, x3);
	// y3 = r * (v - x3)
	curve.mont_mmult(y3, r, y3);
	// t1 = s1 * hhh
	curve.mont_mmult(t1, s1, hhh);
	// y3 = r * (v - x3) - s1 *hhh
	curve.mod_sub_from(y3, t1);

	// z3 = z1 * z2 * h
	// z3 = z1 h
	if ( unlikely(z1_is_one) ) {
		z3 = h;
	} else {
		curve.mont_mmult(z3, z1, h);
	}
#undef	t1
}


template<typename bnT, typename curveT>
forceinline static bool
pointY_recover(const curveT& curve, bnT& y1, const bnT& x1, const bool bOdd)
	noexcept
{
	bnT		t1, t2;
	curve.to_montgomery(t1, x1);
	// t2 = x^2 + a
	curve.mont_msqr(t2, t1);
	curve.mod_add_to(t2, curve.montParamA());
	// t2 = x^3 + ax
	curve.mont_mmult(t2, t2, t1);
	// t1 = t2 reduction
	curve.from_montgomery(t1, t2);
	// t1 = t1 + b = x^3 + ax +b
	if (t1.add_to(curve.paramB()) || t1.cmp(curve.paramP()) >= 0) {
		t1.sub_from(curve.paramP());
	}
	// need mod_sqrt
	// y^2 = x^3 + ax + b
	auto ret = curve.mod_sqrt(y1, t1);
	if ( likely(ret) ) {
		// odd for negative bignum
		if (bOdd ^ y1.is_odd()) y1.sub(curve.paramP(), y1);
	}
	return ret;
}

}	// namespace ecc


/*
 * Computes result = product % mod
 * for special form moduli: p = 2^k-c, for small c (note the minus sign)
 *
 * References:
 * R. Crandall, C. Pomerance. Prime Numbers: A Computational Perspective.
 * 9 Fast Algorithms for Large-Integer Arithmetic. 9.2.3 Moduli of special form
 * Algorithm 9.2.13 (Fast mod operation for special-form moduli).
 */
//__attribute__((optimize("unroll-loops")))
template<const uint N> forceinline
static void
vli_mmod_special(u64 *result, const u64 *product, const u64 *mod) noexcept
{
	u64 c = -mod[0];
	u64 t[N * 2];
	u64 r[N * 2];

	vli_set<N * 2>(r, product);
	while (!vli_is_zero<N>(r + N)) {
		vli_umult<N>(t, r + N, c);
		vli_clear<N>(r + N);
		vli_add_to<N * 2>(r, t);
	}
	vli_set<N>(t, mod);
	vli_clear<N>(t + N);
	while (vli_cmp<N * 2>(r, t) >= 0)
		vli_sub_from<N * 2>(r, t);
	vli_set<N>(result, r);
}

/*
 * Computes result = product % mod
 * for special form moduli: p = 2^{k-1}+c, for small c (note the plus sign)
 * where k-1 does not fit into qword boundary by -1 bit (such as 255).

 * References (loosely based on):
 * A. Menezes, P. van Oorschot, S. Vanstone. Handbook of Applied Cryptography.
 * 14.3.4 Reduction methods for moduli of special form. Algorithm 14.47.
 * URL: http://cacr.uwaterloo.ca/hac/about/chap14.pdf
 *
 * H. Cohen, G. Frey, R. Avanzi, C. Doche, T. Lange, K. Nguyen, F. Vercauteren.
 * Handbook of Elliptic and Hyperelliptic Curve Cryptography.
 * Algorithm 10.25 Fast reduction for special form moduli
 */
template<const uint N> forceinline
static void
vli_mmod_special2(u64 *result, const u64 *product, const u64 *mod) noexcept
{
	u64 c2 = mod[0] * 2;
	u64 q[N];
	u64 r[N * 2];
	u64 m[N * 2]; /* expanded mod */
	bool carry; /* last bit that doesn't fit into q */
	int i;

	vli_set<N>(m, mod);
	vli_clear<N>(m + N);

	vli_set<N>(r, product);
	/* q and carry are top bits */
	vli_set<N>(q, product + N);
	vli_clear<N>(r + N);
	carry = vli_is_negative<N>(r);
	if (carry)
		r[N - 1] &= (1ull << 63) - 1;
	for (i = 1; carry || !vli_is_zero<N>(q); i++) {
		u64 qc[N * 2];

		vli_umult<N>(qc, q, c2);
		if (carry)
			vli_uadd_to<N*2>(qc, mod[0]);
		vli_set<N>(q, qc + N);
		vli_clear<N>(qc + N);
		carry = vli_is_negative<N>(qc);
		if (carry)
			qc[N - 1] &= (1ull << 63) - 1;
		if (i & 1)
			vli_sub_from<N*2>(r, qc);
		else
			vli_add_to<N*2>(r, qc);
		//	vli_add<N*2>(r, r, qc);
	}
	while (vli_is_negative<N*2>(r))
		vli_add_to<N*2>(r, m);
	while (vli_cmp<N*2>(r, m) >= 0)
		vli_sub_from<N*2>(r, m);

	vli_set<N>(result, r);
}

/* Computes result = product % mod using Barrett's reduction with precomputed
 * value mu appended to the mod after ndigits, mu = (2^{2w} / mod) and have
 * length ndigits + 1, where mu * (2^w - 1) should not overflow ndigits
 * boundary.
 *
 * Reference:
 * R. Brent, P. Zimmermann. Modern Computer Arithmetic. 2010.
 * 2.4.1 Barrett's algorithm. Algorithm 2.5.
 */
template<const uint ndigits> static void forceinline
#ifdef	WITH_C2GO
vli_mmod_barrett(u64 *result, const u64 *product, const u64 *mod, u64 *buff) noexcept
#else
vli_mmod_barrett(u64 *result, const u64 *product, const u64 *mod) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*q = buff;
	u64	*r = buff + ndigits * 2;
#else
	u64 q[ndigits * 2];
	u64 r[ndigits];
#endif
	const u64 *mu = mod + ndigits;

	vli_mult<ndigits>(q, product + ndigits, mu);
	if (mu[ndigits])
		vli_add_to<ndigits>(q + ndigits, product + ndigits);
	// add remain * mod
	vli_set<ndigits>(r, q+ndigits);
	vli_umult<ndigits>(q, mu, product[ndigits-1]);
	if (mu[ndigits])
		vli_uadd_to<ndigits>(q + ndigits, product[ndigits-1]);
	vli_rshift1w<ndigits>(q+ndigits);
	vli_add<ndigits>(result, r, q+ndigits);
	vli_mult<ndigits>(q, mod, result);
	vli_sub<ndigits*2>(q, product, q);
	if (!vli_is_zero<ndigits>(q + ndigits) ||
	       vli_cmp<ndigits>(q, mod) >= 0) {
		vli_sub_from<ndigits>(q, mod);
	}
	vli_set<ndigits>(result, q);
}

template<const uint N> forceinline
static void
vli_div_barrett(u64 *result, const u64 *product, const u64 *mu) noexcept
{
	u64 q[N * 2];
	u64 r[N];

	vli_mult<N>(q, product + N, mu);
	if (mu[N])
		vli_add_to<N>(q + N, product + N);
	vli_set<N>(r, q+N);
	vli_umult<N>(q, mu, product[N-1]);
	if (mu[N])
		vli_uadd_to<N>(q + N, product[N-1]);
	vli_rshift1w<N>(q+N);
	vli_add<N>(result, r, q+N);
}

/* Computes result = (1 / p_input) % mod. All VLIs are the same size.
 * See "From Euclid's GCD to Montgomery Multiplication to the Great Divide"
 * https://labs.oracle.com/techrep/2001/smli_tr-2001-95.pdf
 */
template<const uint N> forceinline
static void
#ifdef	WITH_C2GO
vli_mod_inv(u64 *result, const u64 *input, const u64 *mod, u64 *buff) noexcept
#else
vli_mod_inv(u64 *result, const u64 *input, const u64 *mod) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*a=buff;
	u64	*b=a+N;
	u64	*u=b+N;
	u64	*v=u+N;
#else
	u64 a[N], b[N];
	u64 u[N], v[N];
#endif
	int cmp_result;

	if (vli_is_zero<N>(input)) {
		vli_clear<N>(result);
		return;
	}

	vli_set<N>(a, input);
	vli_set<N>(b, mod);
	vli_clear<N>(u);
	u[0] = 1;
	vli_clear<N>(v);

	while ((cmp_result = vli_cmp<N>(a, b)) != 0) {
		bool carry = false;

		if (vli_is_even(a)) {
			vli_rshift1<N>(a);

			if (!vli_is_even(u))
				carry = vli_add_to<N>(u, mod);

			vli_rshift1<N>(u, carry);
		} else if (vli_is_even(b)) {
			vli_rshift1<N>(b);

			if (!vli_is_even(v))
				carry = vli_add_to<N>(v, mod);

			vli_rshift1<N>(v, carry);
		} else if (cmp_result > 0) {
			vli_sub_from<N>(a, b);
			vli_rshift1<N>(a);

			if (vli_cmp<N>(u, v) < 0)
				vli_add_to<N>(u, mod);

			vli_sub_from<N>(u, v);
			if (!vli_is_even(u))
				carry = vli_add_to<N>(u, mod);

			vli_rshift1<N>(u, carry);
		} else {
			vli_sub_from<N>(b, a);
			vli_rshift1<N>(b);

			if (vli_cmp<N>(v, u) < 0)
				vli_add_to<N>(v, mod);

			vli_sub_from<N>(v, u);
			if (!vli_is_even(v))
				carry = vli_add_to<N>(v, mod);

			vli_rshift1<N>(v, carry);
		}
	}

	vli_set<N>(result, u);
}

/* Computes result = (1 / p_input) % mod. All VLIs are the same size.
 * The binary extended gcd algorithm was first described by Knuth
 */
// x is mod, prime > 2
#ifdef	ommit
using namespace vli;
template<const uint N> forceinline
static void
vli_mod_inv_new(u64 *result, const u64 *y, const u64 *x) noexcept
{
	// mod should be prime, >3 MUST BE odd
	bignum<N>	*res = reinterpret_cast<bignum<N> *>(result);
	if ( vli_is_even(x) ) {
		res->clear();
		return;
	}
	bignumz<N>	b(0l), d(1);
	bignum<N>	u(x), v(y);

	while ( !u.is_zero() ) {
		while (u.is_even()) {
			u.rshift1();
			if (b.is_even()) {
				b.rshift1();
			} else {
				b.sub(b, x);
				b.rshift1();
			}
		}

		while (v.is_even()) {
			v.rshift1();
			if (d.is_even()) {
				d.rshift1();
			} else {
				d.sub(d, x);
				d.rshift1();
			}
		}

		if (u >= v) {
			u.sub_from(v);
			b.sub(b, d);
		} else {
			v.sub_from(u);
			d.sub(d, b);
		}
	}
	if (!v.is_one()) {
		res->clear();
	} else {
		if (d.is_negative()) d.add(d, x);
		*res = d.abs();
	}
}
#else
template<const uint N> forceinline
static void
vli_mod_inv_new(u64 *result, const u64 *y, const u64 *x) noexcept
{
	// mod should be prime, >3 MUST BE odd
	if ( vli_is_even(x) ) {
		vli_clear<N>(result);
		return;
	}
	u64		b[N], d[N], u[N], v[N];
	int		bc = 0, dc = 0;
	vli_set<N>(u, x);
	vli_set<N>(v, y);
	vli_clear<N>(b);
	vli_clear<N>(d);
	d[0] = 1;

	while ( !vli_is_zero<N>(u) ) {
		while (vli_is_even(u)) {
			vli_rshift1<N>(u);
			if (vli_is_even(b)) {
				vli_rshift1<N>(b);
				if (bc & 1) b[N-1] |= (1L << 63);
				bc >>= 1;
			} else {
				if (vli_sub_from<N>(b, x)) bc--;
				vli_rshift1<N>(b);
				if (bc & 1) b[N-1] |= (1L << 63);
				bc >>= 1;
			}
		}

		while (vli_is_even(v)) {
			vli_rshift1<N>(v);
			if (vli_is_even(d)) {
				vli_rshift1<N>(d);
				if (dc & 1) d[N-1] |= (1L << 63);
				dc >>= 1;
			} else {
				if (vli_sub_from<N>(d, x)) dc--;
				vli_rshift1<N>(d);
				if (dc & 1) d[N-1] |= (1L << 63);
				dc >>= 1;
			}
		}

		if (vli_cmp<N>(u,v) >= 0) {
			vli_sub_from<N>(u, v);
			bc -= dc;
			if (vli_sub_from<N>(b, d)) bc--;
		} else {
			vli_sub_from<N>(v, u);
			dc -= bc;
			if (vli_sub_from<N>(d, b)) dc--;
		}
	}
	if (!vli_is_one<N>(v)) {
		vli_clear<N>(result);
	} else if ( dc ) {
		vli_add<N>(result, x, d);
	} else {
		vli_set<N>(result, d);
	}
}
#endif

using namespace vli;
template<const uint N> forceinline
static void
mod_inv(bignum<N>& res, const bignum<N>& x, const bignum<N>& mod) noexcept
{
	u64	*rp = const_cast<u64 *>(res.data());
	vli_mod_inv<N>(rp, x.data(), mod.data());
}

/*-
 * This function looks at w+1 scalar bits (5 current, 1 adjacent less
 * significant bit), and recodes them into a signed digit for use in fast point
 * multiplication: the use of signed rather than unsigned digits means that
 * fewer points need to be precomputed, given that point inversion is easy
 * (a precomputed point dP makes -dP available as well).
 *
 * BACKGROUND:
 *
 * Signed digits for multiplication were introduced by Booth ("A signed binary
 * multiplication technique", Quart. Journ. Mech. and Applied Math., vol. IV,
 * pt. 2 (1951), pp. 236-240), in that case for multiplication of integers.
 * Booth's original encoding did not generally improve the density of nonzero
 * digits over the binary representation, and was merely meant to simplify the
 * handling of signed factors given in two's complement; but it has since been
 * shown to be the basis of various signed-digit representations that do have
 * further advantages, including the wNAF, using the following general approach:
 *
 * (1) Given a binary representation
 *
 *       b_k  ...  b_2  b_1  b_0,
 *
 *     of a nonnegative integer (b_k in {0, 1}), rewrite it in digits 0, 1, -1
 *     by using bit-wise subtraction as follows:
 *
 *        b_k b_(k-1)  ...  b_2  b_1  b_0
 *      -     b_k      ...  b_3  b_2  b_1  b_0
 *       -------------------------------------
 *        s_k b_(k-1)  ...  s_3  s_2  s_1  s_0
 *
 *     A left-shift followed by subtraction of the original value yields a new
 *     representation of the same value, using signed bits s_i = b_(i+1) - b_i.
 *     This representation from Booth's paper has since appeared in the
 *     literature under a variety of different names including "reversed binary
 *     form", "alternating greedy expansion", "mutual opposite form", and
 *     "sign-alternating {+-1}-representation".
 *
 *     An interesting property is that among the nonzero bits, values 1 and -1
 *     strictly alternate.
 *
 * (2) Various window schemes can be applied to the Booth representation of
 *     integers: for example, right-to-left sliding windows yield the wNAF
 *     (a signed-digit encoding independently discovered by various researchers
 *     in the 1990s), and left-to-right sliding windows yield a left-to-right
 *     equivalent of the wNAF (independently discovered by various researchers
 *     around 2004).
 *
 * To prevent leaking information through side channels in point multiplication,
 * we need to recode the given integer into a regular pattern: sliding windows
 * as in wNAFs won't do, we need their fixed-window equivalent -- which is a few
 * decades older: we'll be using the so-called "modified Booth encoding" due to
 * MacSorley ("High-speed arithmetic in binary computers", Proc. IRE, vol. 49
 * (1961), pp. 67-91), in a radix-2^5 setting.  That is, we always combine five
 * signed bits into a signed digit:
 *
 *       s_(4j + 4) s_(4j + 3) s_(4j + 2) s_(4j + 1) s_(4j)
 *
 * The sign-alternating property implies that the resulting digit values are
 * integers from -16 to 16.
 *
 * Of course, we don't actually need to compute the signed digits s_i as an
 * intermediate step (that's just a nice way to see how this scheme relates
 * to the wNAF): a direct computation obtains the recoded digit from the
 * six bits b_(4j + 4) ... b_(4j - 1).
 *
 * This function takes those five bits as an integer (0 .. 63), writing the
 * recoded digit to *sign (0 for positive, 1 for negative) and *digit (absolute
 * value, in the range 0 .. 8).  Note that this integer essentially provides the
 * input bits "shifted to the left" by one position: for example, the input to
 * compute the least significant recoded digit, given that there's no bit b_-1,
 * has to be b_4 b_3 b_2 b_1 b_0 0.
 *
 */
// boothW
template<const uint W> forceinline
static int recode_scalar_bits(uint& digit, const uint in)
{
	static_assert(W > 2 && W <= 8, "w of wNAF between 3 .. 8");
    int s, d;

    s = ~((in >> W) - 1);       /* sets all bits to MSB(in), 'in' seen as
                                 * W+1 bit value */
    d = (1 << (W+1)) - in - 1;
    d = (d & s) | (in & ~s);
    d = (d >> 1) + (d & 1);

    digit = d;
	return (s & 1);
}

#endif	//	__ECC_IMPL_H__
