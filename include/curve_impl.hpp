/*
 * Copyright(c) 2020 Jesse Kuang  <jkuang@21cn.com>
 * All rights reserved.
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
#ifndef __CURVE_IMPL_HPP__
#define __CURVE_IMPL_HPP__

#include <string.h>
#include "vli.hpp"
#include "vli_bn.hpp"
#include "ecc_impl.hpp"
#include "mont.hpp"

namespace ecc {

using namespace vli;
constexpr int W = 5;
constexpr int wSize = 1<<(W-1);
constexpr int wNAFsize = 1<<W;

template<const uint N, typename curveT> forceinline
void pre_compute(const curveT& curve, point_t<N> pre_comp[wSize + 1],
		const point_t<N>& p) noexcept
{
	pre_comp[0].clear();
	//memset(&pre_comp[0], 0, sizeof(pre_comp[0]));
	pre_comp[1] = p;
	for (int i = 2; i <= wSize; ++i) {
		if (i & 1) {
			curve.point_add(pre_comp[i], pre_comp[i-1], p);
		} else {
			curve.point_double(pre_comp[i], pre_comp[i >> 1]);
		}
	}
}

// fast scalar Base mult for 256 Bits ECC Curve
/*-
 * Base point pre computation
 * --------------------------
 *
 * Two different sorts of precomputed tables are used in the following code.
 * Each contain various points on the curve, where each point is three field
 * elements (x, y, z).
 *
 * For the base point table, z is usually 1 (0 for the point at infinity).
 * This table has 2 * 16 elements, starting with the following:
 * index | bits    | point
 * ------+---------+------------------------------
 *     0 | 0 0 0 0 | 0G
 *     1 | 0 0 0 1 | 1G
 *     2 | 0 0 1 0 | 2^64G
 *     3 | 0 0 1 1 | (2^64 + 1)G
 *     4 | 0 1 0 0 | 2^128G
 *     5 | 0 1 0 1 | (2^128 + 1)G
 *     6 | 0 1 1 0 | (2^128 + 2^64)G
 *     7 | 0 1 1 1 | (2^128 + 2^64 + 1)G
 *     8 | 1 0 0 0 | 2^192G
 *     9 | 1 0 0 1 | (2^192 + 1)G
 *    10 | 1 0 1 0 | (2^192 + 2^64)G
 *    11 | 1 0 1 1 | (2^192 + 2^64 + 1)G
 *    12 | 1 1 0 0 | (2^192 + 2^128)G
 *    13 | 1 1 0 1 | (2^192 + 2^128 + 1)G
 *    14 | 1 1 1 0 | (2^192 + 2^128 + 2^64)G
 *    15 | 1 1 1 1 | (2^192 + 2^128 + 2^64 + 1)G
 * followed by a copy of this with each element multiplied by 2^32.
 *
 * The reason for this is so that we can clock bits into four different
 * locations when doing simple scalar multiplies against the base point,
 * and then another four locations using the second 16 elements.
 *
 * Tables for other points have table[i] = iG for i in 0 .. 16. */
template<const uint N, typename curveT> forceinline void
pre_compute_base(const curveT& curve, point_t<N>  pre_comps[2][16]) noexcept
{
	int i, j;
	point_t<N>	G;
	curve.to_montgomery(G.x, curve.getGx());
	curve.to_montgomery(G.y, curve.getGy());
	G.z = curve.mont_one();
	/*
	 * compute 2^64*G, 2^128*G, 2^192*G for the first table, 2^32*G, 2^96*G,
	 * 2^160*G, 2^224*G for the second one
	 */
	pre_comps[0][1] = G;
	for (i = 1; i <= 8; i <<= 1) {
		curve.point_double(pre_comps[1][i], pre_comps[0][i]);
		for (j = 0; j < 31; ++j) {
			curve.point_double(pre_comps[1][i], pre_comps[1][i]);
		}
		if (i == 8) break;
		curve.point_double(pre_comps[0][2 * i], pre_comps[1][i]);
		for (j = 0; j < 31; ++j) {
			curve.point_double(pre_comps[0][2 * i], pre_comps[0][2 * i]);
		}
	}
	for (i = 0; i < 2; i++) {
		/* pre_comps[i][0] is the point at infinity */
		//memset(&pre_comps[i][0], 0, sizeof(pre_comps[i][0]));
		pre_comps[i][0].clear();
		/* the remaining multiples */
		/* 2^64*G + 2^128*G resp. 2^96*G + 2^160*G */
		curve.point_add(pre_comps[i][6], pre_comps[i][4], pre_comps[i][2]);
		/* 2^64*G + 2^192*G resp. 2^96*G + 2^224*G */
		curve.point_add(pre_comps[i][10], pre_comps[i][8], pre_comps[i][2]);
		/* 2^128*G + 2^192*G resp. 2^160*G + 2^224*G */
		curve.point_add(pre_comps[i][12], pre_comps[i][8], pre_comps[i][4]);
		/*
		 * 2^64*G + 2^128*G + 2^192*G resp. 2^96*G + 2^160*G + 2^224*G
		 */
		curve.point_add(pre_comps[i][14], pre_comps[i][12], pre_comps[i][2]);
		for (j = 1; j < 8; ++j) {
			/* odd multiples: add G resp. 2^32*G */
			curve.point_add(pre_comps[i][2 * j + 1], pre_comps[i][2 * j],
							pre_comps[i][1]);
		}
	}
	// make_points_affine(31, &(pre->pre_comps[0][1]), tmp_smallfelems);
	for(i=0;i < 2; i++) {
		for (j = 1; j < 16; j++) {
			curve.from_montgomery(pre_comps[i][j].x, pre_comps[i][j].x);
			curve.from_montgomery(pre_comps[i][j].y, pre_comps[i][j].y);
			curve.from_montgomery(pre_comps[i][j].z, pre_comps[i][j].z);
			curve.apply_z_mont(pre_comps[i][j]);
		}
	}
	// convert to affine jacobian
}

/**
 * struct ecc_curve - definition of elliptic curve
 *
 * @name:	Short name of the curve.
 * @g:		Generator point of the curve.
 * @p:		Prime number, if Barrett's reduction is used for this curve
 * @n:		Order of the curve group.
 * @a:		Curve parameter a.
 * @b:		Curve parameter b.
 */
template<const uint N=4>
class ecc_curve {
public:
	using felem_t = bignum<N>;
	ecc_curve(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
			const u64 *_n, const u64 *_a, const u64 *_b, const bool a_n3 = true,
			const bool bBat=false) :
		name(_name), gx(_gx), gy(_gy),
		p(_p), n(_n), a(_a), b(_b), _a_is_neg3(a_n3),
		use_barrett(bBat)
	{
		static_assert(N > 3, "curve only support 256Bits or more");
	}
	ecc_curve(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
			const u64 *_n, const u64 *_a, const u64 *_b, const u64* rrP,
			const u64* rrN, const u64 k0P, const u64 k0N,
			const bool a_n3 = true):
		name(_name), gx(_gx), gy(_gy),
		p(_p), n(_n), a(_a), b(_b), rr_p(rrP), rr_n(rrN), k0_p(k0P), k0_n(k0N),
		_a_is_neg3(a_n3), use_barrett(false)
	{
		static_assert(N > 3, "curve only support 256Bits or more");
	}
	ecc_curve(ecc_curve &&) = default;
	const uint ndigits() const { return N; }
	explicit operator bool() const noexcept {
		//return name != nullptr && name[0] != 0;
		return name != "";
	}
	void getP(u64 *v) const noexcept { p.set(v); }
	void getN(u64 *v) const noexcept { n.set(v); }
	void getA(u64 *v) const noexcept { a.set(v); }
	void getB(u64 *v) const noexcept { b.set(v); }
	void getGx(u64 *v) const noexcept { gx.set(v); }
	void getGy(u64 *v) const noexcept { gy.set(v); }
	const felem_t& getGx() const noexcept { return gx; }
	const felem_t& getGy() const noexcept { return gy; }
	const felem_t& montParamA() const noexcept { return _mont_a; }
	const felem_t& paramP() const noexcept { return p; }
	bool init(const u64 *muP=nullptr, const u64 *muN=nullptr) noexcept {
		if (_inited) return _inited;
#ifdef	WITH_BARRETT
		if (muP == nullptr || muN == nullptr) {
			use_barrett = false;
			return false;
		}
		mu_p = bignum<N+1>(muP);
		mu_n = bignum<N+1>(muN);
		use_barrett = true;
		_inited = true;
		return _inited;
#else
		if (unlikely(k0_p == 0))
		{
			// no calc k0 and rr after instantiation
			return false;
		}
		felem_t	t1;
		t1.clear();
		_mont_one.sub(t1, p);
		if (unlikely( !_a_is_neg3 )) {
			if (likely(a.is_zero())) _a_is_zero = true;
		}
		// should verify calc K0 and RR
		pre_compute_base<N>(*this, g_pre_comp);	// precompute g_pre_comp
		_inited = true;
		return true;
#endif
	}
	const bool a_is_pminus3() const noexcept { return _a_is_neg3; }
	const bool a_is_zero() const noexcept { return _a_is_zero; }
#ifdef	WITH_BARRETT
	void mmod_barrett(felem_t& res, const bn_prod<N>& prod) const noexcept
	{
		if (use_barrett) {
			// res = barrettmod
			prod.mmod_barrett(res, p, mu_p);
		}
		return;
	}
#endif
	const felem_t& mont_one() const noexcept { return _mont_one; }
	void to_montgomery(felem_t& res, const u64 *x) const noexcept
	{
		felem_t   *xx = reinterpret_cast<felem_t *>(const_cast<u64 *>(x));
		res.mont_mult(*xx, rr_p, p, k0_p);
	}
	void to_montgomery(felem_t& res, const felem_t& x) const noexcept
	{
		res.mont_mult(x, rr_p, p, k0_p);
	}
	void from_montgomery(felem_t& res, const felem_t& y) const noexcept
	{
		res.mont_reduction(y, p, k0_p);
	}
	void from_montgomery(u64* result, const felem_t& y) const noexcept
	{
		felem_t   *res = reinterpret_cast<felem_t *>(result);
		res->mont_reduction(y, p, k0_p);
	}
	void
	mont_mmult(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept {
		res.mont_mult(left, right, p, k0_p);
	}
	void mont_msqr(felem_t& res, const felem_t left, const uint nTimes=1)
	const noexcept {
		res.mont_sqr(left, p, k0_p);
		for (uint i=1; i < nTimes; i++) res.mont_sqr(res, p, k0_p);
	}
	// left,right less than p
	void mod_add(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept {
		if (res.add(left, right)) {
			res.sub_from(p);
		}
	}
	// left,right less than p
	void mod_add_to(felem_t& res, const felem_t& right) const noexcept
	{
		if (res.add_to(right)) {
			res.sub_from(p);
		}
	}
	// left,right less than p
	void mod_sub(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept {
		if (res.sub(left, right)) {
			res.add_to(p);
		}
	}
	void mod_sub_from(felem_t& res, const felem_t& right) const noexcept
	{
		if (res.sub_from(right)) {
			res.add_to(p);
		}
	}
	void mont_mult2(felem_t& res, const felem_t& left) const noexcept
	{
#ifdef	ommit
		this->mod_add(res, left, left);
#else
		if (res.lshift1(left) != 0) {
			res.sub_from(p);
		}
#endif
	}
	void mont_mult4(felem_t& res) const noexcept
	{
		this->mont_mult2(res, res);
		this->mont_mult2(res, res);
	}
	void mont_mult8(felem_t& res) const noexcept
	{
		this->mont_mult2(res, res);
		this->mont_mult2(res, res);
		this->mont_mult2(res, res);
	}
/* dbl-1998-cmo-2 algorithm
 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html
 * M ... l1
 * S ... l2
 * T ... x3
 */
	void point_double_jacobian(u64 *x3, u64 *y3, u64 *z3, const u64 *x1,
					const u64 *y1, const u64 *z1 = nullptr) const noexcept
	{
		if ((z1 != nullptr && vli_is_zero<N>(z1)) || vli_is_zero<N>(y1)) {
			/* P_y == 0 || P_z == 0 => [1:1:0] */
			vli_clear<N>(x3);
			vli_clear<N>(y3);
			vli_clear<N>(z3);
			x3[0] = 1;
			y3[0] = 1;
			return;
		}
		bool	z_is_one = (z1 == nullptr || vli_is_one<N>(z1));
		felem_t	xp, yp, zp;
		felem_t	*x3p = reinterpret_cast<felem_t *>(x3);
		felem_t	*y3p = reinterpret_cast<felem_t *>(y3);
		felem_t	*z3p = reinterpret_cast<felem_t *>(z3);
		to_montgomery(xp, x1);
		to_montgomery(yp, y1);
		if (z_is_one) {
			zp = mont_one();
		} else {
			to_montgomery(zp, z1);
		}
#ifdef	ommit
		felem_t	t1, t2, l1, l2;
		if (a_is_pminus3()) {
			/* Use the faster case.  */
			/* L1 = 3(X - Z^2)(X + Z^2) */
			/*						T1: used for Z^2. */
			/*						T2: used for the right term. */
			if (z_is_one) {
				// l1 = X - Z^2
				mod_sub(l1, xp, mont_one());
				// 3(X - Z^2) = 2(X - Z^2) + (X - Z^2)
				// t1 = 2 * l1
				mont_mult2(t1, l1);
				// l1 = 2 * l1 + l1 = 3(X - Z^2)
				mod_add_to(l1, t1);
				// t1 = X + Z^2
				mod_add(t1, xp, mont_one());
				// l1 = 3(X - Z^2)(X + Z^2)
				this->mont_mmult(l1, l1, t1);
			} else {
				// t1 = Z^2
				this->mont_msqr(t1, zp);
				// l1 = X - Z^2
				mod_sub(l1, xp, t1);
				// 3(X - Z^2) = 2(X - Z^2) + (X - Z^2)
				// t2 = 2 * l1
				mont_mult2(t2, l1);
				// l1 = l1 + 2 * l1 = 3(X - Z^2)
				mod_add_to(l1, t2);
				// t2 = X + Z^2
				mod_add(t2, xp, t1);
				// l1 = 3(X - Z^2)(X + Z^2)
				this->mont_mmult(l1, l1, t2);
			}
		} else {
			/* Standard case. */
			/* L1 = 3X^2 + aZ^4 */
			/*					T1: used for aZ^4. */
			// l1 = X^2
			this->mont_msqr(l1, xp);
			mont_mult2(t1, l1);
			// l1 = 3X^2
			mod_add_to(l1, t1);
			if (a_is_zero()) {
				/* Use the faster case.  */
				/* L1 = 3X^2 */
				// do nothing
			} else if (z_is_one) {
				// should be mont_paramA
				mod_add_to(l1, montParamA());
			} else {
				// t1 = Z^4
				this->mont_msqr(t1, zp, 2);
				// t1 = a * Z^4
				this->mont_mmult(t1, t1, montParamA());
				// l1 = 3 X^2 + a Z^4
				mod_add_to(l1, t1);
			}
		}

		/* Z3 = 2YZ */
		if (z_is_one) {
			mont_mult2(*z3p, yp);
		} else {
			// Z3 = YZ
			this->mont_mmult(*z3p, yp, zp);
			// Z3 *= 2
			mont_mult2(*z3p, *z3p);
		}

		/* L2 = 4XY^2 */
		/* t2 = Y^2 */
		this->mont_msqr(t2, yp);
		// t2 = 2 Y^2
		mont_mult2(t2, t2);
		// l2 =  2 XY^2
		this->mont_mmult(l2, t2, xp);
		// l2 = 4 X Y^2
		mont_mult2(l2, l2);

		/* X3 = L1^2 - 2L2 */
		/*						T1: used for 2L2. */
		this->mont_msqr(*x3p, l1);
		mont_mult2(t1, l2);
		mod_sub_from(*x3p, t1);

		/* L3 = 8Y^4 */
		/*   L3 reuse t2, t2: taken from above. */
		this->mont_msqr(t2, t2);		// t2 = t2^2 = 4Y^4
		mont_mult2(t2, t2);	// t2 *= 2, t2 = 8Y^4

		/* Y3 = L1(L2 - X3) - L3 */
		mod_sub(*y3p, l2, *x3p);
		this->mont_mmult(*y3p, l1, *y3p);
		mod_sub_from(*y3p, t2);
#else
		if (z_is_one) point_doublez_jacob(*this, *x3p, *y3p, *z3p, xp, yp); else
			point_double_jacob(*this, *x3p, *y3p, *z3p, xp, yp, zp);
#endif
		// montgomery reduction
		from_montgomery(x3, *x3p);
		from_montgomery(y3, *y3p);
		from_montgomery(z3, *z3p);
	}
	void apply_z(u64 *x1, u64 *y1, u64 *z) const noexcept
	{
		felem_t	t1;

		if (vli_is_one<N>(z)) return;
		if (vli_is_zero<N>(z)) {
			vli_clear<N>(x1);
			vli_clear<N>(y1);
			return;
		}
		vli_mod_inv<N>(z, z, this->p.data());
		felem_t	*xp = reinterpret_cast<felem_t *>(x1);
		felem_t	*yp = reinterpret_cast<felem_t *>(y1);
		felem_t	*zp = reinterpret_cast<felem_t *>(z);
		to_montgomery(*zp, z);
		this->mont_msqr(t1, *zp);	// t1 = z^2
		to_montgomery(*xp, x1);
		to_montgomery(*yp, y1);
		this->mont_mmult(*xp, *xp, t1);	// x1 * z^2
		this->mont_mmult(t1, t1, *zp);	// t1 = z^3
		this->mont_mmult(*yp, *yp, t1);	// y1 * z^3

		// montgomery reduction
		from_montgomery(x1, *xp);
		from_montgomery(y1, *yp);
		// set z to 1
		vli_clear<N>(z);
		z[0] = 1;
	}
	void apply_z(u64* x, u64 *y, const point_t<N>& p) const noexcept
	{
		p.x.set(x);
		p.y.set(y);
		u64	z[N];
		p.z.set(z);
		apply_z(x, y, z);
	}
	void apply_z(point_t<N>& p) const noexcept
	{
		apply_z(p.xd(), p.yd(), p.zd());
	}
	void apply_z_mont(point_t<N>& p) const noexcept
	{
		apply_z(p.xd(), p.yd(), p.zd());
	}
	bool point_eq(const point_t<N>& p, const point_t<N>& q) const noexcept
	{
		if (p == q) return true;
		u64		x1[N], x2[N], y1[N], y2[N];
		this->apply_z(x1, y1, p);
		this->apply_z(x2, y2, q);
		return vli_cmp<N>(x1, x2) == 0 && vli_cmp<N>(y1, y2) == 0;
	}
	void point_neg(point_t<N>& q, const point_t<N>p) const noexcept
	{
		q.x = p.x;
		q.y.sub(this->p, p.y);
		q.z = p.z;
	}
	/* add-1998-cmo-2 algorithm
	 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html
	 */
	void point_add_jacobian(u64 *x3, u64 *y3, u64 *z3, const u64 *x1,
			const u64 *y1, const u64 *z1, const u64 *x2,
			const u64 *y2, const u64 *z2 = nullptr) const noexcept
	{
		felem_t	*x3p = reinterpret_cast<felem_t *>(x3);
		felem_t	*y3p = reinterpret_cast<felem_t *>(y3);
		felem_t	*z3p = reinterpret_cast<felem_t *>(z3);
		if (z1 != nullptr && vli_is_zero<N>(z1)) {
			vli_set<N>(x3, x2);
			vli_set<N>(y3, y2);
			vli_set<N>(z3, z2);
			return;
		} else if (z2 != nullptr && vli_is_zero<N>(z2)) {
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
		bool	z1_is_one = (z1 == nullptr || vli_is_one<N>(z1));
		bool	z2_is_one = (z2 == nullptr || vli_is_one<N>(z2));
#ifdef	WITH_ADD_2007bl
		felem_t	u1, u2, s1, s2, h, i, j, r, v;
		felem_t	t1;
#else
		felem_t	u1, u2, s1, s2, h, hh, hhh, r, v;
#define	t1		z1z1
#endif
		felem_t	z1z1, z1p;
		felem_t	z2z2, z2p;

		/* u1 = x1 z2^2  */
		/* u2 = x2 z1^2  */
		if (z2_is_one) {
			to_montgomery(u1, x1);
		} else {
			to_montgomery(z2p, z2);
			// z2z2 = z2^2
			this->mont_msqr(z2z2, z2p);
			to_montgomery(u1, x1);
			// u1 = x1 z2^2
			this->mont_mmult(u1, u1, z2z2);
		}
		if (z1_is_one) {
			to_montgomery(u2, x2);
		} else {
			to_montgomery(z1p, z1);
			// z1z1 = z1^2
			this->mont_msqr(z1z1, z1p);
			to_montgomery(u2, x2);
			// u2 = x2 z1^2
			this->mont_mmult(u2, u2, z1z1);
		}

		/* h = u2 - u1 */
		mod_sub(h, u2, u1);
		/* s1 = y1 z2^3  */
		// s1 = y1
		to_montgomery(s1, y1);
		if ( ! z2_is_one ) {
			// s1 = y1 z2^3
			this->mont_mmult(s1, s1, z2z2);
			this->mont_mmult(s1, s1, z2p);
		}

		/* s2 = y2 z1^3  */
		// s2 = y2
		to_montgomery(s2, y2);
		if ( !z1_is_one ) {
			// s2 = y2 z1^3
			this->mont_mmult(s2, s2, z1z1);
			this->mont_mmult(s2, s2, z1p);
		}
		/* r = s2 - s1  */
		mod_sub(r, s2, s1);

		if (h.is_zero()) {
			if (r.is_zero()) {
				/* P1 and P2 are the same - use duplicate function. */
				point_double_jacobian(x3, y3, z3, x1, y1, z1);
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
#ifdef	WITH_ADD_2007bl
		// r = 2 * (s2 -s1)
		mont_mult2(r, r);
		// i = (2*h)^2
		mont_mult2(i, h);
		this->mont_msqr(i, i);

		// j = h * i
		this->mont_mmult(j, h, i);
		// v = u1 * i
		this->mont_mmult(v, u1, i);

		// x3 = r^2 - j - 2*v
		this->mont_msqr(*x3p, r);
		mod_sub_from(*x3p, j);
		mod_sub_from(*x3p, v);
		mod_sub_from(*x3p, v);

		// y3 = v - x3
		mod_sub(*y3p, v, *x3p);
		// y3 = r * (v - x3)
		this->mont_mmult(*y3p, r, *y3p);
		// t1 = 2 * s1 * j
		this->mont_mmult(t1, s1, j);
		mont_mult2(t1, t1);
		// y3 = r * (v - x3) - 2 * s1 *j
		mod_sub_from(*y3p, t1);

		if (z2_is_one) {
			// z3 = 2 z1 h
			if (z1_is_one) {
				mont_mult2(*z3p, h);
			} else {
				mont_mult2(t1, z1p);
				this->mont_mmult(*z3p, t1, h);
			}
		} else {
			// z3 = ((z1 + z2)^2 -z1z1 - z2z2) * h
			if (z1_is_one) {
				// t1 = 2 * z2
				// z3 = t1 * h = 2 * z2 * h
				mont_mult2(t1, z2p);
			} else {
				// t1 = (z1 + z2)^2 -z1z1 - z2z2
				mod_add(t1, z1p, z2p);
				this->mont_msqr(t1, t1);
				mod_sub_from(t1, z1z1);
				mod_sub_from(t1, z2z2);
			}
			this->mont_mmult(*z3p, t1, h);
		}
#else
		// hh = h^2
		mont_msqr(hh, h);
		// hhh = h * hh
		mont_mmult(hhh, hh, h);
		// v = u1 * hh
		mont_mmult(v, u1, hh);

		// x3 = r^2 - hhh - 2*v
		this->mont_msqr(*x3p, r);
		mod_sub_from(*x3p, hhh);
		mod_sub_from(*x3p, v);
		mod_sub_from(*x3p, v);

		// y3 = v - x3
		mod_sub(*y3p, v, *x3p);
		// y3 = r * (v - x3)
		this->mont_mmult(*y3p, r, *y3p);
		// t1 = s1 * hhh
		this->mont_mmult(t1, s1, hhh);
		// y3 = r * (v - x3) - s1 *hhh
		mod_sub_from(*y3p, t1);

		// z3 = z1 * z2 * h
		if (z2_is_one) {
			// z3 = z1 h
			if (z1_is_one) {
				*z3p = h;
			} else {
				this->mont_mmult(*z3p, z1p, h);
			}
		} else {
			// z3 = z1 * z2 * h
			if (z1_is_one) {
				// t1 = z2
				t1 =  z2p;
			} else {
				// t1 = z1 * z2
				mont_mmult(t1, z1p, z2p);
			}
			this->mont_mmult(*z3p, t1, h);
		}
#undef	t1
#endif
		// montgomery reduction
		from_montgomery(x3, *x3p);
		from_montgomery(y3, *y3p);
		from_montgomery(z3, *z3p);
	}
	void scalar_mult(point_t<N>& q, const point_t<N>& p,
					const felem_t& scalar) const noexcept
	{
		point_t<N>	tmp;
		point_t<N>	pres[wSize+1];
		uint	nbits = scalar.num_bits();
		q.clear();
		if ( unlikely(nbits == 0) ) return;
		to_montgomery(tmp.x, p.x);
		to_montgomery(tmp.y, p.y);
		to_montgomery(tmp.z, p.z);
		pre_compute<N>(*this, pres, tmp);
		bool	skip = true;
		if (nbits % W != 0) --nbits;
		for (int i = nbits; i >= 0; --i) {
			if (!skip) point_double(q, q);
			if (i % W != 0) continue;
			uint	bits;
#ifdef	ommit
			bits = vli_get_bits<N, W>(scalar.data(), i);
			bool	sign = (bits & (1<<(W-1)));
			int digit = (sign)?(bits - wNAFsize):bits;
			if (i > 0 && scalar.get_bit(i-1)) ++digit;
			if (digit == 0) continue;
			if ( sign ) point_neg(tmp, pres[-digit]); else tmp = pres[digit];
#else
			uint	digit;
			bits = vli_get_bits<N, W+1>(scalar.data(), i-1);
			auto sign = recode_scalar_bits<W>(digit, bits);
			if (digit == 0) continue;
#ifdef	WITH_CONDITIONAL_COPY
			tmp = pres[digit];
			felem_t	ny;
			ny.sub(this->p, tmp.y);
			ny.copy_conditional(tmp.y, sign-1);
			tmp.y = ny;
#else
			if ( sign ) point_neg(tmp, pres[digit]); else tmp = pres[digit];
#endif
#endif
			if (!skip) point_add(q, q, tmp); else {
				q = tmp;
				skip =false;
			}
		}
		// montgomery reduction
		from_montgomery(q.x, q.x);
		from_montgomery(q.y, q.y);
		from_montgomery(q.z, q.z);
		this->apply_z(q);
	}
	void scalar_multNAF2(point_t<N>& q, const point_t<N>& p,
				const felem_t& scalar) const noexcept
	{
		point_t<N>	tmp;
		int	nbits = scalar.num_bits();
		q.clear();
		if (nbits == 0) return;
		point_t<N>	pres[3];
		to_montgomery(tmp.x, p.x);
		to_montgomery(tmp.y, p.y);
		to_montgomery(tmp.z, p.z);
		pres[0].clear();
		pres[1] = tmp;
		point_double(pres[2], tmp);
		bool	skip = true;
		if (nbits & 1) --nbits;
		for (int i = nbits; i >= 0; --i) {
			if (!skip) point_double(q, q);
			if (i & 1) continue;
			uint	bits;
			bits = vli_get_bits<N, 2>(scalar.data(), i);
			bool	sign = (bits & 2);
			int digit = (sign)?(bits - 4):bits;
			if (i > 0 && scalar.get_bit(i-1)) ++digit;
			if (digit == 0) continue;
			if (sign) point_neg(tmp, pres[-digit]); else tmp = pres[digit];
			if (likely(!skip)) point_add(q, q, tmp); else {
				q = tmp;
				skip = false;
			}
		}
		// montgomery reduction
		from_montgomery(q.x, q.x);
		from_montgomery(q.y, q.y);
		from_montgomery(q.z, q.z);
		this->apply_z(q);
	}
	bool select_base_point(point_t<N>& pt, const uint idx) const noexcept {
		if (idx >= 2 * 16) return false;
		pt = g_pre_comp[idx>>4][idx&0xf];
		return true;
	}
	void point_double(point_t<N>& q, const point_t<N>& p) const noexcept {
		point_double_jacob(*this, q.x, q.y, q.z, p.x, p.y, p.z);
	}
	void point_add(point_t<N>& q, const point_t<N>& p1, const point_t<N>& p2)
	const noexcept {
		point_add_jacob(*this, q.x, q.y, q.z, p1.x, p1.y, p1.z,
						p2.x, p2.y, p2.z);
	}
protected:
	const std::string name;
	const felem_t gx;
	const felem_t gy;
	const felem_t p;
	const felem_t n;
	const felem_t a;
	const felem_t b;
	const felem_t rr_p = {};
	const felem_t rr_n = {};
#ifdef	WITH_BARRETT
	bignum<N+1> mu_p;
	bignum<N+1> mu_n;
#endif
	felem_t _mont_one;
	felem_t _mont_a;
	const u64	k0_p = 0;
	const u64	k0_n = 0;
	const uint _ndigits = N;
	point_t<N>	g_pre_comp[2][16];
	const bool _a_is_neg3 = false;
	bool _a_is_zero = false;
	const bool use_barrett = false;
	bool _inited = false;
};

template <const u64 Pk0=1>
class curve256 : public ecc_curve<4> {
public:
	using felem_t = bignum<4>;
	curve256(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
			const u64 *_n, const u64 *_a, const u64 *_b) :
		ecc_curve<4>(_name, _gx, _gy, _p, _n, _a, _b, sm2_p_rr, sm2_n_rr,
					sm2_p_k0, sm2_n_k0)
	{ }
	void to_montgomery(felem_t& res, const u64 *x) const noexcept
	{
		felem_t   *xx = reinterpret_cast<felem_t *>(const_cast<u64 *>(x));
		mont_mult<Pk0>(res, *xx, this->rr_p, this->p);
	}
	void to_montgomery(felem_t& res, const felem_t& x) const noexcept
	{
		mont_mult<Pk0>(res, x, this->rr_p, this->p);
	}
	void from_montgomery(felem_t& res, const felem_t& y) const noexcept
	{
		mont_reduction<Pk0>(res, y, this->p);
	}
	void from_montgomery(u64* result, const felem_t& y) const noexcept
	{
		felem_t   *res = reinterpret_cast<felem_t *>(result);
		mont_reduction<Pk0>(*res, y, this->p);
	}
	void
	mont_mmult(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept {
		mont_mult<Pk0>(res, left, right, this->p);
	}
	void mont_msqr(felem_t& res, const felem_t left, const uint nTimes=1)
	const noexcept {
		mont_sqr<Pk0>(res, left, this->p);
		for (uint i=1; i < nTimes; i++) mont_sqr<Pk0>(res, res, this->p);
	}
	void mont_mult4(felem_t& res) const noexcept
	{
		static_assert(Pk0 == 1, "MUST be sm2");
		u64	carry;
		if ((carry=res.lshift(2)) != 0) {
			this->carry_reduce(res, carry);
		}
	}
	void mont_mult8(felem_t& res) const noexcept
	{
		static_assert(Pk0 == 1, "MUST be sm2");
		u64	carry;
		if ((carry=res.lshift(3)) != 0) {
			this->carry_reduce(res, carry);
		}
	}
	void scalar_mult_base(point_t<4>& q, const felem_t& scalar) noexcept
	{
		bool	skip = true;
		point_t<4>	p;
		q.clear();
		for (int i = 31; i >= 0; --i) {
			/* double */
			if (!skip) point_double(q, q);

			/* add multiples of the generator */
			u64	bits = scalar.get_bit(i + 224) << 3;
			bits |= scalar.get_bit(i + 160) << 2;
			bits |= scalar.get_bit(i + 96) << 1;
			bits |= scalar.get_bit(i + 32);
			if (bits != 0) {
				select_base_point(p, bits | 16);
				if (!skip) point_add(q, q, p.x, p.y); else {
					q = p;
					skip = false;
				}
			}

			/* second, look at the current position */
			bits = scalar.get_bit(i + 192) << 3;
			bits |= scalar.get_bit(i + 128) << 2;
			bits |= scalar.get_bit(i + 64) << 1;
			bits |= scalar.get_bit(i);
			if (bits == 0) continue;
			select_base_point(p, bits);
			if (!skip) point_add(q, q, p.x, p.y); else {
				q = p;
				skip = false;
			}
		}
		this->apply_z(q);
	}
	void point_double_jacobian(u64 *x3, u64 *y3, u64 *z3, const u64 *x1,
					const u64 *y1, const u64 *z1 = nullptr) const noexcept
	{
		if ((z1 != nullptr && vli_is_zero<4>(z1)) || vli_is_zero<4>(y1)) {
			/* P_y == 0 || P_z == 0 => [1:1:0] */
			vli_clear<4>(x3);
			vli_clear<4>(y3);
			vli_clear<4>(z3);
			x3[0] = 1;
			y3[0] = 1;
			return;
		}
		bool	z_is_one = (z1 == nullptr || vli_is_one<4>(z1));
		felem_t	xp, yp, zp;
		felem_t	*x3p = reinterpret_cast<felem_t *>(x3);
		felem_t	*y3p = reinterpret_cast<felem_t *>(y3);
		felem_t	*z3p = reinterpret_cast<felem_t *>(z3);
		to_montgomery(xp, x1);
		to_montgomery(yp, y1);
		if (!z_is_one) {
			to_montgomery(zp, z1);
		}
		if (z_is_one) point_doublez_jacob(*this, *x3p, *y3p, *z3p, xp, yp); else
			point_double_jacob(*this, *x3p, *y3p, *z3p, xp, yp, zp);
		// montgomery reduction
		from_montgomery(x3, *x3p);
		from_montgomery(y3, *y3p);
		from_montgomery(z3, *z3p);
	}
	void point_add_jacobian(u64 *x3, u64 *y3, u64 *z3, const u64 *x1,
			const u64 *y1, const u64 *z1, const u64 *x2,
			const u64 *y2, const u64 *z2 = nullptr) const noexcept
	{
		felem_t	*x3p = reinterpret_cast<felem_t *>(x3);
		felem_t	*y3p = reinterpret_cast<felem_t *>(y3);
		felem_t	*z3p = reinterpret_cast<felem_t *>(z3);
		if (z1 != nullptr && vli_is_zero<4>(z1)) {
			vli_set<4>(x3, x2);
			vli_set<4>(y3, y2);
			vli_set<4>(z3, z2);
			return;
		} else if (z2 != nullptr && vli_is_zero<4>(z2)) {
			vli_set<4>(x3, x1);
			vli_set<4>(y3, y1);
			vli_set<4>(z3, z1);
			return;
		}
		bool	z1_is_one = (z1 == nullptr || vli_is_one<4>(z1));
		bool	z2_is_one = (z2 == nullptr || vli_is_one<4>(z2));
		felem_t	x1p, y1p, z1p;
		felem_t	x2p, y2p, z2p;
		to_montgomery(x1p, x1);
		to_montgomery(y1p, y1);
		if (z1_is_one) {
			z1p = mont_one();
		} else {
			to_montgomery(z1p, z1);
		}
		to_montgomery(x2p, x2);
		to_montgomery(y2p, y2);
		if (z2_is_one) {
			z2p = mont_one();
		} else {
			to_montgomery(z2p, z2);
		}
		if (z2_is_one) {
			point_addz_jacob(*this, *x3p, *y3p, *z3p, x1p, y1p, z1p, x2p, y2p);
		} else {
			point_add_jacob(*this, *x3p, *y3p, *z3p, x1p, y1p, z1p,
							x2p, y2p, z2p);
		}
		// montgomery reduction
		from_montgomery(x3, *x3p);
		from_montgomery(y3, *y3p);
		from_montgomery(z3, *z3p);
	}
	void point_double(point_t<4>& q, const point_t<4>& p) const noexcept {
		point_double_jacobian(q.xd(), q.yd(), q.zd(),
						p.x.data(), p.y.data(), p.z.data());
	}
	void point_add(point_t<4>& q, const point_t<4>& p1, const point_t<4>& p2)
	const noexcept {
		point_add_jacobian(q.xd(), q.yd(), q.zd(),
						p1.x.data(), p1.y.data(), p1.z.data(),
						p2.x.data(), p2.y.data(), p2.z.data());
	}
	void point_add(point_t<4>& q, const point_t<4>& p1, const felem_t& x2,
			const felem_t& y2) const noexcept
	{
		point_add_jacobian(q.xd(), q.yd(), q.zd(), p1.x.data(), p1.y.data(),
						p1.z.data(), x2.data(), y2.data());
	}
private:
	void carry_reduce(felem_t& res, const u64 carry) const noexcept
	{
		// carry < 2^32
		u64		u = carry & ((1L<<32) -1);
		u64		cc[4];
		cc[0] = u;
		cc[1] = (u << 32) - u;
		cc[2] = 0;
		cc[3] = u << 32;
		if (res.add_to(cc)) res.sub_from(this->p);
	}
};


struct slice_t {
	u64	*data;
	int64_t	len;
	int64_t	cap;
	bool isZero() {
		if (len <= 0) return true;
		for (int i=0;i<len;++i) {
			if (data[i] != 0) return false;
		}
		return true;
	}
	explicit operator bool () { return len != 0; }
	template <const uint N>slice_t(u64 vd[N]) noexcept : data(vd),
			len(N), cap(N)
	{
		vli_clear<N>(data);
	}
	void normal() {
		if (len <= 0) return;
		int	i;
		for(i=len-1;i>=0;i--) {
			if (data[i] != 0) break;
		}
		len = i+1;
	}
};


}	// namespace ecc

#endif	//	__CURVE_IMPL_HPP__
