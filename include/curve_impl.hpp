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
void pre_compute(const curveT& curve, point_t<N> pre_comp[wSize],
		const point_t<N>& p) noexcept
{
	pre_comp[0] = p;
	curve.point_double(pre_comp[1], p.x, p.y);
	for (int i = 2; i < wSize; ++i) {
		if ((i & 1) == 0) {
			curve.point_add(pre_comp[i], pre_comp[i-1], p.x, p.y);
			//curve.point_add(pre_comp[i], pre_comp[i-1], p);
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
	static_assert(N == 4, "only 256 bits curve supported");
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
	// convert to affine jacobian, in montgomery bignum
	for(i=0;i < 2; i++) {
		for (j = 1; j < 16; j++) {
			curve.apply_z_mont(pre_comps[i][j]);
		}
	}
}

#ifdef	ommit
template<const uint N, typename curveT> forceinline void
initTable(const curveT& curve, bignum<N>  pre_comps[43][64]) noexcept
{
	int i, j;
	point_t<N>	G;
	static_assert(N == 4, "only 256 bits curve supported");
	curve.to_montgomery(G.x, curve.getGx());
	curve.to_montgomery(G.y, curve.getGy());
	G.z = curve.mont_one();
	p256Precomputed = new([43][32 * 8]uint64)

	basePoint := []uint64{
		0x61328990f418029e, 0x3e7981eddca6c050, 0xd6a1ed99ac24c3c3, 0x91167a5ee1c13b05,
		0xc1354e593c2d0ddd, 0xc1f5e5788d3295fa, 0x8d4cfb066e2a48f8, 0x63cd65d481d735bd,
		0x0000000000000001, 0xffffffff00000000, 0xffffffffffffffff, 0x00000000fffffffe,
	}
	t1 := make([]uint64, 12)
	t2 := make([]uint64, 12)
	copy(t2, basePoint)

	zInv := make([]uint64, 4)
	zInvSq := make([]uint64, 4)
	for j := 0; j < 32; j++ {
		copy(t1, t2)
		for i := 0; i < 43; i++ {
			// The window size is 6 so we need to double 6 times.
			if i != 0 {
				for k := 0; k < 6; k++ {
					p256PointDoubleAsm(t1, t1)
				}
			}
			// Convert the point to affine form. (Its values are
			// still in Montgomery form however.)
			p256Inverse(zInv, t1[8:12])
			p256Sqr(zInvSq, zInv, 1)
			p256Mul(zInv, zInv, zInvSq)

			p256Mul(t1[:4], t1[:4], zInvSq)
			p256Mul(t1[4:8], t1[4:8], zInv)

			copy(t1[8:12], basePoint[8:12])
			// Update the table entry
			copy(p256Precomputed[i][j*8:], t1[:8])
		}
		if j == 0 {
			p256PointDoubleAsm(t2, basePoint)
		} else {
			p256PointAddAsm(t2, t2, basePoint)
		}
	}
}
#endif

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
template<const uint N=4, const bool A_is_n3=true>
class ecc_curve {
public:
	using felem_t = bignum<N>;
	ecc_curve(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
			const u64 *_n, const u64 *_a, const u64 *_b, const u64* rrP,
			const u64* rrN, const u64 k0P, const u64 k0N):
		name(_name), gx(_gx), gy(_gy),
		p(_p), n(_n), a(_a), b(_b), rr_p(rrP), rr_n(rrN), k0_p(k0P), k0_n(k0N),
		_a_is_neg3(A_is_n3)
	{
		init();
		static_assert(N > 3, "curve only support 256Bits or more");
	}
	ecc_curve(ecc_curve &&) = default;
	static const ecc_curve* new_ecc_curve(const char *_name, const u64 *_gx,
			const u64 *_gy, const u64 *_p, const u64 *_n, const u64 *_a,
			const u64 *_b) noexcept
	{
		static_assert(N > 3, "curve only support 256Bits or more");
		bignum<N>	prr, nrr;
		u64		k0P, k0N;
		bignum<N>	prime(_p);
		bignum<N>	n_prime(_n);
		k0P = calcK0<N>(prime);
		calcRR<N>(prr, prime);
		k0N = calcK0<N>(n_prime);
		calcRR<N>(nrr, n_prime);
		return new ecc_curve(_name, _gx, _gy, _p, _n, _a, _b, prr.data(),
						nrr.data(), k0P, k0N);
	}
	const uint ndigits() const { return N; }
	explicit operator bool() const noexcept
	{
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
	const bool a_is_pminus3() const noexcept { return _a_is_neg3; }
	const bool a_is_zero() const noexcept { return _a_is_zero; }
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
	const noexcept
	{
		res.mont_mult(left, right, p, k0_p);
	}
	void mont_msqr(felem_t& res, const felem_t left, const uint nTimes=1)
	const noexcept
	{
		res.mont_sqr(left, p, k0_p);
		for (uint i=1; i < nTimes; i++) res.mont_sqr(res, p, k0_p);
	}
	// left,right less than p
	void mod_add(felem_t& res, const felem_t& left, const felem_t& right)
	const noexcept
	{
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
	const noexcept
	{
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
#ifndef	ommit
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
		if ( z_is_one ) {
			zp = mont_one();
		} else {
			to_montgomery(zp, z1);
		}
		if (z_is_one)
			point_doublez_jacob<A_is_n3>(*this, *x3p, *y3p, *z3p, xp, yp);
		else
			point_double_jacob<A_is_n3>(*this, *x3p, *y3p, *z3p, xp, yp, zp);
		// montgomery reduction
		from_montgomery(x3, *x3p);
		from_montgomery(y3, *y3p);
		from_montgomery(z3, *z3p);
	}
#ifdef	ommit
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
#endif
	void apply_z(bignum<N>& x1, bignum<N>& y1, bignum<N>& z) const noexcept
	{
		felem_t	t1;

		if (z.is_one()) return;
		if (z.is_zero()) {
			x1.clear();
			y1.clear();
			return;
		}
		mod_inv<N>(z, z, this->p);
		to_montgomery(z, z);
		this->mont_msqr(t1, z);	// t1 = z^2
		to_montgomery(x1, x1);
		to_montgomery(y1, y1);
		this->mont_mmult(x1, x1, t1);	// x1 * z^2
		this->mont_mmult(t1, t1, z);	// t1 = z^3
		this->mont_mmult(y1, y1, t1);	// y1 * z^3

		// montgomery reduction
		from_montgomery(x1, x1);
		from_montgomery(y1, y1);
		// set z to 1
		z = felem_t(1);
	}
	void apply_z(felem_t& x, felem_t& y, const point_t<N>& pt) const noexcept
	{
		felem_t	z(pt.z);
		x = pt.x;
		y = pt.y;
		apply_z(x, y, z);
	}
	void apply_z(point_t<N>& p) const noexcept
	{
		apply_z(p.x, p.y, p.z);
	}
	void apply_z_mont(point_t<N>& pt) const noexcept
	{
		if (pt.z.is_zero() || pt.z == mont_one()) return;
		felem_t	t1;
		from_montgomery(pt.z, pt.z);
		mod_inv<N>(pt.z, pt.z, this->p);
		to_montgomery(pt.z, pt.z);		// z = p.z^-1
		this->mont_msqr(t1, pt.z);	// t1 = z^2
		this->mont_mmult(pt.x, pt.x, t1);	// x1 * z^2
		this->mont_mmult(t1, t1, pt.z);	// t1 = z^3
		this->mont_mmult(pt.y, pt.y, t1);	// y1 * z^3
		pt.z = this->mont_one();
	}
	bool point_eq(const point_t<N>& p, const point_t<N>& q) const noexcept
	{
		if (p == q) return true;
		if ( likely(p.z.is_one()) ) {
			felem_t	x2, y2;
			this->apply_z(x2, y2, q);
			return p.x == x2 && p.y == y2;
		} else {
			felem_t	x1, x2, y1, y2;
			this->apply_z(x1, y1, p);
			this->apply_z(x2, y2, q);
			return x1 == x2 && y1 == y2;
		}
	}
	void point_neg(point_t<N>& q, const point_t<N>& p) const noexcept
	{
		q.x = p.x;
		q.y.sub(this->p, p.y);
		q.z = p.z;
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
			point_addz_jacob<A_is_n3>(*this, *x3p, *y3p, *z3p, x1p, y1p, z1p,
							x2p, y2p);
		} else {
			point_add_jacob<A_is_n3>(*this, *x3p, *y3p, *z3p, x1p, y1p, z1p,
							x2p, y2p, z2p);
		}
		// montgomery reduction
		from_montgomery(x3, *x3p);
		from_montgomery(y3, *y3p);
		from_montgomery(z3, *z3p);
	}
	void scalar_mult(point_t<N>& q, const spoint_t<N>& p, const felem_t& scalar)
			const noexcept
	{
		point_t<N>	tmp;
		uint	nbits = scalar.num_bits();
		q.clear();
		if ( unlikely(nbits == 0) ) return;
#ifdef	PRECOMPUTE_INSTACK
		point_t<N>	pres[wSize];
#else
		auto pres = new(point_t<N>[wSize]);
#endif
		to_montgomery(tmp.x, p.x);
		to_montgomery(tmp.y, p.y);
		tmp.z = this->mont_one();
		pre_compute<N>(*this, pres, tmp);
		bool	skip = true;
		if (nbits % W != 0) --nbits;
		for (int i = nbits; i >= 0; --i) {
			if (!skip) point_double(q, q);
			if (i % W != 0) continue;
			uint	bits;
			uint	digit;
			bits = vli_get_bits<N, W+1>(scalar.data(), i-1);
			auto sign = recode_scalar_bits<W>(digit, bits);
			if (digit == 0) continue;
			--digit;
#ifndef	NO_CONDITIONAL_COPY
			tmp = pres[digit];
			felem_t	ny;
			ny.sub(this->p, tmp.y);
			ny.copy_conditional(tmp.y, sign-1);
			tmp.y = ny;
#else
			if ( sign ) point_neg(tmp, pres[digit]); else tmp = pres[digit];
#endif
			if (!skip) point_add(q, q, tmp); else {
				q = tmp;
				skip =false;
			}
		}
#ifndef	PRECOMPUTE_INSTACK
		delete []pres;
#endif
	}
	void scalar_mult(point_t<N>& q, const point_t<N>& p, const felem_t& scalar)
			const noexcept
	{
		spoint_t<N>	pp(p.x, p.y);
		scalar_mult(q, pp, scalar);
		// montgomery reduction
		if ( unlikely(q.z.is_zero()) ) return;
		this->apply_z_mont(q);
		this->from_montgomery(q.x, q.x);
		this->from_montgomery(q.y, q.y);
		q.z = felem_t(1);
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
		if ( unlikely(q.z.is_zero()) ) return;
		this->apply_z_mont(q);
		this->from_montgomery(q.x, q.x);
		this->from_montgomery(q.y, q.y);
		q.z = felem_t(1);
	}
	void point_double(point_t<N>& q, const point_t<N>& p) const noexcept {
		point_double_jacob<A_is_n3>(*this, q.x, q.y, q.z, p.x, p.y, p.z);
	}
	void point_double(point_t<N>& q, const felem_t& x1, const felem_t& y1)
	const noexcept {
		point_doublez_jacob<A_is_n3>(*this, q.x, q.y, q.z, x1, y1);
	}
	void point_add(point_t<N>& q, const point_t<N>& p1, const point_t<N>& p2)
	const noexcept {
		point_add_jacob<A_is_n3>(*this, q.x, q.y, q.z, p1.x, p1.y, p1.z,
						p2.x, p2.y, p2.z);
	}
	void point_add(point_t<N>& q, const point_t<N>& p1, const felem_t& x2,
					const felem_t& y2) const noexcept
	{
		point_addz_jacob<A_is_n3>(*this, q.x, q.y, q.z, p1.x, p1.y, p1.z,
						x2, y2);
	}
protected:
	bool init() noexcept
	{
		if (_inited) return _inited;
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
		_inited = true;
		return true;
	}
	const std::string name;
	const felem_t gx;
	const felem_t gy;
	const felem_t p;
	const felem_t n;
	const felem_t a;
	const felem_t b;
	const felem_t rr_p = {};
	const felem_t rr_n = {};
	felem_t _mont_one;
	felem_t _mont_a;
	const u64	k0_p = 0;
	const u64	k0_n = 0;
	const uint _ndigits = N;
	const bool _a_is_neg3 = false;
	bool _a_is_zero = false;
	const bool use_barrett = false;
	bool _inited = false;
};

template <const u64 Pk0=1, const bool A_is_n3=true>
class curve256 : public ecc_curve<4,A_is_n3> {
public:
	using felem_t = bignum<4>;
	curve256(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
			const u64 *_n, const u64 *_a, const u64 *_b) :
		ecc_curve<4>(_name, _gx, _gy, _p, _n, _a, _b, sm2_p_rr, sm2_n_rr,
					sm2_p_k0, sm2_n_k0)
	{ }
	const felem_t& mont_one() const noexcept { return this->_mont_one; }
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
	void apply_z_mont(point_t<4>& pt) const noexcept
	{
		if (pt.z.is_zero() || pt.z == mont_one()) return;
		felem_t	t1, z;
		from_montgomery(z, pt.z);
		mod_inv<4>(z, z, this->p);
		to_montgomery(z, z);		// z = p.z^-1
		this->mont_msqr(t1, z);	// t1 = z^2
		this->mont_mmult(pt.x, pt.x, t1);	// x1 * z^2
		this->mont_mmult(t1, t1, z);	// t1 = z^3
		this->mont_mmult(pt.y, pt.y, t1);	// y1 * z^3
		pt.z = this->mont_one();
	}
	void mont_mmult(felem_t& res, const felem_t& left, const felem_t& right)
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
	void g_precompute() noexcept
	{
		if (! this->init() ) return;
		pre_compute_base<4>(*this, g_pre_comp);	// precompute g_pre_comp

	}
	bool select_base_point(point_t<4>& pt, const uint idx) const noexcept {
		if (idx >= 2 * 16) return false;
		pt = g_pre_comp[idx>>4][idx&0xf];
		return true;
	}
	// G multiples with 0..31 32 bits, process one bit
	bool base_mult_bit(point_t<4>& q, const felem_t& scalar, const uint bit,
			const bool _skip) const noexcept
	{
		point_t<4>	p;
		bool skip = _skip;
		if ( unlikely(bit >= 32) ) return skip;
		/* add multiples of the generator */
		u64	bits = scalar.get_bit(bit + 224) << 3;
		bits |= scalar.get_bit(bit + 160) << 2;
		bits |= scalar.get_bit(bit + 96) << 1;
		bits |= scalar.get_bit(bit + 32);
		if ( likely(bits != 0) ) {
			this->select_base_point(p, bits | 16);
			if (!skip) point_add(q, q, p.x, p.y); else {
				q = p;
				skip = false;
			}
		}
		/* second, look at the current position */
		bits = scalar.get_bit(bit + 192) << 3;
		bits |= scalar.get_bit(bit + 128) << 2;
		bits |= scalar.get_bit(bit + 64) << 1;
		bits |= scalar.get_bit(bit);
		if (bits == 0) return skip;
		this->select_base_point(p, bits);
		if (!skip) point_add(q, q, p.x, p.y); else {
			q = p;
			skip = false;
		}
		return skip;
	}
	void scalar_mult_base(point_t<4>& q, const felem_t& scalar) const noexcept
	{
		bool	skip = true;
		q.clear();
		for (int i = 31; i >= 0; --i) {
			/* double */
			if (!skip) point_double(q, q);

			skip = this->base_mult_bit(q, scalar, i, skip);
		}
		if ( unlikely(q.z.is_zero()) ) return;
		this->apply_z_mont(q);
		this->from_montgomery(q.x, q.x);
		this->from_montgomery(q.y, q.y);
		q.z = felem_t(1);
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
		zp = this->mont_one();
		if (!z_is_one) {
			to_montgomery(zp, z1);
		}
		if (z_is_one)
			point_doublez_jacob<A_is_n3>(*this, *x3p, *y3p, *z3p, xp, yp);
		else
			point_double_jacob<A_is_n3>(*this, *x3p, *y3p, *z3p, xp, yp, zp);
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
		if (vli_is_zero<4>(z1)) {
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
			point_addz_jacob<A_is_n3>(*this, *x3p, *y3p, *z3p, x1p, y1p, z1p,
							x2p, y2p);
		} else {
			point_add_jacob<A_is_n3>(*this, *x3p, *y3p, *z3p, x1p, y1p, z1p,
							x2p, y2p, z2p);
		}
		// montgomery reduction
		from_montgomery(x3, *x3p);
		from_montgomery(y3, *y3p);
		from_montgomery(z3, *z3p);
	}
	// scalar may be zero, g_scalar may be zero
	void combined_mult(point_t<4>& q, const point_t<4>& p,
				const felem_t& scalar, const felem_t& g_scalar) const noexcept
	{
		point_t<4>	tmp;
		uint	nbits = scalar.num_bits();
		bool	b_gscalar = !(g_scalar.is_zero());
		q.clear();
		if ( unlikely(nbits == 0) && !b_gscalar) return;
		point_t<4>	pres[wSize];
		if (nbits > 0) {
			to_montgomery(tmp.x, p.x);
			to_montgomery(tmp.y, p.y);
			if ( likely(p.z.is_one()) ) tmp.z = this->mont_one(); else
				to_montgomery(tmp.z, p.z);
			pre_compute<4>(*this, pres, tmp);
		}
		bool	skip = true;
		if (nbits % W != 0) --nbits;
		int	cntBits = (b_gscalar && nbits < 31)?31:nbits;
		if (b_gscalar && nbits < 31) nbits = 31;
		for (int i = cntBits; i >= 0; --i) {
			if (!skip) point_double(q, q);

			/* add multiples of the generator */
			skip = this->base_mult_bit(q, g_scalar, i, skip);
			if (nbits == 0) continue;
			if ((i % W) == 0) {
				uint	bits;
				uint	digit;
				bits = vli_get_bits<4, W+1>(scalar.data(), i-1);
				auto sign = recode_scalar_bits<W>(digit, bits);
				if (digit == 0) continue;
				--digit;
#ifndef	NO_CONDITIONAL_COPY
				tmp = pres[digit];
				felem_t	ny;
				ny.sub(this->p, tmp.y);
				ny.copy_conditional(tmp.y, sign-1);
				tmp.y = ny;
#else
				if ( sign ) point_neg(tmp, pres[digit]); else tmp = pres[digit];
#endif
				if (!skip) point_add(q, q, tmp); else {
					q = tmp;
					skip =false;
				}
			}
		}
		// montgomery reduction
		if ( unlikely(q.z.is_zero()) ) return;
		this->apply_z_mont(q);
		this->from_montgomery(q.x, q.x);
		this->from_montgomery(q.y, q.y);
		q.z = felem_t(1);
	}
	void point_double(point_t<4>& q, const point_t<4>& p) const noexcept {
		if ( p.z == this->mont_one() )
			point_doublez_jacob<A_is_n3>(*this, q.x, q.y, q.z, p.x, p.y);
		else
			point_double_jacob<A_is_n3>(*this, q.x, q.y, q.z, p.x, p.y, p.z);
	}
	void point_double(point_t<4>& q, const felem_t& x1, const felem_t& y1)
	const noexcept {
		point_doublez_jacob<A_is_n3>(*this, q.x, q.y, q.z, x1, y1);
	}
	void point_add(point_t<4>& q, const point_t<4>& p1, const point_t<4>& p2)
	const noexcept {
		point_add_jacob<A_is_n3>(*this, q.x, q.y, q.z, p1.x, p1.y, p1.z,
						p2.x, p2.y, p2.z);
	}
	void point_add(point_t<4>& q, const point_t<4>& p1, const felem_t& x2,
			const felem_t& y2) const noexcept
	{
		point_addz_jacob<A_is_n3>(*this, q.x, q.y, q.z, p1.x, p1.y, p1.z, x2, y2);
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
	point_t<4>	g_pre_comp[2][16];
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
