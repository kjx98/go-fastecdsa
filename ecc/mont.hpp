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
#ifndef __MONT_HPP__
#define __MONT_HPP__

#include "vli.hpp"
#include "vli_bn.hpp"
#include "ecc_impl.hpp"


namespace vli {

/**
 * struct ecc_curve - definition of elliptic curve
 *
 * @name:	Short name of the curve.
 * @g:		Generator point of the curve.
 * @p:		Prime number, if Barrett's reduction is used for this curve
 *		pre-calculated value 'mu' is appended to the @p after ndigits.
 *		Use of Barrett's reduction is heuristically determined in
 *		vli_mmod_fast().
 * @n:		Order of the curve group.
 * @a:		Curve parameter a.
 * @b:		Curve parameter b.
 */
template<const uint N=4>
class ecc_curve {
public:
	using mont1_func=std::function<void(bignum<N>&, const bignum<N>&)>;
	using mont2_func=std::function<void(bignum<N>&, const bignum<N>&, const bignum<N>&)>;
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
			const u64* rrN, const u64 k0P, const u64 k0N, const bool a_n3 = true,
			bool is_a0 = false) :
		name(_name), gx(_gx), gy(_gy),
		p(_p), n(_n), a(_a), b(_b), rr_p(rrP), rr_n(rrN), k0_p(k0P), k0_n(k0N),
		_a_is_neg3(a_n3), _a_is_zero(is_a0), use_barrett(false)
	{
		static_assert(N > 3, "curve only support 256Bits or more");
		//static_assert(_name != nullptr && _name[0] != 0, "curve name MUST not empty");
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
	const bignum<N>& montParamA() const noexcept { return _mont_a; }
	const bignum<N>& paramP() const noexcept { return p; }
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
#if	__cplusplus >= 201703L
		if constexpr(k0_p == 0)
#else
		if (k0_p == 0)
#endif
		{
			return false;
		}
		bignum<N>	t1;
		t1.clear();
		_mont_one.sub(t1, p);
		// should be calc K0 and RR
		_inited = true;
		return true;
#endif
	}
	const bool a_is_pminus3() const noexcept { return _a_is_neg3; }
	const bool a_is_zero() const noexcept { return _a_is_zero; }
#ifdef	WITH_BARRETT
	void mmod_barrett(bignum<N>& res, const bn_prod<N>& prod) const noexcept
	{
		if (use_barrett) {
			// res = barrettmod
			prod.mmod_barrett(res, p, mu_p);
		}
		return;
	}
#endif
	const bignum<N>& mont_one() const noexcept { return _mont_one; }
	void to_montgomery(bignum<N>& res, const u64 *x) const noexcept
	{
		bignum<N>   *xx = reinterpret_cast<bignum<N> *>(const_cast<u64 *>(x));
		_mont_mult(res, *xx, rr_p);
	}
	void to_montgomery(bignum<N>& res, const bignum<N>& x) const noexcept
	{
		_mont_mult(res, x, rr_p);
	}
	void from_montgomery(bignum<N>& res, const bignum<N>& y) const noexcept
	{
		_mont_reduction(res, y);
	}
	void from_montgomery(u64* result, const bignum<N>& y) const noexcept
	{
		bignum<N>   *res = reinterpret_cast<bignum<N> *>(result);
		_mont_reduction(*res, y);
	}
	void
	mont_mult(bignum<N>& res, const bignum<N>& left, const bignum<N>& right)
	const noexcept {
		_mont_mult(res, left, right);
	}
	void mont_sqr(bignum<N>& res, const bignum<N> left, const uint nTimes=1)
	const noexcept {
		_mont_sqr(res, left);
		for (uint i=1; i < nTimes; i++) _mont_sqr(res, res);
	}
	// left,right less than p
	void mod_add(bignum<N>& res, const bignum<N>& left, const bignum<N>& right)
	const noexcept {
		if (res.add(left, right)) {
			res.sub_from(p);
		}
	}
	// left,right less than p
	void mod_add_to(bignum<N>& res, const bignum<N>& right) const noexcept
	{
		if (res.add_to(right)) {
			res.sub_from(p);
		}
	}
	// left,right less than p
	void mod_sub(bignum<N>& res, const bignum<N>& left, const bignum<N>& right)
	const noexcept {
		if (res.sub(left, right)) {
			res.add_to(p);
		}
	}
	void mod_sub_from(bignum<N>& res, const bignum<N>& right) const noexcept
	{
		if (res.sub_from(right)) {
			res.add_to(p);
		}
	}
	void mont_mult2(bignum<N>& res, const bignum<N>& left) const noexcept
	{
#ifdef	ommit
		this->mod_add(res, left, left);
#else
		if (res.lshift1(left) != 0 || res.cmp(p) >= 0) {
			res.sub_from(p);
		}
#endif
	}
private:
	mont1_func	_mont_reduction = [this](bignum<N> &res, const bignum<N> &y) {
		res.mont_reduction(y, p, k0_p);
	};
	mont1_func	_mont_sqr = [this](bignum<N> &res, const bignum<N> &y) {
		res.mont_sqr(y, p, k0_p);
	};
	mont2_func	_mont_mult = [this](bignum<N> &res, const bignum<N>& x,
					const bignum<N>& y) {
		res.mont_mult(x, y, p, k0_p);
	};
	const std::string name;
	const bignum<N> gx;
	const bignum<N> gy;
	const bignum<N> p;
	const bignum<N> n;
	const bignum<N> a;
	const bignum<N> b;
	const bignum<N> rr_p;
	const bignum<N> rr_n;
	const bignum<N+1> mu_p;
	const bignum<N+1> mu_n;
	bignum<N> _mont_one;
	bignum<N> _mont_a;
	const u64	k0_p = 0;
	const u64	k0_n = 0;
	const uint _ndigits = N;
	const bool _a_is_neg3 = false;
	const bool _a_is_zero = false;
	const bool use_barrett = false;
	bool _inited = false;
};


template<uint ndigits>
struct point_t {
	u64		x[ndigits];
	u64		y[ndigits];
	u64		z[ndigits];
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
	template <uint ndigits>slice_t(u64 vd[ndigits]) noexcept : data(vd),
			len(ndigits), cap(ndigits)
	{
		vli_clear<ndigits>(data);
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

}	// namespace vli


// SM2 prime optimize
// p is 2^256 - 2^224 - 2^96 + 2^64 -1
forceinline
static void vli_sm2_multP(u64 *result, const u64 u) noexcept
{
	u64	r[4];
	u64	t_low, t_high;
	t_low = u << 32;	// ^192
	t_high = ((u >> 32) & 0xffffffff);
	vli_clear<6>(result);
	result[3] = 0 - t_low;		// ^256 - ^224
	result[4] = u - t_high -1;
	r[0] =0;
	r[1] = t_low;
	r[2] = t_high;
	r[3] =0;
	//vli_sub_from<5>(result, t);	// ^256 - ^224 - ^96
	if (vli_sub_from<4>(result, r)) result[4]--;
	r[2] = 0;
	r[1] = u-1;
	r[0] = -u;		// ^64 -1
	if (vli_add_to<4>(result, r)) result[4]++;
}

#ifdef	ommit
forceinline
static void mont_reductionP(u64 *result, const u64 *y, const u64 *prm) noexcept
{
	u64	s[8];
	u64	r[6];
	vli_clear<6>(r);
#pragma GCC unroll 4
	for (uint i=0; i < 4; i++) {
		u64	u = r[0] + y[i];
#ifdef	WITH_SM2_MULTP
		vli_sm2_multP(s, u);
#else
		vli_umult<4>(s, prm, u);
#endif
		vli_uadd_to<6>(r, y[i]);
		vli_add_to<6>(r, s);
		vli_rshift1w<6>(r);	
	}
	if (r[4] !=0 || vli_cmp<4>(r, prm) >= 0) {
		vli_sub<4>(result, r, prm);
	} else vli_set<4>(result, r);
}

forceinline
static void mont_multP(u64 *result, const u64 *x, const u64 *y,
			const u64 *prime) noexcept
{
	u64	s[8];
	u64	r[6];
	vli_clear<6>(r);
#pragma GCC unroll 4
	for (uint i=0; i < 4;i++) {
		u64	u = r[0] + y[i]*x[0];
#ifdef	WITH_SM2_MULTP
		vli_sm2_multP(s, u);
#else
		vli_umult<4>(s, prime, u);
#endif
		vli_add_to<6>(r, s);
		vli_umult<4>(s, x, y[i]);
		vli_add_to<6>(r, s);
		vli_rshift1w<6>(r);	
	}
	if (r[4] != 0 || vli_cmp<4>(r, prime) >= 0) {
		vli_sub<4>(result, r, prime);
	} else vli_set<4>(result, r);
}
#endif	// NO_BUILTIN_OVERFLOW

template<uint ndigits> forceinline
static void
#ifdef	WITH_C2GO
mont_reduction(u64 *result, const u64 *y, const u64 *prime,
			const u64 k0, u64 *buff) noexcept
#else
mont_reduction(u64 *result, const u64 *y, const u64 *prime,
			const u64 k0) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*s = buff;
	u64	*r = s + ndigits * 2;
#else
	u64	s[ndigits * 2];
	u64	r[ndigits + 2];
#endif
	vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < ndigits; i++) {
		u64	u = (r[0] + y[i]) * k0;
		vli_umult<ndigits>(s, prime, u);
		vli_uadd_to<ndigits + 2>(r, y[i]);
		vli_add_to<ndigits + 2>(r, s);
		vli_rshift1w<ndigits + 2>(r);	
	}
	if (r[ndigits] !=0 || vli_cmp<ndigits>(r, prime) >= 0) {
		vli_sub<ndigits>(result, r, prime);
	} else vli_set<ndigits>(result, r);
}

template<uint ndigits> forceinline
static void
#ifdef	WITH_C2GO
mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 k0, u64 *buff) noexcept
#else
mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 k0) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*s = buff;
	u64	*r = s + ndigits * 2;
#else
	u64	s[ndigits * 2];
	u64	r[ndigits + 2];
#endif
	vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < ndigits;i++) {
		u64	u = (r[0] + y[i]*x[0]) * k0;
		vli_umult<ndigits>(s, prime, u);
		vli_add_to<ndigits + 2>(r, s);
		vli_umult<ndigits>(s, x, y[i]);
		vli_add_to<ndigits + 2>(r, s);
		vli_rshift1w<ndigits + 2>(r);	
	}
	if (r[ndigits] != 0 || vli_cmp<ndigits>(r, prime) >= 0) {
		vli_sub<ndigits>(result, r, prime);
	} else vli_set<ndigits>(result, r);
}

template<uint ndigits> forceinline
static void
#ifdef	WITH_C2GO
mont_sqr(u64 *result, const u64 *x, const u64 *prime, const u64 k0, u64 *buff) noexcept
#else
mont_sqr(u64 *result, const u64 *x, const u64 *prime, const u64 k0) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*s = buff;
	u64	*r = s + ndigits * 2;
#else
	u64	s[ndigits * 2];
	u64	r[ndigits + 2];
#endif
	vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < ndigits;i++) {
		u64	u = (r[0] + x[i]*x[0]) * k0;
		vli_umult<ndigits>(s, prime, u);
		vli_add_to<ndigits + 2>(r, s);
		vli_umult<ndigits>(s, x, x[i]);
		vli_add_to<ndigits + 2>(r, s);
		vli_rshift1w<ndigits + 2>(r);	
	}
	if (r[ndigits] != 0 || vli_cmp<ndigits>(r, prime) >= 0) {
		vli_sub<ndigits>(result, r, prime);
	} else vli_set<ndigits>(result, r);
}

template<const uint N, const u64 k0> forceinline
static void
mont_reduction(u64 *result, const u64 *y, const u64 *prime) noexcept
{
	u64	s[N * 2];
	u64	r[N + 2];
	vli_clear<N + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < N; i++) {
		u64	u = (r[0] + y[i]) * k0;
		vli_umult<N>(s, prime, u);
		vli_uadd_to<N + 2>(r, y[i]);
		vli_add_to<N + 2>(r, s);
		vli_rshift1w<N + 2>(r);	
	}
	if (r[N] !=0 || vli_cmp<N>(r, prime) >= 0) {
		vli_sub<N>(result, r, prime);
	} else vli_set<N>(result, r);
}

template<const uint ndigits, const u64 k0> forceinline
static void
mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime) noexcept
{
	u64	s[ndigits * 2];
	u64	r[ndigits + 2];
	vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < ndigits;i++) {
		u64	u = (r[0] + y[i]*x[0]) * k0;
		vli_umult<ndigits>(s, prime, u);
		vli_add_to<ndigits + 2>(r, s);
		vli_umult<ndigits>(s, x, y[i]);
		vli_add_to<ndigits + 2>(r, s);
		vli_rshift1w<ndigits + 2>(r);	
	}
	if (r[ndigits] != 0 || vli_cmp<ndigits>(r, prime) >= 0) {
		vli_sub<ndigits>(result, r, prime);
	} else vli_set<ndigits>(result, r);
}

template<const uint ndigits, const u64 k0> forceinline
static void
mont_sqr(u64 *result, const u64 *x, const u64 *prime) noexcept
{
	u64	s[ndigits * 2];
	u64	r[ndigits + 2];
	vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < ndigits;i++) {
		u64	u = (r[0] + x[i]*x[0]) * k0;
		vli_umult<ndigits>(s, prime, u);
		vli_add_to<ndigits + 2>(r, s);
		vli_umult<ndigits>(s, x, x[i]);
		vli_add_to<ndigits + 2>(r, s);
		vli_rshift1w<ndigits + 2>(r);	
	}
	if (r[ndigits] != 0 || vli_cmp<ndigits>(r, prime) >= 0) {
		vli_sub<ndigits>(result, r, prime);
	} else vli_set<ndigits>(result, r);
}

#endif	//	__MONT_HPP__
