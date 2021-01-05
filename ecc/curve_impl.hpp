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
			const u64* rrN, const u64 k0P, const u64 k0N,
			const bool a_n3 = true):
		name(_name), gx(_gx), gy(_gy),
		p(_p), n(_n), a(_a), b(_b), rr_p(rrP), rr_n(rrN), k0_p(k0P), k0_n(k0N),
		_a_is_neg3(a_n3), use_barrett(false)
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
		if (unlikely(k0_p == 0))
		{
			// no calc k0 and rr after instantiation
			return false;
		}
		bignum<N>	t1;
		t1.clear();
		_mont_one.sub(t1, p);
		if (unlikely( !_a_is_neg3 )) {
			if (likely(a.is_zero())) _a_is_zero = true;
		}
		// should verify calc K0 and RR
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
	const bignum<N> rr_p = {};
	const bignum<N> rr_n = {};
#ifdef	WITH_BARRETT
	bignum<N+1> mu_p;
	bignum<N+1> mu_n;
#endif
	bignum<N> _mont_one;
	bignum<N> _mont_a;
	const u64	k0_p = 0;
	const u64	k0_n = 0;
	const uint _ndigits = N;
	const bool _a_is_neg3 = false;
	bool _a_is_zero = false;
	const bool use_barrett = false;
	bool _inited = false;
};


template<const uint N>
struct point_t {
	u64		x[N];
	u64		y[N];
	u64		z[N];
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

}	// namespace vli

#endif	//	__CURVE_IMPL_HPP__
