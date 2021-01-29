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
#ifndef	__ECC_KEY_H__
#define __ECC_KEY_H__

#include "cdefs.h"
#include <time.h>
#include <random>
#include "ecc_impl.hpp"


namespace ecc {

using namespace vli;


template<const uint N=4>
class	bn_random {
public:
	bn_random(const bn_random&) = delete;
	static bn_random&	Instance() noexcept {
		static	bn_random	bn_random_(clock());
		return bn_random_;
	}
	bignum<N> get_random() noexcept {
		u64		dd[N];
		for (uint i=0;i<N;++i) dd[i] = _rd();
		return bignum<N>(dd);
	}
private:
	bn_random(u64 _seed) : _rd(_seed) {
	}
	std::mt19937_64	_rd;
};


// definitions for Private/Public key
template<const uint N=4>
class	private_key {
public:
	using felem_t = bignum<N>;
	using public_key = spoint_t<N>;
	private_key() = default;
	template<typename curveT>
	private_key(const curveT& curve)
	{
		auto&	rd = bn_random<N>::Instance();
		felem_t secret;
		do {
			secret = rd.get_random();
		} while ( !init(curve, secret) );
		calc_pa(curve);
	}
	template<typename curveT>
	private_key(const curveT& curve, const felem_t& secret,
				const public_key& pk) : _pa(pk)
	{
		init(curve, secret);
	}
	template<typename curveT>
	bool init(const curveT& curve, const felem_t& secret) noexcept {
		_inited = false;
		curve.modN(_d, secret);
		if (_d.is_zero()) return false;
		_dInv = _d;
		_dInv.uadd_to(1);
		curve.modN(_dInv, _dInv);
		if (_dInv.is_zero()) return false;
		// _dInv = (1 + d)^-1
		mod_inv<N>(_dInv, _dInv, curve.paramN());
		curve.to_montgomeryN(_dInv, _dInv);
#ifdef	WITH_MONT_D
		curve.to_montgomeryN(_mont_d, _d);
#endif
		_inited = true;
		return true;
	}
	explicit operator bool() const noexcept { return _inited; }
	const felem_t&		D() const noexcept { return *(&_d); }
#ifdef	WITH_MONT_D
	const felem_t&		mont_d() const noexcept { return *(&_mont_d); }
#endif
	// Di() return (1 + dA)^-1 in montgomery form modN
	const felem_t&		Di() const noexcept { return *(&_dInv); }
	const public_key&	PubKey() const noexcept { return *(&_pa); }
protected:
	template<typename curveT>
	void calc_pa(const curveT& curve) noexcept {
		if ( unlikely(!_inited) ) return;
		point_t<N>	pt;
		curve.scalar_mult_base(pt, _d);
		_pa.x= pt.x;
		_pa.y = pt.y;
	}
	bool		_inited = false;
	felem_t		_d;			// secret
#ifdef	WITH_MONT_D
	felem_t		_mont_d;	// secret in montgomery form
#endif
	felem_t		_dInv;		// (1+d)^-1  for SM2 sign
	public_key	_pa;
};


template<const uint N=4, typename curveT>
forceinline static
void gen_keypair(const curveT& curve, bignum<N>& secret, bignum<N>& pubX,
		bignum<N>& pubY) noexcept
{
	auto&	rd = bn_random<N>::Instance();
	do {
		secret = rd.get_random();
		curve.modN(secret, secret);
	} while (secret.is_zero());
	point_t<N>	pt;
	curve.scalar_mult_base(pt, secret);
	pubX = pt.x;
	pubY = pt.y;
}

/*
 * SM2 sign/verify
 * 			e = HASH (msg)
 * 			Pa ...	Public key, Px/Py
 * 			dA ...	private key
 * 	sign
 * 			k = random   256bits   Label getk
 * 			x1, y1 = k * G
 * 			r = e + x1
 *			if r == 0 modN  goto getk
 *			s = (1 + dA)^-1 * (k - r * dA)
 *			if s == 0 modN goto getk
 *			if r + s == 0 modN goto getk
 *			return r, s, y1_is_odd
 *
 * 	verify
 * 			t = r + s modN
 * 			if t == 0 fail
 * 			x2, y2 = s*G + t*Pa
 * 			if x2 + e == r modN success else fail
 *
 *	recover public key
 *			Psx = r - e
 *			Psy = get_Py ( Psx, Psy_is_odd )
 *			Ps ...	Psx/Psy
 *			u1 = (r + s)^-1
 *			u2 = s * u1
 *			Pa = (r+s)^1 * ( k*G - s*G) = (r+s)^-1 * Ps - s*(r+s)^-1 * G
 *			Pa = u1 * Ps - u2 * G
 */

// sign return ecInd,  0  for Py(of sign) even, 1 odd
template<const uint N=4, typename curveT>
forceinline static
int ec_sign(const curveT& curve, bignum<N>& r, bignum<N>& s,
		const private_key<N>& priv, const bignum<N>& msg) noexcept
{
	bignum<N>	k, x1;	// y1 reuse tmp
	bignum<N>	tmp;	// r+s
	int		ret=0;
	do {
		do {
			gen_keypair<N>(curve, k, x1, tmp);
			if (r.add(msg, x1)) r.sub_from(curve.paramN());
			curve.modN(r, r);
		} while (r.is_zero());
		ret = tmp.is_odd();
#ifdef	WITH_MONT_D
		bignum<N>	rp;	// r+s
		curve.to_montgomeryN(k, k);
		curve.to_montgomeryN(rp, r);
		// tmp = r * dA
		curve.mont_nmult(tmp, rp, priv.mont_d());
		// k = k - r * dA
		if (k.sub_from(tmp)) k.add_to(curve.paramN());
		// tmp = (1 + dA)^-1 * (k - r * dA)
		curve.mont_nmult(tmp, priv.Di(), k);
		curve.from_montgomeryN(s, tmp);
#else
		if (tmp.add(k, r)) tmp.sub_from(curve.paramN());
		curve.to_montgomeryN(tmp, tmp);
		curve.mont_nmult(tmp, tmp, priv.Di());
		curve.from_montgomeryN(tmp, tmp);
		if (tmp.sub_from(r)) tmp.add_to(curve.paramN());
		curve.modN(s, tmp);
#endif
		if (s.is_zero()) continue;
		if (tmp.add(r, s)) tmp.sub_from(curve.paramN());
		curve.modN(tmp, tmp);
	} while (tmp.is_zero());
	return ret;
}

// verify return bool, true for success(ok)
template<const uint N=4, typename curveT>
forceinline static
bool ec_verify(const curveT& curve, const bignum<N>& r, const bignum<N>& s,
		const spoint_t<N>&	pub, const bignum<N>& msg) noexcept
{
	bignum<N>	t;
	if ( unlikely(r.is_zero()) ) return false;
	if ( unlikely(s.is_zero()) ) return false;
	if (t.add(r, s)) t.sub_from(curve.paramN());
	curve.modN(t, t);
	if ( unlikely(t.is_zero()) ) return false;
	point_t<N>	q, p(pub.x, pub.y);
	// q = s*G + t * Pub
	// x2, y2 = s*G + t * Pub
	curve.combined_mult(q, p, t, s);
	if ( unlikely(q.x.is_zero()) ) return false;
	if (t.add(q.x, msg)) t.sub_from(curve.paramN());
	// t = x2 + msg modN
	curve.modN(t, t);
	if ( unlikely(t.is_zero()) ) return false;
	return r == t;
}

// vecover public key, from r/s/v/msg, v for pubY is odd
// return true for successfully recover public key
template<const uint N=4, typename curveT>
forceinline static
bool ec_recover(const curveT& curve, spoint_t<N>&  pub, const bignum<N>& r,
				const bignum<N>& s, const int v, const bignum<N>& msg) noexcept
{
	point_t<N>	p;
	if (p.x.sub(r, msg)) p.x.add_to(curve.paramN());
	curve.modN(p.x, p.x);
	if ( unlikely(!pointY_recover(curve, p.y, p.x, v)) ) return false;
	bignum<N>	u1, u2;
	if (u1.add(r, s)) u1.sub_from(curve.paramN());
	curve.modN(u1, u1);
	if ( unlikely(u1.is_zero()) ) return false;
	// u1 = (r + s)^1
	mod_inv<N>(u1, u1, curve.paramN());
	bignum<N>	u1p, sp;
	curve.to_montgomeryN(sp, s);
	curve.to_montgomeryN(u1p, u1);
	curve.mont_nmult(u2, sp, u1p);
	// u2 = s * u1
	curve.from_montgomeryN(u2, u2);
	// 0 < u2 < N 
	if ( unlikely(u2.is_zero()) ) return false;
	// u2 = -u2 = - s* u1
	u2.sub(curve.paramN(), u2);
	// Pa = u1 * Ps - u2 * G
	point_t<N>	q;
	curve.combined_mult(q, p, u1, u2);
	pub.x = q.x;
	pub.y = q.y;
	return true;
}

}	// namespace ecc

#endif	//	__ECC_KEY_H__
