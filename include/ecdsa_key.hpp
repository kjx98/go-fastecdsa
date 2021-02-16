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
#ifndef	__ECDSA_KEY_H__
#define __ECDSA_KEY_H__

#include "cdefs.h"
#include <time.h>
#include <unistd.h>
#include <random>
#ifdef	__x86_64__
#include <x86intrin.h>
#endif
#include "ecc_key.hpp"


namespace ecc {

using namespace vli;


/*
 * ECC private keys are generated using the method of extra random bits,
 * equivalent to that described in FIPS 186-4, Appendix B.4.1.
 *
 * d = (c mod(n–1)) + 1    where c is a string of random bits, 64 bits longer
 *                         than requested
 * 0 <= c mod(n-1) <= n-2  and implies that
 * 1 <= d <= n-1
 *
 * This method generates a private key uniformly distributed in the range
 * [1, n-1].
 */
// definitions for Private/Public key
template<const uint N=4>
class	ecdsa_private_key {
public:
	using felem_t = bignum<N>;
	using public_key = spoint_t<N>;
	ecdsa_private_key() = default;
	template<typename curveT>
	ecdsa_private_key(const curveT& curve)
	{
		auto&	rd = bn_random<N>::Instance();
		felem_t secret;
		do {
			secret = rd.get_random();
		} while ( !init(curve, secret) );
		calc_pa(curve);
	}
	template<typename curveT>
	ecdsa_private_key(const curveT& curve, const felem_t& secret,
				const public_key& pk) : _pa(pk)
	{
		init(curve, secret);
	}
	explicit operator bool() const noexcept { return _inited; }
	const felem_t&		D() const noexcept { return *(&_d); }
#ifdef	WITH_MONT_D
	const felem_t&		mont_d() const noexcept { return *(&_mont_d); }
#endif
	const public_key&	PubKey() const noexcept { return *(&_pa); }
protected:
	template<typename curveT>
	bool init(const curveT& curve, const felem_t& secret) noexcept {
		_inited = false;
		if (secret >= curve.paramN()) _d.sub(secret, curve.paramN()); else _d = secret;
		//curve.modN(_d, secret);
		if (_d.is_zero()) return false;
#ifdef	WITH_MONT_D
		curve.to_montgomeryN(_mont_d, _d);
#endif
		_inited = true;
		return true;
	}
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
	public_key	_pa;
};


/*
 * ECDSA sign/verify
 * 			e = HASH (msg)
 * 			Pa ...	Public key, Px/Py
 * 			dA ...	private key
 * 	sign
 * 			k = random   256bits   Label getk
 * 			x1, y1 = k * G
 * 			r = x1
 *			if r == 0 modN  goto getk
 *			s = k^-1 * (e + r * dA)
 *			if s == 0 modN goto getk
 *			return r, s, y1_is_odd	// NOTE r, -s is also true sig
 *	// e + r*dA = s * k
 *	// r*dA  = s * k - e
 *	// dA = r^-1 * s * k - r^-1 * e
 *
 * 	verify
 * 			r, s all in [1, N-1], or fail
 * 			u = s^-1 * e modN
 * 			if u == 0 fail
 * 			v = s^-1 * r modN
 * 			if v == 0 fail
 * 			x2, y2 = u*G + v*Pa
 * 			if x2 == r modN success else fail
 *
 *	recover public key
 *			Psx = r
 *			Psy = get_Py ( Psx, Psy_is_odd )
 *			Ps ...	Psx/Psy
 *			u1 = s * r^-1
 *			u2 = e * r^-1
 *			Pa = r^1 * ( s*k*G - e*G) = s * r^-1 * Ps - e*r^-1 * G
 *			Pa = u1 * Ps - u2 * G
 */

// sign return ecInd,  0  for Py(of sign) even, 1 odd
template<const uint N=4, typename curveT>
forceinline static
int ecdsa_sign(const curveT& curve, bignum<N>& r, bignum<N>& s,
		const private_key<N>& priv, const bignum<N>& msg) noexcept
{
	bignum<N>	k, x1;	// y1 reuse tmp
	bignum<N>	tmp;	// r+s
	int		ret=0;
	do {
		do {
			gen_keypair<N>(curve, k, x1, tmp);
			if (x1 >= curve.paramN()) r.sub(x1, curve.paramN()); else r = x1;
		} while (r.is_zero());
		ret = tmp.is_odd();
		bignum<N>	rp;	// r+s
		curve.to_montgomeryN(rp, r);
		// tmp = r * dA
#ifdef	WITH_MONT_D
		curve.mont_nmult(tmp, rp, priv.mont_d());
#else
		curve.to_montgomeryN(tmp, priv.D());
		curve.mont_nmult(tmp, rp, tmp);
#endif
		curve.from_montgomeryN(tmp, tmp);
		// tmp = e + r * dA
		curve.mod_add_to(tmp, msg);
		// x1 = k^-1
		mod_inv<N>(x1, k, curve.paramN());
		// s = k^-1 * (e + r*dA)
		curve.to_montgomeryN(rp, x1);
		curve.to_montgomeryN(tmp, tmp);
		curve.mont_nmult(tmp, tmp, rp);
		curve.from_montgomeryN(s, tmp);
		if (s.is_zero()) continue;
		// SM2 sign s may above halfN
#ifdef	WITH_HALF_N_ommit
		if (curve.halfN() < s) {
			s.sub(curve.paramN(), s);
		}
#endif
	} while (s.is_zero());
	return ret;
}

// verify return bool, true for success(ok)
template<const uint N=4, typename curveT>
forceinline static
bool ecdsa_verify(const curveT& curve, const bignum<N>& r, const bignum<N>& s,
		const spoint_t<N>&	pub, const bignum<N>& msg) noexcept
{
	bignum<N>	u, v, sInv;
	if ( unlikely(r.is_zero()) ) return false;
	if ( unlikely(s.is_zero()) ) return false;
	mod_inv<N>(sInv, s, curve.paramN());
	curve.to_montgomery(sInv, sInv);
	// u = s^-1 * e
	// v = s^-1 * r
	curve.to_montgomery(u, msg);
	curve.to_montgomery(v, r);
	curve.mont_nmult(u, u, sInv);
	curve.mont_nmult(v, v, sInv);
	if ( unlikely(u.is_zero()) ) return false;
	if ( unlikely(v.is_zero()) ) return false;
	point_t<N>	q, p(pub.x, pub.y);
	// q = u*G + v * Pub
	// x2, y2 = u*G + v * Pub
	curve.combined_mult(q, p, v, u);
	if ( unlikely(q.x.is_zero()) ) return false;
	return r == q.x;
}

// vecover public key, from r/s/v/msg, v for pubY is odd
// return true for successfully recover public key
template<const uint N=4, typename curveT>
forceinline static
bool ecdsa_recover(const curveT& curve, spoint_t<N>&  pub, const bignum<N>& r,
				const bignum<N>& s, const int v, const bignum<N>& msg) noexcept
{
	point_t<N>	p;
	if ( unlikely(r.is_zero()) ) return false;
	if ( unlikely(s.is_zero()) ) return false;
	// p.x = r
	p.x = r;
	if ( unlikely(!pointY_recover(curve, p.y, p.x, v)) ) return false;
	bignum<N>	u1, u2, rInv;
	mod_inv<N>(rInv, r, curve.paramN());
	// u1 = s * r^-1
	// u2 = e * r^-1 = msg * r^-1
	curve.to_montgomeryN(rInv, rInv);
	curve.to_montgomeryN(u1, s);
	curve.to_montgomeryN(u2, msg);
	curve.mont_nmult(u1, u1, rInv);
	curve.mont_nmult(u2, u2, rInv);
	curve.from_montgomeryN(u2, u2);
	curve.from_montgomeryN(u1, u1);
	// 0 < u2 < N 
	if ( unlikely(u1.is_zero()) ) return false;
	if ( unlikely(u2.is_zero()) ) return false;
	// u2 = -u2 = - e * r^-1
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
