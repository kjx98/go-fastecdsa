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
#include "ecc_impl.hpp"


namespace ecc {

using namespace vli;


// definitions for Private/Public key
template<const uint N=4>
class	private_key {
public:
	using felem_t = bignum<N>;
	using public_key = spoint_t<N>;
	private_key() = default;
	template<typename curveT>
	private_key(const curveT& curve, const felem_t& secret,
				const public_key& pk) : _pa(pk)
	{
		init(curve, secret);
	}
	template<typename curveT>
	bool init(const curveT& curve, const felem_t& secret) noexcept {
		curve.modN(_d, secret);
		if (_d.is_zero()) return false;
		_dInv = _d;
		_dInv.uadd_to(1);
		curve.modN(_dInv, _dInv);
		if (_dInv.is_zero()) return false;
		// _dInv = (1 + d)^-1
		mod_inv<N>(_dInv, _dInv, curve.paramN());
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
	const felem_t&		D() const noexcept { return *(&_d); }
	const felem_t&		Di() const noexcept { return *(&_dInv); }
	const public_key&	PA() const noexcept { return *(&_pa); }
protected:
	bool		_inited = false;
	felem_t		_d;	// secret
	felem_t		_dInv;	// (1+d)^-1  for SM2 sign
	public_key	_pa;
};


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
 *	recovery public key
 *			Psx = r - e
 *			Psy = get_Py ( Psx, Psy_is_odd )
 *			Ps ...	Psx/Psy
 *			u1 = (r + s)^-1
 *			u2 = s * u1
 *			Pa = (r+s)^1 * ( k*G - s*G) = (r+s)^-1 * Ps - s*(r+s)^-1 * G
 *			Pa = u1 * Ps - u2 * G
 */

}	// namespace ecc

#endif	//	__ECC_KEY_H__
