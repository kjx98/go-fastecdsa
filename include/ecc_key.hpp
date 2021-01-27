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
using public_key = spoint_t<4>;
class	private_key {
public:
	typedef bignum<4> felem_t;
	private_key() = default;
	template<typename curveT>
	bool init(const curveT& curve, const felem_t& secret) {
		curve.modN(_d, secret);
		if (_d.is_zero()) return false;
		_dInv = _d;
		this->_dInv.uadd_to(1);
		if (this->_dInv.is_zero()) return false;
		// _dInv = (1 + d)^-1
		mod_inv<4>(_dInv, _dInv, curve.paramN());
		// calc public key
		return true;
	}
	const bignum<4>&	D() const noexcept { return *(&_d); }
	const bignum<4>&	Di() const noexcept { return *(&_dInv); }
	const public_key&	PA() const noexcept { return *(&_pa); }
protected:
	felem_t		_d;	// secret
	felem_t		_dInv;	// (1+d)^-1  for SM2 sign
	public_key	_pa;
};


}	// namespace ecc

#endif	//	__ECC_KEY_H__
