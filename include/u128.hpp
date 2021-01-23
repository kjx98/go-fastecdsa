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
#ifndef __U128_HPP__
#define __U128_HPP__

#include <stdint.h>
#include <stdbool.h>
#include <endian.h>
#include <cassert>
#include <type_traits>
#include "cdefs.h"

#if	__cplusplus < 201103L
# error "C++ std MUST at least c++11"
#endif


class alignas(16) uint128_t {
public:
	uint128_t(const __uint128_t vd=0) : _data(vd) { };
	uint128_t(const u64 vl, const u64 vh) :
		_data( (((__uint128_t)vh) << 64) | vl)
	{
	}
	uint128_t(const uint128_t &) = default;
	const __uint128_t data() const { return _data; }
	u64 m_low() const {
		return (u64)_data;
	}
	u64 m_high() const {
		return (u64)(_data >> 64);
	}
	friend uint128_t operator+(uint128_t a, const uint128_t &b) noexcept
	{
		a._data += b._data;
		return a;
	}
	uint128_t& operator+=(const uint128_t& b) noexcept
	{
		_data += b._data;
		return *this;
	}
	uint128_t& mul_64_64(u64 left, u64 right) noexcept
	{
		_data = (__uint128_t)left * right;
		return *this;
	}
private:
	__uint128_t	_data;
};

#endif	//	__U128_HPP__
