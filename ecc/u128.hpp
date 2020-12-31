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
#if defined(__SIZEOF_INT128__)
	uint128_t(const __uint128_t vd=0) : _data(vd) { };
#else
	uint128_t() : _low(0), _high(0) {};
#endif
	uint128_t(const u64 vl, const u64 vh) :
#if defined(__SIZEOF_INT128__)
		_data( (((__uint128_t)vh) << 64) | vl)
#else
		_low( vl ),
		_high( vh )
#endif
	{
	}
	uint128_t(const uint128_t &) = default;
#if defined(__SIZEOF_INT128__)
	const __uint128_t data() const { return _data; }
#endif
	u64 m_low() const {
#if defined(__SIZEOF_INT128__)
		return (u64)_data;
#else
		return _low;
#endif
	}
	u64 m_high() const {
#if defined(__SIZEOF_INT128__)
		return (u64)(_data >> 64);
#else
		return _high;
#endif
	}
	friend uint128_t operator+(uint128_t a, const uint128_t &b) noexcept
	{
#if defined(__SIZEOF_INT128__)
		a._data += b._data;
#else
		a._high += b._high;
		if (__builtin_uaddl_overflow(a._low, b._low, &a._low)) ++a._high;
#endif
		return a;
	}
	uint128_t& operator+=(const uint128_t& b) noexcept
	{
#if defined(__SIZEOF_INT128__)
		_data += b._data;
#else
		if (__builtin_uaddl_overflow(_low, b._low, &_low)) _high++;
		_high += b._high;
#endif
		return *this;
	}
	uint128_t& mul_64_64(u64 left, u64 right) noexcept
	{
#if defined(__SIZEOF_INT128__)
		_data = (__uint128_t)left * right;
#else
		u64 a0 = left & 0xffffffffull;
		u64 a1 = left >> 32;
		u64 b0 = right & 0xffffffffull;
		u64 b1 = right >> 32;
		u64 m0 = a0 * b0;
		u64 m1 = a0 * b1;
		u64 m2 = a1 * b0;
		u64 m3 = a1 * b1;

		m2 += (m0 >> 32);
		m2 += m1;

		/* Overflow */
		if (m2 < m1)
			m3 += 0x100000000ull;

		_low = (m0 & 0xffffffffull) | (m2 << 32);
		_high = m3 + (m2 >> 32);
#endif
		return *this;
	}
private:
#if defined(__SIZEOF_INT128__)
	__uint128_t	_data;
#else
	u64	_low;
	u64	_high;
#endif
};

#endif	//	__U128_HPP__
