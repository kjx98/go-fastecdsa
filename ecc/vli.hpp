// +build ignore

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
#ifndef __VLI_HPP__
#define __VLI_HPP__

#include <stdint.h>
#include <stdbool.h>
#include <endian.h>
#include <cassert>
#include <type_traits>
#include "cdefs.h"

#if	__cplusplus < 201103L
# error "C++ std MUST at least c++11"
#endif


template<typename T> forceinline
static T max(const T&& a, const T&& b) noexcept {
	static_assert(std::is_integral<T>::value, "Integral required.");
	return (a>b)?a:b;
}


//__attribute__((optimize("unroll-loops")))
template<uint ndigits> forceinline
static void vli_clear(u64 *vli) noexcept
{
	uint i;
#pragma GCC unroll 4
	for (i = 0; i < ndigits; i++)
		vli[i] = 0;
}

/**
 * vli_is_zero() - Determine is vli is zero
 *
 * @vli:		vli to check.
 * @ndigits:		length of the @vli
 */
/* Returns true if vli == 0, false otherwise. */
template<uint ndigits> forceinline
static bool vli_is_zero(const u64 *vli) noexcept
{
#pragma GCC unroll 4
	for (uint i = 0; i < ndigits; i++) {
		if (vli[i]) return false;
	}
	return true;
}

template<uint ndigits> forceinline
static bool vli_is_one(const u64 *vli) noexcept
{
	if (vli[0] != 1) return false;
#pragma GCC unroll 4
	for (uint i = 1; i < ndigits; i++) {
		if (vli[i]) return false;
	}
	return true;
}

/* Sets dest = src. */
/* Sets dest = src. */
template<uint ndigits> forceinline
static void vli_set(u64 *dest, const u64 *src) noexcept
{
#pragma GCC unroll 4
	for (uint i = 0; i < ndigits; i++)
		dest[i] = src[i];
}

/* Returns nonzero if bit bit of vli is set. */
static forceinline bool vli_test_bit(const u64 *vli, uint bit)
{
	return (vli[bit / 64] & ((u64)1 << (bit % 64)));
}

/**
 * vli_cmp() - compare left and right vlis
 *
 * @left:		vli
 * @right:		vli
 * @ndigits:		length of both vlis
 *
 * Returns sign of @left - @right, i.e. -1 if @left < @right,
 * 0 if @left == @right, 1 if @left > @right.
 */
template<uint ndigits> forceinline
static int vli_cmp(const u64 *left, const u64 *right) noexcept
{
#pragma GCC unroll 4
	for (int i = ndigits - 1; i >= 0; i--) {
		if (left[i] > right[i]) return 1;
		else if (left[i] < right[i]) return -1;
	}
	return 0;
}

#ifdef	ommit
/* Computes result = in << c, returning carry. Can modify in place
 * (if result == in). 0 < shift < 64.
 */
template<uint ndigits> forceinline static
u64 vli_lshift(u64 *result, const u64 *in, unsigned int shift) noexcept
{
	u64 carry = 0;
#pragma GCC unroll 4
	for (uint i = 0; i < ndigits; i++) {
		u64 temp = in[i];
		result[i] = (temp << shift) | carry;
		carry = temp >> (64 - shift);
	}
	return carry;
}
#endif

template<uint ndigits> forceinline static
u64 vli_lshift1(u64 *result, const u64 *in) noexcept
{
	u64 carry = 0;
#pragma GCC unroll 4
	for (uint i = 0; i < ndigits; i++) {
		u64 temp = in[i];
		result[i] = (temp << 1) | carry;
		carry = temp >> 63;
	}
	return carry;
}

/* Computes vli = vli >> 1. */
template<uint ndigits> forceinline
static void vli_rshift1(u64 *vli) noexcept
{
	u64 carry = 0;
#pragma GCC unroll 4
	for (int i=ndigits - 1; i >= 0; i--) { 
		u64 temp = vli[i];
		vli[i] = (temp >> 1) | carry;
		carry = temp << 63;
	}
}

/* Computes vli = vli >> 1. */
template<uint ndigits> forceinline
static void vli_rshift1w(u64 *vli) noexcept
{
#pragma GCC unroll 4
	for (uint i = 1; i < ndigits; i++) {
		vli[i-1] = vli[i];
	}
	vli[ndigits-1] = 0;
}

/* Computes result = left + right, returning carry. Can modify in place. */
template<uint ndigits> forceinline static
bool vli_add(u64 *result, const u64 *left, const u64 *right) noexcept
{
	bool carry = false;
#pragma GCC unroll 4
	for (uint i = 0; i < ndigits; i++) {
		u64 sum;
#ifdef	NO_BUILTIN_OVERFLOW
		sum = left[i] + right[i] + carry;
		if (sum != left[i])
			carry = (sum < left[i]);
#else
		auto c_carry = __builtin_uaddl_overflow(left[i], right[i], &sum);
		if (unlikely(carry)) {
			carry = c_carry | __builtin_uaddl_overflow(sum, 1, &sum);
		} else carry = c_carry;
#endif
		result[i] = sum;
	}
	return carry;
}

/* Computes result = left + right, returning carry. Can modify in place. */
template<uint ndigits> forceinline
static bool vli_uadd(u64 *result, const u64 *left, u64 right) noexcept
{
#ifdef	NO_BUILTIN_OVERFLOW
	u64 carry = right;
	uint i;
#pragma GCC unroll 4
	for (i = 0; i < ndigits; i++) {
		u64 sum;
		sum = left[i] + carry;
		if (sum != left[i])
			carry = (sum < left[i]);
		else
			carry = !!carry;
		result[i] = sum;
	}
	return carry != 0;
#else
	auto carry = __builtin_uaddl_overflow(left[0], right, result);
#pragma GCC unroll 4
	for (uint i = 1; i < ndigits; i++) {
		if (unlikely(carry)) {
			carry = __builtin_uaddl_overflow(result[i], 1, result+i);
		} else break;
	}
	return carry;
#endif
}

/* Computes result = left + right, returning carry. Can modify in place. */
template<uint ndigits> forceinline static
bool vli_add_to(u64 *result, const u64 *right) noexcept
{
	bool carry = false;
#pragma GCC unroll 4
	for (uint i = 0; i < ndigits; i++) {
		u64 sum;
#ifdef	NO_BUILTIN_OVERFLOW
		sum = result[i] + right[i] + carry;
		if (sum != result[i])
			carry = (sum < result[i]);
#else
		auto c_carry = __builtin_uaddl_overflow(result[i], right[i], &sum);
		if (unlikely(carry)) {
			carry = c_carry | __builtin_uaddl_overflow(sum, 1, &sum);
		} else carry = c_carry;
#endif
		result[i] = sum;
	}
	return carry;
}

/* Computes result = left + right, returning carry. Can modify in place. */
template<uint ndigits> forceinline
static bool vli_uadd_to(u64 *result, u64 right) noexcept
{
#ifdef	NO_BUILTIN_OVERFLOW
	u64 carry = right;
	uint i;
#pragma GCC unroll 4
	for (i = 0; i < ndigits; i++) {
		u64 sum;
		sum = result[i] + carry;
		if (sum != result[i])
			carry = (sum < result[i]);
		else
			carry = !!carry;
		result[i] = sum;
	}
	return carry != 0;
#else
	auto carry = __builtin_uaddl_overflow(result[0], right, result);
#pragma GCC unroll 4
	for (uint i = 1; i < ndigits; i++) {
		if (unlikely(carry)) {
			carry = __builtin_uaddl_overflow(result[i], 1, result+i);
		} else break;
	}
	return carry;
#endif
}

/**
 * vli_sub() - Subtracts right from left
 *
 * @result:		where to write result
 * @left:		vli
 * @right		vli
 * @ndigits:		length of all vlis
 *
 * Note: can modify in-place.
 *
 * Return: carry bit.
 */
template<uint ndigits> forceinline static
bool vli_sub(u64 *result, const u64 *left, const u64 *right) noexcept
{
	bool borrow = false;
#pragma GCC unroll 4
	for (uint i = 0; i < ndigits; i++) {
		u64 diff;
#ifdef	NO_BUILTIN_OVERFLOW
		diff = left[i] - right[i] - borrow;
		if (diff != left[i]) borrow = (diff > left[i]);
#else
		auto c_borrow = __builtin_usubl_overflow(left[i], right[i], &diff);
		if (unlikely(borrow)) {
			borrow = c_borrow | __builtin_usubl_overflow(diff, 1, &diff);
		} else borrow = c_borrow;
#endif
		result[i] = diff;
	}
	return borrow;
}

/* Computes result = left - right, returning borrow. Can modify in place. */
template<uint ndigits> forceinline
static bool vli_usub(u64 *result, const u64 *left, u64 right) noexcept
{
#ifdef	NO_BUILTIN_OVERFLOW
	u64 borrow = right;
	uint i;
#pragma GCC unroll 4
	for (i = 0; i < ndigits; i++) {
		u64 diff;
		diff = left[i] - borrow;
		if (diff != left[i])
			borrow = (diff > left[i]);
		result[i] = diff;
	}
	return borrow != 0;
#else
	auto borrow = __builtin_usubl_overflow(left[0], right, result);
#pragma GCC unroll 4
	for (uint i = 1; i < ndigits; i++) {
		if (unlikely(borrow)) {
			borrow = __builtin_usubl_overflow(result[i], 1, result+i);
		} else break;
	}
	return borrow;
#endif
}

template<uint ndigits> forceinline static
bool vli_sub_from(u64 *result, const u64 *right) noexcept
{
	bool borrow = false;
#pragma GCC unroll 4
	for (uint i = 0; i < ndigits; i++) {
		u64 diff;
#ifdef	NO_BUILTIN_OVERFLOW
		diff = result[i] - right[i] - borrow;
		if (diff != result[i]) borrow = (diff > result[i]);
#else
		auto c_borrow = __builtin_usubl_overflow(result[i], right[i], &diff);
		if (unlikely(borrow)) {
			borrow = c_borrow | __builtin_usubl_overflow(diff, 1, &diff);
		} else borrow = c_borrow;
#endif
		result[i] = diff;
	}
	return borrow;
}

/* Computes result = left - right, returning borrow. Can modify in place. */
template<uint ndigits> forceinline
static bool vli_usub_from(u64 *result, const u64 *left, u64 right) noexcept
{
#ifdef	NO_BUILTIN_OVERFLOW
	u64 borrow = right;
	uint i;
#pragma GCC unroll 4
	for (i = 0; i < ndigits; i++) {
		u64 diff;
		diff = result[i] - borrow;
		if (diff != result[i])
			borrow = (diff > result[i]);
		result[i] = diff;
	}
	return borrow != 0;
#else
	auto borrow = __builtin_usubl_overflow(result[0], right, result);
#pragma GCC unroll 4
	for (uint i = 1; i < ndigits; i++) {
		if (unlikely(borrow)) {
			borrow = __builtin_usubl_overflow(result[i], 1, result+i);
		} else break;
	}
	return borrow;
#endif
}

template<uint ndigits> forceinline
static bool vli_is_negative(const u64 *vli)  noexcept
{
	return vli_test_bit(vli, ndigits * 64 - 1);
}

forceinline static bool vli_is_even(const u64 *vli) noexcept {
	return (vli[0] & 1) == 0;
}

/* Counts the number of 64-bit "digits" in vli. */
template<uint ndigits> forceinline
static uint vli_num_digits(const u64 *vli) noexcept
{
	/* Search from the end until we find a non-zero digit.
	 * We do it in reverse because we expect that most digits will
	 * be nonzero.
	 */
	int i;
	for (i = ndigits - 1; i >= 0 && vli[i] == 0; i--);
	return (i + 1);
}

/* Counts the number of bits required for vli. */
template<uint ndigits> forceinline
static uint vli_num_bits(const u64 *vli) noexcept
{
	auto num_digits = vli_num_digits<ndigits>(vli);
	if (num_digits == 0) return 0;
	auto i = 64 - __builtin_clzl(vli[num_digits - 1]);
	return ((num_digits - 1) * 64 + i);
}

#define ECC_DIGITS_TO_BYTES_SHIFT 3
forceinline static uint vli_bytes(uint ndigits) {
	return ndigits << ECC_DIGITS_TO_BYTES_SHIFT;
}

/**
 * vli_from_be64() - Load vli from big-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 BE values
 * @ndigits:		length of both vli and array
 */
template<uint ndigits> forceinline
static void vli_from_be64(u64 *dest, const void *src) noexcept
{
	const u64 *from = (const u64 *)src;
#pragma GCC unroll 4
	for (uint i = 0; i < ndigits; i++)
		dest[i] = be64toh(from[ndigits - 1 - i]);
}

#ifdef	ommit
/**
 * vli_from_le64() - Load vli from little-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 LE values
 * @ndigits:		length of both vli and array
 */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static void vli_from_le64(u64 *dest, const void *src) noexcept
{
	const u64 *from = (const u64 *)src;
	for (uint i = 0; i < ndigits; i++)
		dest[i] = le64toh(from[ndigits - 1 - i]);
}
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

template<uint ndigits> forceinline static
void vli_mult(u64 *result, const u64 *left, const u64 *right) noexcept
{
	uint128_t r01( 0, 0 );
	u64 r2 = 0;
	unsigned int i, k;
	/* Compute each digit of result in sequence, maintaining the
	 * carries.
	 */
#pragma GCC unroll 8
	for (k = 0; k < ndigits * 2 - 1; k++) {
		unsigned int min;
		if (k < ndigits)
			min = 0;
		else
			min = (k + 1) - ndigits;

		for (i = min; i <= k && i < ndigits; i++) {
			uint128_t product;

			//product = mul_64_64(left[i], right[k - i]);
			product.mul_64_64(left[i], right[k - i]);

			//r01 = add_128_128(r01, product);
			r01 += product;
			r2 += (r01.m_high() < product.m_high());
		}

		result[k] = r01.m_low();
		r01 = uint128_t(r01.m_high(), r2);
		//r01.m_low = r01.m_high;
		//r01.m_high = r2;
		r2 = 0;
	}

	result[ndigits * 2 - 1] = r01.m_low();
}

/* Compute product = left * right, for a small right value. */
template<uint ndigits> forceinline
static void vli_umult(u64 *result, const u64 *left, u64 right) noexcept
{
	uint128_t r01( 0, 0 );
	unsigned int k;

#pragma GCC unroll 4
	for (k = 0; k < ndigits; k++) {
		uint128_t product;

		//if (likely(left[k] != 0))
		{
			//product = mul_64_64(left[k], right);
			product.mul_64_64(left[k], right);
			//r01 = add_128_128(r01, product);
			r01 += product;
		}
		/* no carry */
		result[k] = r01.m_low();
		r01 = uint128_t(r01.m_high(), 0);
		//r01.m_low = r01.m_high;
		//r01.m_high = 0;
	}
	result[ndigits] = r01.m_low();
#pragma GCC unroll 4
	for (k = ndigits+1; k < ndigits * 2; k++)
		result[k] = 0;
}

template<uint ndigits> forceinline
static void vli_square(u64 *result, const u64 *left) noexcept
{
	uint128_t r01( 0, 0 );
	u64 r2 = 0;
	uint i, k;

#pragma GCC unroll 8
	for (k = 0; k < ndigits * 2 - 1; k++) {
		unsigned int min;

		if (k < ndigits)
			min = 0;
		else
			min = (k + 1) - ndigits;

		for (i = min; i <= k && i <= k - i; i++) {
			uint128_t product;

			//product = mul_64_64(left[i], left[k - i]);
			product.mul_64_64(left[i], left[k - i]);

			if (i < k - i) {
				r2 += product.m_high() >> 63;
				u64 _high = (product.m_high() << 1) | (product.m_low() >> 63);
				u64 _low = product.m_low() << 1;
				product = uint128_t(_low, _high);
			}

			//r01 = add_128_128(r01, product);
			r01 += product;
			r2 += (r01.m_high() < product.m_high());
		}

		result[k] = r01.m_low();
		//r01.m_low = r01.m_high;
		//r01.m_high = r2;
		r01 = uint128_t(r01.m_high(), r2);
		r2 = 0;
	}

	result[ndigits * 2 - 1] = r01.m_low();
}

/* Computes result = (left + right) % mod.
 * Assumes that left < mod and right < mod, result != mod.
 */
template<uint ndigits> forceinline
static void vli_mod_add(u64 *result, const u64 *left, const u64 *right,
			const u64 *mod) noexcept
{
	auto carry = vli_add<ndigits>(result, left, right);

	/* result > mod (result = mod + remainder), so subtract mod to
	 * get remainder.
	 */
	if (carry || vli_cmp<ndigits>(result, mod) >= 0)
		vli_sub_from<ndigits>(result, mod);
}

/* Computes result = (left - right) % mod.
 * Assumes that left < mod and right < mod, result != mod.
 */
template<uint ndigits> forceinline
static void vli_mod_sub(u64 *result, const u64 *left, const u64 *right,
			const u64 *mod) noexcept
{
	auto borrow = vli_sub<ndigits>(result, left, right);

	/* In this case, p_result == -diff == (max int) - diff.
	 * Since -x % d == d - x, we can get the correct result from
	 * result + mod (with overflow).
	 */
	if (borrow)
		vli_add_to<ndigits>(result, mod);
}

#endif	//	__VLI_HPP__
