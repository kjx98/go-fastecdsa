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
#ifndef __VLI_BN_HPP__
#define __VLI_BN_HPP__

#include <stdint.h>
#include <stdbool.h>
#include <endian.h>
#include <cassert>
#include <type_traits>
#include "cdefs.h"

#if	__cplusplus < 201103L
# error "C++ std MUST at least c++11"
#endif


namespace vli {

template<uint ndigits>
class bignum {
public:
	using bn_words = u64[ndigits];
	bignum(bn_words init = {}) noexcept : d{} {}
	bignum(const bignum &) = default;
	explicit bignum(const u64 v) noexcept : d{v,0,0,0} {}
	bignum(const u64 *src) {
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++)
			d[i] = src[i];
	}
	void clear() noexcept
	{
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++)
			d[i] = 0;
	}
	const u64* data() const { return d; }
	bool is_zero() noexcept
	{
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			if (d[i]) return false;
		}
		return true;
	}
	bool is_one() noexcept
	{
		if (d[0] != 1) return false;
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (d[i]) return false;
		}
		return true;
	}
	/* Sets dest = this, copyout */
	void set(u64 *dest) noexcept
	{
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++)
			dest[i] = d[i];
	}
	/* Returns nonzero if bit bit of vli is set. */
	bool test_bit(const uint bit)
	{
		if ( bit >= ndigits*64 ) return false;	// out of bound
		return (d[bit / 64] & ((u64)1 << (bit % 64)));
	}
/**
 * cmp() - compare this and right vlis
 *
 * @left:		vli
 * @right:		vli
 * @ndigits:		length of both vlis
 *
 * Returns sign of @left - @right, i.e. -1 if @left < @right,
 * 0 if @left == @right, 1 if @left > @right.
 */
	int cmp(const bignum& right) noexcept
	{
#pragma GCC unroll 4
		for (int i = ndigits - 1; i >= 0; i--) {
			if (d[i] > right.d[i]) return 1;
			else if (d[i] < right.d[i]) return -1;
		}
		return 0;
	}
	u64 lshift1(bignum& result) noexcept
	{
		u64 carry = 0;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			u64 temp = d[i];
			result.d[i] = (temp << 1) | carry;
			carry = temp >> 63;
		}
		return carry;
	}
/* Computes vli = vli >> 1. */
	void rshift1() noexcept
	{
		u64 carry = 0;
#pragma GCC unroll 4
		for (int i=ndigits - 1; i >= 0; i--) { 
			u64 temp = d[i];
			d[i] = (temp >> 1) | carry;
			carry = temp << 63;
		}
	}	
/* Computes vli = vli >> 1 word (64 Bits). */
	void rshift1w() noexcept
	{
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			d[i-1] = d[i];
		}
		d[ndigits-1] = 0;
	}
/* Computes result = this + right, returning carry. Can modify in place. */
	bool add(const bignum& left, const bignum &right) noexcept
	{
		bool carry = false;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			u64 sum;
			auto c_carry = __builtin_uaddl_overflow(left.d[i], right.d[i], &sum);
			if (unlikely(carry)) {
				carry = c_carry | __builtin_uaddl_overflow(sum, 1, &sum);
			} else carry = c_carry;
			d[i] = sum;
		}
		return carry;
	}
/* Computes this = left + right, returning carry. Can modify in place. */
	bool uadd(const bignum& left, const u64 right) noexcept
	{
		auto carry = __builtin_uaddl_overflow(left.d[0], right, this->d);
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (unlikely(carry)) {
				carry = __builtin_uaddl_overflow(left.d[i], 1, d+i);
			} else d[i] = left.d[i];
		}
		return carry;
	}
/* Computes this = this + right, returning carry. Can modify in place. */
	bool add_to(const bignum& right) noexcept
	{
		bool carry = false;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			u64 sum;
			auto c_carry = __builtin_uaddl_overflow(d[i], right.d[i], &sum);
			if (unlikely(carry)) {
				carry = c_carry | __builtin_uaddl_overflow(sum, 1, &sum);
			} else carry = c_carry;
			d[i] = sum;
		}
		return carry;
	}
/* Computes this = this + right, returning carry. Can modify in place. */
	bool uadd_to(const u64 right) noexcept
	{
		auto carry = __builtin_uaddl_overflow(d[0], right, d);
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (unlikely(carry)) {
				carry = __builtin_uaddl_overflow(d[i], 1, d+i);
			} else break;
		}
		return carry;
	}
/**
 * sub() - Subtracts right from this
 *
 * @result:		where to write result(this)
 * @left:		vli
 * @right		vli
 * @ndigits:		length of all vlis
 *
 * Note: can modify in-place.
 *
 * Return: carry bit.
 */
	bool sub(const bignum& left, const bignum& right) noexcept
	{
		bool borrow = false;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			u64 diff;
			auto c_borrow = __builtin_usubl_overflow(left.d[i], right.d[i], &diff);
			if (unlikely(borrow)) {
				borrow = c_borrow | __builtin_usubl_overflow(diff, 1, &diff);
			} else borrow = c_borrow;
			d[i] = diff;
		}
		return borrow;
	}
/* Computes this = left - right, returning borrow. Can modify in place. */
	bool usub(const bignum& left, const u64 right) noexcept
	{
		auto borrow = __builtin_usubl_overflow(left.d[0], right, this->d);
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (unlikely(borrow)) {
				borrow = __builtin_usubl_overflow(left.d[i], 1, d+i);
			} else d[i] = left.d[i];
		}
		return borrow;
	}
	bool sub_from(const bignum& right) noexcept
	{
		bool borrow = false;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			u64 diff;
			auto c_borrow = __builtin_usubl_overflow(d[i], right.d[i], &diff);
			if (unlikely(borrow)) {
				borrow = c_borrow | __builtin_usubl_overflow(diff, 1, &diff);
			} else borrow = c_borrow;
			d[i] = diff;
		}
		return borrow;
	}
/* Computes this = this` - right, returning borrow. Can modify in place. */
	bool usub_from(const u64 right) noexcept
	{
		auto borrow = __builtin_usubl_overflow(d[0], right, d);
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (unlikely(borrow)) {
				borrow = __builtin_usubl_overflow(d[i], 1, d+i);
			} else break;
		}
		return borrow;
	}
	bool is_negative()  noexcept
	{
		return (d[ndigits-1] & ((u64)1 << 63)) != 0;
	}
	bool is_even() noexcept {
		return (d[0] & 1) == 0;
	}
	bool is_odd() noexcept {
		return (d[0] & 1) != 0;
	}
/* Counts the number of 64-bit "digits" in vli. */
	uint num_digits() noexcept
	{
		/* Search from the end until we find a non-zero digit.
		 * We do it in reverse because we expect that most digits will
		 * be nonzero.
		 */
		int i;
		for (i = ndigits - 1; i >= 0 && d[i] == 0; i--);
		return (i + 1);
	}
	uint vli_num_bits() noexcept
	{
		auto num_digits = num_digits();
		if (num_digits == 0) return 0;
		auto i = 64 - __builtin_clzl(d[num_digits - 1]);
		return ((num_digits - 1) * 64 + i);
	}
/**
 * vli_from_be64() - Load vli from big-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 BE values
 * @ndigits:		length of both vli and array
 */
	void from_be64(const void *src) noexcept
	{
		const u64 *from = (const u64 *)src;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++)
			d[i] = be64toh(from[ndigits - 1 - i]);
	}
/* Computes this = (left + right) % mod.
 * Assumes that left < mod and right < mod, result != mod.
 */
	bignum& mod_add(bignum& left, const bignum &right, const bignum &mod) noexcept
	{
		auto carry = this->add(left, right);

	/* result > mod (result = mod + remainder), so subtract mod to
	 * get remainder.
	 */
		if (carry || this->cmp(mod) >= 0)
			this->sub_from(mod);
		return *this;
	}
/* Computes this = (left - right) % mod.
 * Assumes that left < mod and right < mod, result != mod.
 */
	bignum& mod_sub(const bignum& left, const bignum& right, const bignum &mod) noexcept
	{
		auto borrow = this->sub(left, right);

	/* In this case, p_result == -diff == (max int) - diff.
	 * Since -x % d == d - x, we can get the correct result from
	 * result + mod (with overflow).
	 */
		if (borrow)
			this->add_to(mod);
	}

protected:
	bn_words	d;
};


template<uint ndigits>
class bn_prod: public bignum<ndigits*2> {
public:
	// this = left * right
	void mult(const bignum<ndigits>& left, const bignum<ndigits>& right) noexcept
	{
		uint128_t r01( 0, 0 );
	/* Compute each digit of result in sequence, maintaining the
	 * carries.
	 */
#pragma GCC unroll 8
		for (uint k = 0; k < ndigits * 2 - 1; k++) {
			u64 r2 = 0;
			unsigned int min;
			if (k < ndigits) min = 0; else min = (k + 1) - ndigits;
			for (uint i = min; i <= k && i < ndigits; i++) {
				uint128_t product;
				//product = mul_64_64(left[i], right[k - i]);
				product.mul_64_64(left.d[i], right.d[k - i]);
				//r01 = add_128_128(r01, product);
				r01 += product;
				r2 += (r01.m_high() < product.m_high());
			}
			this->d[k] = r01.m_low();
			r01 = uint128_t(r01.m_high(), r2);
			//r01.m_low = r01.m_high;
			//r01.m_high = r2;
			r2 = 0;
		}
		this->d[ndigits * 2 - 1] = r01.m_low();
	}

	/* Compute product = left * right, for a small right value. */
	void umult(const bignum<ndigits>&left, const u64 right) noexcept
	{
		uint128_t r01( 0, 0 );
		unsigned int k;
#pragma GCC unroll 4
		for (k = 0; k < ndigits; k++) {
			uint128_t product;
			//if (likely(left[k] != 0))
			{
				//product = mul_64_64(left[k], right);
				product.mul_64_64(left.d[k], right);
				//r01 = add_128_128(r01, product);
				r01 += product;
			}
			/* no carry */
			this->d[k] = r01.m_low();
			r01 = uint128_t(r01.m_high(), 0);
			//r01.m_low = r01.m_high;
			//r01.m_high = 0;
		}
		this->d[ndigits] = r01.m_low();
#pragma GCC unroll 4
		for (k = ndigits+1; k < ndigits * 2; k++) this->d[k] = 0;
	}

	void square(const bignum<ndigits>& left) noexcept
	{
		uint128_t r01( 0, 0 );
#pragma GCC unroll 8
		for (uint k = 0; k < ndigits * 2 - 1; k++) {
			unsigned int min;
			u64 r2 = 0;
			if (k < ndigits) min = 0; else min = (k + 1) - ndigits;
			for (uint i = min; i <= k && i <= k - i; i++) {
				uint128_t product;
				//product = mul_64_64(left[i], left[k - i]);
				product.mul_64_64(left.d[i], left.d[k - i]);
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
			this->d[k] = r01.m_low();
			//r01.m_low = r01.m_high;
			//r01.m_high = r2;
			r01 = uint128_t(r01.m_high(), r2);
		}
		this->d[ndigits * 2 - 1] = r01.m_low();
	}
};

}
#endif	//	__VLI_BN_HPP__
