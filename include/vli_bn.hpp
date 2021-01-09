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
#if	!defined(__VLI_BN_HPP__)
#define __VLI_BN_HPP__

#include <stdint.h>
#include <stdbool.h>
#include <endian.h>
#include <cassert>
#include <type_traits>
#include <iostream>
//#include <iomanip>
#include "cdefs.h"
#include "vli.hpp"

#if	__cplusplus < 201103L
# error "C++ std MUST at least c++11"
#endif

#include "u128.hpp"		// only for uint128_t

namespace vli {

template<const uint ndigits>
class bignum {
public:
	using bn_words = u64[ndigits];
	bignum() = default;
	bignum(const bignum &) = default;
	explicit bignum(const u64 v) noexcept : d{v,0,0,0}
	{
#if	__cplusplus >= 201703L
		if constexpr(ndigits > 4) {
#pragma GCC unroll 4
			for (uint i = 4; i < ndigits; i++)
				d[i] = 0;
		}
#endif
	}
	explicit bignum(const u64 *src) {
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++)
			d[i] = src[i];
	}
	bignum<4>& bn256() noexcept {
		static_assert(ndigits >= 4, "ndigits >= 4 required.");
		bignum<4>	*res=reinterpret_cast<bignum<4> *>(this->d);
		return *res;
	}
	void clear() noexcept
	{
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++)
			d[i] = 0;
	}
	friend std::ostream& operator<<(std::ostream& os, const bignum& x)
	{
		os << std::hex;
		if (x.is_u64()) {
			os << x.d[0];
		} else {
			for (uint i=0; i < ndigits; i++) {
				os << " " << x.d[i];
			}
		}
		os << std::dec;
		return os;
	}
	bool operator<(const bignum& bn) {
#pragma GCC unroll 4
		for (int i = ndigits - 1; i >= 0; i--) {
			if (d[i] > bn.d[i]) return false;
			else if (d[i] < bn.d[i]) return true;
		}
		return false;
	}
	bool operator>=(const bignum& bn) {
		return !(*this < bn);
	}
	bool operator==(const bignum& bn) {
#pragma GCC unroll 4
		for (int i = ndigits - 1; i >= 0; i--) {
			if (d[i] > bn.d[i]) return false;
			else if (d[i] < bn.d[i]) return false;
		}
		return true;
	}
	const u64* data() const noexcept { return d; }
	//u64* raw_data() noexcept { return d; }
	bool is_zero() const noexcept
	{
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			if (d[i]) return false;
		}
		return true;
	}
	bool is_one() const noexcept
	{
		if (d[0] != 1) return false;
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (d[i]) return false;
		}
		return true;
	}
	bool is_u64() const noexcept
	{
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (d[i]) return false;
		}
		return true;
	}
	/* Sets dest = this, copyout */
	void set(u64 *dest) const noexcept
	{
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++)
			dest[i] = d[i];
	}
	/* Returns nonzero if bit bit of vli is set. */
	bool test_bit(const uint bit) const noexcept
	{
		if ( bit >= ndigits*64 ) return false;	// out of bound
		return (d[bit / 64] & ((u64)1 << (bit % 64)));
	}
	u8 get_bit(const uint bit) const noexcept
	{
		if ( bit >= ndigits*64 ) return 0;	// out of bound
		return (d[bit >> 6] >> (bit & 0x3f)) & 1;
	}
	u8 get_bits(const uint bit, const uint cnt=5) const noexcept
	{
		auto off = bit >> 6;
		auto rem = bit & 0x3f;
		if (off >= ndigits) return 0;
		u8	rr = (d[off] >> rem) & 0xff;
		off++;
		if (off < ndigits && rem > (64 - cnt)) {
			rr |= (d[off] << rem) & 0xff;
		}
		return rr & ((1<<cnt) - 1);
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
	int cmp(const bignum& right) const noexcept
	{
#pragma GCC unroll 4
		for (int i = ndigits - 1; i >= 0; i--) {
			if (d[i] > right.d[i]) return 1;
			else if (d[i] < right.d[i]) return -1;
		}
		return 0;
	}
	int cmp(const u64 *right) const noexcept
	{
#pragma GCC unroll 4
		for (int i = ndigits - 1; i >= 0; i--) {
			if (d[i] > right[i]) return 1;
			else if (d[i] < right[i]) return -1;
		}
		return 0;
	}
	u64 lshift1(const bignum& x) noexcept
	{
		u64 carry = 0;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			u64 temp = x.d[i];
			this->d[i] = (temp << 1) | carry;
			carry = temp >> 63;
		}
		return carry;
	}
	u64 lshift(const uint cnt) noexcept
	{
		u64 carry = 0;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			u64 temp = d[i];
			d[i] = (temp << cnt) | carry;
			carry = temp >> (64 - cnt);
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
#ifdef	NO_BUILTIN_OVERFLOW
		return vli_add<ndigits>(this->d, left.d, right.d);
#else
		bool carry = false;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			u64 sum;
			auto c_carry = __builtin_uaddl_overflow(left.d[i],right.d[i],&sum);
			if (unlikely(carry)) {
				carry = c_carry | __builtin_uaddl_overflow(sum, 1, &sum);
			} else carry = c_carry;
			d[i] = sum;
		}
		return carry;
#endif
	}
/* Computes this = left + right, returning carry. Can modify in place. */
	bool uadd(const bignum& left, const u64 right) noexcept
	{
#ifdef	NO_BUILTIN_OVERFLOW
		return vli_uadd<ndigits>(this->d, left.d, right);
#else
		auto carry = __builtin_uaddl_overflow(left.d[0], right, this->d);
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (unlikely(carry)) {
				carry = __builtin_uaddl_overflow(left.d[i], 1, d+i);
			} else d[i] = left.d[i];
		}
		return carry;
#endif
	}
/* Computes this = this + right, returning carry. Can modify in place. */
	bool add_to(const bignum& right) noexcept
	{
#ifdef	NO_BUILTIN_OVERFLOW
		return vli_add_to<ndigits>(this->d, right.d);
#else
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
#endif
	}
	bool add_to(const u64* right) noexcept
	{
#ifdef	NO_BUILTIN_OVERFLOW
		return vli_add_to<ndigits>(this->d, right);
#else
		bool carry = false;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			u64 sum;
			auto c_carry = __builtin_uaddl_overflow(d[i], right[i], &sum);
			if (unlikely(carry)) {
				carry = c_carry | __builtin_uaddl_overflow(sum, 1, &sum);
			} else carry = c_carry;
			d[i] = sum;
		}
		return carry;
#endif
	}
/* Computes this = this + right, returning carry. Can modify in place. */
	bool uadd_to(const u64 right) noexcept
	{
#ifdef	NO_BUILTIN_OVERFLOW
		return vli_uadd_to<ndigits>(this->d, right);
#else
		auto carry = __builtin_uaddl_overflow(d[0], right, d);
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (unlikely(carry)) {
				carry = __builtin_uaddl_overflow(d[i], 1, d+i);
			} else break;
		}
		return carry;
#endif
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
#ifdef	NO_BUILTIN_OVERFLOW
		return vli_sub<ndigits>(this->d, left.d, right.d);
#else
		bool borrow = false;
#pragma GCC unroll 4
		for (uint i = 0; i < ndigits; i++) {
			u64 diff;
			auto c_borrow=__builtin_usubl_overflow(left.d[i],right.d[i],&diff);
			if (unlikely(borrow)) {
				borrow = c_borrow | __builtin_usubl_overflow(diff, 1, &diff);
			} else borrow = c_borrow;
			d[i] = diff;
		}
		return borrow;
#endif
	}
/* Computes this = left - right, returning borrow. Can modify in place. */
	bool usub(const bignum& left, const u64 right) noexcept
	{
#ifdef	NO_BUILTIN_OVERFLOW
		return vli_usub<ndigits>(this->d, left.d, right);
#else
		auto borrow = __builtin_usubl_overflow(left.d[0], right, this->d);
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (unlikely(borrow)) {
				borrow = __builtin_usubl_overflow(left.d[i], 1, d+i);
			} else d[i] = left.d[i];
		}
		return borrow;
#endif
	}
	bool sub_from(const bignum& right) noexcept
	{
#ifdef	NO_BUILTIN_OVERFLOW
		return vli_sub_from<ndigits>(this->d, right.d);
#else
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
#endif
	}
/* Computes this = this` - right, returning borrow. Can modify in place. */
	bool usub_from(const u64 right) noexcept
	{
#ifdef	NO_BUILTIN_OVERFLOW
		return vli_usub_from<ndigits>(this->d, right);
#else
		auto borrow = __builtin_usubl_overflow(d[0], right, d);
#pragma GCC unroll 4
		for (uint i = 1; i < ndigits; i++) {
			if (unlikely(borrow)) {
				borrow = __builtin_usubl_overflow(d[i], 1, d+i);
			} else break;
		}
		return borrow;
#endif
	}
	bool is_negative() const noexcept
	{
		return (d[ndigits-1] & ((u64)1 << 63)) != 0;
	}
	bool is_even() const noexcept {
		return (d[0] & 1) == 0;
	}
	bool is_odd() const noexcept {
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
		auto _ndigits = num_digits();
		if (_ndigits == 0) return 0;
		auto i = 64 - __builtin_clzl(d[_ndigits - 1]);
		return ((_ndigits - 1) * 64 + i);
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
	bignum& mod_add(const bignum& left, const bignum &right, const bignum &mod)
	noexcept
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
	bignum& mod_sub(const bignum& left, const bignum& right, const bignum &mod)
	noexcept
	{
		auto borrow = this->sub(left, right);
		/* In this case, p_result == -diff == (max int) - diff.
		 * Since -x % d == d - x, we can get the correct result from
		 * result + mod (with overflow).
		 */
		if (borrow)
			this->add_to(mod);
	}

	template<const u64 k0>
	forceinline
	friend void mont_reduction(bignum& res,  const bignum& y,
					const bignum& prime) noexcept
	{
		u64	s[ndigits*2];
		u64	r[ndigits+2];
		vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
		for (uint i=0; i < ndigits; i++) {
			u64	u = (r[0] + y.d[i]) * k0;
			vli_umult<ndigits>(s, prime.d, u);
			vli_uadd_to<ndigits + 2>(r, y.d[i]);
			vli_add_to<ndigits + 2>(r, s);
			vli_rshift1w<ndigits + 2>(r);	
		}
		if (r[ndigits] !=0 || vli_cmp<ndigits>(r, prime.d) >= 0) {
			vli_sub<ndigits>(res.d, r, prime.d);
		} else vli_set<ndigits>(res.d, r);
	}
	void mont_reduction(const bignum& y, const bignum& prime, const u64 k0)
	noexcept
	{
		u64	s[ndigits*2];
		u64	r[ndigits+2];
		vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
		for (uint i=0; i < ndigits; i++) {
			u64	u = (r[0] + y.d[i]) * k0;
			vli_umult<ndigits>(s, prime.d, u);
			vli_uadd_to<ndigits + 2>(r, y.d[i]);
			vli_add_to<ndigits + 2>(r, s);
			vli_rshift1w<ndigits + 2>(r);	
		}
		if (r[ndigits] !=0 || vli_cmp<ndigits>(r, prime.d) >= 0) {
			vli_sub<ndigits>(this->d, r, prime.d);
		} else vli_set<ndigits>(this->d, r);
	}
	template<const u64 k0>
	forceinline
	friend void mont_mult(bignum& res, const bignum& x, const bignum& y,
					const bignum& prime) noexcept
	{
		u64	s[ndigits*2];
		u64	r[ndigits+2];
		vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
		for (uint i=0; i < ndigits;i++) {
			u64	u = (r[0] + y.d[i]*x.d[0]) * k0;
			vli_umult<ndigits>(s, prime.d, u);
			vli_add_to<ndigits + 2>(r, s);
			vli_umult<ndigits>(s, x.d, y.d[i]);
			vli_add_to<ndigits + 2>(r, s);
			vli_rshift1w<ndigits + 2>(r);	
		}
		if (r[ndigits] != 0 || vli_cmp<ndigits>(r, prime.d) >= 0) {
			vli_sub<ndigits>(res.d, r, prime.d);
		} else vli_set<ndigits>(res.d, r);
	}
	template<const u64 k0>
	forceinline
	friend void mont_mult(bignum& res, const u64 *x, const bignum& y,
					const bignum& prime) noexcept
	{
		u64	s[ndigits*2];
		u64	r[ndigits+2];
		vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
		for (uint i=0; i < ndigits;i++) {
			u64	u = (r[0] + y.d[i]*x[0]) * k0;
			vli_umult<ndigits>(s, prime.d, u);
			vli_add_to<ndigits + 2>(r, s);
			vli_umult<ndigits>(s, x, y.d[i]);
			vli_add_to<ndigits + 2>(r, s);
			vli_rshift1w<ndigits + 2>(r);	
		}
		if (r[ndigits] != 0 || vli_cmp<ndigits>(r, prime.d) >= 0) {
			vli_sub<ndigits>(res.d, r, prime.d);
		} else vli_set<ndigits>(res.d, r);
	}
	void mont_mult(const bignum& x, const bignum& y, const bignum& prime,
					const u64 k0) noexcept
	{
		u64	s[ndigits*2];
		u64	r[ndigits+2];
		vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
		for (uint i=0; i < ndigits;i++) {
			u64	u = (r[0] + y.d[i]*x.d[0]) * k0;
			vli_umult<ndigits>(s, prime.d, u);
			vli_add_to<ndigits + 2>(r, s);
			vli_umult<ndigits>(s, x.d, y.d[i]);
			vli_add_to<ndigits + 2>(r, s);
			vli_rshift1w<ndigits + 2>(r);	
		}
		if (r[ndigits] != 0 || vli_cmp<ndigits>(r, prime.d) >= 0) {
			vli_sub<ndigits>(this->d, r, prime.d);
		} else vli_set<ndigits>(this->d, r);
	}
	template<const u64 k0>
	forceinline
	friend void mont_sqr(bignum& res, const bignum& x, const bignum& prime)
	noexcept
	{
		u64	s[ndigits*2];
		u64	r[ndigits+2];
		vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
		for (uint i=0; i < ndigits;i++) {
			u64	u = (r[0] + x.d[i]*x.d[0]) * k0;
			vli_umult<ndigits>(s, prime.d, u);
			vli_add_to<ndigits + 2>(r, s);
			vli_umult<ndigits>(s, x.d, x.d[i]);
			vli_add_to<ndigits + 2>(r, s);
			vli_rshift1w<ndigits + 2>(r);	
		}
		if (r[ndigits] != 0 || vli_cmp<ndigits>(r, prime.d) >= 0) {
			vli_sub<ndigits>(res.d, r, prime.d);
		} else vli_set<ndigits>(res.d, r);
	}
	void mont_sqr(const bignum& x, const bignum& prime, const u64 k0) noexcept
	{
		u64	s[ndigits*2];
		u64	r[ndigits+2];
		vli_clear<ndigits + 2>(r);
#pragma GCC unroll 4
		for (uint i=0; i < ndigits;i++) {
			u64	u = (r[0] + x.d[i]*x.d[0]) * k0;
			vli_umult<ndigits>(s, prime.d, u);
			vli_add_to<ndigits + 2>(r, s);
			vli_umult<ndigits>(s, x.d, x.d[i]);
			vli_add_to<ndigits + 2>(r, s);
			vli_rshift1w<ndigits + 2>(r);	
		}
		if (r[ndigits] != 0 || vli_cmp<ndigits>(r, prime.d) >= 0) {
			vli_sub<ndigits>(this->d, r, prime.d);
		} else vli_set<ndigits>(this->d, r);
	}
protected:
	bn_words	d;
};


template<const uint N>
class bignumz : public bignum<N>{
public:
	bignumz() = default;
	bignumz(const bignumz &) = default;
	bignumz(const bignum<N>& bn) : bignum<N>(bn), carry(0) {}
	explicit bignumz(const s64 v) noexcept : bignum<N>(v), carry(v>=0?0:-1)
	{
		if (v < 0) {
#pragma GCC unroll 4
			for (uint i=1; i<N; i++) this->d[i] = -1;
		}
	}
	explicit bignumz(const u64 *src) noexcept : bignum<N>(src), carry(0) {}
	int is_negative() const noexcept { return carry < 0; }
	bool operator==(const bignumz& bn) {
		return this->carry == bn.carry && vli_cmp<N>(this->d, bn.d) == 0;
	}
	bool operator<(const bignumz& bn) {
		if (this->carry < 0) {
			if (bn.carry >= 0) return true;
			return vli_cmp<N>(this->d, bn.d) < 0;
		} else {
			if (bn.carry < 0) return false;
			return vli_cmp<N>(this->d, bn.d) < 0;
		}
	}
/* Computes result = this + right, returning carry. Can modify in place. */
	void add(const bignumz& left, const bignum<N> &right) noexcept
	{
		this->carry = left.carry;
		if (vli_add<N>(this->d, left.d, right.d)) this->carry++;
	}
	void add(const bignumz& left, const u64* right) noexcept
	{
		this->carry = left.carry;
		if (vli_add<N>(this->d, left.d, right)) this->carry++;
	}
	void sub(const bignumz& left, const bignum<N>& right) noexcept
	{
		this->carry = left.carry;
		if (vli_sub<N>(this->d, left.d, right.d)) this->carry--;
	}
	void sub(const bignumz& left, const u64* right) noexcept
	{
		this->carry = left.carry;
		if (vli_sub<N>(this->d, left.d, right)) this->carry--;
	}
	void add(const bignumz& left, const bignumz& right) noexcept
	{
		this->carry = left.carry + right.carry;
		if (vli_add<N>(this->d, left.d, right.d)) this->carry++;
	}
	void sub(const bignumz& left, const bignumz& right) noexcept
	{
		this->carry = left.carry - right.carry;
		if (vli_sub<N>(this->d, left.d, right.d)) this->carry--;
	}
/* Computes vli = vli >> 1. */
	void rshift1() noexcept
	{
		u64 _rcarry = (carry & 1)?(1L<<63):0;
#pragma GCC unroll 4
		for (int i=N - 1; i >= 0; i--) { 
			u64 temp = this->d[i];
			this->d[i] = (temp >> 1) | _rcarry;
			_rcarry = temp << 63;
		}
		carry >>= 1;
	}	
private:
	int		carry = 0;
};


template<const uint N>
static forceinline
u64 calcK0(const bignum<N>& p) noexcept
{
	u64 t = 1;
	auto n = p.data()[0];
	for (uint i = 1; i < 64; i++) {
		t = t * t * n;		// mod 2^64
	}
	return -t;
}

template<const uint N>
bignum<N>& calcRR(bignum<N>& t, const bignum<N>& p) noexcept
{
	t.clear();
	t.sub_from(p);
	for  (uint i = 256; i<512; i++) {
		if (t.add_to(t) || t.cmp(p) >= 0) {
			t.sub_from(p);
		}
	}
	return t;
}

template<const uint ndigits>
class bn_prod: public bignum<ndigits*2> {
public:
	using bn_words = u64[ndigits*2];
	bn_prod() = default;
	bn_prod(const bn_prod &) = default;
	// this = left * right
	void mult(const bignum<ndigits>& left, const bignum<ndigits>& right)
	noexcept
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
				product.mul_64_64(left.data()[k], right);
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
					u64 _high=(product.m_high() << 1) | (product.m_low() >> 63);
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
	void div_barrett(bignum<ndigits>& result, const bignum<ndigits+1>& mu)
	noexcept
	{
		u64	q[ndigits*2];
		u64	r[ndigits];
		vli_mult<ndigits>(q, this->d + ndigits, mu.data());
		if (mu.data()[ndigits])
			vli_add_to<ndigits>(q + ndigits, this->d + ndigits);
		// add remain * mod
		vli_set<ndigits>(r, q+ndigits);
		vli_umult<ndigits>(q, mu.data(), this->d[ndigits-1]);
		if (mu.data()[ndigits])
			vli_uadd_to<ndigits>(q + ndigits, this->d[ndigits-1]);
		vli_rshift1w<ndigits>(q+ndigits);
		vli_add_to<ndigits>(r, q+ndigits);
		vli_set<ndigits>(const_cast<u64 *>(result.data()), r);
		//result = bignum<ndigits>(r);
	}
	void mmod_barrett(bignum<ndigits>& result, const bignum<ndigits>& prime,
			const bignum<ndigits+1>& mu) noexcept
	{
		u64	q[ndigits*2];
		u64	r[ndigits];
		vli_mult<ndigits>(q, this->d + ndigits, mu.data());
		if (mu.data()[ndigits])
			vli_add_to<ndigits>(q + ndigits, this->d + ndigits);
		// add remain * mod
		vli_set<ndigits>(r, q+ndigits);
		vli_umult<ndigits>(q, mu.data(), this->d[ndigits-1]);
		if (mu.data()[ndigits])
			vli_uadd_to<ndigits>(q + ndigits, this->d[ndigits-1]);
		vli_rshift1w<ndigits>(q+ndigits);
		vli_add_to<ndigits>(r, q+ndigits);
		vli_mult<ndigits>(q, prime.data(), r);
		vli_sub<ndigits*2>(q, this->d, q);
		if (!vli_is_zero<ndigits>(q + ndigits) ||
			vli_cmp<ndigits>(q, prime.data()) >= 0)
		{
			vli_sub_from<ndigits>(q, prime.data());
		}
		vli_set<ndigits>(const_cast<u64 *>(result.data()), q);
		//result = bignum<ndigits>(q);
	}
};

}
#endif	//	__VLI_BN_HPP__
