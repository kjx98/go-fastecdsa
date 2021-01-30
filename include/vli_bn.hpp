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

template<const uint N>
class bn_words {
public:
	bn_words() = default;
	bn_words(const bn_words &) = default;
	explicit bn_words(const u64 v) : d{v, 0, 0, 0}
	{
#if	__cplusplus >= 201703L
		if constexpr(N > 4)
#else
		if ( unlikely(N > 4) )
#endif
		{
			for (uint i = 4; i < N; i++)
				this->d[i] = 0;
		}
	}
protected:
	u64		d[N];
};


template<const uint N>
class bignum : public bn_words<N> {
public:
	bignum() = default;
	bignum(const bignum &) = default;
	constexpr bignum(bn_words<N> init) noexcept : bn_words<N>{init} {}
	explicit bignum(const u64 v) noexcept : bn_words<N>(v) {}
	explicit bignum(const u64 *src) {
		auto	*ss=reinterpret_cast<const bignum<N> *>(src);
		*this = *ss;
	}
	bignum<4>& bn256() noexcept {
		static_assert(N >= 4, "ndigits >= 4 required.");
		bignum<4>	*res=reinterpret_cast<bignum<4> *>(this->d);
		return *res;
	}
	void clear() noexcept
	{
		for (uint i = 0; i < N; i++)
			this->d[i] = 0;
	}
	friend std::ostream& operator<<(std::ostream& os, const bignum& x) noexcept
	{
		os << std::hex;
		if (x.is_u64()) {
			os << x.d[0];
		} else {
			for (uint i=0; i < N; i++) {
				os << " " << x.d[i];
			}
		}
		os << std::dec;
		return os;
	}
	bool operator<(const bignum& bn) const noexcept {
		for (int i = N - 1; i >= 0; i--) {
			if (this->d[i] > bn.d[i]) return false;
			else if (this->d[i] < bn.d[i]) return true;
		}
		return false;
	}
	bool operator>=(const bignum& bn) const noexcept {
		return !(*this < bn);
	}
	bool operator==(const bignum& bn) const noexcept {
#ifdef	NO_U64ZERO
		for (int i = N - 1; i >= 0; i--) {
			if (this->d[i] != bn.d[i]) return false;
		}
		return true;
#else
#if	__cplusplus > 201703L
		if constexpr(N == 4) {
			return u64IsZero(this->d[0] ^ bn.d[0]) &
					u64IsZero(this->d[1] ^ bn.d[1]) &
					u64IsZero(this->d[2] ^ bn.d[2]) &
					u64IsZero(this->d[3] ^ bn.d[3]);
		} else
#endif
		{
			int	ret = u64IsZero(this->d[0] ^ bn.d[0]);
			for (uint i=1; i < N; ++i) ret &= u64IsZero(this->d[i] ^ bn.d[i]);
			return ret;
		}
#endif
	}
	bool operator!=(const bignum& bn) const noexcept {
		return ! (*this == bn);
	}
	const u64* data() const noexcept { return this->d; }
	//u64* raw_data() noexcept { return this->d; }
	int is_zero() const noexcept
	{
#ifndef	NO_U64ZERO
#if	__cplusplus > 201703L
		if constexpr(N == 4) {
			return u64IsZero(this->d[0]) & u64IsZero(this->d[1]) &
					u64IsZero(this->d[2]) & u64IsZero(this->d[3]);
		} else
#endif
		{
			int		ret=u64IsZero(this->d[0]);
			for (uint i = 1; i < N; i++) ret &= u64IsZero(this->d[i]);
			return ret;
		}
#else
		for (uint i = 0; i < N; i++) {
			if (this->d[i]) return 0;
		}
		return true;
#endif
	}
	inline constexpr explicit operator bool() const noexcept {
		return !is_zero();
	}
	int is_one() const noexcept
	{
#ifndef	NO_U64ZERO
#if	__cplusplus > 201703L
		if constexpr(N == 4) {
			return u64IsOne(this->d[0]) & u64IsZero(this->d[1]) &
					u64IsZero(this->d[2]) & u64IsZero(this->d[3]);
		} else
#endif
		{
			int	ret = u64IsOne(this->d[0]);
			for (uint i = 1; i < N; i++) ret &= u64IsZero(this->d[i]);
			return ret;
		}
#else
		if (this->d[0] != 1) return 0;
		for (uint i = 1; i < N; i++) {
			if (this->d[i]) return 0;
		}
		return true;
#endif
	}
	int is_u64() const noexcept
	{
		int	ret = u64IsZero(this->d[1]);
		for (uint i = 2; i < N; i++) {
			ret &= u64IsZero(this->d[i]);
		}
		return ret;
	}
	/* Sets dest = this, copyout */
	void set(u64 *dest) const noexcept
	{
		auto	*res=reinterpret_cast<bignum<N> *>(dest);
		*res = *this;
	}
	/* Returns nonzero if bit bit of vli is set. */
	bool test_bit(const uint bit) const noexcept
	{
		if ( bit >= N*64 ) return false;	// out of bound
		return (this->d[bit >> 6] & ((u64)1 << (bit & 0x3f)));
	}
	uint get_bit(const uint bit) const noexcept
	{
		if (bit >= N*64) return 0;
		auto off = bit >> 6;
		auto rem = bit & 0x3f;
		return (this->d[off] >> rem) & 1;
	}
	/* copy_conditional copies in to out if mask is all ones. */
	void copy_conditional(const bignum& in, u64 mask)
	{
		for (uint i = 0; i < N; ++i) {
			const u64 tmp = mask & (in.d[i] ^ this->d[i]);
			this->d[i] ^= tmp;
		}
	}
/**
 * cmp() - compare this and right vlis
 *
 * @left:		vli
 * @right:		vli
 *
 * Returns sign of @left - @right, i.e. -1 if @left < @right,
 * 0 if @left == @right, 1 if @left > @right.
 */
	int cmp(const bignum& right) const noexcept
	{
		for (int i = N - 1; i >= 0; --i) {
			if (this->d[i] > right.d[i]) return 1;
			else if (this->d[i] < right.d[i]) return -1;
		}
		return 0;
	}
	int cmp(const u64 *right) const noexcept
	{
#ifndef	ommit
		auto	*rt=reinterpret_cast<const bignum<N> *>(right);
		return this->cmp(*rt);
#else
		for (int i = N - 1; i >= 0; --i) {
			if (this->d[i] > right[i]) return 1;
			else if (this->d[i] < right[i]) return -1;
		}
		return 0;
#endif
	}
	u64 lshift1(const bignum& x) noexcept
	{
		u64 carry = 0;
		for (uint i = 0; i < N; i++) {
			u64 temp = x.d[i];
			this->d[i] = (temp << 1) | carry;
			carry = temp >> 63;
		}
		return carry;
	}
	u64 lshift(const uint cnt) noexcept
	{
		u64 carry = 0;
		for (uint i = 0; i < N; i++) {
			u64 temp = this->d[i];
			this->d[i] = (temp << cnt) | carry;
			carry = temp >> (64 - cnt);
		}
		return carry;
	}
/* Computes vli = vli >> 1. */
	void rshift1(u64 carry=0) noexcept
	{
		if (carry) carry = 1L << 63;
		for (int i=N - 1; i >= 0; i--) { 
			u64 temp = this->d[i];
			this->d[i] = (temp >> 1) | carry;
			carry = temp << 63;
		}
	}	
/* Computes vli = vli >> 1 word (64 Bits). */
	void rshift1w() noexcept
	{
		for (uint i = 1; i < N; i++) {
			this->d[i-1] = this->d[i];
		}
		this->d[N-1] = 0;
	}
/* Computes result = this + right, returning carry. Can modify in place. */
	bool add(const bignum& left, const bignum &right) noexcept
	{
		return vli_add<N>(this->d, left.d, right.d);
	}
#ifdef	ommit
/* Computes this = left + right, returning carry. Can modify in place. */
	bool uadd(const bignum& left, const u64 right) noexcept
	{
		return vli_uadd<N>(this->d, left.d, right);
	}
#endif
/* Computes this = this + right, returning carry. Can modify in place. */
	bool add_to(const bignum& right) noexcept
	{
		return vli_add_to<N>(this->d, right.d);
	}
	bool add_to(const u64* right) noexcept
	{
		auto	*rt=reinterpret_cast<const bignum<N> *>(right);
		return this->add_to(*rt);
	}
/* Computes this = this + right, returning carry. Can modify in place. */
	bool uadd_to(const u64 right) noexcept
	{
		return vli_uadd_to<N>(this->d, right);
	}
/* Computes result = this + right, modulo prime. Can modify in place. */
	void mod_add(const bignum& left, const bignum &right, const bignum& prime)
			noexcept
	{
		if (vli_add<N>(this->d, left.d, right.d) ||
			vli_cmp<N>(this->d, prime.d) >= 0)
		{
			vli_sub_from<N>(this->d, prime.d);
		}
	}
#ifdef	ommit
/* Computes this = left + right, modulo prime. Can modify in place. */
	void mod_uadd(const bignum& left, const u64 right, const bignum& prime)
		noexcept
	{
		if (vli_uadd<N>(this->d, left.d, right) ||
			vli_cmp<N>(this->d, prime.d) >= 0)
		{
			vli_sub_from<N>(this->d, prime.d);
		}
	}
#endif
/* Computes this = this + right, modulo prime. Can modify in place. */
	void mod_add_to(const bignum& right, const bignum& prime) noexcept
	{
#if	__cplusplus >= 201703L && defined(__x86_64__)
		if constexpr(N == 4) mod4_add_to(this->d, right.d, prime.d); else
#endif
		if (vli_add_to<N>(this->d, right.d) ||
			vli_cmp<N>(this->d, prime.d) >= 0)
		{
			vli_sub_from<N>(this->d, prime.d);
		}
	}
/* Computes this = this + right, modulo prime. Can modify in place. */
	void mod_uadd_to(const u64 right, const bignum& prime) noexcept
	{
		if (vli_uadd_to<N>(this->d, right) ||
			vli_cmp<N>(this->d, prime.d) >= 0)
		{
			vli_sub_from<N>(this->d, prime.d);
		}
	}
/**
 * sub() - Subtracts right from this
 *
 * @result:		where to write result(this)
 * @left:		vli
 * @right		vli
 *
 * Note: can modify in-place.
 *
 * Return: carry bit.
 */
	bool sub(const bignum& left, const bignum& right) noexcept
	{
		return vli_sub<N>(this->d, left.d, right.d);
	}
#ifdef	ommit
/* Computes this = left - right, returning borrow. Can modify in place. */
	bool usub(const bignum& left, const u64 right) noexcept
	{
		return vli_usub<N>(this->d, left.d, right);
	}
#endif
	bool sub_from(const bignum& right) noexcept
	{
		return vli_sub_from<N>(this->d, right.d);
	}
	bool sub_from(const u64* right) noexcept
	{
		auto	*rt=reinterpret_cast<const bignum<N> *>(right);
		return sub_from(*rt);
	}
/* Computes this = (left - right) % mod.
 * Assumes that left < mod and right < mod, result != mod.
 */
	void mod_sub(const bignum& left, const bignum& right, const bignum& prime) noexcept
	{
		if (vli_sub<N>(this->d, left.d, right.d)) vli_add_to<N>(this->d, prime.d);
	}
	void mod_sub_from(const bignum& right, const bignum& prime) noexcept
	{
		if (vli_sub_from<N>(this->d, right.d)) vli_add_to<N>(this->d, prime.d);
	}
/* Computes this = this` - right, returning borrow. Can modify in place. */
	bool usub_from(const u64 right) noexcept
	{
		return vli_usub_from<N>(this->d, right);
	}
#ifdef	ommit
	bool is_negative() const noexcept
	{
		return (this->d[N-1] & ((u64)1 << 63)) != 0;
	}
#endif
	int is_even() const noexcept {
		return (this->d[0] & 1) ^ 1;
	}
	int is_odd() const noexcept {
		//return (this->d[0] & 1) != 0;
		return this->d[0] & 1;
	}
/* Counts the number of 64-bit "digits" in vli. */
	uint num_digits() const noexcept
	{
		/* Search from the end until we find a non-zero digit.
		 * We do it in reverse because we expect that most digits will
		 * be nonzero.
		 */
		int i;
		for (i = N - 1; i >= 0 && this->d[i] == 0; i--);
		return (i + 1);
	}
	uint num_bits() const noexcept
	{
		auto _ndigits = num_digits();
		if (_ndigits == 0) return 0;
		auto i = 64 - __builtin_clzl(this->d[_ndigits - 1]);
		return ((_ndigits - 1) * 64 + i);
	}
/**
 * vli_from_be64() - Load vli from big-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 BE values
 */
	void from_be64(const void *src) noexcept
	{
		const u64 *from = (const u64 *)src;
		for (uint i = 0; i < N; i++)
			this->d[i] = be64toh(from[N - 1 - i]);
	}
#ifdef	ommit
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
#endif
/* Computes this = (left - right) % mod.
 * Assumes that left < mod and right < mod, result != mod.
 */
#ifdef	ommit
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
#endif
	template<const u64 k0> forceinline
	friend void mont_reduction(bignum& res,  const bignum& y,
					const bignum& prime) noexcept
	{
		u64	s[N*2];
		u64	r[N+2];
		vli_clear<N + 2>(r);
		for (uint i=0; i < N; i++) {
			u64	u = (r[0] + y.d[i]) * k0;
			vli_umult<N>(s, prime.d, u);
			vli_uadd_to<N + 2>(r, y.d[i]);
			vli_add_to<N + 2>(r, s);
			vli_rshift1w<N + 2>(r);	
		}
		if (r[N] !=0 || vli_cmp<N>(r, prime.d) >= 0) {
			vli_sub<N>(res.d, r, prime.d);
		} else vli_set<N>(res.d, r);
	}
	void mont_reduction(const bignum& y, const bignum& prime, const u64 k0)
	noexcept
	{
		u64	s[N*2];
		u64	r[N+2];
		vli_clear<N + 2>(r);
		for (uint i=0; i < N; i++) {
			u64	u = (r[0] + y.d[i]) * k0;
			vli_umult<N>(s, prime.d, u);
			vli_uadd_to<N + 2>(r, y.d[i]);
			vli_add_to<N + 2>(r, s);
			vli_rshift1w<N + 2>(r);	
		}
		if (r[N] !=0 || vli_cmp<N>(r, prime.d) >= 0) {
			vli_sub<N>(this->d, r, prime.d);
		} else vli_set<N>(this->d, r);
	}
	template<const u64 k0> forceinline
	friend void mont_mult(bignum& res, const bignum& x, const bignum& y,
					const bignum& prime) noexcept
	{
		u64	s[N*2];
		u64	r[N+2];
		vli_clear<N + 2>(r);
		for (uint i=0; i < N;i++) {
			u64	u = (r[0] + y.d[i]*x.d[0]) * k0;
			vli_umult<N>(s, prime.d, u);
			vli_add_to<N + 2>(r, s);
			vli_umult<N>(s, x.d, y.d[i]);
			vli_add_to<N + 2>(r, s);
			vli_rshift1w<N + 2>(r);	
		}
		if (r[N] != 0 || vli_cmp<N>(r, prime.d) >= 0) {
			vli_sub<N>(res.d, r, prime.d);
		} else vli_set<N>(res.d, r);
	}
	template<const u64 k0> forceinline
	friend void mont_mult(bignum& res, const u64 *x, const bignum& y,
					const bignum& prime) noexcept
	{
		u64	s[N*2];
		u64	r[N+2];
		vli_clear<N + 2>(r);
		for (uint i=0; i < N;i++) {
			u64	u = (r[0] + y.d[i]*x[0]) * k0;
			vli_umult<N>(s, prime.d, u);
			vli_add_to<N + 2>(r, s);
			vli_umult<N>(s, x, y.d[i]);
			vli_add_to<N + 2>(r, s);
			vli_rshift1w<N + 2>(r);	
		}
		if (r[N] != 0 || vli_cmp<N>(r, prime.d) >= 0) {
			vli_sub<N>(res.d, r, prime.d);
		} else vli_set<N>(res.d, r);
	}
	void mont_mult(const bignum& x, const bignum& y, const bignum& prime,
					const u64 k0) noexcept
	{
		u64	s[N*2];
		u64	r[N+2];
		vli_clear<N + 2>(r);
		for (uint i=0; i < N;i++) {
			u64	u = (r[0] + y.d[i]*x.d[0]) * k0;
			vli_umult<N>(s, prime.d, u);
			vli_add_to<N + 2>(r, s);
			vli_umult<N>(s, x.d, y.d[i]);
			vli_add_to<N + 2>(r, s);
			vli_rshift1w<N + 2>(r);	
		}
		if (r[N] != 0 || vli_cmp<N>(r, prime.d) >= 0) {
			vli_sub<N>(this->d, r, prime.d);
		} else vli_set<N>(this->d, r);
	}
	template<const u64 k0> forceinline friend
	void mont_sqr(bignum& res, const bignum& x, const bignum& prime) noexcept
	{
		u64	s[N*2];
		u64	r[N+2];
		vli_clear<N + 2>(r);
		for (uint i=0; i < N;i++) {
			u64	u = (r[0] + x.d[i]*x.d[0]) * k0;
			vli_umult<N>(s, prime.d, u);
			vli_add_to<N + 2>(r, s);
			vli_umult<N>(s, x.d, x.d[i]);
			vli_add_to<N + 2>(r, s);
			vli_rshift1w<N + 2>(r);	
		}
		if (r[N] != 0 || vli_cmp<N>(r, prime.d) >= 0) {
			vli_sub<N>(res.d, r, prime.d);
		} else vli_set<N>(res.d, r);
	}
	void mont_sqr(const bignum& x, const bignum& prime, const u64 k0) noexcept
	{
		u64	s[N*2];
		u64	r[N+2];
		vli_clear<N + 2>(r);
		for (uint i=0; i < N;i++) {
			u64	u = (r[0] + x.d[i]*x.d[0]) * k0;
			vli_umult<N>(s, prime.d, u);
			vli_add_to<N + 2>(r, s);
			vli_umult<N>(s, x.d, x.d[i]);
			vli_add_to<N + 2>(r, s);
			vli_rshift1w<N + 2>(r);	
		}
		if (r[N] != 0 || vli_cmp<N>(r, prime.d) >= 0) {
			vli_sub<N>(this->d, r, prime.d);
		} else vli_set<N>(this->d, r);
	}
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
			for (uint i=1; i<N; i++) this->d[i] = -1;
		}
	}
	explicit bignumz(const u64 *src) noexcept : bignum<N>(src), carry(0) {}
	int is_negative() const noexcept { return carry < 0; }
	bignum<N>& abs() const noexcept {
		u64	*pp = const_cast<u64 *>(this->d);
		bignum<N> *res = reinterpret_cast<bignum<N> *>(pp);
		return *res;
	}
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
		if (t.add_to(t) || t.cmp(p) >= 0)
		//if (t.lshift1(t) || t.cmp(p) >= 0)
		{
			t.sub_from(p);
		}
	}
	//if (t.cmp(p) >= 0) t.sub_from(p);
	return t;
}

template<const uint N>
class bn_prod: public bignum<N*2> {
public:
	bn_prod() = default;
	bn_prod(const bn_prod &) = default;
	// this = left * right
	void mult(const bignum<N>& left, const bignum<N>& right)
	noexcept
	{
		uint128_t r01( 0, 0 );
	/* Compute each digit of result in sequence, maintaining the
	 * carries.
	 */
		for (uint k = 0; k < N * 2 - 1; k++) {
			u64 r2 = 0;
			unsigned int min;
			if (k < N) min = 0; else min = (k + 1) - N;
			for (uint i = min; i <= k && i < N; i++) {
				uint128_t product;
				product.mul_64_64(left.d[i], right.d[k - i]);
				r01 += product;
				r2 += (r01.m_high() < product.m_high());
			}
			this->d[k] = r01.m_low();
			r01 = uint128_t(r01.m_high(), r2);
			r2 = 0;
		}
		this->d[N * 2 - 1] = r01.m_low();
	}

	/* Compute product = left * right, for a small right value. */
	void umult(const bignum<N>&left, const u64 right) noexcept
	{
		uint128_t r01( 0, 0 );
		unsigned int k;
		for (k = 0; k < N; k++) {
			uint128_t product;
			//if (likely(left[k] != 0))
			{
				product.mul_64_64(left.data()[k], right);
				r01 += product;
			}
			/* no carry */
			this->d[k] = r01.m_low();
			r01 = uint128_t(r01.m_high(), 0);
		}
		this->d[N] = r01.m_low();
		for (k = N+1; k < N * 2; k++) this->d[k] = 0;
	}

	void square(const bignum<N>& left) noexcept
	{
		uint128_t r01( 0, 0 );
		for (uint k = 0; k < N * 2 - 1; k++) {
			unsigned int min;
			u64 r2 = 0;
			if (k < N) min = 0; else min = (k + 1) - N;
			for (uint i = min; i <= k && i <= k - i; i++) {
				uint128_t product;
				product.mul_64_64(left.d[i], left.d[k - i]);
				if (i < k - i) {
					r2 += product.m_high() >> 63;
					u64 _high=(product.m_high() << 1) | (product.m_low() >> 63);
					u64 _low = product.m_low() << 1;
					product = uint128_t(_low, _high);
				}
				r01 += product;
				r2 += (r01.m_high() < product.m_high());
			}
			this->d[k] = r01.m_low();
			r01 = uint128_t(r01.m_high(), r2);
		}
		this->d[N * 2 - 1] = r01.m_low();
	}
	void div_barrett(bignum<N>& result, const bignum<N+1>& mu) noexcept
	{
		u64	q[N*2];
		u64	r[N];
		vli_mult<N>(q, this->d + N, mu.data());
		if (mu.data()[N])
			vli_add_to<N>(q + N, this->d + N);
		// add remain * mod
		vli_set<N>(r, q+N);
		vli_umult<N>(q, mu.data(), this->d[N-1]);
		if (mu.data()[N])
			vli_uadd_to<N>(q + N, this->d[N-1]);
		vli_rshift1w<N>(q+N);
		vli_add_to<N>(r, q+N);
		vli_set<N>(const_cast<u64 *>(result.data()), r);
		//result = bignum<N>(r);
	}
	void mmod_barrett(bignum<N>& result, const bignum<N>& prime,
			const bignum<N+1>& mu) noexcept
	{
		u64	q[N*2];
		u64	r[N];
		vli_mult<N>(q, this->d + N, mu.data());
		if (mu.data()[N])
			vli_add_to<N>(q + N, this->d + N);
		// add remain * mod
		vli_set<N>(r, q+N);
		vli_umult<N>(q, mu.data(), this->d[N-1]);
		if (mu.data()[N])
			vli_uadd_to<N>(q + N, this->d[N-1]);
		vli_rshift1w<N>(q+N);
		vli_add_to<N>(r, q+N);
		vli_mult<N>(q, prime.data(), r);
		vli_sub<N*2>(q, this->d, q);
		if (!vli_is_zero<N>(q + N) ||
			vli_cmp<N>(q, prime.data()) >= 0)
		{
			vli_sub_from<N>(q, prime.data());
		}
		vli_set<N>(const_cast<u64 *>(result.data()), q);
		//result = bignum<N>(q);
	}
};

}
#endif	//	__VLI_BN_HPP__
