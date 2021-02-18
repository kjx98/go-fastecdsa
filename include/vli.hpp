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

#include "u128.hpp"


template<typename T> forceinline
static T max(const T&& a, const T&& b) noexcept
{
	static_assert(std::is_integral<T>::value, "Integral required.");
	return (a>b)?a:b;
}


//__attribute__((optimize("unroll-loops")))
template<const uint N> forceinline
static void vli_clear(u64 *vli) noexcept
{
	for (uint i = 0; i < N; i++)
		vli[i] = 0;
}


/*
 * X86_64 func param:  di/si/dx/cx/r8/r9
 *			retutn register: ax
 *			temp register:	r10/r11
 *			callee-saved register: rbx/rbp,		r12/r13/r14/r15
 *
 * ARM64 func param:	X0-X7
 *			retutn register: X0
 * 			temp register: X9-X15
 * 			callee-saved register: X19-X29
 *
 */
// u64IsZero returns 1 if x is zero and zero otherwise.
static forceinline int u64IsZero(u64 x) noexcept
{
#if	defined(__x86_64__) && __GNUC__ > 5
	bool	ret;
	asm volatile(
		"subq $0, %%rax\n"
			: "=@ccz" (ret)
			: "a" (x)
			: "cc");
	return (int)ret;
#else
	x = ~x;
	x &= x >> 32;
	x &= x >> 16;
	x &= x >> 8;
	x &= x >> 4;
	x &= x >> 2;
	x &= x >> 1;
	return x & 1;
#endif
}

// u64IsOne returns 1 if x is one and zero otherwise.
static forceinline int u64IsOne(u64 x) noexcept
{
#if	defined(__x86_64__) && __GNUC__ > 5
	bool	ret;
	asm volatile(
		"subq $1, %%rax\n"
			: "=@ccz" (ret)
			: "a" (x)
			: "cc");
	return (int)ret;
#else
	return u64IsZero(x ^ 1);
#endif
}


// result = a + b + carry, set new carry
static forceinline u64 u64_addc(const u64 a, const u64 b, u64& carry)
{
#if	defined(__x86_64__) && __GNUC__ > 5
	u64		ret=0;
	asm volatile(
		"movq %1, %0\n"
		"xorq %1, %1\n"
		"addq %2, %0\n"
		"adcq $0, %1\n"
		"addq %3, %0\n"
		"adcq $0, %1\n"
			: "+r" (ret), "+r" (carry)
			: "r" (a),"rm" (b)
			: "cc");
	return ret;
#elif	defined(__aarch64__notwork)
	u64		ret=carry;
	asm volatile(
		"adds	%0, %0, %2\n"
		"adc	%1, xzr, xzr\n"
		"adds	%0, %0, %3\n"
		"adc %1, %1, xzr\n"
			: "+r" (ret), "+r"(carry)
			: "r" (a), "r" (b)
			: "cc");
	return ret;
#elif	__clang__ > 3
	u64		ret, c_out;
	ret = __builtin_addcl(a, b, carry, &c_out);
	carry = c_out;
	return ret;
#else
	u64		ret;
	ret = a + b + carry;
	if (ret != a) carry = ret < a;
	return ret;
#endif
}

static forceinline u64 u64_addcz(const u64 a, u64& carry)
{
#if	defined(__x86_64__) && __GNUC__ > 5
	u64		ret=a;
	asm volatile(
		"addq %1, %0\n"
		"movq $0, %1\n"
		"adcq $0, %1\n"
			: "+r" (ret), "+r" (carry)
			:
			: "cc");
	return ret;
#elif	defined(__aarch64__)
	u64		ret;
	asm volatile(
		"adds	%0, %1, %2\n"
		"adc	%1, xzr, xzr\n"
			: "=r" (ret), "+r"(carry)
			: "r" (a)
			: "cc");
	return ret;
#elif	__clang__ > 3
	u64		ret, c_out;
	ret = __builtin_addcl(a, 0, carry, &c_out);
	carry = c_out;
	return ret;
#else
	u64		ret;
	ret = a + carry;
	if (ret != a) carry = ret < a;
	return ret;
#endif
}

// result = a - b - carry, set new carry
static forceinline u64 u64_subc(const u64 a, const u64 b, u64& carry)
{
#if	defined(__x86_64__) && __GNUC__ > 5
	u64		ret=a;
	asm volatile(
		"subq %1, %0\n"
		"movq $0, %1\n"
		"adcq $0, %1\n"
		"subq %2, %0\n"
		"adcq $0, %1\n"
			: "+r" (ret), "+r" (carry)
			: "rm" (b) 
			: "cc");
	return ret;
#elif	defined(__aarch64__notwork)
	u64		ret=a;
	asm volatile(
		"subs %0, %0, %1\n"
		"sbc %1, xzr, xzr\n"
		"subs %0, %0, %2\n"
		"sbc %1, %1, xzr\n"
			: "+r" (ret), "+r" (carry)
			: "r" (b)
			: "cc");
	return ret;
#elif	__clang__ > 3
	u64		ret, c_out;
	ret = __builtin_subcl(a, b, carry, &c_out);
	carry = c_out;
	return ret;
#else
	u64		ret;
	ret = a - b - carry;
	if (ret != a) carry = ret > a;
	return ret;
#endif
}

static forceinline u64 u64_subcz(const u64 a, u64& carry)
{
#if	defined(__x86_64__) && __GNUC__ > 5
	u64		ret=a;
	asm volatile(
		"subq %1, %0\n"
		"movq $0, %1\n"
		"adcq $0, %1\n"
			: "+r" (ret), "+r" (carry)
			:
			: "cc");
	return ret;
#elif	defined(__aarch64__)
	u64		ret=a;
	asm volatile(
		"subs %0, %0, %1\n"
		"sbc %1, xzr, xzr\n"
			: "+r" (ret), "+r" (carry)
			:
			: "cc");
	return ret;
#elif	__clang__ > 3
	u64		ret, c_out;
	ret = __builtin_subcl(a, 0, carry, &c_out);
	carry = c_out;
	return ret;
#else
	u64		ret;
	ret = a - carry;
	if (ret != a) carry = ret > a;
	return ret;
#endif
}

#ifdef	WITH_ASM
#ifdef	__x86_64__
static forceinline void
vli4_mod(u64 *res, const u64 *left, const u64 *mod, const bool carry) noexcept
{
	/* result > mod (result = mod + remainder), so subtract mod to
	 * get remainder.
	 */
	asm volatile(
		"movq (%%rsi), %%r8\n"	// mov r8/9/10/11, right
		"movq 8(%%rsi), %%r9\n"
		"movq 16(%%rsi), %%r10\n"
		"movq 24(%%rsi), %%r11\n"
		"movq %%r8, %%r12\n"
		"movq %%r9, %%r13\n"
		"movq %%r10, %%r14\n"
		"movq %%r11, %%r15\n"
		"movq %[mod], %%rsi\n"
		"subq (%%rsi), %%r12\n"
		"sbbq 8(%%rsi), %%r13\n"
		"sbbq 16(%%rsi), %%r14\n"
		"sbbq 24(%%rsi), %%r15\n"
		"sbbq $0, %%rax\n"
		"cmovcq %%r8, %%r12\n"
		"cmovcq %%r9, %%r13\n"
		"cmovcq %%r10, %%r14\n"
		"cmovcq %%r11, %%r15\n"
		"movq %%r12, (%%rdi)\n"
		"movq %%r13, 8(%%rdi)\n"
		"movq %%r14, 16(%%rdi)\n"
		"movq %%r15, 24(%%rdi)\n"
			:
			: "S"(left), "D"(res), "a" ((u64)carry), [mod] "m" (mod)
			: "%r8", "%r9", "%r10", "%r11" , "%r12", "%r13", "%r14",
			"%r15", "cc", "memory");
}

static forceinline
void mod4_add_to(u64 *left, const u64 *right, const u64 *mod) noexcept
{
	asm volatile(
		"movq (%%rsi), %%r8\n"	// mov r8/9/10/11, right
		"movq 8(%%rsi), %%r9\n"
		"movq 16(%%rsi), %%r10\n"
		"movq 24(%%rsi), %%r11\n"
		"xorq %%rax, %%rax\n"
		"addq (%%rdi), %%r8\n"		// add
		"adcq 8(%%rdi), %%r9\n"
		"adcq 16(%%rdi), %%r10\n"
		"adcq 24(%%rdi), %%r11\n"
		"adcq $0, %%rax\n"
		"movq %%r8, %%r12\n"
		"movq %%r9, %%r13\n"
		"movq %%r10, %%r14\n"
		"movq %%r11, %%r15\n"
		"movq %[mod], %%rsi\n"
		"subq (%%rsi), %%r12\n"
		"sbbq 8(%%rsi), %%r13\n"
		"sbbq 16(%%rsi), %%r14\n"
		"sbbq 24(%%rsi), %%r15\n"
		"sbbq $0, %%rax\n"
		"cmovcq %%r8, %%r12\n"
		"cmovcq %%r9, %%r13\n"
		"cmovcq %%r10, %%r14\n"
		"cmovcq %%r11, %%r15\n"
		"movq %%r12, (%%rdi)\n"
		"movq %%r13, 8(%%rdi)\n"
		"movq %%r14, 16(%%rdi)\n"
		"movq %%r15, 24(%%rdi)\n"
			:
			: "S"(right), "D"(left), [mod] "m" (mod)
			: "rax", "%r8", "%r9", "%r10", "%r11" , "%r12", "%r13", "%r14",
			"%r15", "cc", "memory");
}

static forceinline void
mod4_add(u64 *res, const u64 *left, const u64 *right, const u64* mod) noexcept
{
	asm volatile(
		"movq (%%rsi), %%r8\n"	// mov r8/9/10/11, right
		"movq 8(%%rsi), %%r9\n"
		"movq 16(%%rsi), %%r10\n"
		"movq 24(%%rsi), %%r11\n"
		"xorq %%rax, %%rax\n"
		"addq (%%rdi), %%r8\n"		// add
		"adcq 8(%%rdi), %%r9\n"
		"adcq 16(%%rdi), %%r10\n"
		"adcq 24(%%rdi), %%r11\n"
		"adcq $0, %%rax\n"
		"movq %[res], %%rdi\n"		// mov rdi, res
		"movq %[mod], %%rsi\n"
		"movq %%r8, %%r12\n"
		"movq %%r9, %%r13\n"
		"movq %%r10, %%r14\n"
		"movq %%r11, %%r15\n"
		"subq (%%rsi), %%r12\n"
		"sbbq 8(%%rsi), %%r13\n"
		"sbbq 16(%%rsi), %%r14\n"
		"sbbq 24(%%rsi), %%r15\n"
		"sbbq $0, %%rax\n"
		"cmovcq %%r8, %%r12\n"
		"cmovcq %%r9, %%r13\n"
		"cmovcq %%r10, %%r14\n"
		"cmovcq %%r11, %%r15\n"
		"movq %%r12, (%%rdi)\n"
		"movq %%r13, 8(%%rdi)\n"
		"movq %%r14, 16(%%rdi)\n"
		"movq %%r15, 24(%%rdi)\n"
			:
			: "S"(right), "D"(left), [res] "m" (res), [mod] "m" (mod)
			: "rax", "%r8", "%r9", "%r10", "%r11" , "%r12", "%r13", "%r14",
			"%r15", "cc", "memory");
}


static forceinline
bool vli4_add_to(u64 *left, const u64 *right) noexcept
{
	bool carry;
	asm volatile(
		"movq (%%rsi), %%r8\n"	// mov r8/9/10/11, right
		"movq 8(%%rsi), %%r9\n"
		"movq 16(%%rsi), %%r10\n"
		"movq 24(%%rsi), %%r11\n"
		"xorq %%rax, %%rax\n"
		"addq (%%rdi), %%r8\n"		// add
		"adcq 8(%%rdi), %%r9\n"
		"adcq 16(%%rdi), %%r10\n"
		"adcq 24(%%rdi), %%r11\n"
		"adcq %%rax, %%rax\n"
		"movq %%r8, (%%rdi)\n"
		"movq %%r9, 8(%%rdi)\n"
		"movq %%r10, 16(%%rdi)\n"
		"movq %%r11, 24(%%rdi)\n"
			: "=a"(carry)
			: "S"(right), "D"(left)
			: "%r8", "%r9", "%r10", "%r11" , "cc", "memory");
	return carry;
}

static forceinline
bool vli4_add(u64 *res, const u64 *left, const u64 *right) noexcept
{
	bool carry;
	asm volatile(
		"movq (%%rsi), %%r8\n"	// mov r8/9/10/11, right
		"movq 8(%%rsi), %%r9\n"
		"movq 16(%%rsi), %%r10\n"
		"movq 24(%%rsi), %%r11\n"
		"xorq %%rax, %%rax\n"
		"addq (%%rdi), %%r8\n"		// add
		"adcq 8(%%rdi), %%r9\n"
		"adcq 16(%%rdi), %%r10\n"
		"adcq 24(%%rdi), %%r11\n"
		"adcq %%rax, %%rax\n"
		"movq %[res], %%rdi\n"
		"movq %%r8, (%%rdi)\n"
		"movq %%r9, 8(%%rdi)\n"
		"movq %%r10, 16(%%rdi)\n"
		"movq %%r11, 24(%%rdi)\n"
			: "=a"(carry)
			: "S"(right), "D"(left), [res] "rm" (res)
			: "%r8", "%r9", "%r10", "%r11" , "cc", "memory");
	return carry;
}

static forceinline
bool vli4_sub_from(u64 *left, const u64 *right) noexcept
{
	bool carry;
	asm volatile(
		"movq (%%rdi), %%r8\n"	// mov r8/9/10/11, left
		"movq 8(%%rdi), %%r9\n"
		"movq 16(%%rdi), %%r10\n"
		"movq 24(%%rdi), %%r11\n"
		"xorq %%rax, %%rax\n"
		"subq (%%rsi), %%r8\n"		// sub
		"sbbq 8(%%rsi), %%r9\n"
		"sbbq 16(%%rsi), %%r10\n"
		"sbbq 24(%%rsi), %%r11\n"
		"adcq %%rax, %%rax\n"
		"movq %%r8, (%%rdi)\n"
		"movq %%r9, 8(%%rdi)\n"
		"movq %%r10, 16(%%rdi)\n"
		"movq %%r11, 24(%%rdi)\n"
			: "=a"(carry)
			: "S"(right), "D"(left)
			: "%r8", "%r9", "%r10", "%r11" , "cc", "memory");
	return carry;
}

static forceinline
bool vli4_sub(u64 *res, const u64 *left, const u64 *right) noexcept
{
	bool carry;
	asm volatile(
		"movq (%%rdi), %%r8\n"	// mov r8/9/10/11, left
		"movq 8(%%rdi), %%r9\n"
		"movq 16(%%rdi), %%r10\n"
		"movq 24(%%rdi), %%r11\n"
		"xorq %%rax, %%rax\n"
		"subq (%%rsi), %%r8\n"		// sub
		"sbbq 8(%%rsi), %%r9\n"
		"sbbq 16(%%rsi), %%r10\n"
		"sbbq 24(%%rsi), %%r11\n"
		"adcq %%rax, %%rax\n"
		"movq %[res], %%rdi\n"		// mov rdi, res
		"movq %%r8, (%%rdi)\n"
		"movq %%r9, 8(%%rdi)\n"
		"movq %%r10, 16(%%rdi)\n"
		"movq %%r11, 24(%%rdi)\n"
			: "=a"(carry)
			: "S"(right), "D"(left), [res] "rm" (res)
			: "%r8", "%r9", "%r10", "%r11" , "cc", "memory");
	return carry;
}
#elif	defined(__aarch64__)
static forceinline void
vli4_mod(u64 *res, const u64 *left, const u64 *mod, const bool carry) noexcept
{
	/* result > mod (result = mod + remainder), so subtract mod to
	 * get remainder.
	 */
	asm volatile(
		"ldp x4, x5, [%2]\n"
		"ldp x6, x7, [%2, 16]\n"
		"ldp x9, x10, [%3]\n"
		"ldp x11, x12, [%3, 16]\n"
		"subs x9, x4, x9\n"
		"sbcs x10, x5, x10\n"
		"sbcs x11, x6, x11\n"
		"sbcs x12, x7, x12\n"
		"sbcs %0, %0, xzr\n"
		"csel x4, x4, x9, cc\n"
		"csel x5, x5, x10, cc\n"
		"csel x6, x6 , x11, cc\n"
		"csel x7, x7, x12, cc\n"
		"stp x4, x5, [%1]\n"
		"stp x6, x7, [%1, 16]\n"
		"adc %0, xzr, xzr\n"
		:
		: "r" ((u64)carry), "r" (res), "r" (left), "r" (mod)
		: "%x4", "%x5", "%x6", "%x7", "%x9", "%x10", "%x11", "%x12",
		"cc", "memory");
}

static forceinline
bool vli4_add_to(u64 *left, const u64 *right) noexcept
{
	u64 carry;
	asm volatile(
		"ldp x4, x5, [%1]\n"
		"ldp x6, x7, [%1, 16]\n"
		"ldp x9, x10, [%2]\n"
		"ldp x11, x12, [%2, 16]\n"
		"adds x4, x4, x9\n"
		"adcs x5, x5, x10\n"
		"adcs x6, x6, x11\n"
		"adcs x7, x7, x12\n"
		"stp x4, x5, [%1]\n"
		"stp x6, x7, [%1, 16]\n"
		"adc %0, xzr, xzr\n"
		: "=r"(carry)
		: "r" (left), "r" (right)
		: "%x4", "%x5", "%x6", "%x7", "%x9", "%x10", "%x11", "%x12",
		"cc", "memory");
	return carry;
}

static forceinline
bool vli4_add(u64 *res, const u64 *left, const u64 *right) noexcept
{
	u64 carry;
	asm volatile(
		"ldp x4, x5, [%2]\n"
		"ldp x6, x7, [%2, 16]\n"
		"ldp x9, x10, [%3]\n"
		"ldp x11, x12, [%3, 16]\n"
		"adds x4, x4, x9\n"
		"adcs x5, x5, x10\n"
		"adcs x6, x6, x11\n"
		"adcs x7, x7, x12\n"
		"stp x4, x5, [%1]\n"
		"stp x6, x7, [%1, 16]\n"
		"adc %0, xzr, xzr\n"
		: "=r"(carry)
		: "r" (res), "r" (left), "r" (right)
		: "%x4", "%x5", "%x6", "%x7", "%x9", "%x10", "%x11", "%x12",
		"cc", "memory");
	return carry;
}

static forceinline
bool vli4_sub_from(u64 *left, const u64 *right) noexcept
{
	u64 carry;
	asm volatile(
		"ldp x4, x5, [%1]\n"
		"ldp x6, x7, [%1, 16]\n"
		"ldp x9, x10, [%2]\n"
		"ldp x11, x12, [%2, 16]\n"
		"subs x4, x4, x9\n"
		"sbcs x5, x5, x10\n"
		"sbcs x6, x6, x11\n"
		"sbcs x7, x7, x12\n"
		"stp x4, x5, [%1]\n"
		"stp x6, x7, [%1, 16]\n"
		"sbc %0, xzr, xzr\n"
		: "=r"(carry)
		: "r" (left), "r" (right)
		: "%x4", "%x5", "%x6", "%x7", "%x9", "%x10", "%x11", "%x12",
		"cc", "memory");
	return carry;
}

static forceinline
bool vli4_sub(u64 *res, const u64 *left, const u64 *right) noexcept
{
	u64 carry;
	asm volatile(
		"ldp x4, x5, [%2]\n"
		"ldp x6, x7, [%2, 16]\n"
		"ldp x9, x10, [%3]\n"
		"ldp x11, x12, [%3, 16]\n"
		"subs x4, x4, x9\n"
		"sbcs x5, x5, x10\n"
		"sbcs x6, x6, x11\n"
		"sbcs x7, x7, x12\n"
		"stp x4, x5, [%1]\n"
		"stp x6, x7, [%1, 16]\n"
		"sbc %0, xzr, xzr\n"
		: "=r"(carry)
		: "r" (res), "r" (left), "r" (right)
		: "%x4", "%x5", "%x6", "%x7", "%x9", "%x10", "%x11", "%x12",
		"cc", "memory");
	return carry;
}
#endif
#endif	// WITH_ASM


forceinline static void
vli4_load(const u64 *vli, u64& r0, u64& r1, u64& r2, u64& r3) noexcept
{
	r0 = vli[0];
	r1 = vli[1];
	r2 = vli[2];
	r3 = vli[3];
}

forceinline static void
vli4_save(u64 *vli, u64 r0, u64 r1, u64 r2, u64 r3) noexcept
{
	vli[0] = r0;
	vli[1] = r1;
	vli[2] = r2;
	vli[3] = r3;
}


/**
 * vli_is_zero() - Determine is vli is zero
 *
 * @vli:		vli to check.
 * @ndigits:		length of the @vli
 */
/* Returns true if vli == 0, false otherwise. */
template<const uint N> forceinline
static int vli_is_zero(const u64 *vli) noexcept
{
	static_assert(N>=4, "N MUST >= 4");
	u64	v=vli[0];
	for (uint i = 1; i < N; i++) {
		v |= vli[i];
	}
#ifdef	NO_U64ZERO
	return v == 0;
#else
	return u64IsZero(v);
#endif
}

template<const uint N> forceinline
static int vli_is_one(const u64 *vli) noexcept
{
	static_assert(N>=4, "N MUST >= 4");
	u64	v=vli[0] ^ 1;
	for (uint i = 1; i < N; i++) {
		v |= vli[i];
	}
#ifdef	NO_U64ZERO
	return v == 0;
#else
	return u64IsZero(v);
#endif
}

/* Sets dest = src. */
template<const uint N> forceinline
static void vli_set(u64 *dest, const u64 *src) noexcept
{
	for (uint i = 0; i < N; i++)
		dest[i] = src[i];
}

/* Returns nonzero if bit bit of vli is set. */
template<const uint N> forceinline
static bool vli_test_bit(const u64 *vli, uint bit) noexcept
{
	if ( bit >= N*64 ) return false;	// out of bound
	return (vli[bit / 64] & ((u64)1 << (bit % 64)));
}


// -1 <= bit < 64*N, cnt <= 8
template<const uint N, const uint cnt=5> forceinline
static uint vli_get_bits(const u64 *vli, const int bit) noexcept
{
	static_assert(cnt <= 8, "w of wNAF MUST no larger than 8");
	uint	rr;
	if (unlikely(bit < 0)) {
		// bit MUST BE -1
		rr = vli[0] << 1;
		return rr & ((1<<cnt) - 1);
	}
	auto off = (uint)bit >> 6;
	if (off >= N) return 0;
	auto rem = (uint)bit & 0x3f;
	rr = (vli[off] >> rem); // & 0xff;
	if ((uint)bit < (N-1)*64 && rem > (64 - cnt)) {
		++off;
		rr |= (vli[off] << (64 - rem));
	}
	return rr & ((1<<cnt) - 1);
}

/**
 * vli_cmp() - compare left and right vlis
 *
 * @left:		vli
 * @right:		vli
 * @ndigits:	N, length of both vlis
 *
 * Returns sign of @left - @right, i.e. -1 if @left < @right,
 * 0 if @left == @right, 1 if @left > @right.
 */
template<const uint N> forceinline
static int vli_cmp(const u64 *left, const u64 *right) noexcept
{
	for (uint i = N - 1; i < N; i--) {
		if (left[i] > right[i]) return 1;
		else if (left[i] < right[i]) return -1;
	}
	return 0;
}

template<const uint N> forceinline
static int vli_less(const u64 *left, const u64 *right) noexcept
{
	u64		carry=0;
	for (uint i = 0; i < N; ++i) {
		u64_subc(left[i], right[i], carry);
	}
	return carry != 0;
}

#ifdef	ommit
/* Computes result = in << c, returning carry. Can modify in place
 * (if result == in). 0 < shift < 32.
 */
template<const uint N, const uint shift> forceinline static
u64 vli_lshift(u64 *result, const u64 *in) noexcept
{
	static_assert(shift < 32, "shift MUST less 32");
	u64 carry = 0;
	for (uint i = 0; i < N; i++) {
		u64 temp = in[i];
		result[i] = (temp << shift) | carry;
		carry = temp >> (64 - shift);
	}
	return carry;
}
#endif

template<const uint N> forceinline static
u64 vli_lshift1(u64 *result, const u64 *in) noexcept
{
	u64 carry = 0;
	for (uint i = 0; i < N; i++) {
		u64 temp = in[i];
		result[i] = (temp << 1) | carry;
		carry = temp >> 63;
	}
	return carry;
}

/* Computes vli = vli >> 1. */
template<const uint N> forceinline
static void vli_rshift1(u64 *vli, u64 carry=0) noexcept
{
	if (carry) carry = 1L << 63;
	for (int i=N - 1; i >= 0; i--) { 
		u64 temp = vli[i];
		vli[i] = (temp >> 1) | carry;
		carry = temp << 63;
	}
}

/* Computes vli = vli >> 1. */
template<const uint N> forceinline
static void vli_rshift1w(u64 *vli, const u64 carry=0) noexcept
{
	for (uint i = 1; i < N; i++) {
		vli[i-1] = vli[i];
	}
	vli[N-1] = carry;
}

/* copy_conditional copies in to out if enabled, set mask to all ones. */
template<const uint N> forceinline
void vli_copy_conditional(u64 *res, const u64 *in, const bool _enabled)
{
	u64	mask = (u64)!_enabled - 1;
	for (uint i = 0; i < N; ++i) {
		const u64 tmp = mask & (in[i] ^ res[i]);
		res[i] ^= tmp;
	}
}

/* Computes result = left + right, returning carry. Can modify in place. */
template<const uint N> forceinline static
bool vli_add(u64 *result, const u64 *left, const u64 *right) noexcept
{
#if	defined(WITH_ASM) && __cplusplus >= 201703L && defined(__clang__)
	if constexpr(N == 4) return vli4_add(result, left, right);
#endif
	u64 carry = 0;
	for (uint i = 0; i < N; ++i) {
		result[i] = u64_addc(left[i], right[i], carry);
	}
	return carry != 0;
}

/* Computes result = left + right, returning carry. Can modify in place. */
template<const uint N> forceinline static
bool vli_add_to(u64 *result, const u64 *right) noexcept
{
#if	defined(WITH_ASM) && __cplusplus >= 201703L && defined(__clang__)
	if constexpr(N == 4) return vli4_add_to(result, right);
#endif
#ifdef	ommit
	if ( unlikely(result == right) ) {
		return vli_lshift1<N>(result, right);
	}
#endif
	u64 carry = 0;
	for (uint i = 0; i < N; ++i) {
		result[i] = u64_addc(result[i], right[i], carry);
	}
	return carry != 0;
}

/* Computes result = left + right, returning carry. Can modify in place. */
template<const uint N> forceinline
static bool vli_uadd_to(u64 *result, u64 right) noexcept
{
	u64 carry=0;
	result[0] = u64_addc(result[0], right, carry);
	for (uint i = 1; i < N; i++) {
		result[i] = u64_addcz(result[i], carry);
	}
	return carry != 0;
}

/**
 * vli_sub() - Subtracts right from left
 *
 * @result:		where to write result
 * @left:		vli
 * @right		vli
 * @ndigits:	N, length of all vlis
 *
 * Note: can modify in-place.
 *
 * Return: carry bit.
 */
template<const uint N> forceinline static
bool vli_sub(u64 *result, const u64 *left, const u64 *right) noexcept
{
#if	defined(WITH_ASM) && __cplusplus >= 201703L && defined(__clang__)
	if constexpr(N == 4) return vli4_sub(result, left, right);
#endif
	u64 borrow = 0;
	for (uint i = 0; i < N; i++) {
		result[i] = u64_subc(left[i], right[i], borrow);
	}
	return borrow != 0;
}

#ifdef	ommit
// carry for left[N]
// result = left[1..N] - right[1..N-1]
// return borrow
template<const uint N> forceinline static bool
vli_subc(u64 *result, const u64 *left, const u64 *right, const u64 carry=0) noexcept
{
	u64 borrow = 0;
	for (uint i = 0; i < N; i++) {
		result[i] = u64_subc(left[i], right[i], borrow);
	}
	u64_subcz(carry, borrow);
	return borrow != 0;
}
#endif

template<const uint N> forceinline static
bool vli_sub_from(u64 *result, const u64 *right) noexcept
{
#if	defined(WITH_ASM) && __cplusplus >= 201703L && defined(__clang__)
	if constexpr(N == 4) return vli4_sub_from(result, right);
#endif
	u64 borrow = 0;
	for (uint i = 0; i < N; i++) {
		result[i] = u64_subc(result[i], right[i], borrow);
	}
	return borrow != 0;
}

/* Computes result = left - right, returning borrow. Can modify in place. */
template<const uint N> forceinline
static bool vli_usub(u64 *result, const u64 *left, u64 right) noexcept
{
	static_assert(N>=4, "N MUST >= 4");
	u64 borrow = 0;
	result[0] = u64_subc(left[0], right, borrow);
	for (uint i = 1; i < N; i++) {
		result[i] = u64_subcz(left[i], borrow);
	}
	return borrow != 0;
}

template<const uint N> forceinline
static bool vli_is_negative(const u64 *vli)  noexcept
{
	return vli_test_bit<N>(vli, N * 64 - 1);
}

forceinline static bool vli_is_even(const u64 *vli) noexcept {
#ifdef	NO_U64ZERO
	return (vli[0] & 1) == 0;
#else
	return u64IsZero(vli[0] & 1);
#endif
}

/* Counts the number of 64-bit "digits" in vli. */
template<const uint N> forceinline
static uint vli_num_digits(const u64 *vli) noexcept
{
	/* Search from the end until we find a non-zero digit.
	 * We do it in reverse because we expect that most digits will
	 * be nonzero.
	 */
	int i;
	for (i = N - 1; i >= 0 && vli[i] == 0; i--);
	return (i + 1);
}

/* Counts the number of bits required for vli. */
template<const uint N> forceinline
static uint vli_num_bits(const u64 *vli) noexcept
{
	auto num_digits = vli_num_digits<N>(vli);
	if (num_digits == 0) return 0;
	auto i = 64 - __builtin_clzl(vli[num_digits - 1]);
	return ((num_digits - 1) * 64 + i);
}

/**
 * vli_from_be64() - Load vli from big-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 BE values
 * @ndigits:		length of both vli and array
 */
template<const uint N> forceinline
static void vli_from_be64(u64 *dest, const void *src) noexcept
{
	const u64 *from = (const u64 *)src;
	for (uint i = 0; i < N; i++)
		dest[i] = be64toh(from[N - 1 - i]);
}

#ifdef	ommit
/**
 * vli_from_le64() - Load vli from little-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 LE values
 * @ndigits:		length of both vli and array
 */
template<uint N> forceinline
__attribute__((optimize("unroll-loops")))
static void vli_from_le64(u64 *dest, const void *src) noexcept
{
	const u64 *from = (const u64 *)src;
	for (uint i = 0; i < N; i++)
		dest[i] = le64toh(from[N - 1 - i]);
}
#endif


template<const uint N> forceinline static
void vli_mult(u64 *result, const u64 *left, const u64 *right) noexcept
{
	uint128_t r01( 0, 0 );
	u64 r2 = 0;
	unsigned int i, k;
	/* Compute each digit of result in sequence, maintaining the
	 * carries.
	 */
	for (k = 0; k < N * 2 - 1; k++) {
		unsigned int min;
		min = (k < N)?0:( (k + 1) - N);

		for (i = min; i <= k && i < N; i++) {
			uint128_t product;
			product.mul_64_64(left[i], right[k - i]);

			r01 += product;
			r2 += (r01.m_high() < product.m_high());
		}

		result[k] = r01.m_low();
		r01 = uint128_t(r01.m_high(), r2);
		r2 = 0;
	}

	result[N * 2 - 1] = r01.m_low();
}

/* Compute product = left * right, for a small right value. */
template<const uint N> forceinline
static void vli_umult(u64 *result, const u64 *left, u64 right) noexcept
{
	uint128_t r01( 0, 0 );
	unsigned int k;

	for (k = 0; k < N; k++) {
		uint128_t product;
		product.mul_64_64(left[k], right);
		// r01 will never overflow
		r01 += product;
		/* no carry */
		result[k] = r01.m_low();
		r01 = uint128_t(r01.m_high(), 0);
	}
	result[N] = r01.m_low();
	for (k = N+1; k < N * 2; k++)
		result[k] = 0;
}

/* Compute product = left * right, reserved for mont_red/mont_mult. */
// result MUST be u64[N+1]
// no clear N+1 .. 2N -1
template<const uint N> forceinline
static void vli_umult2(u64 *result, const u64 *left, u64 right) noexcept
{
	uint128_t r01( 0, 0 );
	unsigned int k;

	for (k = 0; k < N; k++) {
		uint128_t product;
		product.mul_64_64(left[k], right);
		// r01 will never overflow
		r01 += product;
		/* no carry */
		result[k] = r01.m_low();
		r01 = uint128_t(r01.m_high(), 0);
	}
	result[N] = r01.m_low();
}

template<const uint N> forceinline
static void vli_squareOld(u64 *result, const u64 *left) noexcept
{
	uint128_t r01( 0, 0 );
	u64 r2 = 0;
	uint i, k;

	for (k = 0; k < N * 2 - 1; k++) {
		unsigned int min;
		min = (k < N)?0:( (k + 1) - N);

		for (i = min; i <= k && i <= k - i; i++) {
			uint128_t product;
			product.mul_64_64(left[i], left[k - i]);
			if (i < k - i) {
				r2 += product.m_high() >> 63;
				u64 _high = (product.m_high() << 1) | (product.m_low() >> 63);
				u64 _low = product.m_low() << 1;
				product = uint128_t(_low, _high);
			}
			r01 += product;
			r2 += (r01.m_high() < product.m_high());
		}

		result[k] = r01.m_low();
		r01 = uint128_t(r01.m_high(), r2);
		r2 = 0;
	}

	result[N * 2 - 1] = r01.m_low();
}

template<const uint N> forceinline
static void vli_square(u64 *result, const u64 *left) noexcept
{
	vli_clear<N*2>(result);
	for (uint k = 0; k < N - 1; ++k) {
		uint128_t r01( 0, 0 );
		for (uint i = k+1; i < N; ++i) {
			uint128_t product;
			bool		carry, h_carry;
			product.mul_64_64(left[i], left[k]);
			r01 += product;
			h_carry = (r01.m_high() < product.m_high());
			result[i+k] += r01.m_low();
			carry = result[i+k] < r01.m_low();
			r01 = uint128_t(r01.m_high() + carry, h_carry);
		}
		result[k+N] += r01.m_low();
		result[k+N+1] += r01.m_high();
	}
	vli_lshift1<N*2>(result, result);
	u64		dd[N*2];
	for (uint i=0; i < N; ++i) {
		uint128_t product;
		product.mul_64_64(left[i], left[i]);
		dd[i+i] = product.m_low();
		dd[i+i+1] = product.m_high();
	}
	vli_add_to<N*2>(result, dd);
}

/* Computes result = left % mod w/ carry.
 */
template<uint N> forceinline
static void
vli_mod(u64 *result, const u64 *left, const u64 *mod, const bool carry) noexcept
{
	/* result > mod (result = mod + remainder), so subtract mod to
	 * get remainder.
	 */
#if	defined(WITH_ASM) && __cplusplus >= 201703L && defined(__clang__)
	if constexpr(N == 4) {
		vli4_mod(result, left, mod, carry);
	} else
#endif
	//if (carry || vli_cmp<N>(left, mod) >= 0)
	if (carry || !vli_less<N>(left, mod))
	{
		vli_sub<N>(result, left, mod);
	} else vli_set<N>(result, left);
}


#endif	//	__VLI_HPP__
