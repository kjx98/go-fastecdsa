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
#ifndef __MONT_HPP__
#define __MONT_HPP__

#include "vli.hpp"
#include "curve_const.hpp"		// include sm2_p, sm2_p_rr

// SM2 prime optimize
// p is 2^256 - 2^224 - 2^96 + 2^64 -1
forceinline
static void vli_sm2_multP(u64 *result, const u64 u) noexcept
{
	u64	t_low, t_high;
	t_low = u << 32;	// ^192, ^96
	t_high = u >> 32;
	// result = 2^256 + 2^64 - u*2^224 (high 32 bits)
	// r = 2^224 + 2^96 + 1
	u64		carry = 0;
	result[0] = u64_subc(0, u, carry);
	result[1] = u64_subc(u, t_low, carry);
	result[2] = u64_subc(0, t_high, carry);
	result[3] = u64_subc(0, t_low, carry);
	result[4] = u64_subc(u, t_high, carry);
	//result[5] = carry;
}


// p is 2^256 - 2^224 - 2^96 + 2^64 -1
// vli_sm2_multPh return u * (p + 1) >> 64
forceinline
static void vli_sm2_multPh(u64 *result, const u64 u) noexcept
{
	u64	t_low, t_high;
	t_low = u << 32;	// ^192
	t_high = u >> 32;
	u64		carry = 0;
	result[0] = u64_subc(u, t_low, carry);
	result[1] = u64_subc(0, t_high, carry);
	result[2] = u64_subc(0, t_low, carry);
	result[3] = u64_subc(u, t_high, carry);
}

// u * 2^256 mod sm2 prime
// p is 2^256 - 2^224 - 2^96 + 2^64 -1
// R(2^256) - p = 2^224 + 2^96 - 2^64 + 1
// u < 2^32
// used for carry_reduce
forceinline
static void vli_sm2_multR(u64 *result, const u64 uv) noexcept
{
	u64 u = uv & 0xffffffff;
	result[0] = u;
	result[1] = (u << 32) - u;
	result[2] = 0;
	result[3] = u << 32;
	return;
}


// secp256k1 prime optimize
// p is 2^256 - 2^32 - 0x3d1 = 2^256 - 0x1000003d1
forceinline
static void vli_btc_multP(u64 *result, const u64 u) noexcept
{
	uint128_t	pd;
	pd.mul_64_64(0x1000003d1, u);
	u64		carry = 0;
	result[0] = u64_subc(0, pd.m_low(), carry);
	result[1] = u64_subc(0, pd.m_high(), carry);
	result[2] = u64_subc(0, 0, carry);
	result[3] = u64_subc(0, 0, carry);
	result[4] = u64_subc(u, 0, carry);
}

/* result > mod (result = mod + remainder), so subtract mod to
 * get remainder.
 */
static forceinline void
sm2p_mod(u64 *res, const u64 *left, const bool carry) noexcept
{
#ifdef	__x86_64__
	asm volatile("movq (%%rsi), %%r8\n"	// mov r8/9/10/11, right
				"movq 8(%%rsi), %%r9\n"
				"movq 16(%%rsi), %%r10\n"
				"movq 24(%%rsi), %%r11\n"
				"movq %%r8, %%r12\n"
				"movq %%r9, %%r13\n"
				"movq %%r10, %%r14\n"
				"movq %%r11, %%r15\n"
				"subq $-1, %%r12\n"
				"sbbq %[mod0], %%r13\n"
				"sbbq $-1, %%r14\n"
				"sbbq %[mod1], %%r15\n"
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
				: "S"(left), "D"(res), "a" ((u64)carry), [mod0] "m" (sm2_p[1]),
				[mod1] "m" (sm2_p[3])
				: "%r8", "%r9", "%r10", "%r11" , "%r12", "%r13", "%r14", "%r15", "cc", "memory");
#elif	defined(__aarch64__)
	asm volatile( "mov x9, -1\n"
				"mov x10, %3\n"
				"mov x11, -1\n"
				"mov x12, %4\n"
				"ldp x4, x5, [%2]\n"
				"ldp x6, x7, [%2, 16]\n"
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
		: "r" ((u64)carry), "r" (res), "r" (left), "rm" (sm_p[1]),
		"rm" (sm2_p[3])
		: "%x4", "%x5", "%x6", "%x7", "%x9", "%x10", "%x11", "%x12", "cc", "memory");
#else
	vli_mod<4>(res, left, sm2_p, carry);
#endif
}

forceinline static void
sm2p_reduction(u64 *result, const u64 *y, const bool isProd=false) noexcept
{
#ifdef	__x86_64__
	register u64 res0 asm("r12");
	register u64 res1 asm("r13");
	register u64 res2 asm("r8");
	register u64 res3 asm("r9");
	asm volatile("MOVQ (8*0)(%%rsi), %%r8\n"
		"MOVQ (8*1)(%%rsi), %%r9\n"
		"MOVQ (8*2)(%%rsi), %%r10\n"
		"MOVQ (8*3)(%%rsi), %%r11\n"
		"XORQ %%r12, %%r12\n"

	// Only reduce, no multiplications are needed
	// First stage
"MOVQ %%r8, %%rax\n"
"MOVQ %%r8, %%r15\n"
"SHLQ $32, %%r8\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r9\n"
"ADCQ $0, %%r10\n"
"ADCQ $0, %%r11\n"
"ADCQ %%rax, %%r12\n"
"subq %%r8, %%r9\n"
"sbbq %%r15, %%r10\n"
"sbbq %%r8, %%r11\n"
"sbbq %%r15, %%r12\n"
"XORQ %%r13, %%r13\n"
	// Second stage
"MOVQ %%r9, %%rax\n"
"MOVQ %%r9, %%r15\n"
"SHLQ $32, %%r9\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r10\n"
"adcq $0, %%r11\n"
"adcq $0, %%r12\n"
"adcq %%rax, %%r13\n"
"subq %%r9, %%r10\n"
"sbbq %%r15, %%r11\n"
"sbbq %%r9, %%r12\n"
"sbbq %%r15, %%r13\n"
"XORQ %%r8, %%r8\n"
	// Third stage
"MOVQ %%r10, %%rax\n"
"MOVQ %%r10, %%r15\n"
"SHLQ $32, %%r10\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r11\n"
"adcq $0, %%r12\n"
"adcq $0, %%r13\n"
"adcq %%rax, %%r8\n"
"subq %%r10, %%r11\n"
"sbbq %%r15, %%r12\n"
"sbbq %%r10, %%r13\n"
"sbbq %%r15, %%r8\n"
"XORQ %%r9, %%r9\n"
	// Last stage
"MOVQ %%r11, %%rax\n"
"MOVQ %%r11, %%r15\n"
"SHLQ $32, %%r11\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r12\n"
"adcq $0, %%r13\n"
"adcq $0, %%r8\n"
"adcq %%rax, %%r9\n"
"subq %%r11, %%r12\n"
"sbbq %%r15, %%r13\n"
"sbbq %%r11, %%r8\n"
"sbbq %%r15, %%r9\n"
//"movq $0, %%rax\n"
//"sbbq $0, %%rax\n"
				: "=r" (res0), "=r" (res1), "=r" (res2), "=r" (res3)
				: "S" (y)
				: "rax", "r10", "r11", "r14", "r15", "cc", "memory");

	u64	carry = 0;
	//if (carry != 0) carry = 1;
	// add high 256 bits
	if ( unlikely(isProd) )
	{
		u64	cc=0;
		res0 = u64_addc(res0, y[4], cc);
		res1 = u64_addc(res1, y[5], cc);
		res2 = u64_addc(res2, y[6], cc);
		res3 = u64_addc(res3, y[7], cc);
		carry += cc;
	}

#ifdef	__x86_64__
	// mod prime
	asm volatile(
		"MOVQ %%r12, %%r10\n"
		"MOVQ %%r13, %%r11\n"
		"MOVQ %%r8, %%r14\n"
		"MOVQ %%r9, %%r15\n"

		"SUBQ $-1, %%r12\n"
		"SBBQ %[pr1], %%r13\n"
		"SBBQ $-1, %%r8\n"
		"SBBQ %[pr3], %%r9\n"
		"sbbq $0, %%rax\n"

		"CMOVCQ %%r10, %%r12\n"
		"CMOVCQ %%r11, %%r13\n"
		"CMOVCQ %%r14, %%r8\n"
		"CMOVCQ %%r15, %%r9\n"

		"MOVQ %%r12, (8*0)(%%rdi)\n"
		"MOVQ %%r13, (8*1)(%%rdi)\n"
		"MOVQ %%r8, (8*2)(%%rdi)\n"
		"MOVQ %%r9, (8*3)(%%rdi)\n"
			: 			// acc4/5/0/1
			: "D" (result), "a" (carry), [pr1] "m" (sm2_p[1]),
			[pr3] "m" (sm2_p[3]), "r" (res0), "r" (res1), "r" (res2), "r" (res3)
			: "r10", "r11", "r14", "r15", "cc", "memory");
#else
	result[0] = res0;
	result[1] = res1;
	result[2] = res2;
	result[3] = res3;
	// sm2p_mod
	sm2p_mod(result, result, carry != 0);
#endif
#else
	u64	r[4];
	vli_set<4>(r, y);
	u64 carry = 0;
	for (uint i=0; i < 4; i++) {
		u64	u = r[0];
		u64	t_low, t_high;
		t_low = u << 32;	// ^192
		t_high = u >> 32;
		vli_rshift1w<4>(r, carry);	
		u64 cc = 0;
		r[0] = u64_addc(r[0], u, cc);
		r[1] = u64_addcz(r[1], cc);
		r[2] = u64_addcz(r[2], cc);
		r[3] = u64_addc(r[3], u, cc);
		carry = cc;
		cc = 0;
		r[0] = u64_subc(r[0], t_low, cc);
		r[1] = u64_subc(r[1], t_high, cc);
		r[2] = u64_subc(r[2], t_low, cc);
		r[3] = u64_subc(r[3], t_high, cc);
		carry -= cc;
	}
	// add high 256 bits
	if ( unlikely(isProd) )
	{
		u64	cc=0;
		r[0] = u64_addc(r[0], y[4], cc);
		r[1] = u64_addc(r[1], y[5], cc);
		r[2] = u64_addc(r[2], y[6], cc);
		r[3] = u64_addc(r[3], y[7], cc);
		carry += cc;
	}
	sm2p_mod(result, r, carry != 0);
#endif
}


template<const uint N> forceinline
static void
#ifdef	WITH_C2GO_1
vli_mont_reduction(u64 *result, const u64 *y, const u64 *prime,
			const u64 k0, u64 *buff) noexcept
#else
vli_mont_reduction(u64 *result, const u64 *y, const u64 *prime,
			const u64 k0) noexcept
#endif
{
#ifdef	WITH_C2GO_1
	u64	*s = buff;
	u64	*r = s + N * 2;
#else
	u64	s[N + 1];
	u64	r[N + 1];
#endif
	vli_set<N>(r, y);
	r[N] = 0;
	for (uint i=0; i < N; i++) {
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);	
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}

template<const uint N> forceinline
static void
#ifdef	WITH_C2GO_1
vli_mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 k0, u64 *buff) noexcept
#else
vli_mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 k0) noexcept
#endif
{
#ifdef	WITH_C2GO_1
	u64	*s = buff;
	u64	*r = s + N * 2;
#else
	u64	s[N + 1];
	u64	r[N + 1];
#endif
	vli_clear<N + 1>(r);
	for (uint i=0; i < N;i++) {
		vli_umult2<N>(s, x, y[i]);
		vli_add_to<N + 1>(r, s);
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);	
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}

template<const uint N> forceinline
static void
#ifdef	WITH_C2GO_1
vli_mont_sqr(u64 *result, const u64 *x, const u64 *prime, const u64 k0, u64 *buff) noexcept
#else
vli_mont_sqr(u64 *result, const u64 *x, const u64 *prime, const u64 k0) noexcept
#endif
{
#ifdef	WITH_C2GO_1
	u64	*s = buff;
	u64	*r = s + N * 2;
#else
	u64	s[N + 1];
	u64	r[N + 1];
#endif
	vli_clear<N + 1>(r);
	for (uint i=0; i < N;i++) {
		vli_umult2<N>(s, x, x[i]);
		vli_add_to<N + 1>(r, s);
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);	
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}


template<const uint N, const u64 k0> forceinline
static void
mont_reduction(u64 *result, const u64 *y, const u64 *prime) noexcept
{
	u64	s[N + 1];
	u64	r[N + 1];
	vli_set<N>(r, y);
	r[N] = 0;
	for (uint i=0; i < N; i++) {
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);	
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}

template<const uint N, const u64 k0> forceinline
static void
mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime) noexcept
{
	u64	s[N + 1];
	u64	r[N + 1];
	vli_clear<N + 1>(r);
	for (uint i=0; i < N;i++) {
		vli_umult2<N>(s, x, y[i]);
		vli_add_to<N + 1>(r, s);
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);	
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}

forceinline static void
sm2p_mult(u64 *result, const u64 *x, const u64 *y) noexcept
{
#ifdef	__x86_64__
	asm volatile(
	// x * y[0]
"MOVQ (8*0)(%0), %%r14\n"

"MOVQ (8*0)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"MOVQ %%RAX, %%r8\n"
"MOVQ %%RDX, %%r9\n"

"MOVQ (8*1)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%RAX, %%r9\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r10\n"

"MOVQ (8*2)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%RAX, %%r10\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r11\n"

"MOVQ (8*3)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%RAX, %%r11\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r12\n"
"XORQ %%r13, %%r13\n"
	// First reduction step
"MOVQ %%r8, %%RAX\n"
"MOVQ %%r8, %%r15\n"
"SHLQ $32, %%r8\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r9\n"
"ADCQ $0, %%r10\n"
"ADCQ $0, %%r11\n"
"ADCQ %%rax, %%r12\n"
"ADCQ $0, %%r13\n"
"SUBQ %%r8, %%r9\n"
"SBBQ %%r15, %%r10\n"
"SBBQ %%r8, %%r11\n"
"SBBQ %%R15, %%r12\n"
"SBBQ $0, %%r13\n"
"XORQ %%r8, %%r8\n"
	// x * y[1]
"MOVQ (8*1)(%0), %%r14\n"

"MOVQ (8*0)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%RAX, %%r9\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r15\n"

"MOVQ (8*1)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%r15, %%r10\n"
"ADCQ $0, %%RDX\n"
"ADDQ %%RAX, %%r10\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r15\n"

"MOVQ (8*2)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%r15, %%r11\n"
"ADCQ $0, %%RDX\n"
"ADDQ %%RAX, %%r11\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r15\n"

"MOVQ (8*3)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%r15, %%r12\n"
"ADCQ $0, %%RDX\n"
"ADDQ %%RAX, %%r12\n"
"ADCQ %%RDX, %%r13\n"
"ADCQ $0, %%r8\n"
	// Second reduction step
"MOVQ %%r9, %%RAX\n"
"MOVQ %%r9, %%r15\n"
"SHLQ $32, %%r9\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r10\n"
"ADCQ $0, %%r11\n"
"ADCQ $0, %%r12\n"
"ADCQ %%rax, %%r13\n"
"ADCQ $0, %%r8\n"
"SUBQ %%r9, %%r10\n"
"SBBQ %%r15, %%r11\n"
"SBBQ %%r9, %%r12\n"
"SBBQ %%r15, %%r13\n"
"SBBQ $0, %%r8\n"
"XORQ %%r9, %%r9\n"
	// x * y[2]
"MOVQ (8*2)(%0), %%r14\n"

"MOVQ (8*0)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%RAX, %%r10\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r15\n"

"MOVQ (8*1)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%r15, %%r11\n"
"ADCQ $0, %%RDX\n"
"ADDQ %%RAX, %%r11\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r15\n"

"MOVQ (8*2)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%r15, %%r12\n"
"ADCQ $0, %%RDX\n"
"ADDQ %%RAX, %%r12\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r15\n"

"MOVQ (8*3)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%r15, %%r13\n"
"ADCQ $0, %%RDX\n"
"ADDQ %%RAX, %%r13\n"
"ADCQ %%RDX, %%r8\n"
"ADCQ $0, %%r9\n"
	// Third reduction step
"MOVQ %%r10, %%RAX\n"
"MOVQ %%r10, %%r15\n"
"SHLQ $32, %%r10\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r11\n"
"ADCQ $0, %%r12\n"
"ADCQ $0, %%r13\n"
"ADCQ %%rax, %%r8\n"
"ADCQ $0, %%r9\n"
"SUBQ %%r10, %%r11\n"
"SBBQ %%r15, %%r12\n"
"SBBQ %%r10, %%r13\n"
"SBBQ %%r15, %%r8\n"
"SBBQ $0, %%r9\n"
"XORQ %%r10, %%r10\n"
	// x * y[3]
"MOVQ (8*3)(%0), %%r14\n"

"MOVQ (8*0)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%RAX, %%r11\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r15\n"

"MOVQ (8*1)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%r15, %%r12\n"
"ADCQ $0, %%RDX\n"
"ADDQ %%RAX, %%r12\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r15\n"

"MOVQ (8*2)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%r15, %%r13\n"
"ADCQ $0, %%RDX\n"
"ADDQ %%RAX, %%r13\n"
"ADCQ $0, %%RDX\n"
"MOVQ %%RDX, %%r15\n"

"MOVQ (8*3)(%%rsi), %%RAX\n"
"MULQ %%r14\n"
"ADDQ %%r15, %%r8\n"
"ADCQ $0, %%RDX\n"
"ADDQ %%RAX, %%r8\n"
"ADCQ %%RDX, %%r9\n"
"ADCQ $0, %%r10\n"
	// Last reduction step
"MOVQ %%r11, %%RAX\n"
"MOVQ %%r11, %%r15\n"
"SHLQ $32, %%r11\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r12\n"
"ADCQ $0, %%r13\n"
"ADCQ $0, %%r8\n"
"ADCQ %%rax, %%r9\n"
"ADCQ $0, %%r10\n"
"SUBQ %%r11, %%r12\n"
"SBBQ %%r15, %%r13\n"
"SBBQ %%r11, %%r8\n"
"SBBQ %%r15, %%r9\n"
"SBB $0, %%r10\n"
	// Copy result [255:0]
"MOVQ %%r12, %%rax\n"
"MOVQ %%r13, %%r11\n"
"MOVQ %%r8, %%r14\n"
"MOVQ %%r9, %%r15\n"
	// Subtract sm2_p
"SUBQ $-1, %%r12\n"
"SBBQ %[pr1] ,%%r13\n"
"SBBQ $-1, %%r8\n"
"SBBQ %[pr3], %%r9\n"
"SBBQ $0, %%r10\n"

"CMOVCQ %%rax, %%r12\n"
"CMOVCQ %%r11, %%r13\n"
"CMOVCQ %%r14, %%r8\n"
"CMOVCQ %%r15, %%r9\n"

"MOVQ %%r12, (8*0)(%%rdi)\n"
"MOVQ %%r13, (8*1)(%%rdi)\n"
"MOVQ %%r8, (8*2)(%%rdi)\n"
"MOVQ %%r9, (8*3)(%%rdi)\n"
			: 			// acc4/5/0/1
			: "r" (y), "D" (result), "S" (x), [pr1] "m" (sm2_p[1]),
			[pr3] "m" (sm2_p[3])
			: "rax", "rdx", "r8", "r9", "r12", "r13", "r10", "r11", "r14",
			"r15", "cc", "memory");
#else
	u64	s[4+1];
	u64	r[4+1];
	vli_clear<4 + 1>(r);
	for (uint i=0; i < 4;i++) {
		vli_umult2<4>(s, x, y[i]);
		vli_add_to<4 + 1>(r, s);
		vli_sm2_multPh(s, r[0]);
		vli_rshift1w<4>(r, r[4]);
		r[4] = vli_add_to<4>(r, s);
	}
#if	__cplusplus >= 201703L && defined(WITH_ASM)
	sm2p_mod(result, r, r[4] != 0);
#else
	vli_mod<4>(result, r, sm2_p, r[4] != 0);
#endif
#endif
}

template<const uint N, const u64 k0> forceinline
static void
mont_sqr(u64 *result, const u64 *x, const u64 *prime) noexcept
{
	u64	s[N + 1];
	u64	r[N + 1];
	vli_clear<N + 1>(r);
	for (uint i=0; i < N;i++) {
		vli_umult2<N>(s, x, x[i]);
		vli_add_to<N + 1>(r, s);
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);	
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}

template<const uint N, const u64 k0> forceinline
static void
mont_sqrN(u64 *result, const u64 *x, const u64 *prime) noexcept
{
	u64	r[N * 2];
	u64	s[N + 1];
	vli_square<N>(r, x);
	vli_set<N>(result, r + N);
	r[N] = 0;
	for (uint i=0; i < N; i++) {
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);	
	}
#if	__cplusplus >= 201703L && defined(WITH_ASM)
	if constexpr(N==4) {
		r[N] += vli4_add_to(r, result);
		vli4_mod(result, r, prime, r[N] != 0);
	} else
#endif
	{
		r[N] += vli_add_to<N>(r, result);
		vli_mod<N>(result, r, prime, r[N] != 0);
	}
}

forceinline static void
sm2p_sqrN(u64 *result, const u64 *x) noexcept
{
#ifdef	__x86_64__
	asm volatile(
	// y[1:] * y[0]
"MOVQ (8*0)(%%rsi), %%r14\n"

"MOVQ (8*1)(%%rsi), %%rax\n"
"MULQ %%r14\n"
"MOVQ %%rax, %%r9\n"
"MOVQ %%rdx, %%r10\n"

"MOVQ (8*2)(%%rsi), %%rax\n"
"MULQ %%r14\n"
"ADDQ %%rax, %%r10\n"
"ADCQ $0, %%rdx\n"
"MOVQ %%rdx, %%r11\n"

"MOVQ (8*3)(%%rsi), %%rax\n"
"MULQ %%r14\n"
"ADDQ %%rax, %%r11\n"
"ADCQ $0, %%rdx\n"
"MOVQ %%rdx, %%r12\n"
	// y[2:] * y[1]
"MOVQ (8*1)(%%rsi), %%r14\n"

"MOVQ (8*2)(%%rsi), %%rax\n"
"MULQ %%r14\n"
"ADDQ %%rax, %%r11\n"
"ADCQ $0, %%rdx\n"
"MOVQ %%rdx, %%r15\n"

"MOVQ (8*3)(%%rsi), %%rax\n"
"MULQ %%r14\n"
"ADDQ %%r15, %%r12\n"
"ADCQ $0, %%rdx\n"
"ADDQ %%rax, %%r12\n"
"ADCQ $0, %%rdx\n"
"MOVQ %%rdx, %%r13\n"
	// y[3] * y[2]
"MOVQ (8*2)(%%rsi), %%r14\n"

"MOVQ (8*3)(%%rsi), %%rax\n"
"MULQ %%r14\n"
"ADDQ %%rax, %%r13\n"
"ADCQ $0, %%rdx\n"
"MOVQ %%rdx, %%rcx\n"
"XORQ %%r15, %%r15\n"
	// *2
"ADDQ %%r9, %%r9\n"
"ADCQ %%r10, %%r10\n"
"ADCQ %%r11, %%r11\n"
"ADCQ %%r12, %%r12\n"
"ADCQ %%r13, %%r13\n"
"ADCQ %%rcx, %%rcx\n"
"ADCQ $0, %%r15\n"
	// Missing products
"MOVQ (8*0)(%%rsi), %%rax\n"
"MULQ %%rax\n"
"MOVQ %%rax, %%r8\n"
"MOVQ %%rdx, %%r14\n"

"MOVQ (8*1)(%%rsi), %%rax\n"
"MULQ %%rax\n"
"ADDQ %%r14, %%r9\n"
"ADCQ %%rax, %%r10\n"
"ADCQ $0, %%rdx\n"
"MOVQ %%rdx, %%r14\n"

"MOVQ (8*2)(%%rsi), %%rax\n"
"MULQ %%rax\n"
"ADDQ %%r14, %%r11\n"
"ADCQ %%rax, %%r12\n"
"ADCQ $0, %%rdx\n"
"MOVQ %%rdx, %%r14\n"

"MOVQ (8*3)(%%rsi), %%rax\n"
"MULQ %%rax\n"
"ADDQ %%r14, %%r13\n"
"ADCQ %%rax, %%rcx\n"
"ADCQ %%rdx, %%r15\n"
"MOVQ %%r15, %%rbx\n"
	// First reduction step
"MOVQ %%r8, %%rax\n"
"MOVQ %%r8, %%rdx\n"
"MOVQ %%r8, %%r15\n"
"SHLQ $32, %%r8\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r9\n"
"ADCQ $0, %%r10\n"
"ADCQ $0, %%r11\n"
"ADCQ $0, %%rdx\n"
"SUBQ %%r8, %%r9\n"
"SBBQ %%r15, %%r10\n"
"SBBQ %%r8, %%r11\n"
"SBBQ %%R15, %%rdx\n"
"MOVQ %%rdx, %%r8\n"
	// Second reduction step
"MOVQ %%r9, %%rax\n"
"MOVQ %%r9, %%rdx\n"
"MOVQ %%r9, %%r15\n"
"SHLQ $32, %%r9\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r10\n"
"ADCQ $0, %%r11\n"
"ADCQ $0, %%r8\n"
"ADCQ $0, %%rdx\n"
"SUBQ %%r9, %%r10\n"
"SBBQ %%r15, %%r11\n"
"SBBQ %%r9, %%r8\n"
"SBBQ %%r15, %%rdx\n"
"MOVQ %%rdx, %%r9\n"
	// Third reduction step
"MOVQ %%r10, %%rax\n"
"MOVQ %%r10, %%rdx\n"
"MOVQ %%r10, %%r15\n"
"SHLQ $32, %%r10\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r11\n"
"ADCQ $0, %%r8\n"
"ADCQ $0, %%r9\n"
"ADCQ $0, %%rdx\n"
"SUBQ %%r10, %%r11\n"
"SBBQ %%r15, %%r8\n"
"SBBQ %%r10, %%r9\n"
"SBBQ %%r15, %%rdx\n"
"MOVQ %%rdx, %%r10\n"
	// Last reduction step
"XORQ %%r14, %%r14\n"
"MOVQ %%r11, %%rax\n"
"MOVQ %%r11, %%rdx\n"
"MOVQ %%r11, %%r15\n"
"SHLQ $32, %%r11\n"
"SHRQ $32, %%r15\n"
"ADDQ %%rax, %%r8\n"
"ADCQ $0, %%r9\n"
"ADCQ $0, %%r10\n"
"ADCQ $0, %%rdx\n"
"SUBQ %%r11, %%r8\n"
"SBBQ %%r15, %%r9\n"
"SBBQ %%r11, %%r10\n"
"SBBQ %%r15, %%rdx\n"
"MOVQ %%rdx, %%r11\n"
	// Add bits [511:256] of the sqr result
"ADCQ %%r12, %%r8\n"
"ADCQ %%r13, %%r9\n"
"ADCQ %%rcx, %%r10\n"
"ADCQ %%rbx, %%r11\n"
"ADCQ $0, %%r14\n"

"MOVQ %%r8, %%r12\n"
"MOVQ %%r9, %%r13\n"
"MOVQ %%r10, %%rcx\n"
"MOVQ %%r11, %%r15\n"
	// Subtract sm2_p
"SUBQ $-1, %%r8\n"
"SBBQ %[pr1] ,%%r9\n"
"SBBQ $-1, %%r10\n"
"SBBQ %[pr3], %%r11\n"
"SBBQ $0, %%r14\n"

"CMOVCQ %%r12, %%r8\n"
"CMOVCQ %%r13, %%r9\n"
"CMOVCQ %%rcx, %%r10\n"
"CMOVCQ %%r15, %%r11\n"

"MOVQ %%r8, (8*0)(%%rdi)\n"
"MOVQ %%r9, (8*1)(%%rdi)\n"
"MOVQ %%r10, (8*2)(%%rdi)\n"
"MOVQ %%r11, (8*3)(%%rdi)\n"
		: 			// acc4/5/0/1
		: "D" (result), "S" (x), [pr1] "m" (sm2_p[1]), [pr3] "m" (sm2_p[3])
		: "rax", "rbx", "rcx", "rdx", "r8", "r9", "r12", "r13", "r10", "r11",
		"r14",
		"r15", "cc", "memory");
#else
	u64	r[8];
	vli_squareN<4>(r, x);
	sm2p_reduction(result, r, true);
#endif
}

#endif	//	__MONT_HPP__
