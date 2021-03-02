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

// u * 2^256 mod sm2 prime
// p is 2^256 - 2^224 - 2^96 + 2^64 -1
// R(2^256) = 2^224 + 2^96 - 2^64 + 1 Mod p = 2^32 * 2^192 + (2^32-1) * 2^64 + 1
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
	asm volatile(
		"movq (%%rsi), %%r8\n"	// mov r8/9/10/11, right
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
	asm volatile(
		//"mov x9, #-1\n"
		//"mov x10, %3\n"
		"mov x11, #-1\n"
		//"mov x12, %4\n"
		"ldp x4, x5, [%2]\n"
		"ldp x6, x7, [%2, 16]\n"
		"subs x9, x4, #-1\n"
		"sbcs x10, x5, %3\n"
		"sbcs x11, x6, x11\n"
		"sbcs x12, x7, %4\n"
		"sbcs %0, %0, xzr\n"
		"csel x4, x4, x9, cc\n"
		"csel x5, x5, x10, cc\n"
		"csel x6, x6 , x11, cc\n"
		"csel x7, x7, x12, cc\n"
		"stp x4, x5, [%1]\n"
		"stp x6, x7, [%1, 16]\n"
		//"adc %0, xzr, xzr\n"
		:
		: "r" ((u64)carry), "r" (res), "r" (left), "r" (sm2_p[1]), "r" (sm2_p[3])
		: "%x4", "%x5", "%x6", "%x7", "%x9", "%x10", "%x11", "%x12", "cc", "memory");
#else
	vli_mod<4>(res, left, sm2_p, carry);
#endif
}

#ifdef	WITH_SM2P_MOD2
/* result > mod (result = mod + remainder), so subtract mod to
 * get remainder.
 */
static forceinline void
sm2p_mod2(u64 *res, const u64 r0, const u64 r1, const u64 r2, const u64 r3,
		const bool carry) noexcept
{
#ifdef	__x86_64__
	register u64 re0 asm("r8") = r0;
	register u64 re1 asm("r9") = r1;
	register u64 re2 asm("r10") = r2;
	register u64 re3 asm("r11") = r3;
	asm volatile(
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
		: "D"(res), "a" ((u64)carry), "r" (re0), "r" (re1), "r" (re2),
			"r" (re3), [mod0] "m" (sm2_p[1]), [mod1] "m" (sm2_p[3])
		: "%r12", "%r13", "%r14", "%r15", "cc", "memory");
#elif	defined(__aarch64__)
	register u64 re0 asm("x4") = r0;
	register u64 re1 asm("x5") = r1;
	register u64 re2 asm("x6") = r2;
	register u64 re3 asm("x7") = r3;
	asm volatile(
		//"mov x9, #-1\n"
		//"mov x10, %2\n"
		"mov x11, #-1\n"
		//"mov x12, %3\n"
		"subs x9, x4, #-1\n"
		"sbcs x10, x5, %2\n"
		"sbcs x11, x6, x11\n"
		"sbcs x12, x7, %3\n"
		"sbcs %0, %0, xzr\n"
		"csel x4, x4, x9, cc\n"
		"csel x5, x5, x10, cc\n"
		"csel x6, x6 , x11, cc\n"
		"csel x7, x7, x12, cc\n"
		"stp x4, x5, [%1]\n"
		"stp x6, x7, [%1, 16]\n"
		//"adc %0, xzr, xzr\n"
		:
		: "r" ((u64)carry), "r" (res), "r" (sm2_p[1]), "r" (sm2_p[3]),
		"r" (re0), "r" (re1), "r" (re2), "r" (re3)
		: "%x9", "%x10", "%x11", "%x12", "cc", "memory");
#else
	static_assert(false, "only x86_64 and aarch64 supported");
#endif
}
#endif


forceinline static void
sm2p_reductionStep(u64& r0, u64& r1, u64& r2, u64& r3, u64& carry)
{
	u64	u = r0;
	r0 = carry;		// rshift1w
	u64 cc = 0;
	r1 = u64_addc(r1, u, cc);
	r2 = u64_addcz(r2, cc);
	r3 = u64_addcz(r3, cc);
	r0 = u64_addc(r0, u, cc);
	carry = cc;
	u64 t_low = u << 32;	// ^192
	u64 t_high = u >> 32;
	cc = 0;
	r1 = u64_subc(r1, t_low, cc);
	r2 = u64_subc(r2, t_high, cc);
	r3 = u64_subc(r3, t_low, cc);
	r0 = u64_subc(r0, t_high, cc);
	carry -= cc;
}

forceinline static void
sm2p_reductionN(u64 *result, const u64 *y, const bool isProd=false) noexcept
{
	u64	r0, r1, r2, r3;
	vli4_load(y, r0, r1, r2, r3);
	u64 carry = 0;
	sm2p_reductionStep(r0, r1, r2, r3, carry);
	sm2p_reductionStep(r1, r2, r3, r0, carry);
	sm2p_reductionStep(r2, r3, r0, r1, carry);
	sm2p_reductionStep(r3, r0, r1, r2, carry);

	// add high 256 bits
	if ( unlikely(isProd) )
	{
		u64	cc=0;
		r0 = u64_addc(r0, y[4], cc);
		r1 = u64_addc(r1, y[5], cc);
		r2 = u64_addc(r2, y[6], cc);
		r3 = u64_addc(r3, y[7], cc);
		carry += cc;
	}
	// sm2p_mod
#ifdef	WITH_SM2P_MOD2
	sm2p_mod2(result, r0, r1, r2, r3, carry);
#else
	{
		u64	cc=0;
		u64 s0 = u64_subc(r0, sm2_p[0], cc);
		u64 s1 = u64_subc(r1, sm2_p[1], cc);
		u64 s2 = u64_subc(r2, sm2_p[2], cc);
		u64 s3 = u64_subc(r3, sm2_p[3], cc);
		u64_subcz(carry, cc);
		if (cc != 0) vli4_save(result, r0, r1, r2, r3); else
			vli4_save(result, s0, s1, s2, s3);
	}
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
	asm volatile(
		"MOVQ (8*0)(%%rsi), %%r8\n"
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
		: "=r" (res0), "=r" (res1), "=r" (res2), "=r" (res3)
		: "S" (y)
		: "rax", "r10", "r11", "r15", "cc", "memory");

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
#elif	defined(__aarch64__)
	register u64 res0 asm("x4") = y[0];
	register u64 res1 asm("x5") = y[1];
	register u64 res2 asm("x6") = y[2];
	register u64 res3 asm("x7") = y[3];
	register u64 carry asm("x0") = 0;
	asm volatile(
	// Only reduce, no multiplications are needed
	// First reduction step
		"MOV x10, x4\n"
		"MOV x4, xzr\n"
		"ADDS x5, x5, x10\n"
		"ADCS x6, x6, XZR\n"
		"ADCS x7, x7, XZR\n"
		"ADCS x4, x4, x10\n"
		"ADCS x0, XZR, XZR\n"
		"LSR x9, x10, 32\n"
		"LSL x10, x10, 32\n"
		"SUBS x5, x5, x10\n"
		"SBCS x6, x6, x9\n"
		"SBCS x7, x7, x10\n"
		"SBCS x4, x4, x9\n"
		"SBCS x0, x0, xzr\n"
	// Second reduction step
		"MOV x10, x5\n"
		"MOV x5, x0\n"
		"ADDS x6, x6, x10\n"
		"ADCS x7, x7, XZR\n"
		"ADCS x4, x4, XZR\n"
		"ADCS x5, x5, x10\n"
		"ADCS x0, XZR, XZR\n"
		"LSR x9, x10, 32\n"
		"LSL x10, x10, 32\n"
		"SUBS x6, x6, x10\n"
		"SBCS x7, x7, x9\n"
		"SBCS x4, x4, x10\n"
		"SBCS x5, x5, x9\n"
		"SBCS x0, x0, xzr\n"
	// Third reduction step
		"MOV x10, x6\n"
		"MOV x6, x0\n"
		"ADDS x7, x7, x10\n"
		"ADCS x4, x4, XZR\n"
		"ADCS x5, x5, XZR\n"
		"ADCS x6, x6, x10\n"
		"ADCS x0, XZR, XZR\n"
		"LSR x9, x10, 32\n"
		"LSL x10, x10, 32\n"
		"SUBS x7, x7, x10\n"
		"SBCS x4, x4, x9\n"
		"SBCS x5, x5, x10\n"
		"SBCS x6, x6, x9\n"
		"SBCS x0, x0, xzr\n"
	// Last reduction step
		"MOV x10, x7\n"
		"MOV x7, x0\n"
		"ADDS x4, x4, x10\n"
		"ADCS x5, x5, XZR\n"
		"ADCS x6, x6, XZR\n"
		"ADCS x7, x7, x10\n"
		"ADCS x0, XZR, XZR\n"
		"LSR x9, x10, 32\n"
		"LSL x10, x10, 32\n"
		"SUBS x4, x4, x10\n"
		"SBCS x5, x5, x9\n"
		"SBCS x6, x6, x10\n"
		"SBCS x7, x7, x9\n"
		"SBCS x0, x0, xzr\n"
		: "+r" (res0), "+r" (res1), "+r" (res2), "+r" (res3), "+r" (carry)
		:
		: "%x9", "%x10", "%x11", "%x12", "cc", "memory");

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

	// mod prime
	asm volatile(
		//"MOV	x9, #-1\n"
		"MOV	x11, #-1\n"

		"SUBS	x9, x4, #-1\n"
		"SBCS	x10, x5, %1\n"
		"SBCS	x11, x6, x11\n"
		"SBCS	x12, x7, %2\n"
		"SBCS	x0, x0, XZR\n"

		"CSEL	x4, x4, x9, cc\n"
		"CSEL	x5, x5, x10, cc\n"
		"CSEL	x6, x6, x11, cc\n"
		"CSEL	x7, x7, x12, cc\n"

		"STP	x4, x5, [%0]\n"
		"STP	x6, x7, [%0, 16]\n"
		:
		: "r" (result), "r" (sm2_p[1]), "r" (sm2_p[3]), "r" (carry),
			"r" (res0), "r" (res1), "r" (res2), "r" (res3)
		: "%x9", "%x10", "%x11", "%x12", "cc", "memory");
#else
	sm2p_reductionN(result, y, isProd);
#endif
}


template<const uint N> forceinline
static void vli_mont_reduction(u64 *result, const u64 *y, const u64 *prime,
		const u64 k0) noexcept
{
	u64	s[N + 1];
	u64	r[N];
	vli_set<N>(r, y);
	for (uint i=0; i < N; i++) {
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		s[N] += vli_add_to<N>(r, s);
		vli_rshift1w<N>(r, s[N]);
	}
	vli_mod<N>(result, r, prime);
}

template<const uint N> forceinline
static void
vli_mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
		const u64 k0) noexcept
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

template<const uint N> forceinline
static void
vli_mont_sqr(u64 *result, const u64 *x, const u64 *prime, const u64 k0) noexcept
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
mont_reduction(u64 *result, const u64 *y, const u64 *prime) noexcept
{
	u64	s[N + 1];
	u64	r[N];
	vli_set<N>(r, y);
	for (uint i=0; i < N; i++) {
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		s[N] += vli_add_to<N>(r, s);
		vli_rshift1w<N>(r, s[N]);
	}
	vli_mod<N>(result, r, prime);
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


// secp256k1 prime optimize
// p is 2^256 - 2^32 - 0x3d1 = 2^256 - 0x1000003d1
/* result > mod (result = mod + remainder), so subtract mod to
 * get remainder.
 */
static forceinline void
btcp_mod(u64 *res, const u64 *product) noexcept
{
	u64 c = -secp256k1_p[0];
	u64 t[4 + 1];
	u64 r[4 + 2];
	vli_set<4>(r, product);
	r[4] = 0;
	vli_umult2<4>(t, product+4, c);
	r[5] = vli_add_to<5>(r, t);
	vli_clear<4>(t);
	vli_umult2<2>(t, r+4, c);
	r[4] = vli_add_to<4>(r, t);
	if (r[4]) {
		vli_sub<4>(res, r, secp256k1_p);
		return;
	}
	auto carry = vli_sub<4>(res, r, secp256k1_p);
	if (carry) vli_set<4>(res, r);
}


forceinline static void
btc_reduction(u64 *result, const u64 *y) noexcept
{
	u64	s[4 + 1];
	u64	r[4];
	vli_set<4>(r, y);
	for (uint i=0; i < 4; i++) {
		u64	u = r[0] * secp256k1_p_k0;
		vli_btc_multP(s, u);
		s[4] += vli_add_to<4>(r, s);
		vli_rshift1w<4>(r, s[4]);
	}
	vli_mod<4>(result, r, secp256k1_p);
}

forceinline static void
btc_mont_mult(u64 *result, const u64 *x, const u64 *y) noexcept
{
	u64	s[4 + 1];
	u64	r[4 + 1];
	vli_clear<4 + 1>(r);
	for (uint i=0; i < 4;i++) {
		vli_umult2<4>(s, x, y[i]);
		vli_add_to<4 + 1>(r, s);
		u64	u = r[0] * secp256k1_p_k0;
		vli_btc_multP(s, u);
		u = vli_add_to<4 + 1>(r, s);
		vli_rshift1w<4 + 1>(r, u);
	}
	vli_mod<4>(result, r, secp256k1_p, r[4] != 0);
}


#ifdef	WITH_SM2_MULTSTEP
forceinline static void
sm2p_multStep(u64& r0, u64& r1, u64& r2, u64& r3, u64& r4, const u64& x0,
		const u64& x1, const u64& x2, const u64& x3, const u64 yi) noexcept
{
	uint128_t	pd;
	u64			cc, t0;	// t0 never overflow,  0xf * 0xf = 225 + 30 ... 255(ff)
	pd.mul_64_64(x0, yi);
	cc = 0;
	r0 = u64_addc(r0, pd.m_low(), cc);
	t0 = u64_addcz(pd.m_high(), cc);
	pd.mul_64_64(x1, yi);
	//cc = 0;
	r1 = u64_addc(r1, t0, cc);
	t0 = u64_addcz(pd.m_high(), cc);
	//cc = 0;
	r1 = u64_addc(r1, pd.m_low(), cc);
	t0 = u64_addcz(t0, cc);
	pd.mul_64_64(x2, yi);
	//cc = 0;
	r2 = u64_addc(r2, t0, cc);
	t0 = u64_addcz(pd.m_high(), cc);
	//cc = 0;
	r2 = u64_addc(r2, pd.m_low(), cc);
	t0 = u64_addcz(t0, cc);
	pd.mul_64_64(x3, yi);
	//cc = 0;
	r3 = u64_addc(r3, t0, cc);
	t0 = u64_addcz(pd.m_high(), cc);
	//cc = 0;
	r3 = u64_addc(r3, pd.m_low(), cc);
	r4 = u64_addc(r4, t0, cc);
}
#else
// r0..3 = r0..3 + s[0..3]
// return s[4] + carry (of above)
forceinline static u64
vli4_addc_to(u64& r0, u64& r1, u64& r2, u64& r3, const u64* s,
	const u64 carry=0) noexcept
{
	u64	cc = 0;
	r0 = u64_addc(r0, s[0], cc);
	r1 = u64_addc(r1, s[1], cc);
	r2 = u64_addc(r2, s[2], cc);
	r3 = u64_addc(r3, s[3], cc);
	return u64_addc(carry, s[4], cc);
}
#endif

forceinline static void
sm2p_multN(u64 *result, const u64 *x, const u64 *y) noexcept
{
#ifdef	WITH_SM2_MULTSTEP
	u64	r0=0, r1=0, r2=0, r3=0;
	u64	x0=x[0], x1=x[1], x2=x[2], x3=x[3];
	u64	carry=0;
	{
		sm2p_multStep(r0, r1, r2, r3, carry, x0, x1, x2, x3, y[0]);
		sm2p_reductionStep(r0, r1, r2, r3, carry);
		sm2p_multStep(r1, r2, r3, r0, carry, x0, x1, x2, x3, y[1]);
		sm2p_reductionStep(r1, r2, r3, r0, carry);
		sm2p_multStep(r2, r3, r0, r1, carry, x0, x1, x2, x3, y[2]);
		sm2p_reductionStep(r2, r3, r0, r1, carry);
		sm2p_multStep(r3, r0, r1, r2, carry, x0, x1, x2, x3, y[3]);
		sm2p_reductionStep(r3, r0, r1, r2, carry);
	}
#else
	u64	r0, r1, r2, r3;
	u64	carry;
	{
		u64	s[4+1];
		vli_umult2<4>(s, x, y[0]);
		vli4_load(s, r0, r1, r2, r3);
		carry = s[4];
		sm2p_reductionStep(r0, r1, r2, r3, carry);
		vli_umult2<4>(s, x, y[1]);
		carry = vli4_addc_to(r1, r2, r3, r0, s, carry);
		sm2p_reductionStep(r1, r2, r3, r0, carry);
		vli_umult2<4>(s, x, y[2]);
		carry = vli4_addc_to(r2, r3, r0, r1, s, carry);
		sm2p_reductionStep(r2, r3, r0, r1, carry);
		vli_umult2<4>(s, x, y[3]);
		carry = vli4_addc_to(r3, r0, r1, r2, s, carry);
		sm2p_reductionStep(r3, r0, r1, r2, carry);
	}
#endif
	// sm2p_mod
#ifdef	WITH_SM2P_MOD2
	sm2p_mod2(result, r0, r1, r2, r3, carry);
#else
	{
		u64	cc=0;
		u64 s0 = u64_subc(r0, sm2_p[0], cc);
		u64 s1 = u64_subc(r1, sm2_p[1], cc);
		u64 s2 = u64_subc(r2, sm2_p[2], cc);
		u64 s3 = u64_subc(r3, sm2_p[3], cc);
		u64_subcz(carry, cc);
		if (cc != 0) vli4_save(result, r0, r1, r2, r3); else
			vli4_save(result, s0, s1, s2, s3);
	}
#endif
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
#elif	defined(__aarch64__)
	// x0 -- x3   register %%x4 -- %%x7
	register u64 x0 asm("x4") = x[0];
	register u64 x1 asm("x5") = x[1];
	register u64 x2 asm("x6") = x[2];
	register u64 x3 asm("x7") = x[3];
	asm volatile(
	// y[0] * x
	// load y0, y1
		"LDR x3, [%1]\n"
		"MUL	x9, x3, x4\n"
		"UMULH	x10, x3, x4\n"

		"MUL	x14, x3, x5\n"
		"ADDS	x10, x14, x10\n"
		"UMULH	x11, x3, x5\n"

		"MUL	x14, x3, x6\n"
		"ADCS	x11, x14, x11\n"
		"UMULH	x12, x3, x6\n"

		"MUL	x14, x3, x7\n"
		"ADCS	x12, x14, x12\n"
		"UMULH	x13, x3, x7\n"
		"ADC	x13, xzr, x13\n"
	// First reduction step
		"LSR	x14, x9, #32\n"
		"LSL	x15, x9, #32\n"
		"ADDS	x10, x10, x9\n"
		"ADCS	x11, x11, xzr\n"
		"ADCS	x12, x12, xzr\n"
		"ADCS	x13, x13, x9\n"
		"ADC	x9, xzr, xzr\n"
		"SUBS	x10, x10, x15\n"
		"SBCS	x11, x11, x14\n"
		"SBCS	x12, x12, x15\n"
		"SBCS	x13, x13, x14\n"
		"SBC	x9, x9, xzr\n"
	// y[1] * x
		"LDR x3, [%1, 8]\n"
		"MUL	x14, x3, x4\n"
		"ADDS	x10, x14, x10\n"
		"UMULH	x15, x3, x4\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x5\n"
		"ADDS	x11, x15, x11\n"
		"UMULH	x15, x3, x5\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x11, x14, x11\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x6\n"
		"ADDS	x12, x15, x12\n"
		"UMULH	x15, x3, x6\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x12, x14, x12\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x7\n"
		"ADDS	x13, x15, x13\n"
		"UMULH	x15, x3, x7\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x13, x14, x13\n"
		"ADC	x9, x15, x9\n"

	// Second reduction step
		"LSR	x14, x10, 32\n"
		"LSL	x15, x10, 32\n"
		"ADDS	x11, x11, x10\n"
		"ADCS	x12, x12, xzr\n"
		"ADCS	x13, x13, xzr\n"
		"ADCS	x9, x9, x10\n"
		"ADC	x10, xzr, xzr\n"
		"SUBS	x11, x11, x15\n"
		"SBCS	x12, x12, x14\n"
		"SBCS	x13, x13, x15\n"
		"SBCS	x9, x9, x14\n"
		"SBC	x10, x10, xzr\n"
	// y[2] * x
	// load y2, y3
		"LDR x3, [%1, 16]\n"
		"MUL	x14, x3, x4\n"
		"ADDS	x11, x14, x11\n"
		"UMULH	x15, x3, x4\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x5\n"
		"ADDS	x12, x15, x12\n"
		"UMULH	x15, x3, x5\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x12, x14, x12\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x6\n"
		"ADDS	x13, x15, x13\n"
		"UMULH	x15, x3, x6\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x13, x14, x13\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x7\n"
		"ADDS	x9, x15, x9\n"
		"UMULH	x15, x3, x7\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x9, x14, x9\n"
		"ADC	x10, x15, x10\n"

	// Third reduction step
		"LSR	x14, x11, 32\n"
		"LSL	x15, x11, 32\n"
		"ADDS	x12, x12, x11\n"
		"ADCS	x13, x13, xzr\n"
		"ADCS	x9, x9, xzr\n"
		"ADCS	x10, x10, x11\n"
		"ADC	x11, xzr, xzr\n"
		"SUBS	x12, x12, x15\n"
		"SBCS	x13, x13, x14\n"
		"SBCS	x9, x9, x15\n"
		"SBCS	x10, x10, x14\n"
		"SBC	x11, x11, xzr\n"
	// y[3] * x
		"LDR x3, [%1, 24]\n"
		"MUL	x14, x3, x4\n"
		"ADDS	x12, x14, x12\n"
		"UMULH	x15, x3, x4\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x5\n"
		"ADDS	x13, x15, x13\n"
		"UMULH	x15, x3, x5\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x13, x14, x13\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x6\n"
		"ADDS	x9, x15, x9\n"
		"UMULH	x15, x3, x6\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x9, x14, x9\n"
		"ADC	x15, x15, xzr\n"

		"MUL	x14, x3, x7\n"
		"ADDS	x10, x15, x10\n"
		"UMULH	x15, x3, x7\n"
		"ADC	x15, x15, xzr\n"
		"ADDS	x10, x14, x10\n"
		"ADC	x11, x15, x11\n"

	// Last reduction step
		"LSR	x14, x12, 32\n"
		"LSL	x15, x12, 32\n"
		"ADDS	x13, x13, x12\n"
		"ADCS	x9, x9, xzr\n"
		"ADCS	x10, x10, xzr\n"
		"ADCS	x11, x11, x12\n"
		"ADC	x12, xzr, xzr\n"
		"SUBS	x13, x13, x15\n"
		"SBCS	x9, x9, x14\n"
		"SBCS	x10, x10, x15\n"
		"SBCS	x11, x11, x14\n"
		"SBC	x12, x12, xzr\n"

		"ldp	x4, x5, [%2]\n"
		"ldp	x6, x7, [%2, 16]\n"
		"SUBS	x4, x13, x4\n"
		"SBCS	x5, x9, x5\n"
		"SBCS	x6, x10, x6\n"
		"SBCS	x7, x11, x7\n"
		"SBCS	x12, x12, xzr\n"

		"CSEL	x4, x4, x13, cs\n"
		"CSEL	x5, x5, x9, cs\n"
		"CSEL	x6, x6, x10, cs\n"
		"CSEL	x7, x7, x11, cs\n"
		"stp	x4, x5, [%0]\n"
		"stp	x6, x7, [%0, 16]\n"
		:
		: "r" (result), "r" (y), "r" (sm2_p), "r" (x0), "r" (x1), "r" (x2),
		"r" (x3)
		: "%x3", "%x9", "%x10", "%x11", "%x12", "%x13", "%x14", "%x15",
		"cc", "memory");
#else
	sm2p_multN(result, x, y);
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
	// reduction r[0..N-1]
	for (uint i=0; i < N; i++) {
		u64	u = r[0] * k0;
		vli_umult2<N>(s, prime, u);
		s[N] += vli_add_to<N>(r, s);
		vli_rshift1w<N>(r, s[N]);
	}
	// add r[N...2N-1], then mod P
#if	__cplusplus >= 201703L && defined(WITH_ASM)
	if constexpr(N==4) {
		s[N] = vli4_add_to(r, r + N);
		vli4_mod(result, r, prime, s[N] != 0);
	} else
#endif
	{
		s[N] = vli_add_to<N>(r, r + N);
		vli_mod<N>(result, r, prime, s[N] != 0);
	}
}

forceinline static void
btc_sqrN(u64 *result, const u64 *x) noexcept
{
	u64	r[4 * 2];
	vli_square<4>(r, x);
	btcp_mod(result, r);
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
		"r14", "r15", "cc", "memory");
#elif	defined(__aarch64__)
#ifdef	ommit_notwork
	sm2p_mult(result, x, x);
	// x0 -- x3   register %%x4 -- %%x7
#else
	register u64 x0 asm("x4") = x[0];
	register u64 x1 asm("x5") = x[1];
	register u64 x2 asm("x6") = x[2];
	register u64 x3 asm("x7") = x[3];
	asm volatile(
	// x[1:] * x[0]
		"MUL	x9, x4, x5\n"
		"UMULH	x10, x4, x5\n"

		"MUL	x21, x4, x6\n"
		"ADDS	x10, x21, x10\n"
		"UMULH	x11, x4, x6\n"

		"MUL	x21, x4, x7\n"
		"ADCS	x11, x21, x11\n"
		"UMULH	x12, x4, x7\n"
		"ADC	x12, xzr, x12\n"
	// x[2:] * x[1]
		"MUL	x21, x5, x6\n"
		"ADDS	x11, x21, x11\n"
		"UMULH	x21, x5, x6\n"
		"ADCS	x12, x21, x12\n"
		"ADC	x13, xzr, xzr\n"

		"MUL	x21, x5, x7\n"
		"ADDS	x12, x21, x12\n"
		"UMULH	x21, x5, x7\n"
		"ADC	x13, x21, x13\n"
	// x[3] * x[2]
		"MUL	x21, x6, x7\n"
		"ADDS	x13, x21, x13\n"
		"UMULH	x14, x6, x7\n"
		"ADC	x14, xzr, x14\n"

		//"MOV	x15, xzr\n"
	// *2
		"ADDS	x9, x9, x9\n"
		"ADCS	x10, x10, x10\n"
		"ADCS	x11, x11, x11\n"
		"ADCS	x12, x12, x12\n"
		"ADCS	x13, x13, x13\n"
		"ADCS	x14, x14, x14\n"
		"ADC	x15, xzr, xzr\n"
		//"ADC	x15, xzr, x15\n"
	// Missing products
		"MUL	x3, x4, x4\n"
		"UMULH	x21, x4, x4\n"
		"ADDS	x9, x21, x9\n"

		"MUL	x21, x5, x5\n"
		"ADCS	x10, x21, x10\n"
		"UMULH	x21, x5, x5\n"
		"ADCS	x11, x21, x11\n"

		"MUL	x21, x6, x6\n"
		"ADCS	x12, x21, x12\n"
		"UMULH	x21, x6, x6\n"
		"ADCS	x13, x21, x13\n"

		"MUL	x21, x7, x7\n"
		"ADCS	x14, x21, x14\n"
		"UMULH	x21, x7, x7\n"
		"ADCS	x15, x21, x15\n"
	// First reduction step
		"MOV	x19, x3\n"
		"MOV	x3, xzr\n"
		"ADDS	x9, x9, x19\n"
		"ADCS	x10, x10, xzr\n"
		"ADCS	x11, x11, xzr\n"
		"ADCS	x3, xzr, x19\n"
		"ADC	x20, xzr, xzr\n"
		"LSL	x21, x19, 32\n"
		"LSR	x19, x19, 32\n"
		"SUBS	x9, x9, x21\n"
		"SBCS	x10, x10, x19\n"
		"SBCS	x11, x11, x21\n"
		"SBCS	x3, x3, x19\n"
		"SBCS	x20, x20, xzr\n"
	// Second reduction step
		"MOV	x19, x9\n"
		"MOV	x9, xzr\n"
		"ADDS	x10, x10, x19\n"
		"ADCS	x11, x11, xzr\n"
		"ADCS	x3, x3, xzr\n"
		"ADCS	x9, x20, x19\n"
		"ADC	x20, xzr, xzr\n"
		"LSL	x21, x19, 32\n"
		"LSR	x19, x19, 32\n"
		"SUBS	x10, x10, x21\n"
		"SBCS	x11, x11, x19\n"
		"SBCS	x3, x3, x21\n"
		"SBCS	x9, x9, x19\n"
		"SBCS	x20, x20, xzr\n"
	// Third reduction step
		"MOV	x19, x10\n"
		"MOV	x10, xzr\n"
		"ADDS	x11, x11, x19\n"
		"ADCS	x3, x3, xzr\n"
		"ADCS	x9, x9, xzr\n"
		"ADCS	x10, x20, x19\n"
		"ADC	x20, xzr, xzr\n"
		"LSL	x21, x19, 32\n"
		"LSR	x19, x19, 32\n"
		"SUBS	x11, x11, x21\n"
		"SBCS	x3, x3, x19\n"
		"SBCS	x9, x11, x21\n"
		"SBCS	x10, x10, x19\n"
		"SBCS	x20, x20, xzr\n"
	// Last reduction step
		"MOV	x19, x11\n"
		"MOV	x11, xzr\n"
		"ADDS	x3, x3, x19\n"
		"ADCS	x9, x9, xzr\n"
		"ADCS	x10, x10, xzr\n"
		"ADCS	x11, x20, x19\n"
		"ADC	x20, xzr, xzr\n"
		"LSL	x21, x19, 32\n"
		"LSR	x19, x19, 32\n"
		"SUBS	x3, x3, x21\n"
		"SBCS	x9, x9, x19\n"
		"SBCS	x10, x10, x21\n"
		"SBCS	x11, x11, x19\n"
		"SBCS	x20, x20, xzr\n"
	// Add bits [511:256] of the sqr result
		"ADDS	x3, x12, x3\n"
		"ADCS	x9, x13, x9\n"
		"ADCS	x10, x14, x10\n"
		"ADCS	x11, x15, x11\n"
		"ADC	x20, xzr, xzr\n"

		"LDP	x4, x5, [%1]\n"
		"LDP	x6, x7, [%1, 16]\n"
		"SUBS	x4, x3, x4\n"
		"SBCS	x5, x9, x5\n"
		"SBCS	x6, x10, x6\n"
		"SBCS	x7, x11, x7\n"
		"SBCS	x20, x20, xzr\n"

		"CSEL	x4, x4, x3, cc\n"
		"CSEL	x5, x5, x9, cc\n"
		"CSEL	x6, x6, x10, cc\n"
		"CSEL	x7, x7, x11, cc\n"
		"STP	x4, x5, [%0]\n"
		"STP	x6, x7, [%0, 16]\n"
		:
		: "r" (result), "r" (sm2_p), "r" (x0), "r" (x1), "r" (x2),
		"r" (x3)
		: "%x3", "%x9", "%x10", "%x11", "%x12", "%x13", "%x14", "%x15",
		"%x21", "%x19", "%x20", "cc", "memory");
#endif
#else
	u64	r[8];
	vli_square<4>(r, x);
	sm2p_reduction(result, r, true);
#endif
}

#endif	//	__MONT_HPP__
