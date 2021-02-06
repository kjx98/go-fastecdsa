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
	result[5] = carry;
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



forceinline static void
sm2p_reduction(u64 *result, const u64 *y, const bool isProd=false) noexcept
{
#ifdef	__x86_64__1
	asm volatile("MOVQ (8*0)(%%rsi), acc0\n"
			MOVQ (8*1)(%%rsi), acc1
			MOVQ (8*2)(%%rsi), acc2
			MOVQ (8*3)(%%rsi), acc3
			XORQ acc4, acc4

	// Only reduce, no multiplications are needed
	// First stage
			MOVQ acc0, AX
			MOVQ acc0, t1
			SHLQ $32, acc0
			SHRQ $32, t1
			ADDQ %%rax, acc1
			ADCQ $0, acc2
			ADCQ $0, acc3
			ADCQ %%rax, acc4
			subq acc0, acc1
			sbbq t1, acc2
			sbbq acc0, acc3
			sbbq t1, acc4
	XORQ acc5, acc5
	// Second stage
			MOVQ acc1, AX
			MOVQ acc1, t1
			SHLQ $32, acc1
			SHRQ $32, t1
			ADDQ %%rax, acc2
			adcq $0, acc3
			adcq $0, acc4
			adcq %%rax, acc5
			subq acc1, acc2
			sbbq t1, acc3
			sbbq acc1, acc4
			sbbq t1, acc5
	XORQ acc0, acc0
	// Third stage
	MOVQ acc2, AX
	MOVQ acc2, t1
	SHLQ $32, acc2
	SHRQ $32, t1
	ADDQ %%rax, acc3
	adcq $0, acc4
	adcq $0, acc5
	adcq %%rax, acc0
	subq acc2, acc3
	sbbq t1, acc4
	sbbq acc2, acc5
	sbbq t1, acc0
	XORQ acc1, acc1
	// Last stage
	MOVQ acc3, AX
	MOVQ acc3, t1
	SHLQ $32, acc3
	SHRQ $32, t1
	ADDQ %%rax, acc4
	adcq $0, acc5
	adcq $0, acc0
	adcq %%rax, acc1
	subq acc3, acc4
	sbbq t1, acc5
	sbbq acc3, acc0
	sbbq t1, acc1

	MOVQ acc4, x_ptr
	MOVQ acc5, acc3
	MOVQ acc0, t0
	MOVQ acc1, t1

	SUBQ $-1, acc4
	SBBQ p256const0<>(SB), acc5
	SBBQ $-1, acc0
	SBBQ p256const1<>(SB), acc1

	CMOVQCS x_ptr, acc4
	CMOVQCS acc3, acc5
	CMOVQCS t0, acc0
	CMOVQCS t1, acc1

	MOVQ acc4, (8*0)(%%rdi)
	MOVQ acc5, (8*1)(%%rdi)
	MOVQ acc0, (8*2)(%%rdi)
	MOVQ acc1, (8*3)(%%rdi)
				: 			// acc4/5/0/1
				: "S" (y), "D" (result)
				: "cc", "memory");
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
	sm2p_mod(result, r, sm2_p, carry != 0);
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
		u64	u = (r[0] + y[i]*x[0]) * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_umult2<N>(s, x, y[i]);
		u += vli_add_to<N + 1>(r, s);
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
		u64	u = (r[0] + x[i]*x[0]) * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_umult2<N>(s, x, x[i]);
		u += vli_add_to<N + 1>(r, s);
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
		u64	u = (r[0] + y[i]*x[0]) * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_umult2<N>(s, x, y[i]);
		u += vli_add_to<N + 1>(r, s);
		vli_rshift1w<N + 1>(r, u);	
	}
	vli_mod<N>(result, r, prime, r[N] != 0);
}

template<const uint N, const u64 k0> forceinline
static void
mont_sqr(u64 *result, const u64 *x, const u64 *prime) noexcept
{
	u64	s[N + 1];
	u64	r[N + 1];
	vli_clear<N + 1>(r);
	for (uint i=0; i < N;i++) {
		u64	u = (r[0] + x[i]*x[0]) * k0;
		vli_umult2<N>(s, prime, u);
		u = vli_add_to<N + 1>(r, s);
		vli_umult2<N>(s, x, x[i]);
		u += vli_add_to<N + 1>(r, s);
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
#if	__cplusplus >= 201703L
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
sm2p_sqrN(u64 *result, const u64 *x, const int nTimes=1) noexcept
{
	u64	r[8];
	vli_set<4>(result, x);
	for (int i=0; i < nTimes; ++i) {
		vli_squareN<4>(r, result);
		sm2p_reduction(result, r, true);
	}
}

#endif	//	__MONT_HPP__
