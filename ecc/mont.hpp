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
#include "vli_bn.hpp"
#include "ecc_impl.hpp"


// SM2 prime optimize
// p is 2^256 - 2^224 - 2^96 + 2^64 -1
forceinline
static void vli_sm2_multP(u64 *result, const u64 u) noexcept
{
	u64	r[4];
	u64	t_low, t_high;
	t_low = u << 32;	// ^192
	t_high = ((u >> 32) & 0xffffffff);
	vli_clear<6>(result);
	result[3] = 0 - t_low;		// ^256 - ^224
	result[4] = u - t_high -1;
	r[0] =0;
	r[1] = t_low;
	r[2] = t_high;
	r[3] =0;
	//vli_sub_from<5>(result, t);	// ^256 - ^224 - ^96
	if (vli_sub_from<4>(result, r)) result[4]--;
	r[2] = 0;
	r[1] = u-1;
	r[0] = -u;		// ^64 -1
	if (vli_add_to<4>(result, r)) result[4]++;
}


template<const uint N> forceinline
static void
#ifdef	WITH_C2GO
mont_reduction(u64 *result, const u64 *y, const u64 *prime,
			const u64 k0, u64 *buff) noexcept
#else
mont_reduction(u64 *result, const u64 *y, const u64 *prime,
			const u64 k0) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*s = buff;
	u64	*r = s + N * 2;
#else
	u64	s[N * 2];
	u64	r[N + 2];
#endif
	vli_clear<N + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < N; i++) {
		u64	u = (r[0] + y[i]) * k0;
		vli_umult<N>(s, prime, u);
		vli_uadd_to<N + 2>(r, y[i]);
		vli_add_to<N + 2>(r, s);
		vli_rshift1w<N + 2>(r);	
	}
	if (r[N] !=0 || vli_cmp<N>(r, prime) >= 0) {
		vli_sub<N>(result, r, prime);
	} else vli_set<N>(result, r);
}

template<const uint N> forceinline
static void
#ifdef	WITH_C2GO
mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 k0, u64 *buff) noexcept
#else
mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 k0) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*s = buff;
	u64	*r = s + N * 2;
#else
	u64	s[N * 2];
	u64	r[N + 2];
#endif
	vli_clear<N + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < N;i++) {
		u64	u = (r[0] + y[i]*x[0]) * k0;
		vli_umult<N>(s, prime, u);
		vli_add_to<N + 2>(r, s);
		vli_umult<N>(s, x, y[i]);
		vli_add_to<N + 2>(r, s);
		vli_rshift1w<N + 2>(r);	
	}
	if (r[N] != 0 || vli_cmp<N>(r, prime) >= 0) {
		vli_sub<N>(result, r, prime);
	} else vli_set<N>(result, r);
}

template<const uint N> forceinline
static void
#ifdef	WITH_C2GO
mont_sqr(u64 *result, const u64 *x, const u64 *prime, const u64 k0, u64 *buff) noexcept
#else
mont_sqr(u64 *result, const u64 *x, const u64 *prime, const u64 k0) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*s = buff;
	u64	*r = s + N * 2;
#else
	u64	s[N * 2];
	u64	r[N + 2];
#endif
	vli_clear<N + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < N;i++) {
		u64	u = (r[0] + x[i]*x[0]) * k0;
		vli_umult<N>(s, prime, u);
		vli_add_to<N + 2>(r, s);
		vli_umult<N>(s, x, x[i]);
		vli_add_to<N + 2>(r, s);
		vli_rshift1w<N + 2>(r);	
	}
	if (r[N] != 0 || vli_cmp<N>(r, prime) >= 0) {
		vli_sub<N>(result, r, prime);
	} else vli_set<N>(result, r);
}

template<const uint N, const u64 k0> forceinline
static void
mont_reduction(u64 *result, const u64 *y, const u64 *prime) noexcept
{
	u64	s[N * 2];
	u64	r[N + 2];
	vli_clear<N + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < N; i++) {
		u64	u = (r[0] + y[i]) * k0;
		vli_umult<N>(s, prime, u);
		vli_uadd_to<N + 2>(r, y[i]);
		vli_add_to<N + 2>(r, s);
		vli_rshift1w<N + 2>(r);	
	}
	if (r[N] !=0 || vli_cmp<N>(r, prime) >= 0) {
		vli_sub<N>(result, r, prime);
	} else vli_set<N>(result, r);
}

template<const uint N, const u64 k0> forceinline
static void
mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime) noexcept
{
	u64	s[N * 2];
	u64	r[N + 2];
	vli_clear<N + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < N;i++) {
		u64	u = (r[0] + y[i]*x[0]) * k0;
		vli_umult<N>(s, prime, u);
		vli_add_to<N + 2>(r, s);
		vli_umult<N>(s, x, y[i]);
		vli_add_to<N + 2>(r, s);
		vli_rshift1w<N + 2>(r);	
	}
	if (r[N] != 0 || vli_cmp<N>(r, prime) >= 0) {
		vli_sub<N>(result, r, prime);
	} else vli_set<N>(result, r);
}

template<const uint N, const u64 k0> forceinline
static void
mont_sqr(u64 *result, const u64 *x, const u64 *prime) noexcept
{
	u64	s[N * 2];
	u64	r[N + 2];
	vli_clear<N + 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < N;i++) {
		u64	u = (r[0] + x[i]*x[0]) * k0;
		vli_umult<N>(s, prime, u);
		vli_add_to<N + 2>(r, s);
		vli_umult<N>(s, x, x[i]);
		vli_add_to<N + 2>(r, s);
		vli_rshift1w<N + 2>(r);	
	}
	if (r[N] != 0 || vli_cmp<N>(r, prime) >= 0) {
		vli_sub<N>(result, r, prime);
	} else vli_set<N>(result, r);
}

#endif	//	__MONT_HPP__
