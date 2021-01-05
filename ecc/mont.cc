/*
 * Copyright(c) 2020 Jesse Kuang  <jkuang@21cn.com>
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

#include <errno.h>
#ifdef	__x86_64__
#include <cpuid.h>
#endif
#include "ecc.h"
#include "curve_const.hpp"
#include "mont.hpp"

//#pragma GCC push_options
//#pragma GCC optimize ("unroll-loops")
//#pragma GCC pop_options


using namespace vli;

u64 vli_asm_acc()
{
#ifdef	__x86_64__
	uint	_eax, _ebx, _ecx, _edx;
	__cpuid(1, _eax, _ebx, _ecx, _edx);
	return _ecx & bit_FMA;
#else
	return 0;
#endif
}

/* Computes result = product % mod using Barrett's reduction with precomputed
 * value mu appended to the mod after ndigits, mu = (2^{2w} / mod) and have
 * length ndigits + 1, where mu * (2^w - 1) should not overflow ndigits
 * boundary.
 *
 * Reference:
 * R. Brent, P. Zimmermann. Modern Computer Arithmetic. 2010.
 * 2.4.1 Barrett's algorithm. Algorithm 2.5.
 */
void vli_mmod_barrett(u64 *result, const u64 *product, const u64 *mod)
{
#ifdef	WITH_C2GO
	vli_mmod_barrett<4>(result, product, mod);
#else
	bn_prod<4>	*prod = reinterpret_cast<bn_prod<4> *>(const_cast<u64 *>(product));
	bignum<4>	*res = reinterpret_cast<bignum<4> *>(result);
	bignum<4>	*p = reinterpret_cast<bignum<4> *>(const_cast<u64 *>(mod));
	bignum<5>	*mu = reinterpret_cast<bignum<5> *>(const_cast<u64 *>(mod+4));
	prod->mmod_barrett(*res, *p, *mu);
#endif
}

#ifndef  WITH_C2GO
void vli_div_barrett(u64 *result, const u64 *product, const u64 *mu)
{
#ifdef	ommit
	vli_div_barrett<4>(result, product, mu);
#else
	bn_prod<4>	*prod = reinterpret_cast<bn_prod<4> *>(const_cast<u64 *>(product));
	bignum<4>	*res = reinterpret_cast<bignum<4> *>(result);
	bignum<5>	*muPtr = reinterpret_cast<bignum<5> *>(const_cast<u64 *>(mu));
	prod->div_barrett(*res, *muPtr);
#endif
}
#endif

void to_montgomery(u64 *res, const u64 *x, const montParams *pa)
{
	const u64	*rr = pa->rr;
	const u64	*prime = pa->p;
	const u64	k0 = pa->k0;
	mont_mult<4>(res, x, rr, prime, k0);
}

void from_montgomery(u64 *res, const u64 *y, const montParams *pa)
{
	const u64	*prime = pa->p;
	const u64	k0 = pa->k0;
	mont_reduction<4>(res, y, prime, k0);
}

//static u64 montOne[]={1, 0, 0, 0};
void mont_mod_mult(u64 *res, const u64 *x, const u64 *y, const montParams *pa)
{
	const u64	*prime = pa->p;
	const u64	k0 = pa->k0;
	mont_mult<4>(res, x, y, prime, k0);
}

void mont_mod_sqr(u64 *res, const u64 *x, const montParams *pa, const u64 n)
{
	const u64	*prime = pa->p;
	const u64	k0 = pa->k0;
	mont_sqr<4>(res, x, prime, k0);
	for (uint i=1; i < n; i++) {
		mont_sqr<4>(res, res, prime, k0);
	}
}

#ifndef  WITH_C2GO
void mont_mod_exp(u64 *result, const u64 *x, const u64 *y, const montParams *pa)
{
	const u64	*rr = pa->rr;
	const u64	*prime = pa->p;
	const u64	k0 = pa->k0;
	u64	xp[4];
	u64	t[4];
	int	num_bits = vli_num_bits<4>(y);
	mont_mult<4>(xp, x, rr, prime, k0);
	mont_reduction<4>(t, rr, prime, k0);
	for (int i = num_bits - 1;i >= 0; i--) {
		mont_sqr<4>(t, t, prime, k0);
		if (vli_test_bit(y, i)) mont_mult<4>(t, t, xp, prime, k0);
	}
	//mont_mult<4>(result, montOne, t, prime, k0);
	mont_reduction<4>(result, t, prime, k0);
}

void mont_sm2_mod_mult_p(u64 *result, const u64 *x, const u64 *y)
{
#ifdef	ommit
	u64	xp[4];
	u64	yp[4];
	u64	r[4];
	mont_multP(xp, x, sm2_p_rr, sm2_p);
	mont_multP(yp, y, sm2_p_rr, sm2_p);
	mont_multP(r, xp, yp, sm2_p);
	mont_reductionP(result, r, sm2_p);
#else
	bignum<4> xp;
	bignum<4> yp;
	bignum<4> r;
	bignum<4>	*res = reinterpret_cast<bignum<4> *>(result);
	bignum<4>	*rr = reinterpret_cast<bignum<4> *>(const_cast<u64 *>(sm2_p_rr));
	bignum<4>	*p = reinterpret_cast<bignum<4> *>(const_cast<u64 *>(sm2_p));
	mont_mult<sm2_p_k0>(xp, x, *rr, *p);
	mont_mult<sm2_p_k0>(yp, y, *rr, *p);
	mont_mult<sm2_p_k0>(r, xp, yp, *p);
	mont_reduction<sm2_p_k0>(*res, r, *p);
#endif
}

void mont_sm2_mod_mult_n(u64 *result, const u64 *x, const u64 *y)
{
#ifdef	ommit
	u64	xp[4];
	u64	yp[4];
	u64	r[4];
	mont_mult<4>(xp, x, sm2_n_rr, sm2_n, sm2_n_k0);
	mont_mult<4>(yp, y, sm2_n_rr, sm2_n, sm2_n_k0);
	mont_mult<4>(r, xp, yp, sm2_n, sm2_n_k0);
	mont_reduction<4>(result, r, sm2_n, sm2_n_k0);
#else
	bignum<4> xp;
	bignum<4> yp;
	bignum<4> r;
	bignum<4>	*res=reinterpret_cast<bignum<4> *>(result);
	bignum<4>	*rr = reinterpret_cast<bignum<4> *>(const_cast<u64 *>(sm2_n_rr));
	bignum<4>	*p = reinterpret_cast<bignum<4> *>(const_cast<u64 *>(sm2_n));
	mont_mult<sm2_n_k0>(xp, x, *rr, *p);
	mont_mult<sm2_n_k0>(yp, y, *rr, *p);
	mont_mult<sm2_n_k0>(r, xp, yp, *p);
	mont_reduction<sm2_n_k0>(*res, r, *p);
#endif
}
#endif

#ifdef	WITH_C2GO
void vli_sm2_mult_p(GoSlice *result, const u64 u)
#else
void vli_sm2_mult_p(u64 *result, const u64 rLen, const u64 u)
#endif
{
#ifdef	WITH_C2GO
	if (result->len < 5) return;
	vli_sm2_multP((u64 *)result->data, u);
#else
	if (rLen < 5) return;
	vli_sm2_multP(result, u);
#endif
}

#ifdef	WITH_C2GO
void vli_mod_inv(u64 *result, const u64 *input, const u64 *mod, u64 *buff)
{
	vli_mod_inv<4>(result, input, mod, buff);
}
#else
void vli_mod_inv(u64 *result, const u64 *input, const u64 *mod)
{
	vli_mod_inv<4>(result, input, mod);
}
#endif
