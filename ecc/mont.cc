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
#include "vli.hpp"
#include "ecc.h"

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#pragma GCC pop_options

#ifdef	__cplusplus
extern "C" {
#endif
#ifdef	__cplusplus
}
#endif

static void mont_reduction(u64 *result, const u64 *prod, const u64 *prime,
			const u64 k0) noexcept
{
	u64	t[ECC_MAX_DIGITS * 2];
	u64	s[ECC_MAX_DIGITS * 2];
	u64	r[ECC_MAX_DIGITS * 2];
	vli_clear<8>(r);
	for (uint i=0; i<4; i++) {
		u64	u = (r[0] + prod[i]) * k0;
		//vli_umult<4>(s, prime, u);
		vli_uadd<8>(t, s, prod[i]);
		vli_add<8>(r, r, t);
		//vli_rshift1w<8>(r);
	}
	if (!vli_is_zero<4>(r+4) || vli_cmp<4>(r, prime) >= 0) {
		vli_sub<4>(r, r, prime);
	}
	vli_set<4>(result, r);
}


static void mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 k0) noexcept
{
	u64	t[ECC_MAX_DIGITS * 2];
	u64	s[ECC_MAX_DIGITS * 2];
	u64	r[ECC_MAX_DIGITS * 2];
	vli_clear<8>(r);
	for (uint i=0; i<4;i++) {
		u64	u = (r[0] + y[i]*x[0]) * k0;
		//vli_umult<4>(s, prime, u);
		//vli_umult<4>(t, x, y[i]);
		vli_add<8>(t, t, s);
		vli_add<8>(r, r, t);
		//vli_rshift1w<8>(r);	
	}
	if (!vli_is_zero<4>(r+4) || vli_cmp<4>(r, prime) >= 0) {
		vli_sub<4>(r, r, prime);
	}
	vli_set<4>(result, r);
}

void mont_MulMod(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 *rr, const u64 k0) noexcept
{
	u64	xp[ECC_MAX_DIGITS];
	u64	yp[ECC_MAX_DIGITS];
	u64	r[ECC_MAX_DIGITS];
	mont_mult(xp, x, rr, prime, k0);
	mont_mult(yp, y, rr, prime, k0);
	mont_mult(r, xp, yp, prime, k0);
	mont_reduction(result, r, prime, k0);
}
