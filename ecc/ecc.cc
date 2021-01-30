/*
 * Copyright (c) 2013, 2014 Kenneth MacKay. All rights reserved.
 * Copyright (c) 2019 Vitaly Chikunov <vt@altlinux.org>
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
#include "ecc.h"
#include "ecc_impl.hpp"

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#pragma GCC pop_options

#ifndef	ARRAY_SIZE
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))
#endif


/* Computes p_result = p_product % curve_p.
 * See algorithm 5 and 6 from
 * http://www.isys.uni-klu.ac.at/PDF/2001-0126-MT.pdf
 */
static void forceinline vli_mmod_fast_192(u64 *result, const u64 *product,
			      const u64 *curve_prime, u64 *tmp)
{
	int carry;

	vli_set<3>(result, product);

	vli_set<3>(tmp, &product[3]);
	carry = vli_add_to<3>(result, tmp);

	tmp[0] = 0;
	tmp[1] = product[3];
	tmp[2] = product[4];
	carry += vli_add_to<3>(result, tmp);

	tmp[0] = tmp[1] = product[5];
	tmp[2] = 0;
	carry += vli_add_to<3>(result, tmp);

	while (carry || vli_cmp<3>(curve_prime, result) != 1)
		carry -= vli_sub_from<3>(result, curve_prime);
}

/* Computes result = product % curve_prime
 * from http://www.nsa.gov/ia/_files/nist-routines.pdf
 */
static void forceinline vli_mmod_fast_256(u64 *result, const u64 *product,
			      const u64 *curve_prime, u64 *tmp)
{
	int carry;

	/* t */
	vli_set<4>(result, product);

	/* s1 */
	tmp[0] = 0;
	tmp[1] = product[5] & 0xffffffff00000000ull;
	tmp[2] = product[6];
	tmp[3] = product[7];
	carry = vli_lshift1<4>(tmp, tmp);
	carry += vli_add_to<4>(result, tmp);

	/* s2 */
	tmp[1] = product[6] << 32;
	tmp[2] = (product[6] >> 32) | (product[7] << 32);
	tmp[3] = product[7] >> 32;
	carry += vli_lshift1<4>(tmp, tmp);
	carry += vli_add_to<4>(result, tmp);

	/* s3 */
	tmp[0] = product[4];
	tmp[1] = product[5] & 0xffffffff;
	tmp[2] = 0;
	tmp[3] = product[7];
	carry += vli_add_to<4>(result, tmp);

	/* s4 */
	tmp[0] = (product[4] >> 32) | (product[5] << 32);
	tmp[1] = (product[5] >> 32) | (product[6] & 0xffffffff00000000ull);
	tmp[2] = product[7];
	tmp[3] = (product[6] >> 32) | (product[4] << 32);
	carry += vli_add_to<4>(result, tmp);

	/* d1 */
	tmp[0] = (product[5] >> 32) | (product[6] << 32);
	tmp[1] = (product[6] >> 32);
	tmp[2] = 0;
	tmp[3] = (product[4] & 0xffffffff) | (product[5] << 32);
	carry -= vli_sub_from<4>(result, tmp);

	/* d2 */
	tmp[0] = product[6];
	tmp[1] = product[7];
	tmp[2] = 0;
	tmp[3] = (product[4] >> 32) | (product[5] & 0xffffffff00000000ull);
	carry -= vli_sub_from<4>(result, tmp);

	/* d3 */
	tmp[0] = (product[6] >> 32) | (product[7] << 32);
	tmp[1] = (product[7] >> 32) | (product[4] << 32);
	tmp[2] = (product[4] >> 32) | (product[5] << 32);
	tmp[3] = (product[6] << 32);
	carry -= vli_sub_from<4>(result, tmp);

	/* d4 */
	tmp[0] = product[7];
	tmp[1] = product[4] & 0xffffffff00000000ull;
	tmp[2] = product[5];
	tmp[3] = product[6] & 0xffffffff00000000ull;
	carry -= vli_sub_from<4>(result, tmp);

	if (carry < 0) {
		do {
			carry += vli_add_to<4>(result, curve_prime);
		} while (carry < 0);
	} else {
		while (carry || vli_cmp<4>(curve_prime, result) != 1)
			carry -= vli_sub_from<4>(result, curve_prime);
	}
}

/* Computes result = product % curve_prime for different curve_primes.
 *
 * Note that curve_primes are distinguished just by heuristic check and
 * not by complete conformance check.
 */
static forceinline bool vli_mmod_fast(u64 *result, u64 *product,
			  const u64 *curve_prime, unsigned int ndigits)
{
	u64 tmp[2 * 4];

	/* Currently, both NIST primes have -1 in lowest qword. */
	if (curve_prime[0] != -1ull) {
		/* Try to handle Pseudo-Marsenne primes. */
		if (ndigits != 4) return false;
		if (curve_prime[ndigits - 1] == -1ull) {
			vli_mmod_special<4>(result, product, curve_prime);
			return true;
		} else if (curve_prime[ndigits - 1] == 1ull << 63 &&
			   curve_prime[ndigits - 2] == 0) {
			vli_mmod_special2<4>(result, product, curve_prime);
			return true;
		}
		vli_mmod_barrett<4>(result, product, curve_prime);
		return true;
	}
	if ((curve_prime[1] >> 32) == 0) {
		// is SM2, curve_prime MUST following with mu
		switch (ndigits) {
		case 4:
			vli_mmod_barrett<4>(result, product, curve_prime);
			break;
		default:
			return false;
		}
		return true;
	}

	switch (ndigits) {
	case 3:
		vli_mmod_fast_192(result, product, curve_prime, tmp);
		break;
	case 4:
		vli_mmod_fast_256(result, product, curve_prime, tmp);
		break;
	default:
		//pr_err_ratelimited("ecc: unsupported digits size!\n");
		return false;
	}

	return true;
}

/* Computes result = left^2 % curve_prime. */
template<uint ndigits> forceinline
static void vli_mod_square_fast(u64 *result, const u64 *left,
				const u64 *curve_prime) noexcept
{
	u64 product[2 * ndigits];

	vli_square<ndigits>(product, left);
	vli_mmod_fast(result, product, curve_prime, ndigits);
}

/* Computes result = (left * right) % mod.
 * Assumes that mod is big enough curve order.
 */
/* Computes result = (left * right) % curve_prime. */
void vli_mod_mult_fast(u64 *result, const u64 *left, const u64 *right,
			      const u64 *curve_prime, unsigned int ndigits)
{
	u64 product[2 * 4];

	switch (ndigits) {
	case 3:
		vli_mult<3>(product, left, right);
		break;
	case 4:
		vli_mult<4>(product, left, right);
		break;
	default:	// error, no proc
		return;
	}
	vli_mmod_fast(result, product, curve_prime, ndigits);
}


void vli_mult(u64 *result, const u64 *left, const u64 *right)
{
	vli_mult<4>(result, left, right);
}


#ifdef	ommit
void vli_from_be64(u64 *dest, const void *src, uint ndigits)
{
	vli_from_be64<4>(dest, src);
}
#endif

/* ------ Point operations ------ */

forceinline static
void ecc_swap_digits(const u64 *in, u64 *out, uint ndigits)
{
	const be64 *src = (be64 *)in;
	uint i;
	for (i = 0; i < ndigits; i++)
		out[i] = be64toh(src[ndigits - 1 - i]);
}

#ifdef	ommit
template<uint ndigits> forceinline
static int __ecc_is_key_valid(const struct ecc_curve *curve,
			      const u64 *private_key) noexcept
{
	u64 one[ndigits] = { 1, };
	u64 res[ndigits];

	if (!private_key)
		return -EINVAL;

	if (curve->ndigits != ndigits)
		return -EINVAL;

	/* Make sure the private key is in the range [2, n-3]. */
	if (vli_cmp<ndigits>(one, private_key) != -1)
		return -EINVAL;
	vli_sub<ndigits>(res, curve->n, one);
	vli_sub_from<ndigits>(res, one);
	if (vli_cmp<ndigits>(res, private_key) != 1)
		return -EINVAL;

	return 0;
}


int ecc_make_pub_key(unsigned int curve_id, unsigned int ndigits,
		     const u64 *priv_key, u64 *public_key)
{
	int ret = 0;
	const ecc_curve *curve = ecc_get_curve(curve_id);

	if (!priv_key || !curve || curve->ndigits != ndigits) {
		ret = -EINVAL;
		return ret;
	}
	bignum<4>	secr(priv_key);
	point_t<4>	pt;
	curve->scalar_mult_base(pt, secr);

	pt.x.set(public_key);
	pt.y.set(&public_key[4]);

	return ret;
}

/* SP800-56A section 5.6.2.3.4 partial verification: ephemeral keys only */
int ecc_is_pubkey_valid_partial(const uint curve_id,
				const u64 *pk_x, const u64 *pk_y)
{
	const ecc_curve *curve = ecc_get_curve(curve_id);

	if (curve == nullptr || curve->ndigits != 4) return -EINVAL;
	bignum<4>	x1(pk_x), y1(pk_y);
	uint	ndigits = curve->ndigits;
	/* Check 1: Verify key is not the zero point. */
	if (x1.is_zero() || y1.is_zero()) return -EINVAL;

	/* Check 2: Verify key is in the range [1, p-1]. */
	if (curve->p < x1) return -EINVAL;
	if (curve->p < y1) return -EINVAL;

	/* Check 3: Verify that y^2 == (x^3 + a·x + b) mod p */
	if (! curve->is_on_curve(x1, y1) ) /* Equation */
		return -EINVAL;

	return 0;
}
#endif

