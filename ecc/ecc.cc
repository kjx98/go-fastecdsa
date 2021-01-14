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
#ifdef	WITH_SYS_RANDOM
#include <sys/random.h>
#endif
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

#ifdef	WITH_SHAMIR
/* Computes R = u1P + u2Q mod p using Shamir's trick.
 * Based on: Kenneth MacKay's micro-ecc (2014).
 */
void ecc_point_mult_shamir(const u64 *result_x, const u64 *result_y,
			   const u64 *u1, const u64 *p_x, const u64 *p_y,
			   const u64 *u2, const u64 *q_x, const u64 *q_y,
			   const struct ecc_curve *curve)
{
	u64 z[ndigits];
	u64 *rx = (u64 *)result_x;
	u64 *ry = (u64 *)result_y;
	unsigned int ndigits = curve->ndigits; //curve->g.ndigits;
	unsigned int num_bits;
	struct ecc_point sum;//  = ECC_POINT_INIT({}, {}, ndigits);
	const struct ecc_point *points[4];
	const struct ecc_point *point;
	unsigned int idx;
	int i;

	ecc_point_add(&sum, p, q, curve);
	points[0] = nullptr;
	points[1] = p;
	points[2] = q;
	points[3] = &sum;

	num_bits = max<uint>(vli_num_bits(u1, ndigits), vli_num_bits(u2, ndigits));
	i = num_bits - 1;
	idx = (!!vli_test_bit(u1, i)) | ((!!vli_test_bit(u2, i)) << 1);
	point = points[idx];

	vli_set(rx, point->x, ndigits);
	vli_set(ry, point->y, ndigits);
	vli_clear(z + 1, ndigits - 1);
	z[0] = 1;

	for (--i; i >= 0; i--) {
		ecc_point_double_jacobian(rx, ry, z, curve->p, ndigits);
		idx = (!!vli_test_bit(u1, i)) | ((!!vli_test_bit(u2, i)) << 1);
		point = points[idx];
		if (point) {
			u64 tx[ndigits];
			u64 ty[ndigits];
			u64 tz[ndigits];

			vli_set(tx, point->x, ndigits);
			vli_set(ty, point->y, ndigits);
			apply_z(tx, ty, z, curve->p, ndigits);
			vli_mod_sub(tz, rx, tx, curve->p, ndigits);
			xycz_add(tx, ty, rx, ry, curve->p, ndigits);
			vli_mod_mult_fast(z, z, tz, curve->p, ndigits);
		}
	}
	apply_z(rx, ry, z, curve->p, ndigits);
}
#endif

forceinline static
void ecc_swap_digits(const u64 *in, u64 *out, uint ndigits)
{
	const be64 *src = (be64 *)in;
	uint i;
//#pragma GCC unroll 4
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

#ifdef	WITH_ECCKEY
int ecc_is_key_valid(unsigned int curve_id, unsigned int ndigits,
		     const u64 *private_key, unsigned int private_key_len)
{
	uint nbytes;
	const struct ecc_curve *curve = ecc_get_curve(curve_id);

	nbytes = vli_bytes(ndigit);

	if (private_key_len != nbytes || ndigits != 4)
		return -EINVAL;

	return __ecc_is_key_valid<4>(curve, private_key);
}
#endif

/*
 * ECC private keys are generated using the method of extra random bits,
 * equivalent to that described in FIPS 186-4, Appendix B.4.1.
 *
 * d = (c mod(n–1)) + 1    where c is a string of random bits, 64 bits longer
 *                         than requested
 * 0 <= c mod(n-1) <= n-2  and implies that
 * 1 <= d <= n-1
 *
 * This method generates a private key uniformly distributed in the range
 * [1, n-1].
 */
#ifdef	WITH_SYS_RANDOM
int ecc_gen_privkey(unsigned int curve_id, unsigned int ndigits, u64 *privkey)
{
	const ecc_curve *curve = ecc_get_curve(curve_id);
	u64 priv[4];
	unsigned int nbytes = vli_bytes(ndigits);
	unsigned int nbits = vli_num_bits(curve->n, ndigits);
	int err;

	/* Check that N is included in Table 1 of FIPS 186-4, section 6.1.1 */
	if (nbits < 160 || ndigits > ARRAY_SIZE(priv))
		return -EINVAL;

	/*
	 * FIPS 186-4 recommends that the private key should be obtained from a
	 * RBG with a security strength equal to or greater than the security
	 * strength associated with N.
	 *
	 * The maximum security strength identified by NIST SP800-57pt1r4 for
	 * ECC is 256 (N >= 512).
	 *
	 * This condition is met by the default RNG because it selects a favored
	 * DRBG with a security strength of 256.
	 */

#ifdef	WITH_SYS_RANDOM
	err = getrandom(priv, nbytes, 0);
	if (err)
		return err;
#endif

	/* Make sure the private key is in the valid range. */
	if (__ecc_is_key_valid(curve, priv, ndigits))
		return -EINVAL;

	ecc_swap_digits(priv, privkey, ndigits);

	return 0;
}
#endif

int ecc_make_pub_key(unsigned int curve_id, unsigned int ndigits,
		     const u64 *private_key, u64 *public_key)
{
	int ret = 0;
	u64	pk_x[4];
	u64	pk_y[4];
	u64 priv[4];
	const ecc_curve *curve = ecc_get_curve(curve_id);

	if (!private_key || !curve || curve->ndigits != ndigits
	|| ndigits > ARRAY_SIZE(priv)) {
		ret = -EINVAL;
		return ret;
	}

	ecc_swap_digits(private_key, priv, ndigits);

	ecc_point_mult<4>(pk_x, pk_y, curve->gx, curve->gy, priv, nullptr, curve);
	if (ecc_point_is_zero<4>(pk_x, pk_y)) {
		ret = -EAGAIN;
		return ret;
	}

	ecc_swap_digits(pk_x, public_key, ndigits);
	ecc_swap_digits(pk_y, &public_key[ndigits], ndigits);

	return ret;
}

/* SP800-56A section 5.6.2.3.4 partial verification: ephemeral keys only */
int ecc_is_pubkey_valid_partial(const uint curve_id,
				const u64 *pk_x, const u64 *pk_y)
{
	const ecc_curve *curve = ecc_get_curve(curve_id);
	u64 yy[4], xxx[4], w[4];

	if (curve == nullptr || curve->ndigits != 4) return -EINVAL;
	uint	ndigits = curve->ndigits;
	/* Check 1: Verify key is not the zero point. */
	if (ecc_point_is_zero<4>(pk_x, pk_y)) return -EINVAL;

	/* Check 2: Verify key is in the range [1, p-1]. */
	if (vli_cmp<4>(curve->p, pk_x) != 1)
		return -EINVAL;
	if (vli_cmp<4>(curve->p, pk_y) != 1)
		return -EINVAL;

	/* Check 3: Verify that y^2 == (x^3 + a·x + b) mod p */
	vli_mod_square_fast<4>(yy, pk_y, curve->p); /* y^2 */
	vli_mod_square_fast<4>(xxx, pk_x, curve->p); /* x^2 */
	vli_mod_mult_fast(xxx, xxx, pk_x, curve->p, ndigits); /* x^3 */
	vli_mod_mult_fast(w, curve->a, pk_x, curve->p, ndigits); /* a·x */
	vli_mod_add<4>(w, w, curve->b, curve->p); /* a·x + b */
	vli_mod_add<4>(w, w, xxx, curve->p); /* x^3 + a·x + b */
	if (vli_cmp<4>(yy, w) != 0) /* Equation */
		return -EINVAL;

	return 0;
}
#endif

#ifdef	WITH_SYS_RANDOM
int crypto_ecdh_shared_secret(unsigned int curve_id, unsigned int ndigits,
			      const u64 *private_key, const u64 *public_key,
			      u64 *secret)
{
	int ret = 0;
	struct ecc_point *product, *pk;
	u64 priv[4];
	u64 rand_z[4];
	unsigned int nbytes;
	const struct ecc_curve *curve = ecc_get_curve(curve_id);

	if (!private_key || !public_key || !curve || curve->ndigits != 4 ||
	    ndigits > ARRAY_SIZE(priv) || ndigits > ARRAY_SIZE(rand_z)) {
		ret = -EINVAL;
		goto out;
	}

	nbytes = vli_bytes(ndigits);

#ifdef	WITH_SYS_RANDOM
	// TODO: rng, should check return for error
	if (getrandom(rand_z, nbytes, 0) < 0) {
		ret = -ENOMEM;
		goto out;
	}
#endif

	pk = ecc_alloc_point(ndigits);
	if (!pk) {
		ret = -ENOMEM;
		goto out;
	}

	ecc_swap_digits(public_key, pk->x, ndigits);
	ecc_swap_digits(&public_key[ndigits], pk->y, ndigits);
	ret = ecc_is_pubkey_valid_partial(curve, pk);
	if (ret)
		goto err_alloc_product;

	ecc_swap_digits(private_key, priv, ndigits);

	product = ecc_alloc_point(ndigits);
	if (!product) {
		ret = -ENOMEM;
		goto err_alloc_product;
	}

	ecc_point_mult(product, pk, priv, rand_z, curve, ndigits);

	ecc_swap_digits(product->x, secret, ndigits);

	if (ecc_point_is_zero(product, ndigits))
		ret = -EFAULT;

	ecc_free_point(product);
err_alloc_product:
	ecc_free_point(pk);
out:
	return ret;
}
#endif
