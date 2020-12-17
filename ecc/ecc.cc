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
#include "vli.hpp"
#include "ecc.h"
#include "ecc_impl.hpp"
#include "ecc_curve_defs.h"

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#pragma GCC pop_options

#define ECC_DIGITS_TO_BYTES_SHIFT 3

static inline const ecc_curve *ecc_get_curve(uint curve_id)
{
	switch (curve_id) {
	/* In FIPS mode only allow P256 and higher */
	case ECC_CURVE_SECP256K1:
		return &secp256k1;
	case ECC_CURVE_NIST_P256:
		return &nist_p256;
	case ECC_CURVE_SM2:
		return &sm2_p256;
	default:
		return nullptr;
	}
}

#ifdef	ommit
__attribute__((optimize("unroll-loops")))
static inline void vli_clear(u64 *vli, uint ndigits)
{
	uint i;
	for (i = 0; i < ndigits; i++)
		vli[i] = 0;
}

/* Returns true if vli == 0, false otherwise. */
bool vli_is_zero(const u64 *vli, uint ndigits)
{
	uint i;
//#pragma unroll 4
	for (i = 0; i < ndigits; i++) {
		if (vli[i])
			return false;
	}

	return true;
}
#endif

#ifdef	ommit
static bool inline vli_is_negative(const u64 *vli, uint ndigits)
{
	return vli_test_bit(vli, ndigits * 64 - 1);
}

/* Counts the number of 64-bit "digits" in vli. */
static uint inline vli_num_digits(const u64 *vli, uint ndigits)
{
	int i;

	/* Search from the end until we find a non-zero digit.
	 * We do it in reverse because we expect that most digits will
	 * be nonzero.
	 */
	for (i = ndigits - 1; i >= 0 && vli[i] == 0; i--);

	return (i + 1);
}

/* Counts the number of bits required for vli. */
static uint inline vli_num_bits(const u64 *vli, uint ndigits)
{
	uint i, num_digits;
	u64 digit;

	num_digits = vli_num_digits(vli, ndigits);
	if (num_digits == 0)
		return 0;

	digit = vli[num_digits - 1];
	i = 64 - __builtin_clzl(digit);

	return ((num_digits - 1) * 64 + i);
}

/* Set dest from unaligned bit string src. */
void vli_from_be64(u64 *dest, const void *src, uint ndigits)
{
	uint i;
	const u64 *from = (u64 *)src;

	for (i = 0; i < ndigits; i++)
		dest[i] = be64toh(from[ndigits - 1 - i]);
}

void vli_from_le64(u64 *dest, const void *src, uint ndigits)
{
	uint i;
	const u64 *from = (u64 *)src;

	for (i = 0; i < ndigits; i++)
		dest[i] = le64toh(from[i]);
}

/* Sets dest = src. */
static inline void vli_set(u64 *dest, const u64 *src, unsigned int ndigits)
{
	uint i;

	for (i = 0; i < ndigits; i++)
		dest[i] = src[i];
}

/* Returns sign of left - right. */
int vli_cmp(const u64 *left, const u64 *right, unsigned int ndigits)
{
	int i;

	for (i = ndigits - 1; i >= 0; i--) {
		if (left[i] > right[i])
			return 1;
		else if (left[i] < right[i])
			return -1;
	}

	return 0;
}

/* Computes result = in << c, returning carry. Can modify in place
 * (if result == in). 0 < shift < 64.
 */
static u64 vli_lshift(u64 *result, const u64 *in, unsigned int shift,
		      unsigned int ndigits)
{
	u64 carry = 0;
	uint i;

	for (i = 0; i < ndigits; i++) {
		u64 temp = in[i];

		result[i] = (temp << shift) | carry;
		carry = temp >> (64 - shift);
	}

	return carry;
}

/* Computes vli = vli >> 1. */
static void vli_rshift1(u64 *vli, unsigned int ndigits)
{
	u64 *end = vli;
	u64 carry = 0;

	vli += ndigits;

	while (vli-- > end) {
		u64 temp = *vli;
		*vli = (temp >> 1) | carry;
		carry = temp << 63;
	}
}

/* Computes result = left + right, returning carry. Can modify in place. */
static u64 vli_add(u64 *result, const u64 *left, const u64 *right,
		   unsigned int ndigits)
{
	u64 carry = 0;
	uint i;

	for (i = 0; i < ndigits; i++) {
		u64 sum;

		sum = left[i] + right[i] + carry;
		if (sum != left[i])
			carry = (sum < left[i]);

		result[i] = sum;
	}

	return carry;
}

/* Computes result = left + right, returning carry. Can modify in place. */
static u64 vli_uadd(u64 *result, const u64 *left, u64 right,
		    unsigned int ndigits)
{
	u64 carry = right;
	uint i;

	for (i = 0; i < ndigits; i++) {
		u64 sum;

		sum = left[i] + carry;
		if (sum != left[i])
			carry = (sum < left[i]);
		else
			carry = !!carry;

		result[i] = sum;
	}

	return carry;
}

/* Computes result = left - right, returning borrow. Can modify in place. */
u64 vli_sub(u64 *result, const u64 *left, const u64 *right,
		   unsigned int ndigits)
{
	u64 borrow = 0;
	uint i;

	for (i = 0; i < ndigits; i++) {
		u64 diff;

		diff = left[i] - right[i] - borrow;
		if (diff != left[i])
			borrow = (diff > left[i]);

		result[i] = diff;
	}

	return borrow;
}

/* Computes result = left - right, returning borrow. Can modify in place. */
static u64 vli_usub(u64 *result, const u64 *left, u64 right,
	     unsigned int ndigits)
{
	u64 borrow = right;
	uint i;

	for (i = 0; i < ndigits; i++) {
		u64 diff;

		diff = left[i] - borrow;
		if (diff != left[i])
			borrow = (diff > left[i]);

		result[i] = diff;
	}

	return borrow;
}

static forceinline uint128_t mul_64_64(u64 left, u64 right)
{
#if defined(__SIZEOF_INT128__)
	uint128_t result( (__uint128_t)left * right );
#else
	u64 a0 = left & 0xffffffffull;
	u64 a1 = left >> 32;
	u64 b0 = right & 0xffffffffull;
	u64 b1 = right >> 32;
	u64 m0 = a0 * b0;
	u64 m1 = a0 * b1;
	u64 m2 = a1 * b0;
	u64 m3 = a1 * b1;

	m2 += (m0 >> 32);
	m2 += m1;

	/* Overflow */
	if (m2 < m1)
		m3 += 0x100000000ull;

	u64 m_low = (m0 & 0xffffffffull) | (m2 << 32);
	u64 m_high = m3 + (m2 >> 32);
	uint128_t	result(m_low, m_high);
#endif
	return result;
}

static forceinline uint128_t add_128_128(uint128_t a, uint128_t b)
{
#if defined(__SIZEOF_INT128__)
	uint128_t result(a.data() + b.data());
#else
	u64 m_low = a.m_low + b.m_low;
	u64 m_high = a.m_high + b.m_high + (result.m_low < a.m_low);
	uint128_t result(m_low, m_high);
#endif
	return result;
}

void vli_mult(u64 *result, const u64 *left, const u64 *right,
		     unsigned int ndigits)
{
	uint128_t r01( 0, 0 );
	u64 r2 = 0;
	unsigned int i, k;

	/* Compute each digit of result in sequence, maintaining the
	 * carries.
	 */
	for (k = 0; k < ndigits * 2 - 1; k++) {
		unsigned int min;

		if (k < ndigits)
			min = 0;
		else
			min = (k + 1) - ndigits;

		for (i = min; i <= k && i < ndigits; i++) {
			uint128_t product;

			//product = mul_64_64(left[i], right[k - i]);
			product.mul_64_64(left[i], right[k - i]);

			//r01 = add_128_128(r01, product);
			r01 += product;
			r2 += (r01.m_high() < product.m_high());
		}

		result[k] = r01.m_low();
		r01 = uint128_t(r01.m_high(), r2);
		//r01.m_low = r01.m_high;
		//r01.m_high = r2;
		r2 = 0;
	}

	result[ndigits * 2 - 1] = r01.m_low();
}

/* Compute product = left * right, for a small right value. */
static void vli_umult(u64 *result, const u64 *left, u64 right,
		      unsigned int ndigits)
{
	uint128_t r01( 0, 0 );
	unsigned int k;

	for (k = 0; k < ndigits; k++) {
		uint128_t product;

		//product = mul_64_64(left[k], right);
		product.mul_64_64(left[k], right);
		//r01 = add_128_128(r01, product);
		r01 += product;
		/* no carry */
		result[k] = r01.m_low();
		r01 = uint128_t(r01.m_high(), 0);
		//r01.m_low = r01.m_high;
		//r01.m_high = 0;
	}
	result[k] = r01.m_low();
	for (++k; k < ndigits * 2; k++)
		result[k] = 0;
}

static void vli_square(u64 *result, const u64 *left, unsigned int ndigits)
{
	uint128_t r01( 0, 0 );
	u64 r2 = 0;
	uint i, k;

	for (k = 0; k < ndigits * 2 - 1; k++) {
		unsigned int min;

		if (k < ndigits)
			min = 0;
		else
			min = (k + 1) - ndigits;

		for (i = min; i <= k && i <= k - i; i++) {
			uint128_t product;

			//product = mul_64_64(left[i], left[k - i]);
			product.mul_64_64(left[i], left[k - i]);

			if (i < k - i) {
				r2 += product.m_high() >> 63;
#ifdef	ommit
				product.m_high = (product.m_high << 1) |
						 (product.m_low >> 63);
				product.m_low <<= 1;
#endif
				u64 _high = (product.m_high() << 1) | (product.m_low() >> 63);
				u64 _low = product.m_low() << 1;
				product = uint128_t(_low, _high);
			}

			//r01 = add_128_128(r01, product);
			r01 += product;
			r2 += (r01.m_high() < product.m_high());
		}

		result[k] = r01.m_low();
		//r01.m_low = r01.m_high;
		//r01.m_high = r2;
		r01 = uint128_t(r01.m_high(), r2);
		r2 = 0;
	}

	result[ndigits * 2 - 1] = r01.m_low();
}

/* Computes result = (left + right) % mod.
 * Assumes that left < mod and right < mod, result != mod.
 */
static void vli_mod_add(u64 *result, const u64 *left, const u64 *right,
			const u64 *mod, unsigned int ndigits)
{
	u64 carry;

	carry = vli_add(result, left, right, ndigits);

	/* result > mod (result = mod + remainder), so subtract mod to
	 * get remainder.
	 */
	if (carry || vli_cmp(result, mod, ndigits) >= 0)
		vli_sub(result, result, mod, ndigits);
}

/* Computes result = (left - right) % mod.
 * Assumes that left < mod and right < mod, result != mod.
 */
static void vli_mod_sub(u64 *result, const u64 *left, const u64 *right,
			const u64 *mod, unsigned int ndigits)
{
	u64 borrow = vli_sub(result, left, right, ndigits);

	/* In this case, p_result == -diff == (max int) - diff.
	 * Since -x % d == d - x, we can get the correct result from
	 * result + mod (with overflow).
	 */
	if (borrow)
		vli_add(result, result, mod, ndigits);
}

/*
 * Computes result = product % mod
 * for special form moduli: p = 2^k-c, for small c (note the minus sign)
 *
 * References:
 * R. Crandall, C. Pomerance. Prime Numbers: A Computational Perspective.
 * 9 Fast Algorithms for Large-Integer Arithmetic. 9.2.3 Moduli of special form
 * Algorithm 9.2.13 (Fast mod operation for special-form moduli).
 */
static void vli_mmod_special(u64 *result, const u64 *product,
			      const u64 *mod, unsigned int ndigits)
{
	u64 c = -mod[0];
	u64 t[ECC_MAX_DIGITS * 2];
	u64 r[ECC_MAX_DIGITS * 2];

	vli_set(r, product, ndigits * 2);
	while (!vli_is_zero(r + ndigits, ndigits)) {
		vli_umult(t, r + ndigits, c, ndigits);
		vli_clear(r + ndigits, ndigits);
		vli_add(r, r, t, ndigits * 2);
	}
	vli_set(t, mod, ndigits);
	vli_clear(t + ndigits, ndigits);
	while (vli_cmp(r, t, ndigits * 2) >= 0)
		vli_sub(r, r, t, ndigits * 2);
	vli_set(result, r, ndigits);
}

/*
 * Computes result = product % mod
 * for special form moduli: p = 2^{k-1}+c, for small c (note the plus sign)
 * where k-1 does not fit into qword boundary by -1 bit (such as 255).

 * References (loosely based on):
 * A. Menezes, P. van Oorschot, S. Vanstone. Handbook of Applied Cryptography.
 * 14.3.4 Reduction methods for moduli of special form. Algorithm 14.47.
 * URL: http://cacr.uwaterloo.ca/hac/about/chap14.pdf
 *
 * H. Cohen, G. Frey, R. Avanzi, C. Doche, T. Lange, K. Nguyen, F. Vercauteren.
 * Handbook of Elliptic and Hyperelliptic Curve Cryptography.
 * Algorithm 10.25 Fast reduction for special form moduli
 */
static void vli_mmod_special2(u64 *result, const u64 *product,
			       const u64 *mod, unsigned int ndigits)
{
	u64 c2 = mod[0] * 2;
	u64 q[ECC_MAX_DIGITS];
	u64 r[ECC_MAX_DIGITS * 2];
	u64 m[ECC_MAX_DIGITS * 2]; /* expanded mod */
	int carry; /* last bit that doesn't fit into q */
	int i;

	vli_set(m, mod, ndigits);
	vli_clear(m + ndigits, ndigits);

	vli_set(r, product, ndigits);
	/* q and carry are top bits */
	vli_set(q, product + ndigits, ndigits);
	vli_clear(r + ndigits, ndigits);
	carry = vli_is_negative(r, ndigits);
	if (carry)
		r[ndigits - 1] &= (1ull << 63) - 1;
	for (i = 1; carry || !vli_is_zero(q, ndigits); i++) {
		u64 qc[ECC_MAX_DIGITS * 2];

		vli_umult(qc, q, c2, ndigits);
		if (carry)
			vli_uadd(qc, qc, mod[0], ndigits * 2);
		vli_set(q, qc + ndigits, ndigits);
		vli_clear(qc + ndigits, ndigits);
		carry = vli_is_negative(qc, ndigits);
		if (carry)
			qc[ndigits - 1] &= (1ull << 63) - 1;
		if (i & 1)
			vli_sub(r, r, qc, ndigits * 2);
		else
			vli_add(r, r, qc, ndigits * 2);
	}
	while (vli_is_negative(r, ndigits * 2))
		vli_add(r, r, m, ndigits * 2);
	while (vli_cmp(r, m, ndigits * 2) >= 0)
		vli_sub(r, r, m, ndigits * 2);

	vli_set(result, r, ndigits);
}

/*
 * Computes result = product % mod, where product is 2N words long.
 * Reference: Ken MacKay's micro-ecc.
 * Currently only designed to work for curve_p or curve_n.
 */
static void vli_mmod_slow(u64 *result, u64 *product, const u64 *mod,
			  unsigned int ndigits)
{
	u64 mod_m[2 * ECC_MAX_DIGITS];
	u64 tmp[2 * ECC_MAX_DIGITS];
	u64 *v[2] = { tmp, product };
	u64 carry = 0;
	unsigned int i;
	/* Shift mod so its highest set bit is at the maximum position. */
	int shift = (ndigits * 2 * 64) - vli_num_bits(mod, ndigits);
	int word_shift = shift / 64;
	int bit_shift = shift % 64;

	vli_clear(mod_m, word_shift);
	if (bit_shift > 0) {
		for (i = 0; i < ndigits; ++i) {
			mod_m[word_shift + i] = (mod[i] << bit_shift) | carry;
			carry = mod[i] >> (64 - bit_shift);
		}
	} else
		vli_set(mod_m + word_shift, mod, ndigits);

	for (i = 1; shift >= 0; --shift) {
		u64 borrow = 0;
		unsigned int j;

		for (j = 0; j < ndigits * 2; ++j) {
			u64 diff = v[i][j] - mod_m[j] - borrow;

			if (diff != v[i][j])
				borrow = (diff > v[i][j]);
			v[1 - i][j] = diff;
		}
		i = !(i ^ borrow); /* Swap the index if there was no borrow */
		vli_rshift1(mod_m, ndigits);
		mod_m[ndigits - 1] |= mod_m[ndigits] << (64 - 1);
		vli_rshift1(mod_m + ndigits, ndigits);
	}
	vli_set(result, v[i], ndigits);
}
#endif

/* Computes result = product % mod using Barrett's reduction with precomputed
 * value mu appended to the mod after ndigits, mu = (2^{2w} / mod) and have
 * length ndigits + 1, where mu * (2^w - 1) should not overflow ndigits
 * boundary.
 *
 * Reference:
 * R. Brent, P. Zimmermann. Modern Computer Arithmetic. 2010.
 * 2.4.1 Barrett's algorithm. Algorithm 2.5.
 */
void vli_mmod_barrett(u64 *result, u64 *product, const u64 *mod,
			     unsigned int ndigits)
{
	vli_mmod_barrett<4>(result, product, mod);
#ifdef	ommit
	u64 q[ECC_MAX_DIGITS * 2];
	u64 r[ECC_MAX_DIGITS * 2];
	const u64 *mu = mod + ndigits;

	vli_mult(q, product + ndigits, mu, ndigits);
	if (mu[ndigits])
		vli_add(q + ndigits, q + ndigits, product + ndigits, ndigits);
#ifdef	ommit
	vli_mult(r, mod, q + ndigits, ndigits);
#else
	// add remain * mod
	vli_set(r, q+ndigits, ndigits);
	q[2*ndigits] = 0;
	vli_umult(q, mu, product[ndigits-1], ndigits);
	if (mu[ndigits])
		vli_uadd(q + ndigits, q + ndigits, product[ndigits-1], ndigits);
	vli_add(result, r, q+ndigits+1, ndigits);
	vli_mult(r, mod, result, ndigits);
#endif
	vli_sub(r, product, r, ndigits * 2);
	while (!vli_is_zero(r + ndigits, ndigits) ||
	       vli_cmp(r, mod, ndigits) != -1) {
		u64 carry;

		carry = vli_sub(r, r, mod, ndigits);
		vli_usub(r + ndigits, r + ndigits, carry, ndigits);
	}
	vli_set(result, r, ndigits);
#endif
}

void vli_div_barrett(u64 *result, u64 *product, const u64 *mu,
			     unsigned int ndigits)
{
	vli_div_barrett<4>(result, product, mu);
#ifdef	ommit
	u64 q[ECC_MAX_DIGITS * 2];
	u64 r[ECC_MAX_DIGITS * 2];

	vli_mult(q, product + ndigits, mu, ndigits);
	if (mu[ndigits])
		vli_add(q + ndigits, q + ndigits, product + ndigits, ndigits);
	vli_set(r, q+ndigits, ndigits);
	q[2*ndigits] = 0;
	vli_umult(q, mu, product[ndigits-1], ndigits);
	if (mu[ndigits])
		vli_uadd(q + ndigits, q + ndigits, product[ndigits-1], ndigits);
	vli_add(result, r, q+ndigits+1, ndigits);
#endif
}

/* Computes p_result = p_product % curve_p.
 * See algorithm 5 and 6 from
 * http://www.isys.uni-klu.ac.at/PDF/2001-0126-MT.pdf
 */
static void vli_mmod_fast_192(u64 *result, const u64 *product,
			      const u64 *curve_prime, u64 *tmp)
{
	int carry;

	vli_set<3>(result, product);

	vli_set<3>(tmp, &product[3]);
	carry = vli_add<3>(result, result, tmp);

	tmp[0] = 0;
	tmp[1] = product[3];
	tmp[2] = product[4];
	carry += vli_add<3>(result, result, tmp);

	tmp[0] = tmp[1] = product[5];
	tmp[2] = 0;
	carry += vli_add<3>(result, result, tmp);

	while (carry || vli_cmp<3>(curve_prime, result) != 1)
		carry -= vli_sub<3>(result, result, curve_prime);
}

/* Computes result = product % curve_prime
 * from http://www.nsa.gov/ia/_files/nist-routines.pdf
 */
static void vli_mmod_fast_256(u64 *result, const u64 *product,
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
	carry = vli_lshift<4>(tmp, tmp, 1);
	carry += vli_add<4>(result, result, tmp);

	/* s2 */
	tmp[1] = product[6] << 32;
	tmp[2] = (product[6] >> 32) | (product[7] << 32);
	tmp[3] = product[7] >> 32;
	carry += vli_lshift<4>(tmp, tmp, 1);
	carry += vli_add<4>(result, result, tmp);

	/* s3 */
	tmp[0] = product[4];
	tmp[1] = product[5] & 0xffffffff;
	tmp[2] = 0;
	tmp[3] = product[7];
	carry += vli_add<4>(result, result, tmp);

	/* s4 */
	tmp[0] = (product[4] >> 32) | (product[5] << 32);
	tmp[1] = (product[5] >> 32) | (product[6] & 0xffffffff00000000ull);
	tmp[2] = product[7];
	tmp[3] = (product[6] >> 32) | (product[4] << 32);
	carry += vli_add<4>(result, result, tmp);

	/* d1 */
	tmp[0] = (product[5] >> 32) | (product[6] << 32);
	tmp[1] = (product[6] >> 32);
	tmp[2] = 0;
	tmp[3] = (product[4] & 0xffffffff) | (product[5] << 32);
	carry -= vli_sub<4>(result, result, tmp);

	/* d2 */
	tmp[0] = product[6];
	tmp[1] = product[7];
	tmp[2] = 0;
	tmp[3] = (product[4] >> 32) | (product[5] & 0xffffffff00000000ull);
	carry -= vli_sub<4>(result, result, tmp);

	/* d3 */
	tmp[0] = (product[6] >> 32) | (product[7] << 32);
	tmp[1] = (product[7] >> 32) | (product[4] << 32);
	tmp[2] = (product[4] >> 32) | (product[5] << 32);
	tmp[3] = (product[6] << 32);
	carry -= vli_sub<4>(result, result, tmp);

	/* d4 */
	tmp[0] = product[7];
	tmp[1] = product[4] & 0xffffffff00000000ull;
	tmp[2] = product[5];
	tmp[3] = product[6] & 0xffffffff00000000ull;
	carry -= vli_sub<4>(result, result, tmp);

	if (carry < 0) {
		do {
			carry += vli_add<4>(result, result, curve_prime);
		} while (carry < 0);
	} else {
		while (carry || vli_cmp<4>(curve_prime, result) != 1)
			carry -= vli_sub<4>(result, result, curve_prime);
	}
}

/* Computes result = product % curve_prime for different curve_primes.
 *
 * Note that curve_primes are distinguished just by heuristic check and
 * not by complete conformance check.
 */
static bool vli_mmod_fast(u64 *result, u64 *product,
			  const u64 *curve_prime, unsigned int ndigits)
{
	u64 tmp[2 * ECC_MAX_DIGITS];

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
		case 8:
			vli_mmod_barrett<8>(result, product, curve_prime);
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

/* Computes result = (left * right) % mod.
 * Assumes that mod is big enough curve order.
 */
#ifdef	ommit
void vli_mod_mult_slow(u64 *result, const u64 *left, const u64 *right,
		       const u64 *mod, unsigned int ndigits)
{
	u64 product[ECC_MAX_DIGITS * 2];

	vli_mult<4>(product, left, right);
	vli_mmod_slow<4>(result, product, mod);
}
#endif

/* Computes result = (left * right) % curve_prime. */
void vli_mod_mult_fast(u64 *result, const u64 *left, const u64 *right,
			      const u64 *curve_prime, unsigned int ndigits)
{
	u64 product[2 * ECC_MAX_DIGITS];

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

#ifdef	ommit
/* Computes result = left^2 % curve_prime. */
static void vli_mod_square_fast(u64 *result, const u64 *left,
				const u64 *curve_prime, unsigned int ndigits)
{
	u64 product[2 * ECC_MAX_DIGITS];

	vli_square(product, left, ndigits);
	vli_mmod_fast(result, product, curve_prime, ndigits);
}

#define EVEN(vli) (!(vli[0] & 1))
/* Computes result = (1 / p_input) % mod. All VLIs are the same size.
 * See "From Euclid's GCD to Montgomery Multiplication to the Great Divide"
 * https://labs.oracle.com/techrep/2001/smli_tr-2001-95.pdf
 */
void vli_mod_inv(u64 *result, const u64 *input, const u64 *mod,
			unsigned int ndigits)
{
	u64 a[ECC_MAX_DIGITS], b[ECC_MAX_DIGITS];
	u64 u[ECC_MAX_DIGITS], v[ECC_MAX_DIGITS];
	u64 carry;
	int cmp_result;

	if (vli_is_zero(input, ndigits)) {
		vli_clear(result, ndigits);
		return;
	}

	vli_set(a, input, ndigits);
	vli_set(b, mod, ndigits);
	vli_clear(u, ndigits);
	u[0] = 1;
	vli_clear(v, ndigits);

	while ((cmp_result = vli_cmp(a, b, ndigits)) != 0) {
		carry = 0;

		if (EVEN(a)) {
			vli_rshift1(a, ndigits);

			if (!EVEN(u))
				carry = vli_add(u, u, mod, ndigits);

			vli_rshift1(u, ndigits);
			if (carry)
				u[ndigits - 1] |= 0x8000000000000000ull;
		} else if (EVEN(b)) {
			vli_rshift1(b, ndigits);

			if (!EVEN(v))
				carry = vli_add(v, v, mod, ndigits);

			vli_rshift1(v, ndigits);
			if (carry)
				v[ndigits - 1] |= 0x8000000000000000ull;
		} else if (cmp_result > 0) {
			vli_sub(a, a, b, ndigits);
			vli_rshift1(a, ndigits);

			if (vli_cmp(u, v, ndigits) < 0)
				vli_add(u, u, mod, ndigits);

			vli_sub(u, u, v, ndigits);
			if (!EVEN(u))
				carry = vli_add(u, u, mod, ndigits);

			vli_rshift1(u, ndigits);
			if (carry)
				u[ndigits - 1] |= 0x8000000000000000ull;
		} else {
			vli_sub(b, b, a, ndigits);
			vli_rshift1(b, ndigits);

			if (vli_cmp(v, u, ndigits) < 0)
				vli_add(v, v, mod, ndigits);

			vli_sub(v, v, u, ndigits);
			if (!EVEN(v))
				carry = vli_add(v, v, mod, ndigits);

			vli_rshift1(v, ndigits);
			if (carry)
				v[ndigits - 1] |= 0x8000000000000000ull;
		}
	}

	vli_set(result, u, ndigits);
}
#endif

void vli_mod_inv(u64 *result, const u64 *input, const u64 *mod,
			unsigned int ndigits)
{
	vli_mod_inv<4>(result, input, mod);
}

void vli_mult(u64 *result, const u64 *left, const u64 *right,
		     unsigned int ndigits)
{
	vli_mult<4>(result, left, right);
}

u64 vli_sub(u64 *result, const u64 *left, const u64 *right,
	    unsigned int ndigits)
{
	return vli_sub<4>(result, left, right);
}

void vli_from_be64(u64 *dest, const void *src, uint ndigits)
{
	vli_from_be64<4>(dest, src);
}

/* ------ Point operations ------ */

/* Returns true if p_point is the point at infinity, false otherwise. */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static bool ecc_point_is_zero(const u64 *p_x, const u64 *p_y)
{
	return (vli_is_zero<ndigits>(p_x) &&
		vli_is_zero<ndigits>(p_y));
}

/* Point multiplication algorithm using Montgomery's ladder with co-Z
 * coordinates. From http://eprint.iacr.org/2011/338.pdf
 */

/* Double in place */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static void ecc_point_double_jacobian(u64 *x1, u64 *y1, u64 *z1,
				      const u64 *curve_prime)
{
	/* t1 = x, t2 = y, t3 = z */
	u64 t4[ECC_MAX_DIGITS];
	u64 t5[ECC_MAX_DIGITS];

	if (vli_is_zero<ndigits>(z1))
		return;

	/* t4 = y1^2 */
	vli_mod_square_fast<ndigits>(t4, y1, curve_prime);
	/* t5 = x1*y1^2 = A */
	vli_mod_mult_fast(t5, x1, t4, curve_prime, ndigits);
	/* t4 = y1^4 */
	vli_mod_square_fast<ndigits>(t4, t4, curve_prime);
	/* t2 = y1*z1 = z3 */
	vli_mod_mult_fast(y1, y1, z1, curve_prime, ndigits);
	/* t3 = z1^2 */
	vli_mod_square_fast<ndigits>(z1, z1, curve_prime);

	/* t1 = x1 + z1^2 */
	vli_mod_add<ndigits>(x1, x1, z1, curve_prime);
	/* t3 = 2*z1^2 */
	vli_mod_add<ndigits>(z1, z1, z1, curve_prime);
	/* t3 = x1 - z1^2 */
	vli_mod_sub<ndigits>(z1, x1, z1, curve_prime);
	/* t1 = x1^2 - z1^4 */
	vli_mod_mult_fast(x1, x1, z1, curve_prime, ndigits);

	/* t3 = 2*(x1^2 - z1^4) */
	vli_mod_add<ndigits>(z1, x1, x1, curve_prime);
	/* t1 = 3*(x1^2 - z1^4) */
	vli_mod_add<ndigits>(x1, x1, z1, curve_prime);
	if (vli_test_bit(x1, 0)) {
		u64 carry = vli_add<ndigits>(x1, x1, curve_prime);

		vli_rshift1<ndigits>(x1);
		x1[ndigits - 1] |= carry << 63;
	} else {
		vli_rshift1<ndigits>(x1);
	}
	/* t1 = 3/2*(x1^2 - z1^4) = B */

	/* t3 = B^2 */
	vli_mod_square_fast<ndigits>(z1, x1, curve_prime);
	/* t3 = B^2 - A */
	vli_mod_sub<ndigits>(z1, z1, t5, curve_prime);
	/* t3 = B^2 - 2A = x3 */
	vli_mod_sub<ndigits>(z1, z1, t5, curve_prime);
	/* t5 = A - x3 */
	vli_mod_sub<ndigits>(t5, t5, z1, curve_prime);
	/* t1 = B * (A - x3) */
	vli_mod_mult_fast(x1, x1, t5, curve_prime, ndigits);
	/* t4 = B * (A - x3) - y1^4 = y3 */
	vli_mod_sub<ndigits>(t4, x1, t4, curve_prime);

	vli_set<ndigits>(x1, z1);
	vli_set<ndigits>(z1, y1);
	vli_set<ndigits>(y1, t4);
}

/* Modify (x1, y1) => (x1 * z^2, y1 * z^3) */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static void apply_z(u64 *x1, u64 *y1, u64 *z, const u64 *curve_prime) noexcept
{
	u64 t1[ECC_MAX_DIGITS];

	vli_mod_square_fast<ndigits>(t1, z, curve_prime);    /* z^2 */
	vli_mod_mult_fast(x1, x1, t1, curve_prime, ndigits); /* x1 * z^2 */
	vli_mod_mult_fast(t1, t1, z, curve_prime, ndigits);  /* z^3 */
	vli_mod_mult_fast(y1, y1, t1, curve_prime, ndigits); /* y1 * z^3 */
}

/* P = (x1, y1) => 2P, (x2, y2) => P' */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static void xycz_initial_double(u64 *x1, u64 *y1, u64 *x2, u64 *y2,
				u64 *p_initial_z, const u64 *curve_prime) noexcept
{
	u64 z[ECC_MAX_DIGITS];

	vli_set<ndigits>(x2, x1);
	vli_set<ndigits>(y2, y1);

	vli_clear<ndigits>(z);
	z[0] = 1;

	if (p_initial_z) {
		vli_set<ndigits>(z, p_initial_z);
		apply_z<ndigits>(x1, y1, z, curve_prime);
	}

	ecc_point_double_jacobian<ndigits>(x1, y1, z, curve_prime);

	apply_z<ndigits>(x2, y2, z, curve_prime);
}

/* Input P = (x1, y1, Z), Q = (x2, y2, Z)
 * Output P' = (x1', y1', Z3), P + Q = (x3, y3, Z3)
 * or P => P', Q => P + Q
 */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static void
xycz_add(u64 *x1, u64 *y1, u64 *x2, u64 *y2, const u64 *curve_prime) noexcept
{
	/* t1 = X1, t2 = Y1, t3 = X2, t4 = Y2 */
	u64 t5[ECC_MAX_DIGITS];

	/* t5 = x2 - x1 */
	vli_mod_sub<ndigits>(t5, x2, x1, curve_prime);
	/* t5 = (x2 - x1)^2 = A */
	vli_mod_square_fast<ndigits>(t5, t5, curve_prime);
	/* t1 = x1*A = B */
	vli_mod_mult_fast(x1, x1, t5, curve_prime, ndigits);
	/* t3 = x2*A = C */
	vli_mod_mult_fast(x2, x2, t5, curve_prime, ndigits);
	/* t4 = y2 - y1 */
	vli_mod_sub<ndigits>(y2, y2, y1, curve_prime);
	/* t5 = (y2 - y1)^2 = D */
	vli_mod_square_fast<ndigits>(t5, y2, curve_prime);

	/* t5 = D - B */
	vli_mod_sub<ndigits>(t5, t5, x1, curve_prime);
	/* t5 = D - B - C = x3 */
	vli_mod_sub<ndigits>(t5, t5, x2, curve_prime);
	/* t3 = C - B */
	vli_mod_sub<ndigits>(x2, x2, x1, curve_prime);
	/* t2 = y1*(C - B) */
	vli_mod_mult_fast(y1, y1, x2, curve_prime, ndigits);
	/* t3 = B - x3 */
	vli_mod_sub<ndigits>(x2, x1, t5, curve_prime);
	/* t4 = (y2 - y1)*(B - x3) */
	vli_mod_mult_fast(y2, y2, x2, curve_prime, ndigits);
	/* t4 = y3 */
	vli_mod_sub<ndigits>(y2, y2, y1, curve_prime);

	vli_set<ndigits>(x2, t5);
}

/* Input P = (x1, y1, Z), Q = (x2, y2, Z)
 * Output P + Q = (x3, y3, Z3), P - Q = (x3', y3', Z3)
 * or P => P - Q, Q => P + Q
 */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static void xycz_add_c(u64 *x1, u64 *y1, u64 *x2, u64 *y2,
				const u64 *curve_prime) noexcept
{
	/* t1 = X1, t2 = Y1, t3 = X2, t4 = Y2 */
	u64 t5[ECC_MAX_DIGITS];
	u64 t6[ECC_MAX_DIGITS];
	u64 t7[ECC_MAX_DIGITS];

	/* t5 = x2 - x1 */
	vli_mod_sub<ndigits>(t5, x2, x1, curve_prime);
	/* t5 = (x2 - x1)^2 = A */
	vli_mod_square_fast<ndigits>(t5, t5, curve_prime);
	/* t1 = x1*A = B */
	vli_mod_mult_fast(x1, x1, t5, curve_prime, ndigits);
	/* t3 = x2*A = C */
	vli_mod_mult_fast(x2, x2, t5, curve_prime, ndigits);
	/* t4 = y2 + y1 */
	vli_mod_add<ndigits>(t5, y2, y1, curve_prime);
	/* t4 = y2 - y1 */
	vli_mod_sub<ndigits>(y2, y2, y1, curve_prime);

	/* t6 = C - B */
	vli_mod_sub<ndigits>(t6, x2, x1, curve_prime);
	/* t2 = y1 * (C - B) */
	vli_mod_mult_fast(y1, y1, t6, curve_prime, ndigits);
	/* t6 = B + C */
	vli_mod_add<ndigits>(t6, x1, x2, curve_prime);
	/* t3 = (y2 - y1)^2 */
	vli_mod_square_fast<ndigits>(x2, y2, curve_prime);
	/* t3 = x3 */
	vli_mod_sub<ndigits>(x2, x2, t6, curve_prime);

	/* t7 = B - x3 */
	vli_mod_sub<ndigits>(t7, x1, x2, curve_prime);
	/* t4 = (y2 - y1)*(B - x3) */
	vli_mod_mult_fast(y2, y2, t7, curve_prime, ndigits);
	/* t4 = y3 */
	vli_mod_sub<ndigits>(y2, y2, y1, curve_prime);

	/* t7 = (y2 + y1)^2 = F */
	vli_mod_square_fast<ndigits>(t7, t5, curve_prime);
	/* t7 = x3' */
	vli_mod_sub<ndigits>(t7, t7, t6, curve_prime);
	/* t6 = x3' - B */
	vli_mod_sub<ndigits>(t6, t7, x1, curve_prime);
	/* t6 = (y2 + y1)*(x3' - B) */
	vli_mod_mult_fast(t6, t6, t5, curve_prime, ndigits);
	/* t2 = y3' */
	vli_mod_sub<ndigits>(y1, t6, y1, curve_prime);

	vli_set<ndigits>(x1, t7);
}

template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static void ecc_point_mult(u64 *result_x, u64 *result_y,
			   const u64 *point_x, const u64 *point_y, const u64 *scalar,
			   u64 *initial_z, const struct ecc_curve *curve) noexcept
{
	/* R0 and R1 */
	u64 rx[2][ECC_MAX_DIGITS];
	u64 ry[2][ECC_MAX_DIGITS];
	u64 z[ECC_MAX_DIGITS];
	u64 sk[2][ECC_MAX_DIGITS];
	const u64 *curve_prime = curve->p;
	int i, nb;
	int num_bits;
	int carry;

	carry = vli_add<ndigits>(sk[0], scalar, curve->n);
	vli_add<ndigits>(sk[1], sk[0], curve->n);
	scalar = sk[!carry];
	num_bits = sizeof(u64) * ndigits * 8 + 1;

	vli_set<ndigits>(rx[1], point_x);
	vli_set<ndigits>(ry[1], point_y);

	xycz_initial_double<ndigits>(rx[1], ry[1], rx[0], ry[0], initial_z,
					curve_prime);

	for (i = num_bits - 2; i > 0; i--) {
		nb = !vli_test_bit(scalar, i);
		xycz_add_c<ndigits>(rx[1 - nb], ry[1 - nb], rx[nb], ry[nb],
						curve_prime);
		xycz_add<ndigits>(rx[nb], ry[nb], rx[1 - nb], ry[1 - nb], curve_prime);
	}

	nb = !vli_test_bit(scalar, 0);
	xycz_add_c<ndigits>(rx[1 - nb], ry[1 - nb], rx[nb], ry[nb], curve_prime);

	/* Find final 1/Z value. */
	/* X1 - X0 */
	vli_mod_sub<ndigits>(z, rx[1], rx[0], curve_prime);
	/* Yb * (X1 - X0) */
	vli_mod_mult_fast(z, z, ry[1 - nb], curve_prime, ndigits);
	/* xP * Yb * (X1 - X0) */
	vli_mod_mult_fast(z, z, point_x, curve_prime, ndigits);

	/* 1 / (xP * Yb * (X1 - X0)) */
	vli_mod_inv<ndigits>(z, z, curve_prime);

	/* yP / (xP * Yb * (X1 - X0)) */
	vli_mod_mult_fast(z, z, point_y, curve_prime, ndigits);
	/* Xb * yP / (xP * Yb * (X1 - X0)) */
	vli_mod_mult_fast(z, z, rx[1 - nb], curve_prime, ndigits);
	/* End 1/Z calculation */

	xycz_add<ndigits>(rx[nb], ry[nb], rx[1 - nb], ry[1 - nb], curve_prime);

	apply_z<ndigits>(rx[0], ry[0], z, curve_prime);

	vli_set<ndigits>(result_x, rx[0]);
	vli_set<ndigits>(result_y, ry[0]);
}

/* Computes R = P + Q mod p */
static void ecc_point_add(u64 *result_x, u64 *result_y,
		   const u64 *p_x, const u64 *p_y, const u64 *q_x, const u64 *q_y,
		   const struct ecc_curve *curve)
{
	u64 z[ECC_MAX_DIGITS];
	u64 px[ECC_MAX_DIGITS];
	u64 py[ECC_MAX_DIGITS];
	unsigned int ndigits = curve->ndigits;

	if (ndigits != 4) return;	// NOOOO
	vli_set<4>((u64 *)result_x, (u64 *)q_x);
	vli_set<4>((u64 *)result_y, (u64 *)q_y);
	vli_mod_sub<4>(z, result_x, p_x, curve->p);
	vli_set<4>(px, p_x);
	vli_set<4>(py, p_y);
	xycz_add<4>(px, py, result_x, result_y, curve->p);
	vli_mod_inv<4>(z, z, curve->p);
	apply_z<4>(result_x, result_y, z, curve->p);
}

static void mont_reduction(u64 *result, const u64 *prod, const u64 *prime,
			const u64 k0) noexcept
{
	u64	t[ECC_MAX_DIGITS * 2];
	u64	s[ECC_MAX_DIGITS * 2];
	u64	r[ECC_MAX_DIGITS * 2];
	vli_clear<8>(r);
	for (uint i=0; i<4; i++) {
		u64	u = (r[0] + prod[i]) * k0;
		vli_umult<4>(s, prime, u);
		vli_uadd<8>(t, s, prod[i]);
		vli_add<8>(r, r, t);
		vli_rshift1w<8>(r);
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
		vli_umult<4>(s, prime, u);
		vli_umult<4>(t, x, y[i]);
		vli_add<8>(t, t, s);
		vli_add<8>(r, r, t);
		vli_rshift1w<8>(r);	
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

#ifdef	ommit
/* Computes R = u1P + u2Q mod p using Shamir's trick.
 * Based on: Kenneth MacKay's micro-ecc (2014).
 */
void ecc_point_mult_shamir(const u64 *result_x, const u64 *result_y,
			   const u64 *u1, const u64 *p_x, const u64 *p_y,
			   const u64 *u2, const u64 *q_x, const u64 *q_y,
			   const struct ecc_curve *curve)
{
	u64 z[ECC_MAX_DIGITS];
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
			u64 tx[ECC_MAX_DIGITS];
			u64 ty[ECC_MAX_DIGITS];
			u64 tz[ECC_MAX_DIGITS];

			vli_set(tx, point->x, ndigits);
			vli_set(ty, point->y, ndigits);
			apply_z(tx, ty, z, curve->p, ndigits);
			vli_mod_sub(tz, rx, tx, curve->p, ndigits);
			xycz_add(tx, ty, rx, ry, curve->p, ndigits);
			vli_mod_mult_fast(z, z, tz, curve->p, ndigits);
		}
	}
	vli_mod_inv(z, z, curve->p, ndigits);
	apply_z(rx, ry, z, curve->p, ndigits);
}
#endif

static inline void ecc_swap_digits(const u64 *in, u64 *out,
				   unsigned int ndigits)
{
	const be64 *src = (be64 *)in;
	uint i;

	for (i = 0; i < ndigits; i++)
		out[i] = be64toh(src[ndigits - 1 - i]);
}

template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static int __ecc_is_key_valid(const struct ecc_curve *curve,
			      const u64 *private_key) noexcept
{
	u64 one[ECC_MAX_DIGITS] = { 1, };
	u64 res[ECC_MAX_DIGITS];

	if (!private_key)
		return -EINVAL;

	if (curve->ndigits != ndigits)
		return -EINVAL;

	/* Make sure the private key is in the range [2, n-3]. */
	if (vli_cmp<ndigits>(one, private_key) != -1)
		return -EINVAL;
	vli_sub<ndigits>(res, curve->n, one);
	vli_sub<ndigits>(res, res, one);
	if (vli_cmp<ndigits>(res, private_key) != 1)
		return -EINVAL;

	return 0;
}

int ecc_is_key_valid(unsigned int curve_id, unsigned int ndigits,
		     const u64 *private_key, unsigned int private_key_len)
{
	uint nbytes;
	const struct ecc_curve *curve = ecc_get_curve(curve_id);

	nbytes = ndigits << ECC_DIGITS_TO_BYTES_SHIFT;

	if (private_key_len != nbytes || ndigits != 4)
		return -EINVAL;

	return __ecc_is_key_valid<4>(curve, private_key);
}

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
	u64 priv[ECC_MAX_DIGITS];
	unsigned int nbytes = ndigits << ECC_DIGITS_TO_BYTES_SHIFT;
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
	u64	pk_x[ECC_MAX_DIGITS];
	u64	pk_y[ECC_MAX_DIGITS];
	u64 priv[ECC_MAX_DIGITS];
	const ecc_curve *curve = ecc_get_curve(curve_id);

	if (!private_key || !curve || ndigits > ARRAY_SIZE(priv)) {
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
	u64 yy[ECC_MAX_DIGITS], xxx[ECC_MAX_DIGITS], w[ECC_MAX_DIGITS];

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

#ifdef	WITH_SYS_RANDOM
int crypto_ecdh_shared_secret(unsigned int curve_id, unsigned int ndigits,
			      const u64 *private_key, const u64 *public_key,
			      u64 *secret)
{
	int ret = 0;
	struct ecc_point *product, *pk;
	u64 priv[ECC_MAX_DIGITS];
	u64 rand_z[ECC_MAX_DIGITS];
	unsigned int nbytes;
	const struct ecc_curve *curve = ecc_get_curve(curve_id);

	if (!private_key || !public_key || !curve || curve->ndigits != 4 ||
	    ndigits > ARRAY_SIZE(priv) || ndigits > ARRAY_SIZE(rand_z)) {
		ret = -EINVAL;
		goto out;
	}

	nbytes = ndigits << ECC_DIGITS_TO_BYTES_SHIFT;

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
