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
	u64 product[2 * ECC_MAX_DIGITS];

	vli_square<ndigits>(product, left);
	vli_mmod_fast(result, product, curve_prime, ndigits);
}

/* Computes result = (left * right) % mod.
 * Assumes that mod is big enough curve order.
 */
#ifdef	ommit
void forceinline
vli_mod_mult_slow(u64 *result, const u64 *left, const u64 *right,
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

/* Returns true if p_point is the point at infinity, false otherwise. */
template<uint ndigits> forceinline
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
		u64 carry = vli_add_to<ndigits>(x1, curve_prime);

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
template<uint ndigits> forceinline
static void ecc_point_add(u64 *result_x, u64 *result_y,
		   const u64 *p_x, const u64 *p_y, const u64 *q_x, const u64 *q_y,
		   const struct ecc_curve *curve)
{
	u64 z[ECC_MAX_DIGITS];
	u64 px[ECC_MAX_DIGITS];
	u64 py[ECC_MAX_DIGITS];

	if (ndigits != curve->ndigits) return;	// NOOOO
	vli_set<ndigits>((u64 *)result_x, (u64 *)q_x);
	vli_set<ndigits>((u64 *)result_y, (u64 *)q_y);
	vli_mod_sub<ndigits>(z, result_x, p_x, curve->p);
	vli_set<ndigits>(px, p_x);
	vli_set<ndigits>(py, p_y);
	xycz_add<ndigits>(px, py, result_x, result_y, curve->p);
	vli_mod_inv<ndigits>(z, z, curve->p);
	apply_z<ndigits>(result_x, result_y, z, curve->p);
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
#ifndef  WITH_C2GO
void vli_div_barrett(u64 *result, u64 *product, const u64 *mu)
{
	vli_div_barrett<4>(result, product, mu);
}
#endif

#ifdef	WITH_SHAMIR
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

forceinline static
void ecc_swap_digits(const u64 *in, u64 *out, uint ndigits)
{
	const be64 *src = (be64 *)in;
	uint i;
#pragma GCC unroll 4
	for (i = 0; i < ndigits; i++)
		out[i] = be64toh(src[ndigits - 1 - i]);
}

template<uint ndigits> forceinline
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
	u64	pk_x[ECC_MAX_DIGITS];
	u64	pk_y[ECC_MAX_DIGITS];
	u64 priv[ECC_MAX_DIGITS];
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
#endif

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
