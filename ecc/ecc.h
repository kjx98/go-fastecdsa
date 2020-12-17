/*
 * Copyright (c) 2013, Kenneth MacKay
 * All rights reserved.
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
#pragma once
#ifndef __ECC_H__
#define __ECC_H__

#include "cdefs.h"

/* One digit is u64 qword. */

#ifndef	ECC_CURVE_NIST_P256
/* Curves IDs */
#define ECC_CURVE_NIST_P192	0x0001
#define ECC_CURVE_NIST_P256	0x0002
#define ECC_CURVE_SM2		0x0003
#define	ECC_CURVE_SECP256K1	0x0004
#endif

#ifdef	__cplusplus
extern "C" {
#endif

/**
 * ecc_is_key_valid() - Validate a given ECDH private key
 *
 * @curve_id:		id representing the curve to use
 * @ndigits:		curve's number of digits
 * @private_key:	private key to be used for the given curve
 * @private_key_len:	private key length
 *
 * Returns 0 if the key is acceptable, a negative value otherwise
 */
int ecc_is_key_valid(uint curve_id, unsigned int ndigits,
		     const u64 *private_key, unsigned int private_key_len);

/**
 * ecc_gen_privkey() -  Generates an ECC private key.
 * The private key is a random integer in the range 0 < random < n, where n is a
 * prime that is the order of the cyclic subgroup generated by the distinguished
 * point G.
 * @curve_id:		id representing the curve to use
 * @ndigits:		curve number of digits
 * @private_key:	buffer for storing the generated private key
 *
 * Returns 0 if the private key was generated successfully, a negative value
 * if an error occurred.
 */
int ecc_gen_privkey(uint curve_id, unsigned int ndigits, u64 *privkey);

/**
 * ecc_make_pub_key() - Compute an ECC public key
 *
 * @curve_id:		id representing the curve to use
 * @ndigits:		curve's number of digits
 * @private_key:	pregenerated private key for the given curve
 * @public_key:		buffer for storing the generated public key
 *
 * Returns 0 if the public key was generated successfully, a negative value
 * if an error occurred.
 */
int ecc_make_pub_key(const unsigned int curve_id, unsigned int ndigits,
		     const u64 *private_key, u64 *public_key);

/**
 * crypto_ecdh_shared_secret() - Compute a shared secret
 *
 * @curve_id:		id representing the curve to use
 * @ndigits:		curve's number of digits
 * @private_key:	private key of part A
 * @public_key:		public key of counterpart B
 * @secret:		buffer for storing the calculated shared secret
 *
 * Note: It is recommended that you hash the result of crypto_ecdh_shared_secret
 * before using it for symmetric encryption or HMAC.
 *
 * Returns 0 if the shared secret was generated successfully, a negative value
 * if an error occurred.
 */
int crypto_ecdh_shared_secret(unsigned int curve_id, unsigned int ndigits,
			      const u64 *private_key, const u64 *public_key,
			      u64 *secret);

/**
 * ecc_is_pubkey_valid_partial() - Partial public key validation
 *
 * @curve:		elliptic curve domain parameters
 * @pk:			public key as a point
 *
 * Valdiate public key according to SP800-56A section 5.6.2.3.4 ECC Partial
 * Public-Key Validation Routine.
 *
 * Note: There is no check that the public key is in the correct elliptic curve
 * subgroup.
 *
 * Return: 0 if validation is successful, -EINVAL if validation is failed.
 */
int ecc_is_pubkey_valid_partial(const uint curve_id,
				const u64 *px, const u64 *py);

/**
 * vli_from_be64() - Load vli from big-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 BE values
 * @ndigits:		length of both vli and array
 */
void vli_from_be64(u64 *dest, const void *src, uint ndigits);

/**
 * vli_from_le64() - Load vli from little-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 LE values
 * @ndigits:		length of both vli and array
 */
void vli_from_le64(u64 *dest, const void *src, uint ndigits);

/**
 * vli_mod_inv() - Modular inversion
 *
 * @result:		where to write vli number
 * @input:		vli value to operate on
 * @mod:		modulus
 * @ndigits:		length of all vlis
 */
void vli_mod_inv(u64 *result, const u64 *input, const u64 *mod, uint ndigits);
void vli_mult(u64 *result, const u64 *left, const u64 *right, uint ndigits);

/* Computes result = product % mod using Barrett's reduction with precomputed
 * value mu appended to the mod after ndigits, mu = (2^{2w} / mod) and have
 * length ndigits + 1, where mu * (2^w - 1) should not overflow ndigits
 * boundary.
 *
 * Reference:
 * R. Brent, P. Zimmermann. Modern Computer Arithmetic. 2010.
 * 2.4.1 Barrett's algorithm. Algorithm 2.5.
 */
void vli_mmod_barrett(u64 *result, u64 *product, const u64 *mod, uint ndigits);
void vli_div_barrett(u64 *result, u64 *product, const u64 *mod, uint ndigits);
void mont_MulMod(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 *rr, const u64 k0);
void mont_ExpMod(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 *rr, const u64 k0);

/**
 * vli_mod_mult_slow() - Modular multiplication
 *
 * @result:		where to write result value
 * @left:		vli number to multiply with @right
 * @right:		vli number to multiply with @left
 * @mod:		modulus
 * @ndigits:		length of all vlis
 *
 * Note: Assumes that mod is big enough curve order.
 */
void vli_mod_mult_slow(u64 *result, const u64 *left, const u64 *right,
		       const u64 *mod, unsigned int ndigits);
/* Computes result = (left * right) % curve_prime. */
// SM2 use barret reduction
void vli_mod_mult_fast(u64 *result, const u64 *left, const u64 *right,
			      const u64 *curve_prime, unsigned int ndigits);

/**
 * ecc_point_mult_shamir() - Add two points multiplied by scalars
 *
 * @result:		resulting point
 * @x:			scalar to multiply with @p
 * @p:			point to multiply with @x
 * @y:			scalar to multiply with @q
 * @q:			point to multiply with @y
 * @curve:		curve
 *
 * Returns result = x * p + x * q over the curve.
 * This works faster than two multiplications and addition.
 */
void ecc_point_mult_shamir(const u64 *result_x, const u64 *result_y,
			   const u64 *x, const u64 *px, const u64 *py,
			   const u64 *y, const u64 *qx, const u64 *qy,
			   const uint curve_id);
#ifdef	__cplusplus
}
#endif

#endif	//	__ECC_H__
