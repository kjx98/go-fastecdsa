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


typedef struct {
	void	*data;
	int64_t	len;
	int64_t	cap;
}	GoSlice;


typedef	struct {
	u64	p[4];
	u64	rr[4];
	u64	k0;
} montParams;

typedef	void*	CURVE_HND;


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

#ifdef	ommit
/**
 * vli_from_le64() - Load vli from little-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 LE values
 * @ndigits:		length of both vli and array
 */
void vli_from_le64(u64 *dest, const void *src, uint ndigits);
#endif

/**
 * vli_mod_inv() - Modular inversion
 *
 * @result:		where to write vli number, 256 Bits
 * @input:		vli value to operate on
 * @mod:		modulus
 * @ndigits:		length of all vlis
 */
#ifdef	WITH_C2GO
void vli_mod_inv(u64 *result, const u64 *input, const u64 *mod, u64 *buff);
#else
void vli_mod_inv(u64 *result, const u64 *input, const u64 *mod);
#endif
void vli_mult(u64 *result, const u64 *left, const u64 *right);

/*
 * using Mongtgomey reduction/multiply with precomputed RR and K0
 * RR		b^2n mod P
 * K0		- P^(-1) mod b
 */
void to_montgomery(u64 *res, const u64 *x, const montParams *pa);
void from_montgomery(u64 *res, const u64 *y, const montParams *pa);
void mont_mod_mult(u64 *res, const u64 *x, const u64 *y, const montParams *pa);
void mont_mod_sqr(u64 *res, const u64 *x, const montParams *pa, const u64 n);
void mont_mod_exp(u64 *re, const u64 *x, const u64 *y, const montParams *pa);

/*
 * vli_sm2_mult_p
 * vli_sm2_mult_n
 * mont_sm2_mod_mult_p
 * mont_sm2_mod_mult_n
 */
#ifdef	WITH_C2GO
void vli_sm2_mult_p(GoSlice *result, const u64 u);
#else
void vli_sm2_mult_p(u64 *result, const u64 rLen, const u64 u);
#endif
void mont_sm2_mod_mult_p(u64 *result, const u64 *x, const u64 *y);
void mont_sm2_mod_mult_n(u64 *result, const u64 *x, const u64 *y);
u64 vli_asm_acc();
u64 vli_asm_cmov();
bool bn256_add_to(u64 *left, const u64 *right);
bool bn256_sub_from(u64 *left, const u64 *right);


/**
 * get_curve()		--	get curve defines
 *
 */
CURVE_HND   get_curve(uint curve_id);
/**
 * get_curve_params	--	get curve params
 * p, n, b, gx, gy	--	bn_t 256 Bits
 */
void	get_curve_params(u64 *p, u64 *n, u64 *b, u64 *gx, u64 *gy,
				CURVE_HND curveH);
void	point_double_jacobian(Point *pt, const Point *p, CURVE_HND curveH);
void	point_add_jacobian(Point *pt, const Point *p, const Point *q,
				CURVE_HND curveH);
void	point_double(Point *pt, const Point *p, CURVE_HND curveH);
void	point_add(Point *pt, const Point *p, const Point *q, CURVE_HND curvH);
void	point_mult(Point *pt, const Point *p, const u64 *scalar,
				CURVE_HND curveH, void *scratchBuff);
void	point_cmult(Point *pt, const Point *p, const u64 *scalar,
				const u64 *gscalar, CURVE_HND curveH, void *scratchBuff);
void	affine_from_jacobian(u64 *x, u64 *y, const Point *pt, CURVE_HND curveH);
int		ecc_verify(const u64 *rP, const u64 *sP, const u64 *msgP,
				const Point *pubKey, CURVE_HND curveH);
int		ecc_sign(u64 *rP, u64 *sP, const u64 *msgP, const Point *privKey,
				CURVE_HND curveH);
int		ecc_recover(const u64 *rP, const u64 *sP, const u64 *msgP, const uint v,
				Point *pubKey, CURVE_HND curveH);

#ifdef	__cplusplus
}
#endif

#endif	//	__ECC_H__
