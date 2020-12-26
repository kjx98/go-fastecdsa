/* SPDX-License-Identifier: GPL-2.0 */
#pragma once
#ifndef __CURVE_DEFS_HPP__
#define __CURVE_DEFS_HPP__
#include "curve_const.hpp"
#include "ecc_impl.hpp"

#ifdef	ommit
/* NIST P-192: a = p - 3 */
static ecc_curve nist_p192( "nist_192",
	// .gx
	{ 0xF4FF0AFD82FF1012ull, 0x7CBF20EB43A18800ull,
			0x188DA80EB03090F6ull },
	//.gy
	{ 0x73F977A11E794811ull, 0x631011ED6B24CDD5ull,
			0x07192B95FFC8DA78ull },
	//.p
	nist_p192_p,
	//.n
	nist_p192_n,
	//.a
	nist_p192_a,
	//.b
	nist_p192_b,
	//.ndigits
	3);
#endif

/* NIST P-256: a = p - 3 */
static ecc_curve nist_p256( //.name
	"nist_256",
	//.gx
	nist_p256_gx,
	//.gy
	nist_p256_gy,
	//.p
	nist_p256_p,
	//.n
	nist_p256_n,
	//.a
	nist_p256_a,
	//.b
	nist_p256_b
	//.ndigits = 4
	);

/* GM/T 0003.5-2012 SM2: a = p - 3 */
/* prime following mu for Barrett's reduction */
static ecc_curve sm2_p256( //.name
	"sm2p256",
	// .gx
	sm2_p256_gx,
	//.gy
	sm2_p256_gy,
	//.p
	sm2_p256_p,
	//.n
	sm2_p256_n,
	//.a
	sm2_p256_a,
	//.b
	sm2_p256_b,
	// rr_p, rr_n
	sm2_p256_p_rr,
	sm2_p256_n_rr,
	// k0_p, k0_n
	1, 0x327f9e8872350975,
	//.ndigits
	4,
	//.use_barrett
	true);

static ecc_curve secp256k1( //.name
	"secp256k1",
	//.gx
	secp256k1_gx,
	//.gy
	secp256k1_gy,
	//.p
	secp256k1_p,
	//.n
	secp256k1_n,
	//.a
	secp256k1_a,
	//.b
	secp256k1_b
	//.ndigits = 4,
	);

#endif	//__CURVE_DEFS_HPP__
