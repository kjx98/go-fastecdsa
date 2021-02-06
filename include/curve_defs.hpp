/* SPDX-License-Identifier: GPL-2.0 */
#pragma once
#ifndef __CURVE_DEFS_HPP__
#define __CURVE_DEFS_HPP__
#include "curve_const.hpp"
#include "curve_impl.hpp"


/* NIST P-256: a = p - 3 */
const auto nist_p256 = ecc::ecc_curve<4>::new_ecc_curve( //.name
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
	);

/* GM/T 0003.5-2012 SM2: a = p - 3 */
/* prime following mu for Barrett's reduction */
//			ecc::build_curve<4, 1, 0x327f9e8872350975>( //.name
#ifdef	ommit
const auto sm2_p256p = ecc::ecc_curve<4>::new_ecc_curve(//.name
	"sm2p256",
	// .gx
	sm2_gx,
	//.gy
	sm2_gy,
	//.p
	sm2_p,
	//.n
	sm2_n,
	//.a
	sm2_a,
	//.b
	sm2_b);
#else
static ecc::ecc_curve<4>	sm2_p256(//.name
	"sm2p256",
	// .gx
	sm2_gx,
	//.gy
	sm2_gy,
	//.p
	sm2_p,
	//.n
	sm2_n,
	//.a
	sm2_a,
	//.b
	sm2_b,
	// rr_p, rr_n
	sm2_p_rr,
	sm2_n_rr,
	// k0_p, k0_n
	1, 0x327f9e8872350975 //
	);
const auto sm2_p256p = &sm2_p256;
#endif

static ecc::curve256 sm2_k256(//.name
	"sm2k256",
	// .gx
	sm2_gx,
	//.gy
	sm2_gy,
	//.p
	sm2_p,
	//.n
	sm2_n,
	//.a
	sm2_a,
	//.b
	sm2_b
	);

const auto secp256k1 = ecc::ecc_curve<4,false>::new_ecc_curve( //.name
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
	//.b, a_is_pminus3
	secp256k1_b
	);

#endif	//__CURVE_DEFS_HPP__
