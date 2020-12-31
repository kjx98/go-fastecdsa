#pragma once
#ifndef __CURVE_CONST_HPP__
#define __CURVE_CONST_HPP__
#include "cdefs.h"


/* NIST P-192: a = p - 3 */
/*
constexpr u64 nist_p192_p[] = { 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFEull,
				0xFFFFFFFFFFFFFFFFull };
constexpr u64 nist_p192_n[] = { 0x146BC9B1B4D22831ull, 0xFFFFFFFF99DEF836ull,
				0xFFFFFFFFFFFFFFFFull };
constexpr u64 nist_p192_a[] = { 0xFFFFFFFFFFFFFFFCull, 0xFFFFFFFFFFFFFFFEull,
				0xFFFFFFFFFFFFFFFFull };
constexpr u64 nist_p192_b[] = { 0xFEB8DEECC146B9B1ull, 0x0FA7E9AB72243049ull,
				0x64210519E59C80E7ull };
*/

/* NIST P-256: a = p - 3 */
constexpr u64 nist_p256_p[] = { 0xFFFFFFFFFFFFFFFFull, 0x00000000FFFFFFFFull,
				0x0000000000000000ull, 0xFFFFFFFF00000001ull };
constexpr u64 nist_p256_n[] = { 0xF3B9CAC2FC632551ull, 0xBCE6FAADA7179E84ull,
				0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFF00000000ull };
constexpr u64 nist_p256_a[] = { 0xFFFFFFFFFFFFFFFCull, 0x00000000FFFFFFFFull,
				0x0000000000000000ull, 0xFFFFFFFF00000001ull };
constexpr u64 nist_p256_b[] = { 0x3BCE3C3E27D2604Bull, 0x651D06B0CC53B0F6ull,
				0xB3EBBD55769886BCull, 0x5AC635D8AA3A93E7ull };
constexpr u64 nist_p256_gx[]= { 0xF4A13945D898C296ull, 0x77037D812DEB33A0ull,
				0xF8BCE6E563A440F2ull, 0x6B17D1F2E12C4247ull };
constexpr u64 nist_p256_gy[]= { 0xCBB6406837BF51F5ull, 0x2BCE33576B315ECEull,
				0x8EE7EB4A7C0F9E16ull, 0x4FE342E2FE1A7F9Bull };

/* GM/T 0003.5-2012 SM2: a = p - 3 */
/* prime following mu for Barrett's reduction */
constexpr u64 sm2_p[] = { 0xFFFFFFFFFFFFFFFFull, 0xffffffff00000000ull,
				0xffffffffffffffffull, 0xfffffffeffffffffull };
constexpr u64 sm2_p_mu[] = {0x200000003, 0x200000002, 0x100000001,
				0x100000001, 1};
constexpr u64 sm2_n[] = { 0x53bbf40939d54123ull, 0x7203df6b21c6052bull,
				0xFFFFFFFFFFFFFFFFull, 0xfffffffeffffffffull };
constexpr u64 sm2_n_mu[] = {0x12AC6361F15149A0ull, 0x8DFC2096FA323C01ull,
				0x100000001ull, 0x100000001ull, 1 };
constexpr bn_words sm2_p_rr = { 0x200000003ull, 0x2ffffffffull,
				0x100000001ull, 0x400000002ull };
constexpr bn_words sm2_n_rr = {0x901192AF7C114F20ull, 0x3464504ADE6FA2FAull,
	   			0x620FC84C3AFFE0D4ull, 0x1EB5E412A22B3D3Bull };
constexpr u64 sm2_p_k0 = 1;
constexpr u64 sm2_n_k0 = 0x327f9e8872350975;
constexpr u64 sm2_a[] = { 0xFFFFFFFFFFFFFFFCull, 0xffffffff00000000ull,
				0xffffffffffffffffull, 0xfffffffeffffffffull };
constexpr u64 sm2_b[] = { 0xDDBCBD414D940E93ull, 0xF39789F515AB8F92ull,
				0x4D5A9E4BCF6509A7ull, 0x28E9FA9E9D9F5E34ull };
constexpr u64 sm2_gx[]= { 0x715A4589334C74C7ull, 0x8FE30BBFF2660BE1ull,
				0x5F9904466A39C994ull, 0x32C4AE2C1F198119ull };
constexpr u64 sm2_gy[]= { 0x2DF32E52139F0A0ull, 0xD0A9877CC62A4740ull,
				0x59BDCEE36B692153ull, 0xBC3736A2F4F6779Cull };

// ECC const for secp256k1, used by BTC/ETH...
constexpr u64 secp256k1_p[] = { 0xFFFFFFFEFFFFFC2Full, 0xFFFFFFFFFFFFFFFFull,
				0xffffffffffffffffull, 0xFFFFFFFFFFFFFFFFull };
constexpr u64 secp256k1_n[] = { 0xbfd25e8cd0364141ull, 0xbaaedce6af48a03bull,
				0xFFFFFFFFFFFFFFFEull, 0xffffffffffffffffull };
constexpr u64 secp256k1_a[] = { 0, 0, 0, 0 };
constexpr u64 secp256k1_b[] = { 0x7, 0, 0, 0};
constexpr u64 secp256k1_gx[]= { 0x59F2815B16F81798ull, 0x29BFCDB2DCE28D9ull,
				0x55A06295CE870B07ull, 0x79BE667EF9DCBBACull };
constexpr u64 secp256k1_gy[]= { 0x9C47D08FFB10D4B8ull, 0xFD17B448A6855419ull,
				0x5DA4FBFC0E1108A8ull, 0x483ADA7726A3C465ull };

#endif	// __CURVE_CONST_HPP__
