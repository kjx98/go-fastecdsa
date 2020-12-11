/* SPDX-License-Identifier: GPL-2.0 */
#ifndef _CRYTO_ECC_CURVE_DEFS_H
#define _CRYTO_ECC_CURVE_DEFS_H

/* NIST P-192: a = p - 3 */
static u64 nist_p192_p[] = { 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFEull,
				0xFFFFFFFFFFFFFFFFull };
static u64 nist_p192_n[] = { 0x146BC9B1B4D22831ull, 0xFFFFFFFF99DEF836ull,
				0xFFFFFFFFFFFFFFFFull };
static u64 nist_p192_a[] = { 0xFFFFFFFFFFFFFFFCull, 0xFFFFFFFFFFFFFFFEull,
				0xFFFFFFFFFFFFFFFFull };
static u64 nist_p192_b[] = { 0xFEB8DEECC146B9B1ull, 0x0FA7E9AB72243049ull,
				0x64210519E59C80E7ull };
static struct ecc_curve nist_p192 = {
	.name = "nist_192",
	.g = {
		.x = { 0xF4FF0AFD82FF1012ull, 0x7CBF20EB43A18800ull,
			0x188DA80EB03090F6ull },
		.y = { 0x73F977A11E794811ull, 0x631011ED6B24CDD5ull,
			0x07192B95FFC8DA78ull },
	},
	.ndigits = 3,
	.p = nist_p192_p,
	.n = nist_p192_n,
	.a = nist_p192_a,
	.b = nist_p192_b
};

/* NIST P-256: a = p - 3 */
static u64 nist_p256_p[] = { 0xFFFFFFFFFFFFFFFFull, 0x00000000FFFFFFFFull,
				0x0000000000000000ull, 0xFFFFFFFF00000001ull };
static u64 nist_p256_n[] = { 0xF3B9CAC2FC632551ull, 0xBCE6FAADA7179E84ull,
				0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFF00000000ull };
static u64 nist_p256_a[] = { 0xFFFFFFFFFFFFFFFCull, 0x00000000FFFFFFFFull,
				0x0000000000000000ull, 0xFFFFFFFF00000001ull };
static u64 nist_p256_b[] = { 0x3BCE3C3E27D2604Bull, 0x651D06B0CC53B0F6ull,
				0xB3EBBD55769886BCull, 0x5AC635D8AA3A93E7ull };
static struct ecc_curve nist_p256 = {
	.name = "nist_256",
	.g = {
		.x = { 0xF4A13945D898C296ull, 0x77037D812DEB33A0ull,
			0xF8BCE6E563A440F2ull, 0x6B17D1F2E12C4247ull },
		.y = { 0xCBB6406837BF51F5ull, 0x2BCE33576B315ECEull,
			0x8EE7EB4A7C0F9E16ull, 0x4FE342E2FE1A7F9Bull },
	},
	.ndigits = 4,
	.p = nist_p256_p,
	.n = nist_p256_n,
	.a = nist_p256_a,
	.b = nist_p256_b
};

/* GM/T 0003.5-2012 SM2: a = p - 3 */
/* prime following mu for Barrett's reduction */
static u64 sm2_p256_p[] = { 0xFFFFFFFFFFFFFFFFull, 0xffffffff00000000ull,
				0xffffffffffffffffull, 0xfffffffeffffffffull,
				8589934595, 8589934594, 4294967297, 4294967297, 1};
static u64 sm2_p256_n[] = { 0x53bbf40939d54123ull, 0x7203df6b21c6052bull,
				0xFFFFFFFFFFFFFFFFull, 0xfffffffeffffffffull };
static u64 sm2_p256_a[] = { 0xFFFFFFFFFFFFFFFCull, 0xffffffff00000000ull,
				0xffffffffffffffffull, 0xfffffffeffffffffull };
static u64 sm2_p256_b[] = { 0xDDBCBD414D940E93ull, 0xF39789F515AB8F92ull,
				0x4D5A9E4BCF6509A7ull, 0x28E9FA9E9D9F5E34ull };
static struct ecc_curve sm2_p256 = {
	.name = "sm2p256",
	.g = {
		.x = { 0x715A4589334C74C7ull, 0x8FE30BBFF2660BE1ull,
			0x5F9904466A39C994ull, 0x32C4AE2C1F198119ull },
		.y = { 0x2DF32E52139F0A0ull, 0xD0A9877CC62A4740ull,
			0x59BDCEE36B692153ull, 0xBC3736A2F4F6779Cull },
	},
	.ndigits = 4,
	.p = sm2_p256_p,
	.n = sm2_p256_n,
	.a = sm2_p256_a,
	.b = sm2_p256_b
};

#endif
