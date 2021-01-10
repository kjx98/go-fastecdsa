// +build curve sm2p

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
#include "curve_defs.hpp"
#include "mont.hpp"

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#pragma GCC pop_options

#ifndef	ARRAY_SIZE
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))
#endif

using namespace vli;
using curve_t = vli::ecc_curve<4>;

static forceinline const curve_t *ecc_get_curve(uint curve_id) noexcept
{
	switch (curve_id) {
#ifdef	ommit
	/* In FIPS mode only allow P256 and higher */
	case ECC_CURVE_SECP256K1:
		return &secp256k1;
	case ECC_CURVE_NIST_P256:
		return &nist_p256;
#endif
	case ECC_CURVE_SM2:
		if (sm2_p256.init()) return &sm2_p256;
		return nullptr;
	default:
		return nullptr;
	}
}


CURVE_HND	get_curve(uint curve_id)
{
	return (CURVE_HND)ecc_get_curve(curve_id);
}


/**
 * get_curve_params	--	get curve params
 * p, n, b, gx, gy	--	bn_t 256 Bits
 */
void	get_curve_params(u64 *p, u64 *n, u64 *b, u64 *gx, u64 *gy,
				CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if ( !(*curve) || curve->ndigits() != 4) return;
	curve->getP(p);
	curve->getN(n);
	curve->getB(b);
	curve->getGx(gx);
	curve->getGy(gy);
}


/* ------ Point operations ------ */

/* Point multiplication algorithm using Montgomery's ladder with co-Z
 * coordinates. From http://eprint.iacr.org/2011/338.pdf
 */

// point_add_jacobian/point_double_jacobian move under ecc_curve class


void    point_double_jacobian(Point *pt, const Point *p, CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	curve->point_double_jacobian(pt->x, pt->y, pt->z, p->x, p->y, p->z);
}

// p->z MUST be one
void    point_double(Point *pt, const Point *p, CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	curve->point_double_jacobian(pt->x, pt->y, pt->z, p->x, p->y, p->z);
	u64 z[4];
	vli_mod_inv<4>(z, pt->z, curve->paramP().data());
	curve->apply_z(pt->x, pt->y, z);
}

void    point_add_jacobian(Point *pt, const Point *p, const Point *q,
                CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	curve->point_add_jacobian(pt->x, pt->y, pt->z, p->x, p->y, p->z,
				q->x, q->y, q->z);
}

// p->z and q->z MUST be one
void    point_add(Point *pt, const Point *p, const Point *q,
                CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	curve->point_add_jacobian(pt->x, pt->y, pt->z, p->x, p->y, p->z,
				q->x, q->y, q->z);
	u64	z[4];
	vli_mod_inv<4>(z, pt->z, curve->paramP().data());
	curve->apply_z(pt->x, pt->y, z);
}

void    affine_from_jacobian(u64 *x, u64 *y, const Point *pt, CURVE_HND curveH)
{
	u64 z[4];
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	vli_mod_inv<4>(z, pt->z, curve->paramP().data());
	vli_set<4>(x, pt->x);
	vli_set<4>(y, pt->y);
	curve->apply_z(x, y, z);
}

#ifdef	ommit
void	point_mult(Point *pt, const Point *p, const u64 *scalar,
				CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	ecc_curve	*curve=(ecc_curve *)curveH;
	if (curve->name[0] == 0 || curve->ndigits != 4) return;
	ecc_point_mult<4>(pt->x, pt->y, p->x, p->y, scalar, nullptr, *curve);
}
#endif


using felem = bn_words_t;

/*-
 * Base point pre computation
 * --------------------------
 *
 * Two different sorts of precomputed tables are used in the following code.
 * Each contain various points on the curve, where each point is three field
 * elements (x, y, z).
 *
 * For the base point table, z is usually 1 (0 for the point at infinity).
 * This table has 2 * 16 elements, starting with the following:
 * index | bits    | point
 * ------+---------+------------------------------
 *     0 | 0 0 0 0 | 0G
 *     1 | 0 0 0 1 | 1G
 *     2 | 0 0 1 0 | 2^64G
 *     3 | 0 0 1 1 | (2^64 + 1)G
 *     4 | 0 1 0 0 | 2^128G
 *     5 | 0 1 0 1 | (2^128 + 1)G
 *     6 | 0 1 1 0 | (2^128 + 2^64)G
 *     7 | 0 1 1 1 | (2^128 + 2^64 + 1)G
 *     8 | 1 0 0 0 | 2^192G
 *     9 | 1 0 0 1 | (2^192 + 1)G
 *    10 | 1 0 1 0 | (2^192 + 2^64)G
 *    11 | 1 0 1 1 | (2^192 + 2^64 + 1)G
 *    12 | 1 1 0 0 | (2^192 + 2^128)G
 *    13 | 1 1 0 1 | (2^192 + 2^128 + 1)G
 *    14 | 1 1 1 0 | (2^192 + 2^128 + 2^64)G
 *    15 | 1 1 1 1 | (2^192 + 2^128 + 2^64 + 1)G
 * followed by a copy of this with each element multiplied by 2^32.
 *
 * The reason for this is so that we can clock bits into four different
 * locations when doing simple scalar multiplies against the base point,
 * and then another four locations using the second 16 elements.
 *
 * Tables for other points have table[i] = iG for i in 0 .. 16. */

/* gmul is the table of precomputed base points */
static const felem gmul[2][16][3] = {
    {{{0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}},
     {{0x715a4589334c74c7, 0x8fe30bbff2660be1, 0x5f9904466a39c994,
       0x32c4ae2c1f198119},
      {0x2df32e52139f0a0, 0xd0a9877cc62a4740, 0x59bdcee36b692153,
       0xbc3736a2f4f6779c},
      {1, 0, 0, 0}},
     {{0xe18bd546b5824517, 0x673891d791caa486, 0xba220b99df9f9a14,
       0x95afbd1155c1da54},
      {0x8e4450eb334acdcb, 0xc3c7d1898a53f20d, 0x2eee750f4053017c,
       0xe8a6d82c517388c2},
      {1, 0, 0, 0}},
     {{0xf81c8da9b99fba55, 0x137f6c6149feef6e, 0xcb129aa494da9ad4,
       0x82a0f5407d123db6},
      {0xfdeca00772c4dbc9, 0xa961b58f0cf58373, 0xecacab94e973f9c3,
       0xf12fa4696a22ca3f},
      {1, 0, 0, 0}},
     {{0xeae3d9a9d13a42ed, 0x2b2308f6484e1b38, 0x3db7b24888c21f3a,
       0xb692e5b574d55da9},
      {0xd186469de295e5ab, 0xdb61ac1773438e6d, 0x5a924f85544926f9,
       0xa175051b0f3fb613},
      {1, 0, 0, 0}},
     {{0xa72d084f62c8d58b, 0xe3d6467deaf48fd7, 0x8fe75e5a128a56a7,
       0xc0023fe7ff2b68bd},
      {0x64f67782316815f9, 0xb52b6d9b19a69cd2, 0x5d1ed6fa89cbbade,
       0x796c910ee7f4ccdb},
      {1, 0, 0, 0}},
     {{0x1b2150c1c5f13015, 0xdaaba91b5d952c9b, 0xe8cc24c3f546142,
       0x75a34b243705f260},
      {0x77d195421cef1339, 0x636644aa0c3a0623, 0x4683df176eeb2444,
       0x642ce3bd3535e74d},
      {1, 0, 0, 0}},
     {{0x4a59ac2c6e7ecc08, 0xaf2b71164f191d63, 0x3622a87fb284554f,
       0xd9eb397b441e9cd0},
      {0xa66b8a4893b6a54d, 0x26fb89a40b4a663a, 0xafa87501eedfc9f4,
       0xf3f000bc66f98108},
      {1, 0, 0, 0}},
     {{0xad8bc68ce031d616, 0x16888d8ee4003187, 0x44c0757f3bb8b600,
       0x793fae7af0164245},
      {0x210cd042973f333b, 0x8666ff52dbd25f9, 0x65c5b129f5f7ad5d,
       0xe03d7a8d19b3219a},
      {1, 0, 0, 0}},
     {{0xd68bfbace0e00392, 0x261014f7d3445dc7, 0xd9f46b2714a071ee,
       0x1b200af30810b682},
      {0xd91d8b12ae69bcd, 0x74a08f17bf8cd981, 0xd822913cf0d2b82d,
       0x248b7af0b05bfad2},
      {1, 0, 0, 0}},
     {{0xba119a049e62f2e2, 0xf278e8a34df05ae5, 0xd269f3564eb5d180,
       0x8e74ad0f4f957cb1},
      {0x112ff4dabd76e2dd, 0x91373f20630fdb7f, 0xf43eab474992904c,
       0x55a5ccc7af3b6db4},
      {1, 0, 0, 0}},
     {{0x5ad104a8bdd23de9, 0xf5a9e515eb71c2c1, 0x390542a0ba95c174,
       0x4c55fb20426491bf},
      {0x91525735ef626289, 0xd2ed977f88f09635, 0xfd48731b7a8a8521,
       0x8f89a03b8fdebea},
      {1, 0, 0, 0}},
     {{0x7e8e61ea35eb8e2e, 0x1bb2700db98a762c, 0xd81ea23b7738c17c,
       0xf9def2a46dba26a3},
      {0x183a7912d05e329f, 0x34664a0896ccde0e, 0x56c22652614283bb,
       0x91692899d5ff0513},
      {1, 0, 0, 0}},
     {{0x449d48d8f3bdbe19, 0xab95de03cc8510cb, 0xaef159463f8bfb25,
       0xda72c379dae3ca8b},
      {0xcba9315ce82cc3ea, 0x4e524bac38a58020, 0x36ba2752538e348c,
       0xb170d0da75ed450f},
      {1, 0, 0, 0}},
     {{0x947af0f52b4f8da6, 0x7eda17d917827976, 0x5ba79a0c705853a0,
       0xa5d9873b3fb2ddc7},
      {0xc2a48162a5fd9ce9, 0x80ee8ae526f25f02, 0xf60c8ef6633be6a9,
       0xe2e23f0229a84a35},
      {1, 0, 0, 0}},
     {{0xbc4945bd86bb6afb, 0x237eb711eba46fee, 0x7c1db58b7b86eb33,
       0xd94eb728273b3ac7},
      {0xbe1717e59568d0a4, 0x4a6067cc45f70212, 0x19b32eb5afc2fb17,
       0xbe3c1e7ac3ac9d3c},
      {1, 0, 0, 0}}},
    {{{0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}},
     {{0x68a88405ae53c1e9, 0x51e46707fd558656, 0x71e834cf86896c10,
       0x3d251b54e10d581f},
      {0x1884d5b0eeb19032, 0xeeaf729853e526fe, 0x5931f6831a8d8c11,
       0x87891d33fb98b4d8},
      {1, 0, 0, 0}},
     {{0x9047673fcac14893, 0xf5df5d83bfb58659, 0xa6230c81642e71a,
       0xef14b33800777791},
      {0xcf1e99afa3386fca, 0x7ace937791313d53, 0x36fe159b6dcd01bb,
       0xc9bc50d02e2b960a},
      {1, 0, 0, 0}},
     {{0x716e5a7ee12e162d, 0xbbf9bb2c62dd5a00, 0xca235ccb4144dd05,
       0xbcb7de0f8f70520e},
      {0x981e8964947cb8eb, 0x53c7102ea04de08d, 0xe9076332afc6a10d,
       0x93d90f776b58c35d},
      {1, 0, 0, 0}},
     {{0x834dbff6678337ee, 0xc607e811fef0785a, 0xaaefc62be30a298b,
       0xeb5ca335326afad3},
      {0x9774fe1384af54a8, 0xca4b6ef5785388b4, 0x1346c82d66f6c642,
       0xedcc0c2aaa2d53ce},
      {1, 0, 0, 0}},
     {{0xb896b3f764b9e6f4, 0x47e4018c736fb3d0, 0xfc2fc86707413920,
       0x1a8526428e1aeae7},
      {0x1386802650e2ae60, 0x7474dedc995384d0, 0x2c4cc396dd43b011,
       0x63b0e9c7141de1b0},
      {1, 0, 0, 0}},
     {{0xeb5fb3b369d17771, 0x1fe07b18933ed257, 0xdfc4c81ce3673912,
       0x913614c66a91a647},
      {0x18aee853c0ba877f, 0x3109c2deceff091, 0x8532307e7e4ee08c,
       0xcef0791a6e6ce0bb},
      {1, 0, 0, 0}},
     {{0xf0e9f5d8057a4a0f, 0xbbf7f8b49f125aa9, 0x51e8fdd6283187c2,
       0xe0997d4759d36298},
      {0x67ec3c5c6f4221c3, 0x3ea275dbc860722f, 0x152d01e23859f5e2,
       0xfb57404312680f44},
      {1, 0, 0, 0}},
     {{0x21ac3df849be2a1f, 0x11006e9fc51d112f, 0x9151aa584775c857,
       0x5159d218ba04a8d9},
      {0x98b7d1a925fd1866, 0x8f4753cafc2ad9d8, 0x8eb91ec1569c05a9,
       0x4abbd1ae27e13f11},
      {1, 0, 0, 0}},
     {{0x616f6644b2c11f4c, 0x251cd7140e540758, 0xf927a40110f02017,
       0x92ff3cc3c1c941b6},
      {0x3249906213f565fe, 0x4633e3ddeb9dbd4e, 0xea9a9d1ec402e6c2,
       0xdc84ce34b14bb7cf},
      {1, 0, 0, 0}},
     {{0xa93e23e5436ff69a, 0x52dcb0a79b63efce, 0x34f6538a9e90cb41,
       0x9cac08f200234bc0},
      {0x6661825b5174a02d, 0x7d4d06de036be57, 0x589d74610ae6bd27,
       0xa296f5577fc91a93},
      {1, 0, 0, 0}},
     {{0x10acefa9d29721d0, 0x8b0f6b8bb5bcd340, 0x921d318c3d86785c,
       0xd6916f3bc16aa378},
      {0x2a0d646a7ad84a0e, 0x7b93256c2fe7e97a, 0x5765e27626479e41,
       0xae9da2272daaced3},
      {1, 0, 0, 0}},
     {{0x56fdc215f7f34ac5, 0xebcb4ff2da3877d3, 0x1eb96792aba6b832,
       0x807ce6bea24741aa},
      {0xff1c10109c721fb4, 0xd187d4bc796353a7, 0x7639ae749af2d303,
       0xaff6d783d56c9286},
      {1, 0, 0, 0}},
     {{0x6002d51b6290dd01, 0xcba3ab0099a836a5, 0x71776611e00d2528,
       0xfaf2cb8c87fce119},
      {0xd445228bdf6882ae, 0xcbbfade17cbce919, 0x837b6335a2eb2453,
       0x11ad7c4b8597f6b6},
      {1, 0, 0, 0}},
     {{0x48de8f368cf2e399, 0x7ae3d25630a74277, 0xdef1a9a6c505323f,
       0xe55f203b4b8d9672},
      {0xc58d8f0d9a1e6e97, 0xe160e6d4b2737a76, 0xd60bd087d47cbdd8,
       0x687d41364d5fef53},
      {1, 0, 0, 0}},
     {{0x83f21bbe056bbf9b, 0x4c2a9d120b4ba5ab, 0xff383d1845b64e4f,
       0x8f13cc8d06dd7867},
      {0xf3a292d8424f0995, 0xfd2546eae7cbe44b, 0x67d14dee6c1e75a3,
       0x53b49e6cc93fb5a8},
      {1, 0, 0, 0}}}
};

#ifdef	ommit
/*
 * select_point selects the |idx|th point from a precomputation table and
 * copies it to out.
 */
static void select_point(const u64 idx, unsigned int size,
                         const felem pre_comp[16][3], felem out[3])
{
    unsigned j;
    u64 *outlimbs = &out[0][0];

#ifdef SM2_NO_CONST_TIME
    const u64 *inlimbs = (u64 *)&pre_comp[idx][0][0];
    for (j = 0; j < NLIMBS * 3; j++) {
        outlimbs[j] = inlimbs[j];
    }
#else
    int i;
    memset(out, 0, sizeof(*out) * 3);

    for (i = 0; i < size; i++) {
        const u64 *inlimbs = (u64 *)&pre_comp[i][0][0];
        u64 mask = i ^ idx;
        mask |= mask >> 4;
        mask |= mask >> 2;
        mask |= mask >> 1;
        mask &= 1;
        mask--;
        for (j = 0; j < NLIMBS * 3; j++)
            outlimbs[j] |= inlimbs[j] & mask;
    }
#endif
}

/*-
 * This function looks at 5+1 scalar bits (5 current, 1 adjacent less
 * significant bit), and recodes them into a signed digit for use in fast point
 * multiplication: the use of signed rather than unsigned digits means that
 * fewer points need to be precomputed, given that point inversion is easy
 * (a precomputed point dP makes -dP available as well).
 *
 * BACKGROUND:
 *
 * Signed digits for multiplication were introduced by Booth ("A signed binary
 * multiplication technique", Quart. Journ. Mech. and Applied Math., vol. IV,
 * pt. 2 (1951), pp. 236-240), in that case for multiplication of integers.
 * Booth's original encoding did not generally improve the density of nonzero
 * digits over the binary representation, and was merely meant to simplify the
 * handling of signed factors given in two's complement; but it has since been
 * shown to be the basis of various signed-digit representations that do have
 * further advantages, including the wNAF, using the following general approach:
 *
 * (1) Given a binary representation
 *
 *       b_k  ...  b_2  b_1  b_0,
 *
 *     of a nonnegative integer (b_k in {0, 1}), rewrite it in digits 0, 1, -1
 *     by using bit-wise subtraction as follows:
 *
 *        b_k b_(k-1)  ...  b_2  b_1  b_0
 *      -     b_k      ...  b_3  b_2  b_1  b_0
 *       -------------------------------------
 *        s_k b_(k-1)  ...  s_3  s_2  s_1  s_0
 *
 *     A left-shift followed by subtraction of the original value yields a new
 *     representation of the same value, using signed bits s_i = b_(i+1) - b_i.
 *     This representation from Booth's paper has since appeared in the
 *     literature under a variety of different names including "reversed binary
 *     form", "alternating greedy expansion", "mutual opposite form", and
 *     "sign-alternating {+-1}-representation".
 *
 *     An interesting property is that among the nonzero bits, values 1 and -1
 *     strictly alternate.
 *
 * (2) Various window schemes can be applied to the Booth representation of
 *     integers: for example, right-to-left sliding windows yield the wNAF
 *     (a signed-digit encoding independently discovered by various researchers
 *     in the 1990s), and left-to-right sliding windows yield a left-to-right
 *     equivalent of the wNAF (independently discovered by various researchers
 *     around 2004).
 *
 * To prevent leaking information through side channels in point multiplication,
 * we need to recode the given integer into a regular pattern: sliding windows
 * as in wNAFs won't do, we need their fixed-window equivalent -- which is a few
 * decades older: we'll be using the so-called "modified Booth encoding" due to
 * MacSorley ("High-speed arithmetic in binary computers", Proc. IRE, vol. 49
 * (1961), pp. 67-91), in a radix-2^5 setting.  That is, we always combine five
 * signed bits into a signed digit:
 *
 *       s_(4j + 4) s_(4j + 3) s_(4j + 2) s_(4j + 1) s_(4j)
 *
 * The sign-alternating property implies that the resulting digit values are
 * integers from -16 to 16.
 *
 * Of course, we don't actually need to compute the signed digits s_i as an
 * intermediate step (that's just a nice way to see how this scheme relates
 * to the wNAF): a direct computation obtains the recoded digit from the
 * six bits b_(4j + 4) ... b_(4j - 1).
 *
 * This function takes those five bits as an integer (0 .. 63), writing the
 * recoded digit to *sign (0 for positive, 1 for negative) and *digit (absolute
 * value, in the range 0 .. 8).  Note that this integer essentially provides the
 * input bits "shifted to the left" by one position: for example, the input to
 * compute the least significant recoded digit, given that there's no bit b_-1,
 * has to be b_4 b_3 b_2 b_1 b_0 0.
 *
 */
static void recode_scalar_bits(unsigned char *sign,
                                     unsigned char *digit, unsigned char in)
{
    unsigned char s, d;

    s = ~((in >> 5) - 1);       /* sets all bits to MSB(in), 'in' seen as
                                 * 6-bit value */
    d = (1 << 6) - in - 1;
    d = (d & s) | (in & ~s);
    d = (d >> 1) + (d & 1);

    *sign = s & 1;
    *digit = d;
}

/*
 * Interleaved point multiplication using precomputed point multiples: The
 * small point multiples 0*P, 1*P, ..., 17*P are in pre_comp[], the scalars
 * in scalars[]. If g_scalar is non-NULL, we also add this multiple of the
 * generator, using certain (large) precomputed multiples in g_pre_comp.
 * Output point (X, Y, Z) is stored in x_out, y_out, z_out
 */
static void batch_mul(felem x_out, felem y_out, felem z_out,
                      const felem *scalar,
                      const unsigned num_points, const felem *g_scalar,
                      const int mixed, const felem pre_comp[][17][3],
                      const felem g_pre_comp[2][16][3])
{
    int i, skip;
    unsigned num, gen_mul = (g_scalar != NULL);
    felem nq[3], ftmp;
    felem tmp[3];
    u64 bits;
    u8 sign, digit;

    /* set nq to the point at infinity */
    memset(nq, 0, sizeof(nq));


    /*
     * Loop over all scalars msb-to-lsb, interleaving additions of multiples
     * of the generator (two in each of the last 32 rounds) and additions of
     * other points multiples (every 5th round).
     */
    skip = 1;                   /* save two point operations in the first
                                 * round */
    for (i = (num_points ? 255 : 31); i >= 0; --i) {
        /* double */
        if (!skip)
            point_double(nq[0], nq[1], nq[2], nq[0], nq[1], nq[2]);

        /* add multiples of the generator */
        if (gen_mul && (i <= 31)) {
            /* first, look 32 bits upwards */
            bits = g_scalar->get_bit(i + 224) << 3;
            bits |= g_scalar->get_bit(i + 160) << 2;
            bits |= g_scalar->get_bit(i + 96) << 1;
            bits |= g_scalar->get_bit(i + 32);
            /* select the point to add, in constant time */
            select_point(bits, 16, g_pre_comp[1], tmp);

            if (!skip) {
                /* Arg 1 below is for "mixed" */
                point_add(nq[0], nq[1], nq[2],
                          nq[0], nq[1], nq[2], 1, tmp[0], tmp[1], tmp[2]);
            } else {
                nq[0] = tmp[0];
                nq[1] = tmp[1];
                nq[2] = tmp[2];
                skip = 0;
            }

            /* second, look at the current position */
            bits = g_scalar->get_bit(i + 192) << 3;
            bits |= g_scalar->get_bit(i + 128) << 2;
            bits |= g_scalar->get_bit(i + 64) << 1;
            bits |= g_scalar->get_bit(i);
            /* select the point to add, in constant time */
            select_point(bits, 16, g_pre_comp[0], tmp);
            /* Arg 1 below is for "mixed" */
            point_add(nq[0], nq[1], nq[2],
                      nq[0], nq[1], nq[2], 1, tmp[0], tmp[1], tmp[2]);
        }

        /* do other additions every 5 doublings */
        if (num_points && (i % 5 == 0)) {
            /* loop over all scalars */
			bits = scalar->get_bits(i-1);
            recode_scalar_bits(&sign, &digit, bits);
                /*
                 * select the point to add or subtract, in constant time
                 */
            select_point(digit, 17, pre_comp[num], tmp);
#ifdef	ommit
            smallfelem_neg(ftmp, tmp[1]); /* (X, -Y, Z) is the negative
                                               * point */
            copy_small_conditional(ftmp, tmp[1], (((limb) sign) - 1));
            felem_contract(tmp[1], ftmp);
#endif

            if (!skip) {
                point_add(nq[0], nq[1], nq[2],
                          nq[0], nq[1], nq[2],
                          mixed, tmp[0], tmp[1], tmp[2]);
            } else {
                nq[0] = tmp[0];
                nq[1] = tmp[1];
                nq[2] = tmp[2];
                skip = 0;
            }
        }
    }
    x_out = nq[0];
    y_out = nq[1];
    z_out = nq[2];
}

/*
 * Computes scalar*generator + \sum scalars[i]*points[i], ignoring NULL
 * values Result is stored in r (r can equal one of the inputs).
 */
int points_mul(Point *r, const felem *scalar, size_t num,
						const Point *points, const u64 (*scalars)[4])
{
    int ret = 0;
    int j;
    int mixed = 0;
    bn_words_t *x, *y, *z, *tmp_scalar;
    felem_bytearray g_secret;
    felem_bytearray *secrets = NULL;
    smallfelem (*pre_comp)[17][3] = NULL;
    smallfelem *tmp_smallfelems = NULL;
    felem_bytearray tmp;
    unsigned i, num_bytes;
    int have_pre_comp = 0;
    size_t num_points = num;
    smallfelem x_in, y_in, z_in;
    felem x_out, y_out, z_out;
    const smallfelem(*g_pre_comp)[16][3] = gmul;
    Point *generator = NULL;
    const Point *p = NULL;
    const bn_words_t *p_scalar = NULL;

    if (scalar != NULL) {
        //if (0 == Point_cmp(generator, group->generator, ctx))
            /* precomputation matches generator */
            have_pre_comp = 1;
    }
    if (num_points > 0) {
        if (num_points >= 3) {
            /*
             * unless we precompute multiples for just one or two points,
             * converting those into affine form is time well spent
             */
			// MUST never HERE
            mixed = 1;
        }
        secrets = malloc(sizeof(*secrets) * num_points);
        pre_comp = malloc(sizeof(*pre_comp) * num_points);
        if (mixed)
            tmp_smallfelems =
              malloc(sizeof(*tmp_smallfelems) * (num_points * 17 + 1));
        if ((secrets == NULL) || (pre_comp == NULL)
            || (mixed && (tmp_smallfelems == NULL))) {
            goto err;
        }

        /*
         * we treat NULL scalars as 0, and NULL points as points at infinity,
         * i.e., they contribute nothing to the linear combination
         */
        memset(secrets, 0, sizeof(*secrets) * num_points);
        memset(pre_comp, 0, sizeof(*pre_comp) * num_points);
        for (i = 0; i < num_points; ++i) {
                /* the i^th point */
            {
                p = points+i;
                p_scalar = scalars + i;
            }
            //flip_endian(secrets[i], tmp, num_bytes);
			memcpy(secrets[i], p_scalar, 32);
            /* precompute multiples */
            smallfelem_expand(x_out, p->x);
            smallfelem_expand(y_out, p->y);
            smallfelem_expand(z_out, p->z);
            felem_shrink(pre_comp[i][1][0], x_out);
            felem_shrink(pre_comp[i][1][1], y_out);
            felem_shrink(pre_comp[i][1][2], z_out);
            for (j = 2; j <= 16; ++j) {
                if (j & 1) {
                    point_add_small(pre_comp[i][j][0], pre_comp[i][j][1],
                                    pre_comp[i][j][2], pre_comp[i][1][0],
                                    pre_comp[i][1][1], pre_comp[i][1][2],
                                    pre_comp[i][j - 1][0],
                                    pre_comp[i][j - 1][1],
                                    pre_comp[i][j - 1][2]);
                } else {
                    point_double_small(pre_comp[i][j][0],
                                       pre_comp[i][j][1],
                                       pre_comp[i][j][2],
                                       pre_comp[i][j / 2][0],
                                       pre_comp[i][j / 2][1],
                                       pre_comp[i][j / 2][2]);
                }
            }
        }
#ifdef	ommit
        if (mixed)
            make_points_affine(num_points * 17, pre_comp[0], tmp_smallfelems);
#endif
    }

    /* the scalar for the generator */
    if ((scalar != NULL) && (have_pre_comp)) {
        memcpy(g_secret, scalar, sizeof(g_secret));
        /* reduce scalar to 0 <= scalar < 2^256 */
        //flip_endian(g_secret, tmp, num_bytes);
        /* do the multiplication with generator precomputation */
        batch_mul(x_out, y_out, z_out,
                  (const felem_bytearray(*))secrets, num_points,
                  g_secret,
                  mixed, (const smallfelem(*)[17][3])pre_comp, g_pre_comp);
    } else
        /* do the multiplication without generator precomputation */
        batch_mul(x_out, y_out, z_out,
                  (const felem_bytearray(*))secrets, num_points,
                  NULL, mixed, (const smallfelem(*)[17][3])pre_comp, NULL);
    /* reduce the output to its unique minimal representation */
	point_get_affine_jacobian(r->x, r->y, x_out, y_out, z_out);
    //felem_contract(r->x, x_out);
    //felem_contract(r->y, y_out);
    //felem_contract(r->z, z_out);
    //ret = EC_POINT_set_Jprojective_coordinates_GFp(group, r, x, y, z, ctx);

 err:
    free(secrets);
    free(pre_comp);
    free(tmp_smallfelems);
    return ret;
}
#endif
