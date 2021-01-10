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
	curve->apply_z(pt->x, pt->y, pt->z);
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
	curve->apply_z(pt->x, pt->y, pt->z);
}

void    affine_from_jacobian(u64 *x, u64 *y, const Point *pt, CURVE_HND curveH)
{
	u64 z[4];
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	vli_set<4>(x, pt->x);
	vli_set<4>(y, pt->y);
	vli_set<4>(z, pt->z);
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


#ifdef	ommit
using felem = bn_words_t;
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
