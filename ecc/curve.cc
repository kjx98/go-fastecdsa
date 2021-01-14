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
using namespace ecc;
using curve_t = ecc::ecc_curve<4>;

static forceinline const curve_t *ecc_get_curve(uint curve_id) noexcept
{
	curve_t	*ret=nullptr;
	switch (curve_id) {
	/* In FIPS mode only allow P256 and higher */
#ifdef	ommit
	case ECC_CURVE_SECP256K1:
		ret =  &secp256k1;
		break;
#endif
	case ECC_CURVE_NIST_P256:
		ret =  &nist_p256;
		break;
	case ECC_CURVE_SM2:
		ret =  &sm2_p256;
		break;
	}
	if (ret != nullptr && !ret->init()) ret = nullptr;
	return ret;
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
