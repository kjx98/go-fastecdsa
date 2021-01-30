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
#include "ecc.h"
#include "curve_defs.hpp"
#include "mont.hpp"
#include "ecc_key.hpp"

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#pragma GCC pop_options

#ifndef	ARRAY_SIZE
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))
#endif

using namespace vli;
using namespace ecc;
#ifndef	WITH_SM2K256
using curve_t = ecc::ecc_curve<4>;
#else
using curve_t = ecc::curve256<1>;
#endif

static forceinline const curve_t *ecc_get_curve(uint curve_id) noexcept
{
	const curve_t	*ret=nullptr;
	switch (curve_id) {
	/* In FIPS mode only allow P256 and higher */
#ifndef	WITH_SM2K256
#ifdef	ommit
	case ECC_CURVE_SECP256K1:
		ret =  secp256k1;
		break;
#endif
	case ECC_CURVE_NIST_P256:
		ret =  nist_p256;
		break;
	case ECC_CURVE_SM2:
		ret =  sm2_p256p;
		break;
#else
	case ECC_CURVE_SM2:
		sm2_k256.g_precompute();
		ret =  &sm2_k256;
		break;
#endif
	}
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
	bignum<4>	*xx = reinterpret_cast<bignum<4> *>(pt->x);
	bignum<4>	*yy = reinterpret_cast<bignum<4> *>(pt->y);
	bignum<4>	*zz = reinterpret_cast<bignum<4> *>(pt->z);
	curve->apply_z(*xx, *yy, *zz);
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
	bignum<4>	*xx = reinterpret_cast<bignum<4> *>(pt->x);
	bignum<4>	*yy = reinterpret_cast<bignum<4> *>(pt->y);
	bignum<4>	*zz = reinterpret_cast<bignum<4> *>(pt->z);
	curve->apply_z(*xx, *yy, *zz);
}

void    affine_from_jacobian(u64 *x, u64 *y, const Point *pt, CURVE_HND curveH)
{
	u64 z[4];
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	bignum<4>	*xx = reinterpret_cast<bignum<4> *>(const_cast<u64 *>(pt->x));
	bignum<4>	*yy = reinterpret_cast<bignum<4> *>(const_cast<u64 *>(pt->y));
	bignum<4>	*zz = reinterpret_cast<bignum<4> *>(const_cast<u64 *>(pt->z));
	curve->apply_z(*xx, *yy, *zz);
	vli_set<4>(x, pt->x);
	vli_set<4>(y, pt->y);
	vli_set<4>(z, pt->z);
}

void	point_mult(Point *pt, const Point *p, const u64 *scalar,
				CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	point_t<4> *q=reinterpret_cast<point_t<4> *>(pt);
	const point_t<4> *pp = reinterpret_cast<const point_t<4> *>(p);
	const bignum<4>	*sp = reinterpret_cast<const bignum<4> *>(scalar);
	curve->scalar_mult(*q, *pp, *sp);
}

void	point_cmult(Point *pt, const Point *p, const u64 *scalar,
				const u64 *gscalar, CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	if (gscalar == nullptr && scalar == nullptr) return;
	if (scalar == nullptr) {
		point_t<4> *q=reinterpret_cast<point_t<4> *>(pt);
		const bignum<4>	*sp = reinterpret_cast<const bignum<4> *>(gscalar);
		curve->scalar_mult_base(*q, *sp);
	} else {
		point_t<4> *q=reinterpret_cast<point_t<4> *>(pt);
		const point_t<4> *pp = reinterpret_cast<const point_t<4> *>(p);
		const bignum<4>	*sp = reinterpret_cast<const bignum<4> *>(scalar);
		const bignum<4>	*gsp = reinterpret_cast<const bignum<4> *>(gscalar);
		curve->combined_mult(*q, *pp, *sp, *gsp);
	}
}


int	ecc_verify(const u64 *rP, const u64 *sP, const u64 *msgP,
				const Point *pubKey, CURVE_HND curveH)
{
	if (curveH == nullptr) return false;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return false;
	if (rP == nullptr || sP == nullptr || msgP == nullptr || pubKey == nullptr)
		return false;
	bignum<4>	r(rP), s(sP), msg(msgP);
	spoint_t<4>	pk(pubKey->x, pubKey->y);
	return ec_verify(*curve, r, s, pk, msg);
}

void	ecc_sign(u64 *rP, u64 *sP, const u64 *msgP,
				const Point *privKey, CURVE_HND curveH)
{
	if (curveH == nullptr) return;
	auto	*curve=(curve_t *)curveH;
	if (!(*curve) || curve->ndigits() != 4) return;
	if (rP == nullptr || sP == nullptr || msgP == nullptr || privKey == nullptr)
		return;
	bignum<4>	r(rP), s(sP), msg(msgP);
	spoint_t<4>	pk(privKey->x, privKey->y);
	bignum<4>	secr(privKey->z);
	private_key<4>	priv(*curve, secr, pk);
	ec_sign(*curve, r, s, priv, msg);
}
