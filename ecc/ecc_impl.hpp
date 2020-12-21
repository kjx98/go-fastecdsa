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
#ifndef __ECC_IMPL_H__
#define __ECC_IMPL_H__

#include "cdefs.h"

/* One digit is u64 qword. */
#define ECC_MAX_DIGITS             (256 / 64)

/**
 * struct ecc_curve - definition of elliptic curve
 *
 * @name:	Short name of the curve.
 * @g:		Generator point of the curve.
 * @p:		Prime number, if Barrett's reduction is used for this curve
 *		pre-calculated value 'mu' is appended to the @p after ndigits.
 *		Use of Barrett's reduction is heuristically determined in
 *		vli_mmod_fast().
 * @n:		Order of the curve group.
 * @a:		Curve parameter a.
 * @b:		Curve parameter b.
 */
struct ecc_curve {
	ecc_curve(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
				const u64 *_n, const u64 *_a, const u64 *_b, const uint ndig=4,
				const bool bBat=false) : name(_name), gx(_gx), gy(_gy),
			p(_p), n(_n), a(_a), b(_b),
			ndigits(ndig), use_barrett(bBat) {}
	const char *name;
	const u64 *gx;
	const u64 *gy;
	const u64 *p;
	const u64 *n;
	const u64 *a;
	const u64 *b;
	const u64 *rr_pi = nullptr;
	const u64 *rr_n = nullptr;
	const u64 k0_p = 0;
	const u64 k0_n = 0;
	const uint ndigits = 4;
	const bool use_barrett = false;
};

/*
 * Computes result = product % mod
 * for special form moduli: p = 2^k-c, for small c (note the minus sign)
 *
 * References:
 * R. Crandall, C. Pomerance. Prime Numbers: A Computational Perspective.
 * 9 Fast Algorithms for Large-Integer Arithmetic. 9.2.3 Moduli of special form
 * Algorithm 9.2.13 (Fast mod operation for special-form moduli).
 */
//__attribute__((optimize("unroll-loops")))
template<uint ndigits> forceinline
static void
vli_mmod_special(u64 *result, const u64 *product, const u64 *mod) noexcept
{
	u64 c = -mod[0];
	u64 t[ECC_MAX_DIGITS * 2];
	u64 r[ECC_MAX_DIGITS * 2];

	vli_set<ndigits * 2>(r, product);
	while (!vli_is_zero<ndigits>(r + ndigits)) {
		vli_umult<ndigits>(t, r + ndigits, c);
		vli_clear<ndigits>(r + ndigits);
		vli_add_to<ndigits * 2>(r, t);
	}
	vli_set<ndigits>(t, mod);
	vli_clear<ndigits>(t + ndigits);
	while (vli_cmp<ndigits * 2>(r, t) >= 0)
		vli_sub_from<ndigits * 2>(r, t);
	vli_set<ndigits>(result, r);
}

/*
 * Computes result = product % mod
 * for special form moduli: p = 2^{k-1}+c, for small c (note the plus sign)
 * where k-1 does not fit into qword boundary by -1 bit (such as 255).

 * References (loosely based on):
 * A. Menezes, P. van Oorschot, S. Vanstone. Handbook of Applied Cryptography.
 * 14.3.4 Reduction methods for moduli of special form. Algorithm 14.47.
 * URL: http://cacr.uwaterloo.ca/hac/about/chap14.pdf
 *
 * H. Cohen, G. Frey, R. Avanzi, C. Doche, T. Lange, K. Nguyen, F. Vercauteren.
 * Handbook of Elliptic and Hyperelliptic Curve Cryptography.
 * Algorithm 10.25 Fast reduction for special form moduli
 */
template<uint ndigits> forceinline
static void
vli_mmod_special2(u64 *result, const u64 *product, const u64 *mod) noexcept
{
	u64 c2 = mod[0] * 2;
	u64 q[ECC_MAX_DIGITS];
	u64 r[ECC_MAX_DIGITS * 2];
	u64 m[ECC_MAX_DIGITS * 2]; /* expanded mod */
	bool carry; /* last bit that doesn't fit into q */
	int i;

	vli_set<ndigits>(m, mod);
	vli_clear<ndigits>(m + ndigits);

	vli_set<ndigits>(r, product);
	/* q and carry are top bits */
	vli_set<ndigits>(q, product + ndigits);
	vli_clear<ndigits>(r + ndigits);
	carry = vli_is_negative<ndigits>(r);
	if (carry)
		r[ndigits - 1] &= (1ull << 63) - 1;
	for (i = 1; carry || !vli_is_zero<ndigits>(q); i++) {
		u64 qc[ECC_MAX_DIGITS * 2];

		vli_umult<ndigits>(qc, q, c2);
		if (carry)
			vli_uadd_to<ndigits*2>(qc, mod[0]);
		vli_set<ndigits>(q, qc + ndigits);
		vli_clear<ndigits>(qc + ndigits);
		carry = vli_is_negative<ndigits>(qc);
		if (carry)
			qc[ndigits - 1] &= (1ull << 63) - 1;
		if (i & 1)
			vli_sub_from<ndigits*2>(r, qc);
		else
			vli_add_to<ndigits*2>(r, qc);
		//	vli_add<ndigits*2>(r, r, qc);
	}
	while (vli_is_negative<ndigits*2>(r))
		vli_add_to<ndigits*2>(r, m);
	while (vli_cmp<ndigits*2>(r, m) >= 0)
		vli_sub_from<ndigits*2>(r, m);

	vli_set<ndigits>(result, r);
}

/* Computes result = product % mod using Barrett's reduction with precomputed
 * value mu appended to the mod after ndigits, mu = (2^{2w} / mod) and have
 * length ndigits + 1, where mu * (2^w - 1) should not overflow ndigits
 * boundary.
 *
 * Reference:
 * R. Brent, P. Zimmermann. Modern Computer Arithmetic. 2010.
 * 2.4.1 Barrett's algorithm. Algorithm 2.5.
 */
template<uint ndigits> static void forceinline
#ifdef  WITH_C2GO
vli_mmod_barrett(u64 *result, u64 *product, const u64 *mod, u64 *buff) noexcept
#else
vli_mmod_barrett(u64 *result, u64 *product, const u64 *mod) noexcept
#endif
{
#ifdef  WITH_C2GO
	u64	*q = buff;
	u64	*r = buff + ECC_MAX_DIGITS * 2;
#else
	u64 q[ECC_MAX_DIGITS * 2];
	u64 r[ECC_MAX_DIGITS * 2];
#endif
	const u64 *mu = mod + ndigits;

	vli_mult<ndigits>(q, product + ndigits, mu);
	if (mu[ndigits])
		vli_add_to<ndigits>(q + ndigits, product + ndigits);
	// add remain * mod
	vli_set<ndigits>(r, q+ndigits);
	vli_umult<ndigits>(q, mu, product[ndigits-1]);
	if (mu[ndigits])
		vli_uadd_to<ndigits>(q + ndigits, product[ndigits-1]);
	vli_rshift1w<ndigits>(q+ndigits);
	vli_add<ndigits>(result, r, q+ndigits);
	vli_mult<ndigits>(r, mod, result);
	vli_sub<ndigits*2>(r, product, r);
	if (!vli_is_zero<ndigits>(r + ndigits) ||
	       vli_cmp<ndigits>(r, mod) >= 0) {
		vli_sub_from<ndigits>(r, mod);
	}
	vli_set<ndigits>(result, r);
}

template<uint ndigits> forceinline
static void vli_div_barrett(u64 *result, u64 *product, const u64 *mu) noexcept
{
	u64 q[ECC_MAX_DIGITS * 2];
	u64 r[ECC_MAX_DIGITS * 2];

	vli_mult<ndigits>(q, product + ndigits, mu);
	if (mu[ndigits])
		vli_add_to<ndigits>(q + ndigits, product + ndigits);
	vli_set<ndigits>(r, q+ndigits);
	vli_umult<ndigits>(q, mu, product[ndigits-1]);
	if (mu[ndigits])
		vli_uadd_to<ndigits>(q + ndigits, product[ndigits-1]);
	vli_rshift1w<ndigits>(q+ndigits);
	vli_add<ndigits>(result, r, q+ndigits);
}


// SM2 prime optimize
// p is 2^256 - 2^224 - 2^96 + 2^64 -1
forceinline
static void vli_sm2_multP(u64 *result, const u64 u) noexcept
{
	u64	r[ECC_MAX_DIGITS];
	u64	t_low, t_high;
	t_low = u << 32;	// ^192
	t_high = ((u >> 32) & 0xffffffff);
	result[5] = 0;
	result[4] = u - t_high;
	result[4]--;
	result[3] = 0 - t_low;		// ^256 - ^224
	result[2] = 0;
	result[1] = 0;
	result[0] = 0;
	r[0] =0;
	r[3] =0;
	r[1] = t_low;
	r[2] = t_high;
	//vli_sub_from<5>(result, t);	// ^256 - ^224 - ^96
	if (vli_sub_from<4>(result, r)) result[4]--;
	r[2] = 0;
	r[1] = u-1;
	r[0] = -u;		// ^64 -1
	if (vli_add_to<4>(result, r)) result[4]++;
}

forceinline
static void mont_reductionP(u64 *result, const u64 *y, const u64 *prm) noexcept
{
	u64	s[ECC_MAX_DIGITS * 2];
	u64	r[ECC_MAX_DIGITS * 2];
	vli_clear<8>(r);
#pragma GCC unroll 4
	for (uint i=0; i < 4; i++) {
		u64	u = r[0] + y[i];
#ifdef	WITH_SM2_MULTP
		vli_sm2_multP(s, u);
#else
		vli_umult<4>(s, prm, u);
#endif
		vli_uadd_to<8>(r, y[i]);
		vli_add_to<8>(r, s);
		vli_rshift1w<6>(r);	
	}
	if (r[4] !=0 || vli_cmp<4>(r, prm) >= 0) {
		vli_sub<4>(result, r, prm);
	} else vli_set<4>(result, r);
}

forceinline
static void mont_multP(u64 *result, const u64 *x, const u64 *y,
			const u64 *prime) noexcept
{
	u64	t[ECC_MAX_DIGITS * 2];
	u64	s[ECC_MAX_DIGITS * 2];
	u64	r[ECC_MAX_DIGITS * 2];
	vli_clear<8>(r);
#pragma GCC unroll 4
	for (uint i=0; i < 4;i++) {
		u64	u = r[0] + y[i]*x[0];
#ifdef	WITH_SM2_MULTP
		vli_sm2_multP(s, u);
#else
		vli_umult<4>(s, prime, u);
#endif
		vli_umult<4>(t, x, y[i]);
		vli_add_to<8>(r, s);
		vli_add_to<8>(r, t);
		vli_rshift1w<6>(r);	
	}
	if (r[4] != 0 || vli_cmp<4>(r, prime) >= 0) {
		vli_sub<4>(result, r, prime);
	} else vli_set<4>(result, r);
}

template<uint ndigits> forceinline
static void mont_reduction(u64 *result, const u64 *y, const u64 *prime,
			const u64 k0) noexcept
{
	u64	s[ECC_MAX_DIGITS * 2];
	u64	r[ECC_MAX_DIGITS * 2];
	vli_clear<ndigits * 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < ndigits; i++) {
		u64	u = (r[0] + y[i]) * k0;
		vli_umult<ndigits>(s, prime, u);
		vli_uadd_to<ndigits * 2>(r, y[i]);
		vli_add_to<ndigits * 2>(r, s);
		vli_rshift1w<ndigits + 2>(r);	
	}
	if (r[ndigits] !=0 || vli_cmp<ndigits>(r, prime) >= 0) {
		vli_sub<ndigits>(result, r, prime);
	} else vli_set<ndigits>(result, r);
}

template<uint ndigits> forceinline
static void mont_mult(u64 *result, const u64 *x, const u64 *y, const u64 *prime,
				const u64 k0) noexcept
{
	u64	t[ECC_MAX_DIGITS * 2];
	u64	s[ECC_MAX_DIGITS * 2];
	u64	r[ECC_MAX_DIGITS * 2];
	vli_clear<ndigits * 2>(r);
#pragma GCC unroll 4
	for (uint i=0; i < ndigits;i++) {
		u64	u = (r[0] + y[i]*x[0]) * k0;
		vli_umult<ndigits>(s, prime, u);
		vli_umult<ndigits>(t, x, y[i]);
		vli_add_to<ndigits * 2>(r, s);
		vli_add_to<ndigits * 2>(r, t);
		vli_rshift1w<ndigits + 2>(r);	
	}
	if (r[ndigits] != 0 || vli_cmp<ndigits>(r, prime) >= 0) {
		vli_sub<ndigits>(result, r, prime);
	} else vli_set<ndigits>(result, r);
}


/* Computes result = (1 / p_input) % mod. All VLIs are the same size.
 * See "From Euclid's GCD to Montgomery Multiplication to the Great Divide"
 * https://labs.oracle.com/techrep/2001/smli_tr-2001-95.pdf
 */
template<uint ndigits> forceinline
static void vli_mod_inv(u64 *result, const u64 *input, const u64 *mod) noexcept
{
	u64 a[ECC_MAX_DIGITS], b[ECC_MAX_DIGITS];
	u64 u[ECC_MAX_DIGITS], v[ECC_MAX_DIGITS];
	int cmp_result;

	if (vli_is_zero<ndigits>(input)) {
		vli_clear<ndigits>(result);
		return;
	}

	vli_set<ndigits>(a, input);
	vli_set<ndigits>(b, mod);
	vli_clear<ndigits>(u);
	u[0] = 1;
	vli_clear<ndigits>(v);

	while ((cmp_result = vli_cmp<ndigits>(a, b)) != 0) {
		bool carry = false;

		if (vli_is_even(a)) {
			vli_rshift1<ndigits>(a);

			if (!vli_is_even(u))
				carry = vli_add_to<ndigits>(u, mod);
				//carry = vli_add<ndigits>(u, u, mod);

			vli_rshift1<ndigits>(u);
			if (carry)
				u[ndigits - 1] |= 0x8000000000000000ull;
		} else if (vli_is_even(b)) {
			vli_rshift1<ndigits>(b);

			if (!vli_is_even(v))
				carry = vli_add_to<ndigits>(v, mod);
				//carry = vli_add<ndigits>(v, v, mod);

			vli_rshift1<ndigits>(v);
			if (carry)
				v[ndigits - 1] |= 0x8000000000000000ull;
		} else if (cmp_result > 0) {
			vli_sub_from<ndigits>(a, b);
			vli_rshift1<ndigits>(a);

			if (vli_cmp<ndigits>(u, v) < 0)
				vli_add_to<ndigits>(u, mod);

			vli_sub_from<ndigits>(u, v);
			if (!vli_is_even(u))
				carry = vli_add_to<ndigits>(u, mod);
				//carry = vli_add<ndigits>(u, u, mod);

			vli_rshift1<ndigits>(u);
			if (carry)
				u[ndigits - 1] |= 0x8000000000000000ull;
		} else {
			vli_sub_from<ndigits>(b, a);
			vli_rshift1<ndigits>(b);

			if (vli_cmp<ndigits>(v, u) < 0)
				vli_add_to<ndigits>(v, mod);

			vli_sub_from<ndigits>(v, u);
			if (!vli_is_even(v))
				carry = vli_add_to<ndigits>(v, mod);
				//carry = vli_add<ndigits>(v, v, mod);

			vli_rshift1<ndigits>(v);
			if (carry)
				v[ndigits - 1] |= 0x8000000000000000ull;
		}
	}

	vli_set<ndigits>(result, u);
}

#endif	//	__ECC_IMPL_H__
