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
#ifndef	__ECC_IMPL_H__
#define __ECC_IMPL_H__

#include "cdefs.h"
#include <string>
#include <functional>
#include "vli.hpp"
#include "vli_bn.hpp"


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
template<const uint N> forceinline
static void
vli_mmod_special(u64 *result, const u64 *product, const u64 *mod) noexcept
{
	u64 c = -mod[0];
	u64 t[N * 2];
	u64 r[N * 2];

	vli_set<N * 2>(r, product);
	while (!vli_is_zero<N>(r + N)) {
		vli_umult<N>(t, r + N, c);
		vli_clear<N>(r + N);
		vli_add_to<N * 2>(r, t);
	}
	vli_set<N>(t, mod);
	vli_clear<N>(t + N);
	while (vli_cmp<N * 2>(r, t) >= 0)
		vli_sub_from<N * 2>(r, t);
	vli_set<N>(result, r);
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
template<const uint N> forceinline
static void
vli_mmod_special2(u64 *result, const u64 *product, const u64 *mod) noexcept
{
	u64 c2 = mod[0] * 2;
	u64 q[N];
	u64 r[N * 2];
	u64 m[N * 2]; /* expanded mod */
	bool carry; /* last bit that doesn't fit into q */
	int i;

	vli_set<N>(m, mod);
	vli_clear<N>(m + N);

	vli_set<N>(r, product);
	/* q and carry are top bits */
	vli_set<N>(q, product + N);
	vli_clear<N>(r + N);
	carry = vli_is_negative<N>(r);
	if (carry)
		r[N - 1] &= (1ull << 63) - 1;
	for (i = 1; carry || !vli_is_zero<N>(q); i++) {
		u64 qc[N * 2];

		vli_umult<N>(qc, q, c2);
		if (carry)
			vli_uadd_to<N*2>(qc, mod[0]);
		vli_set<N>(q, qc + N);
		vli_clear<N>(qc + N);
		carry = vli_is_negative<N>(qc);
		if (carry)
			qc[N - 1] &= (1ull << 63) - 1;
		if (i & 1)
			vli_sub_from<N*2>(r, qc);
		else
			vli_add_to<N*2>(r, qc);
		//	vli_add<N*2>(r, r, qc);
	}
	while (vli_is_negative<N*2>(r))
		vli_add_to<N*2>(r, m);
	while (vli_cmp<N*2>(r, m) >= 0)
		vli_sub_from<N*2>(r, m);

	vli_set<N>(result, r);
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
template<const uint ndigits> static void forceinline
#ifdef  WITH_C2GO
vli_mmod_barrett(u64 *result, const u64 *product, const u64 *mod, u64 *buff) noexcept
#else
vli_mmod_barrett(u64 *result, const u64 *product, const u64 *mod) noexcept
#endif
{
#ifdef  WITH_C2GO
	u64	*q = buff;
	u64	*r = buff + ndigits * 2;
#else
	u64 q[ndigits * 2];
	u64 r[ndigits];
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
	vli_mult<ndigits>(q, mod, result);
	vli_sub<ndigits*2>(q, product, q);
	if (!vli_is_zero<ndigits>(q + ndigits) ||
	       vli_cmp<ndigits>(q, mod) >= 0) {
		vli_sub_from<ndigits>(q, mod);
	}
	vli_set<ndigits>(result, q);
}

template<const uint N> forceinline
static void
vli_div_barrett(u64 *result, const u64 *product, const u64 *mu) noexcept
{
	u64 q[N * 2];
	u64 r[N];

	vli_mult<N>(q, product + N, mu);
	if (mu[N])
		vli_add_to<N>(q + N, product + N);
	vli_set<N>(r, q+N);
	vli_umult<N>(q, mu, product[N-1]);
	if (mu[N])
		vli_uadd_to<N>(q + N, product[N-1]);
	vli_rshift1w<N>(q+N);
	vli_add<N>(result, r, q+N);
}

/* Computes result = (1 / p_input) % mod. All VLIs are the same size.
 * See "From Euclid's GCD to Montgomery Multiplication to the Great Divide"
 * https://labs.oracle.com/techrep/2001/smli_tr-2001-95.pdf
 */
template<const uint N> forceinline
static void
#ifdef	WITH_C2GO
vli_mod_inv(u64 *result, const u64 *input, const u64 *mod, u64 *buff) noexcept
#else
vli_mod_inv(u64 *result, const u64 *input, const u64 *mod) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*a=buff;
	u64 *b=a+N;
	u64	*u=b+N;
	u64	*v=u+N;
#else
	u64 a[N], b[N];
	u64 u[N], v[N];
#endif
	int cmp_result;

	if (vli_is_zero<N>(input)) {
		vli_clear<N>(result);
		return;
	}

	vli_set<N>(a, input);
	vli_set<N>(b, mod);
	vli_clear<N>(u);
	u[0] = 1;
	vli_clear<N>(v);

	while ((cmp_result = vli_cmp<N>(a, b)) != 0) {
		bool carry = false;

		if (vli_is_even(a)) {
			vli_rshift1<N>(a);

			if (!vli_is_even(u))
				carry = vli_add_to<N>(u, mod);

			vli_rshift1<N>(u);
			if (carry)
				u[N - 1] |= 0x8000000000000000ull;
		} else if (vli_is_even(b)) {
			vli_rshift1<N>(b);

			if (!vli_is_even(v))
				carry = vli_add_to<N>(v, mod);

			vli_rshift1<N>(v);
			if (carry)
				v[N - 1] |= 0x8000000000000000ull;
		} else if (cmp_result > 0) {
			vli_sub_from<N>(a, b);
			vli_rshift1<N>(a);

			if (vli_cmp<N>(u, v) < 0)
				vli_add_to<N>(u, mod);

			vli_sub_from<N>(u, v);
			if (!vli_is_even(u))
				carry = vli_add_to<N>(u, mod);

			vli_rshift1<N>(u);
			if (carry)
				u[N - 1] |= 0x8000000000000000ull;
		} else {
			vli_sub_from<N>(b, a);
			vli_rshift1<N>(b);

			if (vli_cmp<N>(v, u) < 0)
				vli_add_to<N>(v, mod);

			vli_sub_from<N>(v, u);
			if (!vli_is_even(v))
				carry = vli_add_to<N>(v, mod);

			vli_rshift1<N>(v);
			if (carry)
				v[N - 1] |= 0x8000000000000000ull;
		}
	}

	vli_set<N>(result, u);
}

/* Computes result = (1 / p_input) % mod. All VLIs are the same size.
 * The binary extended gcd algorithm was first described by Knuth
 */
// not work yet
#ifndef	ommit
template<const uint N> forceinline
static void
#ifdef	WITH_C2GO
vli_mod_inv_new(u64 *result, const u64 *n, const u64 *mod, u64 *buff) noexcept
#else
vli_mod_inv_new(u64 *result, const u64 *n, const u64 *mod) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*a=buff;
	u64 *b=a+N;
	u64	*x=b+N;
	u64	*y=x+N;
#else
	u64 a[N], b[N];
	u64 x[N], y[N];
#endif
	// mod should be prime, >3 MUST BE odd
	if ( vli_is_even(mod) ) {
		vli_clear<N>(result);
		return;
	}

	vli_set<N>(a, mod);
	vli_set<N>(b, n);
	vli_clear<N>(x);
	vli_clear<N>(y);
	x[0] = 1;
	int cmp_result;
	int	xc=0, yc=0;
	int bc = 0, ac = 0;

	while ( !vli_is_zero<N>(b) ) {

		while (vli_is_even(b)) {
			vli_rshift1<N>(b);
			if ( !vli_is_even(x) ) {
				if (vli_add_to<N>(x, mod)) xc++;
			}
			vli_rshift1<N>(x);
			if (xc & 1) x[N-1] |= (1L << 63);
			xc >>= 1;
		}

		while (vli_is_even(a)) {
			vli_rshift1<N>(a);
			if ( !vli_is_even(y) ) {
				if (vli_add_to<N>(y, mod)) yc++;
			}
			vli_rshift1<N>(y);
			if (yc & 1) y[N-1] |= (1L << 63);
			yc >>= 1;
		}

		// a,b may be neg
		if (bc) {
			cmp_result = (ac)?vli_cmp<N>(a, b):-1;
		} else {
			cmp_result = (ac)?1:vli_cmp<N>(b, a);
		}
		if (cmp_result >= 0) {
			xc += yc;
			if (vli_add_to<N>(x, y)) xc++;
			bc -= ac;
			if (vli_sub_from<N>(b, a)) bc--;
		} else {
			yc += xc;
			if (vli_add_to<N>(y, x)) yc++;
			ac -= bc;
			if (vli_sub_from<N>(a, b)) ac--;
		}
	}
	if (likely(vli_is_one<N>(a))) {
		if (yc) vli_add_to<N>(y, mod);
		if (vli_cmp<N>(y, mod) >= 0) {
			vli_sub<N>(result, y, mod);
		} else vli_set<N>(result, y);
	} else {
		// no inverse
		vli_clear<N>(result);
	}
}
#else
// x is mod, prime
template<const uint N> forceinline
static void
vli_mod_inv_new(u64 *result, const u64 *y, const u64 *x) noexcept
{
	u64 a[N], b[N];
	u64 c[N], d[N];
	u64 u[N], v[N];
	// mod should be prime, >3 MUST BE odd
	if ( vli_is_even(x) ) {
		vli_clear<N>(result);
		return;
	}

	vli_set<N>(u, x);
	vli_set<N>(v, y);
	vli_clear<N>(a);
	vli_clear<N>(b);
	vli_clear<N>(c);
	vli_clear<N>(d);
	a[0] = 1;
	d[0] = 1;
	//int cmp_result;
	int	b_carry=0, d_carry = 0;

	while ( !vli_is_zero<N>(u) ) {
		bool carry;

		while (vli_is_even(u)) {
			vli_rshift1<N>(u);
			if (vli_is_even(a) && vli_is_even(b)) {
				vli_rshift1<N>(a);
				vli_rshift1<N>(b);
			} else {
				carry = vli_add_to<N>(a, y);
				vli_rshift1<N>(a);
				if (carry) a[N-1] |= (1L << 63);
				if (vli_sub_from<N>(b, x)) b_carry--;
				vli_rshift1<N>(b);
				b[N-1] |= ((u64)(b_carry & 1) << 63);
				b_carry >>= 1;
			}
		}

		while (vli_is_even(v)) {
			vli_rshift1<N>(v);
			if (vli_is_even(c) && vli_is_even(d)) {
				vli_rshift1<N>(c);
				vli_rshift1<N>(d);
			} else {
				carry = vli_add_to<N>(c, y);
				vli_rshift1<N>(c);
				if (carry) c[N-1] |= (1L << 63);
				if (vli_sub_from<N>(d, x)) d_carry--;
				vli_rshift1<N>(d);
				d[N-1] |= ((u64)(d_carry & 1) << 63);
				d_carry >>= 1;
			}
		}

		if (vli_cmp<N>(u, v) >= 0) {
			vli_sub_from<N>(u, v);
			vli_sub_from<N>(a, c);
			b_carry -= d_carry;
			if (vli_sub_from<N>(b, d)) b_carry--;
		} else {
			vli_sub_from<N>(v, u);
			vli_sub_from<N>(c, a);
			d_carry -= b_carry;
			if (vli_sub_from<N>(d, b)) d_carry--;
		}
	}
	vli_set<N>(a, c);
	vli_set<N>(b, d);
	if (likely(vli_is_one<N>(v))) {
		if (d_carry) vli_add<N>(result, d, x);
		if (unlikely(vli_cmp<N>(d, x) >= 0)) {
			vli_sub<N>(result, d, x);
		} else vli_set<N>(result, d);
	} else {
		// no inverse
		vli_clear<N>(result);
	}
}
#endif


#endif	//	__ECC_IMPL_H__
