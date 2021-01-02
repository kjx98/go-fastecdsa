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


// vli_bn requires builtin_usubl_overflow...
namespace vli {

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
template<const uint N=4, const u64 k0_p=0, const u64 k0_n=0>
class ecc_curve {
public:
	ecc_curve(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
			const u64 *_n, const u64 *_a, const u64 *_b, const bool a_n3 = true,
			const bool bBat=false) :
		name(_name), gx(_gx), gy(_gy),
		p(_p), n(_n), a(_a), b(_b), _a_is_neg3(a_n3),
		use_barrett(bBat)
	{
		static_assert(k0_p == 0 && k0_n == 0, "No use montgomery mult/reduction");
	}
	ecc_curve(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
			const u64 *_n, const u64 *_a, const u64 *_b, const u64* rrP,
			const u64* rrN, const bool a_n3 = true) :
		name(_name), gx(_gx), gy(_gy),
		p(_p), n(_n), a(_a), b(_b), rr_p(rrP), rr_n(rrN), _a_is_neg3(a_n3),
		use_barrett(false)
	{
		static_assert(k0_p != 0 && k0_n != 0, "Do use montgomery mult/reduction");
	}
	const uint ndigits() const { return N; }
	explicit operator bool() const noexcept {
		//return name != nullptr && name[0] != 0;
		return name != "";
	}
	void getP(u64 *v) { p.set(v); }
	void getN(u64 *v) { n.set(v); }
	void getA(u64 *v) { a.set(v); }
	void getB(u64 *v) { b.set(v); }
	void getGx(u64 *v) { gx.set(v); }
	void getGy(u64 *v) { gy.set(v); }
	const bignum<N>& paramA() const noexcept { return a; }
	const bignum<N>& paramP() const noexcept { return p; }
	bool init(const u64 *muP=nullptr, const u64 *muN=nullptr) noexcept {
		if (_inited) return _inited;
#ifdef	WITH_BARRETT
		if (muP == nullptr || muN == nullptr) {
			use_barrett = false;
			return false;
		}
		mu_p = bignum<N+1>(muP);
		mu_n = bignum<N+1>(muN);
		use_barrett = true;
		_inited = true;
		return _inited;
#else
#if	__cplusplus >= 201703L
		if constexpr(k0_p == 0)
#else
		if (k0_p == 0)
#endif
		{
			return false;
		}
		// should be calc K0 and RR
		return true;
#endif
	}
	const bool a_is_pminus3() const noexcept { return _a_is_neg3; }
#ifdef	WITH_BARRETT
	void mmod_barrett(bignum<N>& res, const bn_prod<N>& prod) noexcept
	{
		if (use_barrett) {
			// res = barrettmod
			prod.mmod_barrett(res, p, mu_p);
		}
		return;
	}
#endif
	const bignum<N>& mont_one() const noexcept { return _mont_one; }
	const bignum<N>& mont_inv2() const noexcept { return _mont_inv2; }
	void to_montgomery(bignum<N>& res, const u64 *x) noexcept
	{
		res.mont_mult<k0_p>(x, rr_p, p);
	}
	void to_montgomery(bignum<N>& res, const bignum<N>& x) noexcept
	{
		res.mont_mult<k0_p>(x, rr_p, p);
	}
	void from_montgomery(bignum<N>& res, const bignum<N>& y) noexcept
	{
		res.mont_reduction<k0_p>(y, p);
	}
	void from_montgomery(u64* result, const bignum<N>& y) noexcept
	{
		bignum<N>   *res = reinterpret_cast<bignum<N> *>(result);
		res.mont_reduction<k0_p>(y, p);
	}
	void
	mont_mult(bignum<N>& res, const bignum<N>& left, const bignum<N>& right)
	noexcept {
#if	__cplusplus >= 201703L
		if constexpr(k0_p != 0)
#else
		if (k0_p != 0)
#endif
		{
			bignum<N> xp;
			bignum<N> yp;
			bignum<N> r;
			xp.mont_mult<k0_p>(left, rr_p, p);
			yp.mont_mult<k0_p>(right, rr_p, p);
			r.mont_mult<k0_p>(xp, yp, p);
			res.mont_reduction<k0_p>(r, p);
			return;
		}
	}
	void mont_sqr(bignum<N>& res, const bignum<N> left, const uint nTimes=1)
	noexcept {
#if	__cplusplus >= 201703L
		if constexpr(k0_p != 0)
#else
		if (k0_p != 0)
#endif
		{
			bignum<N> xp;
			bignum<N> r;
			xp.mont_mult<k0_p>(left, rr_p, p);
			for (uint i=0; i < nTimes; i++) r.mont_sqr<k0_p>(xp, p);
			res.mont_reduction<k0_p>(r, p);
			return;
		}
	}
	// left,right less than p
	void mod_add(bignum<N>& res, const bignum<N>& left, const bignum<N>& right)
	noexcept {
		if (res.add(left, right)) {
			res.sub_from(p);
		}
	}
	// left,right less than p
	void mod_add_to(bignum<N>& res, const bignum<N>& right) noexcept
	{
		if (res.add_to(right)) {
			res.sub_from(p);
		}
	}
	// left,right less than p
	void mod_sub(bignum<N>& res, const bignum<N>& left, const bignum<N>& right)
	noexcept {
		if (res.sub(left, right)) {
			res.add_to(p);
		}
	}
	void mod_sub_from(bignum<N>& res, const bignum<N>& right) noexcept
	{
		if (res.sub_from(right)) {
			res.add_to(p);
		}
	}
	void mont_mult2(bignum<N>& res, const bignum<N>& left) noexcept
	{
#ifdef	ommit
		this->mod_add(res, left, left);
#else
		if (res.lshift1(left) != 0 || res.cmp(p) >= 0) {
			res.sub_from(p);
		}
#endif
	}
private:
	const std::string name;
	const bignum<N> gx;
	const bignum<N> gy;
	const bignum<N> p;
	const bignum<N> n;
	const bignum<N> a;
	const bignum<N> b;
	const bignum<N> rr_p;
	const bignum<N> rr_n;
	const bignum<N+1> mu_p;
	const bignum<N+1> mu_n;
	const bignum<N> _mont_one;
	const bignum<N> _mont_inv2;
	const uint _ndigits = N;
	const bool _a_is_neg3;
	const bool use_barrett = false;
	bool _inited = false;
};

template<uint ndigits>
struct point_t {
	u64		x[ndigits];
	u64		y[ndigits];
	u64		z[ndigits];
};


struct slice_t {
	u64	*data;
	int64_t	len;
	int64_t	cap;
	bool isZero() {
		if (len <= 0) return true;
		for (int i=0;i<len;++i) {
			if (data[i] != 0) return false;
		}
		return true;
	}
	explicit operator bool () { return len != 0; }
	template <uint ndigits>slice_t(u64 vd[ndigits]) noexcept : data(vd),
			len(ndigits), cap(ndigits)
	{
		vli_clear<ndigits>(data);
	}
	void normal() {
		if (len <= 0) return;
		int	i;
		for(i=len-1;i>=0;i--) {
			if (data[i] != 0) break;
		}
		len = i+1;
	}
};

}	// namespace vli


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
	u64 t[ndigits * 2];
	u64 r[ndigits * 2];

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
	u64 q[ndigits];
	u64 r[ndigits * 2];
	u64 m[ndigits * 2]; /* expanded mod */
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
		u64 qc[ndigits * 2];

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

template<uint ndigits> forceinline
static void
vli_div_barrett(u64 *result, const u64 *product, const u64 *mu) noexcept
{
	u64 q[ndigits * 2];
	u64 r[ndigits];

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

/* Computes result = (1 / p_input) % mod. All VLIs are the same size.
 * See "From Euclid's GCD to Montgomery Multiplication to the Great Divide"
 * https://labs.oracle.com/techrep/2001/smli_tr-2001-95.pdf
 */
template<uint ndigits> forceinline
static void
#ifdef	WITH_C2GO
vli_mod_inv(u64 *result, const u64 *input, const u64 *mod, u64 *buff) noexcept
#else
vli_mod_inv(u64 *result, const u64 *input, const u64 *mod) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*a=buff;
	u64 *b=a+ndigits;
	u64	*u=b+ndigits;
	u64	*v=u+ndigits;
#else
	u64 a[ndigits], b[ndigits];
	u64 u[ndigits], v[ndigits];
#endif
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

			vli_rshift1<ndigits>(u);
			if (carry)
				u[ndigits - 1] |= 0x8000000000000000ull;
		} else if (vli_is_even(b)) {
			vli_rshift1<ndigits>(b);

			if (!vli_is_even(v))
				carry = vli_add_to<ndigits>(v, mod);

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

			vli_rshift1<ndigits>(v);
			if (carry)
				v[ndigits - 1] |= 0x8000000000000000ull;
		}
	}

	vli_set<ndigits>(result, u);
}

/* Computes result = (1 / p_input) % mod. All VLIs are the same size.
 * The binary extended gcd algorithm was first described by Knuth
 */
// not work yet
#ifndef	ommit
template<uint ndigits> forceinline
static void
#ifdef	WITH_C2GO
vli_mod_inv_new(u64 *result, const u64 *n, const u64 *mod, u64 *buff) noexcept
#else
vli_mod_inv_new(u64 *result, const u64 *n, const u64 *mod) noexcept
#endif
{
#ifdef	WITH_C2GO
	u64	*a=buff;
	u64 *b=a+ndigits;
	u64	*x=b+ndigits;
	u64	*y=x+ndigits;
#else
	u64 a[ndigits], b[ndigits];
	u64 x[ndigits], y[ndigits];
#endif
	// mod should be prime, >3 MUST BE odd
	if ( vli_is_even(mod) ) {
		vli_clear<ndigits>(result);
		return;
	}

	vli_set<ndigits>(a, mod);
	vli_set<ndigits>(b, n);
	vli_clear<ndigits>(x);
	vli_clear<ndigits>(y);
	x[0] = 1;
	int cmp_result;
	int	xc=0, yc=0;
	int bc = 0, ac = 0;

	while ( !vli_is_zero<ndigits>(b) ) {

		while (vli_is_even(b)) {
			vli_rshift1<ndigits>(b);
			if ( !vli_is_even(x) ) {
				if (vli_add_to<ndigits>(x, mod)) xc++;
			}
			vli_rshift1<ndigits>(x);
			if (xc & 1) x[ndigits-1] |= (1L << 63);
			xc >>= 1;
		}

		while (vli_is_even(a)) {
			vli_rshift1<ndigits>(a);
			if ( !vli_is_even(y) ) {
				if (vli_add_to<ndigits>(y, mod)) yc++;
			}
			vli_rshift1<ndigits>(y);
			if (yc & 1) y[ndigits-1] |= (1L << 63);
			yc >>= 1;
		}

		// a,b may be neg
		if (bc) {
			cmp_result = (ac)?vli_cmp<ndigits>(a, b):-1;
		} else {
			cmp_result = (ac)?1:vli_cmp<ndigits>(b, a);
		}
		if (cmp_result >= 0) {
			xc += yc;
			if (vli_add_to<ndigits>(x, y)) xc++;
			bc -= ac;
			if (vli_sub_from<ndigits>(b, a)) bc--;
		} else {
			yc += xc;
			if (vli_add_to<ndigits>(y, x)) yc++;
			ac -= bc;
			if (vli_sub_from<ndigits>(a, b)) ac--;
		}
	}
	if (likely(vli_is_one<ndigits>(a))) {
		if (yc) vli_add_to<ndigits>(y, mod);
		if (vli_cmp<ndigits>(y, mod) >= 0) {
			vli_sub<ndigits>(result, y, mod);
		} else vli_set<ndigits>(result, y);
	} else {
		// no inverse
		vli_clear<ndigits>(result);
	}
}
#else
// x is mod, prime
template<uint ndigits> forceinline
static void
vli_mod_inv_new(u64 *result, const u64 *y, const u64 *x) noexcept
{
	u64 a[ndigits], b[ndigits];
	u64 c[ndigits], d[ndigits];
	u64 u[ndigits], v[ndigits];
	// mod should be prime, >3 MUST BE odd
	if ( vli_is_even(x) ) {
		vli_clear<ndigits>(result);
		return;
	}

	vli_set<ndigits>(u, x);
	vli_set<ndigits>(v, y);
	vli_clear<ndigits>(a);
	vli_clear<ndigits>(b);
	vli_clear<ndigits>(c);
	vli_clear<ndigits>(d);
	a[0] = 1;
	d[0] = 1;
	//int cmp_result;
	int	b_carry=0, d_carry = 0;

	while ( !vli_is_zero<ndigits>(u) ) {
		bool carry;

		while (vli_is_even(u)) {
			vli_rshift1<ndigits>(u);
			if (vli_is_even(a) && vli_is_even(b)) {
				vli_rshift1<ndigits>(a);
				vli_rshift1<ndigits>(b);
			} else {
				carry = vli_add_to<ndigits>(a, y);
				vli_rshift1<ndigits>(a);
				if (carry) a[ndigits-1] |= (1L << 63);
				if (vli_sub_from<ndigits>(b, x)) b_carry--;
				vli_rshift1<ndigits>(b);
				b[ndigits-1] |= ((u64)(b_carry & 1) << 63);
				b_carry >>= 1;
			}
		}

		while (vli_is_even(v)) {
			vli_rshift1<ndigits>(v);
			if (vli_is_even(c) && vli_is_even(d)) {
				vli_rshift1<ndigits>(c);
				vli_rshift1<ndigits>(d);
			} else {
				carry = vli_add_to<ndigits>(c, y);
				vli_rshift1<ndigits>(c);
				if (carry) c[ndigits-1] |= (1L << 63);
				if (vli_sub_from<ndigits>(d, x)) d_carry--;
				vli_rshift1<ndigits>(d);
				d[ndigits-1] |= ((u64)(d_carry & 1) << 63);
				d_carry >>= 1;
			}
		}

		if (vli_cmp<ndigits>(u, v) >= 0) {
			vli_sub_from<ndigits>(u, v);
			vli_sub_from<ndigits>(a, c);
			b_carry -= d_carry;
			if (vli_sub_from<ndigits>(b, d)) b_carry--;
		} else {
			vli_sub_from<ndigits>(v, u);
			vli_sub_from<ndigits>(c, a);
			d_carry -= b_carry;
			if (vli_sub_from<ndigits>(d, b)) d_carry--;
		}
	}
	vli_set<ndigits>(a, c);
	vli_set<ndigits>(b, d);
	if (likely(vli_is_one<ndigits>(v))) {
		if (d_carry) vli_add<ndigits>(result, d, x);
		if (unlikely(vli_cmp<ndigits>(d, x) >= 0)) {
			vli_sub<ndigits>(result, d, x);
		} else vli_set<ndigits>(result, d);
	} else {
		// no inverse
		vli_clear<ndigits>(result);
	}
}
#endif


#endif	//	__ECC_IMPL_H__
