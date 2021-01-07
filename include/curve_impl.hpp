/*
 * Copyright(c) 2020 Jesse Kuang  <jkuang@21cn.com>
 * All rights reserved.
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
#ifndef __CURVE_IMPL_HPP__
#define __CURVE_IMPL_HPP__

#include "vli.hpp"
#include "vli_bn.hpp"
#include "ecc_impl.hpp"


namespace vli {

/**
 * struct ecc_curve - definition of elliptic curve
 *
 * @name:	Short name of the curve.
 * @g:		Generator point of the curve.
 * @p:		Prime number, if Barrett's reduction is used for this curve
 * @n:		Order of the curve group.
 * @a:		Curve parameter a.
 * @b:		Curve parameter b.
 */
template<const uint N=4, const u64 Pk0=0>
class ecc_curve {
public:
	using mont1_func=std::function<void(bignum<N>&, const bignum<N>&)>;
	using mont2_func=std::function<void(bignum<N>&, const bignum<N>&, const bignum<N>&)>;
	ecc_curve(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
			const u64 *_n, const u64 *_a, const u64 *_b, const bool a_n3 = true,
			const bool bBat=false) :
		name(_name), gx(_gx), gy(_gy),
		p(_p), n(_n), a(_a), b(_b), _a_is_neg3(a_n3),
		use_barrett(bBat)
	{
		static_assert(N > 3, "curve only support 256Bits or more");
	}
	ecc_curve(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
			const u64 *_n, const u64 *_a, const u64 *_b, const u64* rrP,
			const u64* rrN, const u64 k0P, const u64 k0N,
			const bool a_n3 = true):
		name(_name), gx(_gx), gy(_gy),
		p(_p), n(_n), a(_a), b(_b), rr_p(rrP), rr_n(rrN), k0_p(k0P), k0_n(k0N),
		_a_is_neg3(a_n3), use_barrett(false)
	{
		static_assert(N > 3, "curve only support 256Bits or more");
	}
	ecc_curve(ecc_curve &&) = default;
	const uint ndigits() const { return N; }
	explicit operator bool() const noexcept {
		//return name != nullptr && name[0] != 0;
		return name != "";
	}
	void getP(u64 *v) const noexcept { p.set(v); }
	void getN(u64 *v) const noexcept { n.set(v); }
	void getA(u64 *v) const noexcept { a.set(v); }
	void getB(u64 *v) const noexcept { b.set(v); }
	void getGx(u64 *v) const noexcept { gx.set(v); }
	void getGy(u64 *v) const noexcept { gy.set(v); }
	const bignum<N>& montParamA() const noexcept { return _mont_a; }
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
		if (unlikely(k0_p == 0))
		{
			// no calc k0 and rr after instantiation
			return false;
		}
		bignum<N>	t1;
		t1.clear();
		_mont_one.sub(t1, p);
		if (unlikely( !_a_is_neg3 )) {
			if (likely(a.is_zero())) _a_is_zero = true;
		}
		// should verify calc K0 and RR
		_inited = true;
		return true;
#endif
	}
	const bool a_is_pminus3() const noexcept { return _a_is_neg3; }
	const bool a_is_zero() const noexcept { return _a_is_zero; }
#ifdef	WITH_BARRETT
	void mmod_barrett(bignum<N>& res, const bn_prod<N>& prod) const noexcept
	{
		if (use_barrett) {
			// res = barrettmod
			prod.mmod_barrett(res, p, mu_p);
		}
		return;
	}
#endif
	const bignum<N>& mont_one() const noexcept { return _mont_one; }
	void to_montgomery(bignum<N>& res, const u64 *x) const noexcept
	{
		bignum<N>   *xx = reinterpret_cast<bignum<N> *>(const_cast<u64 *>(x));
#if	__cplusplus >= 201703L
		if constexpr(Pk0 != 0) {
			mont_mult<Pk0>(res, *xx, rr_p, p);
		} else _mont_mult(res, *xx, rr_p);
#else
		res.mont_mult(*xx, rr_p, p, k0_p);
#endif
	}
	void to_montgomery(bignum<N>& res, const bignum<N>& x) const noexcept
	{
#if	__cplusplus >= 201703L
		if constexpr(Pk0 != 0) {
			mont_mult<Pk0>(res, x, rr_p, p);
		} else  _mont_mult(res, x, rr_p);
#else
		res.mont_mult(x, rr_p, p, k0_p);
#endif
	}
	void from_montgomery(bignum<N>& res, const bignum<N>& y) const noexcept
	{
#if	__cplusplus >= 201703L
		if constexpr(Pk0 != 0) {
			mont_reduction<PK0>(res, y, p);
		} else _mont_reduction(res, y);
#else
		res.mont_reduction(y, p, k0_p);
#endif
	}
	void from_montgomery(u64* result, const bignum<N>& y) const noexcept
	{
		bignum<N>   *res = reinterpret_cast<bignum<N> *>(result);
#if	__cplusplus >= 201703L
		if constexpr(Pk0 != 0) {
			mont_reduction<Pk0>(*res, y, p);
		} else _mont_reduction(*res, y);
#else
		res->mont_reduction(y, p, k0_p);
#endif
	}
	void
	mont_mmult(bignum<N>& res, const bignum<N>& left, const bignum<N>& right)
	const noexcept {
#if	__cplusplus >= 201703L
		if constexpr(Pk0 != 0) {
			mont_mult<Pk0>(res, left, right, p);
		} else _mont_mult(res, left, right);
#else
		res.mont_mult(left, right, p, k0_p);
#endif
	}
	void mont_msqr(bignum<N>& res, const bignum<N> left, const uint nTimes=1)
	const noexcept {
#if	__cplusplus >= 201703L
		if constexpr(Pk0 != 0) {
			mont_sqr<Pk0>(res, left, p);
			for (uint i=1; i < nTimes; i++) mont_sqr<Pk0>(res, res, p);
		} else {
			_mont_sqr(res, left);
			for (uint i=1; i < nTimes; i++) _mont_sqr(res, res);
		}
#else 
		res.mont_sqr(left, p, k0_p);
		for (uint i=1; i < nTimes; i++) res.mont_sqr(res, p, k0_p);
#endif
	}
	// left,right less than p
	void mod_add(bignum<N>& res, const bignum<N>& left, const bignum<N>& right)
	const noexcept {
		if (res.add(left, right)) {
			res.sub_from(p);
		}
	}
	// left,right less than p
	void mod_add_to(bignum<N>& res, const bignum<N>& right) const noexcept
	{
		if (res.add_to(right)) {
			res.sub_from(p);
		}
	}
	// left,right less than p
	void mod_sub(bignum<N>& res, const bignum<N>& left, const bignum<N>& right)
	const noexcept {
		if (res.sub(left, right)) {
			res.add_to(p);
		}
	}
	void mod_sub_from(bignum<N>& res, const bignum<N>& right) const noexcept
	{
		if (res.sub_from(right)) {
			res.add_to(p);
		}
	}
	void mont_mult2(bignum<N>& res, const bignum<N>& left) const noexcept
	{
#ifdef	ommit
		this->mod_add(res, left, left);
#else
		if (res.lshift1(left) != 0 || res.cmp(p) >= 0) {
			res.sub_from(p);
		}
#endif
	}
/* dbl-1998-cmo-2 algorithm
 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html
 * M ... l1
 * S ... l2
 * T ... x3
 */
	void point_double_jacobian(u64 *x3, u64 *y3, u64 *z3, const u64 *x1,
					const u64 *y1, const u64 *z1 = nullptr) const noexcept
	{
		if ((z1 != nullptr && vli_is_zero<N>(z1)) || vli_is_zero<N>(y1)) {
			/* P_y == 0 || P_z == 0 => [1:1:0] */
			vli_clear<N>(x3);
			vli_clear<N>(y3);
			vli_clear<N>(z3);
			x3[0] = 1;
			y3[0] = 1;
			return;
		}
		bool	z_is_one = (z1 == nullptr || vli_is_one<N>(z1));
		bignum<N>	t1, t2, l1, l2;
		bignum<N>	xp, yp, zp;
		bignum<N>	*x3p = reinterpret_cast<bignum<N> *>(x3);
		bignum<N>	*y3p = reinterpret_cast<bignum<N> *>(y3);
		bignum<N>	*z3p = reinterpret_cast<bignum<N> *>(z3);
		to_montgomery(xp, x1);
		to_montgomery(yp, y1);
		if (z_is_one) {
			zp = mont_one();
		} else {
			to_montgomery(zp, z1);
		}
		if (a_is_pminus3()) {
			/* Use the faster case.  */
			/* L1 = 3(X - Z^2)(X + Z^2) */
			/*                          T1: used for Z^2. */
			/*                          T2: used for the right term. */
			if (z_is_one) {
				// l1 = X - Z^2
				mod_sub(l1, xp, mont_one());
				// 3(X - Z^2) = 2(X - Z^2) + (X - Z^2)
				// t1 = 2 * l1
				mont_mult2(t1, l1);
				// l1 = 2 * l1 + l1 = 3(X - Z^2)
				mod_add_to(l1, t1);
				// t1 = X + Z^2
				mod_add(t1, xp, mont_one());
				// l1 = 3(X - Z^2)(X + Z^2)
				this->mont_mmult(l1, l1, t1);
			} else {
				// t1 = Z^2
				this->mont_msqr(t1, zp);
				// l1 = X - Z^2
				mod_sub(l1, xp, t1);
				// 3(X - Z^2) = 2(X - Z^2) + (X - Z^2)
				// t2 = 2 * l1
				mont_mult2(t2, l1);
				// l1 = l1 + 2 * l1 = 3(X - Z^2)
				mod_add_to(l1, t2);
				// t2 = X + Z^2
				mod_add(t2, xp, t1);
				// l1 = 3(X - Z^2)(X + Z^2)
				this->mont_mmult(l1, l1, t2);
			}
		} else {
			/* Standard case. */
			/* L1 = 3X^2 + aZ^4 */
			/*                          T1: used for aZ^4. */
			// l1 = X^2
			this->mont_msqr(l1, xp);
			mont_mult2(t1, l1);
			// l1 = 3X^2
			mod_add_to(l1, t1);
			if (a_is_zero()) {
				/* Use the faster case.  */
				/* L1 = 3X^2 */
				// do nothing
			} else if (z_is_one) {
				// should be mont_paramA
				mod_add_to(l1, montParamA());
			} else {
				// t1 = Z^4
				this->mont_msqr(t1, zp, 2);
				// t1 = a * Z^4
				this->mont_mmult(t1, t1, montParamA());
				// l1 = 3 X^2 + a Z^4
				mod_add_to(l1, t1);
			}
		}

		/* Z3 = 2YZ */
		if (z_is_one) {
			mont_mult2(*z3p, yp);
		} else {
			// Z3 = YZ
			this->mont_mmult(*z3p, yp, zp);
			// Z3 *= 2
			mont_mult2(*z3p, *z3p);
		}

		/* L2 = 4XY^2 */
		/* t2 = Y^2 */
		this->mont_msqr(t2, yp);
		// t2 = 2 Y^2
		mont_mult2(t2, t2);
		// l2 =  2 XY^2
		this->mont_mmult(l2, t2, xp);
		// l2 = 4 X Y^2
		mont_mult2(l2, l2);

		/* X3 = L1^2 - 2L2 */
		/*                              T1: used for 2L2. */
		this->mont_msqr(*x3p, l1);
		mont_mult2(t1, l2);
		mod_sub_from(*x3p, t1);

		/* L3 = 8Y^4 */
		/*   L3 reuse t2, t2: taken from above. */
		this->mont_msqr(t2, t2);		// t2 = t2^2 = 4Y^4
		mont_mult2(t2, t2);	// t2 *= 2, t2 = 8Y^4

		/* Y3 = L1(L2 - X3) - L3 */
		mod_sub(*y3p, l2, *x3p);
		this->mont_mmult(*y3p, l1, *y3p);
		mod_sub_from(*y3p, t2);

		// montgomery reduction
		from_montgomery(x3, *x3p);
		from_montgomery(y3, *y3p);
		from_montgomery(z3, *z3p);
	}
	void apply_z(u64 *x1, u64 *y1, u64 *z) const noexcept
	{
		bignum<N>	t1;

		if (vli_is_one<N>(z) || vli_is_zero<N>(z)) return; 
		bignum<N>	*xp = reinterpret_cast<bignum<N> *>(x1);
		bignum<N>	*yp = reinterpret_cast<bignum<N> *>(y1);
		bignum<N>	*zp = reinterpret_cast<bignum<N> *>(z);
		to_montgomery(*zp, z);
		this->mont_msqr(t1, *zp);	// t1 = z^2
		to_montgomery(*xp, x1);
		to_montgomery(*yp, y1);
		this->mont_mmult(*xp, *xp, t1);	// x1 * z^2
		this->mont_mmult(t1, t1, *zp);	// t1 = z^3
		this->mont_mmult(*yp, *yp, t1);	// y1 * z^3

		// montgomery reduction
		from_montgomery(x1, *xp);
		from_montgomery(y1, *yp);
	}
	/* add-1998-cmo-2 algorithm
	 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html
	 */
	void point_add_jacobian(u64 *x3, u64 *y3, u64 *z3, const u64 *x1,
			const u64 *y1, const u64 *z1, const u64 *x2,
			const u64 *y2, const u64 *z2 = nullptr) const noexcept
	{
		bignum<N>	*x3p = reinterpret_cast<bignum<N> *>(x3);
		bignum<N>	*y3p = reinterpret_cast<bignum<N> *>(y3);
		bignum<N>	*z3p = reinterpret_cast<bignum<N> *>(z3);
		if (vli_is_zero<N>(z1)) {
			vli_set<N>(x3, x2);
			vli_set<N>(y3, y2);
			vli_set<N>(z3, z2);
			return;
		} else if (z2 != nullptr && vli_is_zero<N>(z2)) {
			vli_set<N>(x3, x1);
			vli_set<N>(y3, y1);
			vli_set<N>(z3, z1);
			return;
		}
		// add-2007-bl
		/* add-2007-bl algorithm
		 * http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html
		 */
		// U1 ... l1
		// U2 ... l2
		// S1 ... l4
		// S2 ... l5
		bool	z1_is_one = vli_is_one<N>(z1);
		bool	z2_is_one = (z2 == nullptr || vli_is_one<N>(z2));
#ifdef	WITH_ADD_2007bl
		bignum<N>	u1, u2, s1, s2, h, i, j, r, v;
		bignum<N>	t1;
#else
		bignum<N>	u1, u2, s1, s2, h, hh, hhh, r, v;
#define	t1		z1z1
#endif
		bignum<N>	z1z1, z1p;
		bignum<N>	z2z2, z2p;

		/* u1 = x1 z2^2  */
		/* u2 = x2 z1^2  */
		if (z2_is_one) {
			to_montgomery(u1, x1);
		} else {
			to_montgomery(z2p, z2);
			// z2z2 = z2^2
			this->mont_msqr(z2z2, z2p);
			to_montgomery(u1, x1);
			// u1 = x1 z2^2
			this->mont_mmult(u1, u1, z2z2);
		}
		if (z1_is_one) {
			to_montgomery(u2, x2);
		} else {
			to_montgomery(z1p, z1);
			// z1z1 = z1^2
			this->mont_msqr(z1z1, z1p);
			to_montgomery(u2, x2);
			// u2 = x2 z1^2
			this->mont_mmult(u2, u2, z1z1);
		}

		/* h = u2 - u1 */
		mod_sub(h, u2, u1);
		/* s1 = y1 z2^3  */
		// s1 = y1
		to_montgomery(s1, y1);
		if ( ! z2_is_one ) {
			// s1 = y1 z2^3
			this->mont_mmult(s1, s1, z2z2);
			this->mont_mmult(s1, s1, z2p);
		}

		/* s2 = y2 z1^3  */
		// s2 = y2
		to_montgomery(s2, y2);
		if ( !z1_is_one ) {
			// s2 = y2 z1^3
			this->mont_mmult(s2, s2, z1z1);
			this->mont_mmult(s2, s2, z1p);
		}
		/* r = s2 - s1  */
		mod_sub(r, s2, s1);

		if (h.is_zero()) {
			if (r.is_zero()) {
				/* P1 and P2 are the same - use duplicate function. */
				point_double_jacobian(x3, y3, z3, x1, y1, z1);
				return;
			}
			/* P1 is the inverse of P2.  */
			vli_clear<N>(x3);
			vli_clear<N>(y3);
			vli_clear<N>(z3);
			x3[0] = 1;
			y3[0] = 1;
			return;
		}
#ifdef	WITH_ADD_2007bl
		// r = 2 * (s2 -s1)
		mont_mult2(r, r);
		// i = (2*h)^2
		mont_mult2(i, h);
		this->mont_msqr(i, i);

		// j = h * i
		this->mont_mmult(j, h, i);
		// v = u1 * i
		this->mont_mmult(v, u1, i);

		// x3 = r^2 - j - 2*v
		this->mont_msqr(*x3p, r);
		mod_sub_from(*x3p, j);
		mod_sub_from(*x3p, v);
		mod_sub_from(*x3p, v);

		// y3 = v - x3
		mod_sub(*y3p, v, *x3p);
		// y3 = r * (v - x3)
		this->mont_mmult(*y3p, r, *y3p);
		// t1 = 2 * s1 * j
		this->mont_mmult(t1, s1, j);
		mont_mult2(t1, t1);
		// y3 = r * (v - x3) - 2 * s1 *j
		mod_sub_from(*y3p, t1);

		if (z2_is_one) {
			// z3 = 2 z1 h
			if (z1_is_one) {
				mont_mult2(*z3p, h);
			} else {
				mont_mult2(t1, z1p);
				this->mont_mmult(*z3p, t1, h);
			}
		} else {
			// z3 = ((z1 + z2)^2 -z1z1 - z2z2) * h
			if (z1_is_one) {
				// t1 = 2 * z2
				// z3 = t1 * h = 2 * z2 * h
				mont_mult2(t1, z2p);
			} else {
				// t1 = (z1 + z2)^2 -z1z1 - z2z2
				mod_add(t1, z1p, z2p);
				this->mont_msqr(t1, t1);
				mod_sub_from(t1, z1z1);
				mod_sub_from(t1, z2z2);
			}
			this->mont_mmult(*z3p, t1, h);
		}
#else
		// hh = h^2
		mont_msqr(hh, h);
		// hhh = h * hh
		mont_mmult(hhh, hh, h);
		// v = u1 * hh
		mont_mmult(v, u1, hh);

		// x3 = r^2 - hhh - 2*v
		this->mont_msqr(*x3p, r);
		mod_sub_from(*x3p, hhh);
		mod_sub_from(*x3p, v);
		mod_sub_from(*x3p, v);

		// y3 = v - x3
		mod_sub(*y3p, v, *x3p);
		// y3 = r * (v - x3)
		this->mont_mmult(*y3p, r, *y3p);
		// t1 = s1 * hhh
		this->mont_mmult(t1, s1, hhh);
		// y3 = r * (v - x3) - s1 *hhh
		mod_sub_from(*y3p, t1);

		// z3 = z1 * z2 * h
		if (z2_is_one) {
			// z3 = z1 h
			if (z1_is_one) {
				*z3p = h;
			} else {
				this->mont_mmult(*z3p, z1p, h);
			}
		} else {
			// z3 = z1 * z2 * h
			if (z1_is_one) {
				// t1 = z2
				t1 =  z2p;
			} else {
				// t1 = z1 * z2
				mont_mmult(t1, z1p, z2p);
			}
			this->mont_mmult(*z3p, t1, h);
		}
#undef	t1
#endif
		// montgomery reduction
		from_montgomery(x3, *x3p);
		from_montgomery(y3, *y3p);
		from_montgomery(z3, *z3p);
	}
private:
#if	__cplusplus >= 201703L
	mont1_func	_mont_reduction = [this](bignum<N> &res, const bignum<N> &y) {
		res.mont_reduction(y, p, k0_p);
	};
	mont1_func	_mont_sqr = [this](bignum<N> &res, const bignum<N> &y) {
		res.mont_sqr(y, p, k0_p);
	};
	mont2_func	_mont_mult = [this](bignum<N> &res, const bignum<N>& x,
					const bignum<N>& y) {
		res.mont_mult(x, y, p, k0_p);
	};
#endif
	const std::string name;
	const bignum<N> gx;
	const bignum<N> gy;
	const bignum<N> p;
	const bignum<N> n;
	const bignum<N> a;
	const bignum<N> b;
	const bignum<N> rr_p = {};
	const bignum<N> rr_n = {};
#ifdef	WITH_BARRETT
	bignum<N+1> mu_p;
	bignum<N+1> mu_n;
#endif
	bignum<N> _mont_one;
	bignum<N> _mont_a;
	const u64	k0_p = 0;
	const u64	k0_n = 0;
	const uint _ndigits = N;
	const bool _a_is_neg3 = false;
	bool _a_is_zero = false;
	const bool use_barrett = false;
	bool _inited = false;
};


template<const uint N>
struct point_t {
	u64		x[N];
	u64		y[N];
	u64		z[N];
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
	template <const uint N>slice_t(u64 vd[N]) noexcept : data(vd),
			len(N), cap(N)
	{
		vli_clear<N>(data);
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

#endif	//	__CURVE_IMPL_HPP__
