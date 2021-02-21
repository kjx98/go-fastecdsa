// Copyright 2010 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package elliptic implements several standard elliptic curves over prime
// fields.
package sm2

// This package operates, internally, on Jacobian coordinates. For a given
// (x, y) position on the curve, the Jacobian coordinates are (x1, y1, z1)
// where x = x1/z1² and y = y1/z1³. The greatest speedups come when the whole
// calculation can be performed within the transform (as in ScalarMult and
// ScalarBaseMult). But even for Add and Double, it's faster to apply and
// reverse the transform than to operate in affine coordinates.

import (
	//"log"
	"math/big"
)

// sm2Curve contains the parameters of an elliptic curve and also provides
// a generic, non-constant time implementation of Curve.
type sm2Curve struct {
	*CurveParams
	mu *big.Int
	k0 uint64
	rr *big.Int
}

var sm2g sm2Curve
var sm2Base *big.Int
var montOne *big.Int

func initSM2go() {
	// Use pure Go implementation.
	sm2g.CurveParams = sm2Params
	sm2g.k0 = 1 //0x327f9e8872350975
	n512 := new(big.Int).SetUint64(1)
	sm2Base = new(big.Int).Lsh(n512, 256)
	n512.Lsh(n512, 512)
	sm2g.mu = new(big.Int).Div(n512, sm2g.P)
	rrBits := []big.Word{0x200000003, 0x2ffffffff, 0x100000001, 0x400000002}
	sm2g.rr = new(big.Int).SetBits(rrBits)
	//rrBits = []big.Word{0x00000001, 0xffffffff, 0x00000000, 0x100000000}
	//montOne = new(big.Int).SetBits(rrBits)
	montOne = new(big.Int).SetUint64(1)
	rr := new(big.Int).Mul(sm2g.mu, sm2g.P)
	if rr.Cmp(n512) >= 0 {
		panic("mu large not floor")
	}
	rem := new(big.Int).Sub(n512, rr)
	if rem.Cmp(sm2g.P) >= 0 {
		panic("mu small not floor")
	}
}

func (curve sm2Curve) GetMu() *big.Int {
	return curve.mu
}

func (curve sm2Curve) BarrettMod(prod *big.Int) *big.Int {
	qBits := prod.Bits()
	if len(qBits) < 5 {
		if prod.Cmp(curve.P) >= 0 {
			return new(big.Int).Sub(prod, curve.P)
		}
		return prod
	}
	q := BarrettDiv(prod, curve.mu)
	q3 := new(big.Int).Mul(q, curve.P)
	rr := new(big.Int).Sub(prod, q3)
	if rr.Cmp(curve.P) >= 0 {
		rr.Sub(rr, curve.P)
		//log.Printf("r : %s", r.Text(16))
	}
	return rr
}

func (curve sm2Curve) montK0() uint64 {
	t := uint64(1)
	N := uint64(curve.P.Bits()[0])
	for i := 1; i < 4; i++ {
		t = t * t * N
	}
	return -t
}

//go: noinline
func getWordAt(x *big.Int, i int) big.Word {
	xx := x.Bits()
	if len(xx) <= i || i < 0 {
		return big.Word(0)
	}
	return xx[i]
}

// use poly, shift/add/sub
func (curve sm2Curve) multP(u uint64) *big.Int {
	var uBN big.Int
	if u == 0 {
		return &uBN
	}
	var uBits [5]big.Word
	uBits[4] = big.Word(u)
	uBits[1] = big.Word(u)
	res := new(big.Int).SetBits(uBits[:])
	var uBits2 [5]big.Word
	uBits2[0] = big.Word(u)
	uBits2[1] = big.Word(u << 32)
	uBits2[2] = big.Word(u >> 32)
	uBits2[3] = big.Word(u << 32)
	uBits2[4] = big.Word(u >> 32)
	n19296 := new(big.Int).SetBits(uBits2[:])
	res.Sub(res, n19296)
	return res
}

func (curve sm2Curve) multPh(u uint64) *big.Int {
	var uBN big.Int
	if u == 0 {
		return &uBN
	}
	var uBits [4]big.Word
	uBits[3] = big.Word(u)
	uBits[0] = big.Word(u)
	res := new(big.Int).SetBits(uBits[:])
	var uBits2 [4]big.Word
	uBits2[0] = big.Word(u << 32)
	uBits2[1] = big.Word(u >> 32)
	uBits2[2] = big.Word(u << 32)
	uBits2[3] = big.Word(u >> 32)
	n19296 := new(big.Int).SetBits(uBits2[:])
	res.Sub(res, n19296)
	return res
}

func (curve sm2Curve) montRed(y *big.Int) *big.Int {
	r := new(big.Int)
	r.SetBytes(y.Bytes())
	for i := 0; i < 4; i++ {
		r0 := getWordAt(r, 0)
		u := uint64(r0) * curve.k0 // uint64 same as modula 2^64
		s := curve.multPh(u)
		r.Add(r, s)
		ww := r.Bits()
		r.SetBits(ww[1:])
	}
	if r.Cmp(curve.P) >= 0 {
		r.Sub(r, curve.P)
	}
	return r
}

func (curve sm2Curve) montMul(x, y *big.Int) *big.Int {
	r := new(big.Int)
	yy := y.Bits()
	var t big.Int
	for i := 0; i < 4; i++ {
		t.SetUint64(uint64(yy[i]))
		t.Mul(&t, x)
		r.Add(r, &t)
		r0 := getWordAt(r, 0)
		ww := r.Bits()
		r.SetBits(ww[1:])
		u := uint64(r0) * curve.k0 // uint64 same as modula 2^64
		s := curve.multPh(u)
		r.Add(r, s)
	}
	if r.Cmp(curve.P) >= 0 {
		r.Sub(r, curve.P)
	}
	return r
}

func (curve sm2Curve) montModMul(x, y *big.Int) *big.Int {
	xp := curve.montMul(x, curve.rr)
	yp := curve.montMul(y, curve.rr)
	res := curve.montMul(xp, yp)
	//return curve.montRed(res)
	return curve.montMul(montOne, res)
}

func (curve sm2Curve) IsOnCurve(x, y *big.Int) bool {
	// y² = x³ - 3x + b
	y2 := new(big.Int).Mul(y, y)
	y2.Mod(y2, curve.P)

	x3 := new(big.Int).Mul(x, x)
	x3.Mul(x3, x)

	threeX := new(big.Int).Lsh(x, 1)
	threeX.Add(threeX, x)

	x3.Sub(x3, threeX)
	x3.Add(x3, curve.B)
	x3.Mod(x3, curve.P)

	return x3.Cmp(y2) == 0
}

// affineFromJacobian reverses the Jacobian transform. See the comment at the
// top of the file. If the point is ∞ it returns 0, 0.
func (curve sm2Curve) AffineFromJacobian(x, y, z *big.Int) (xOut, yOut *big.Int) {
	if z.Sign() == 0 {
		return new(big.Int), new(big.Int)
	}

	zinv := new(big.Int).ModInverse(z, curve.P)
	zinvsq := new(big.Int).Mul(zinv, zinv)

	xOut = new(big.Int).Mul(x, zinvsq)
	xOut.Mod(xOut, curve.P)
	zinvsq.Mul(zinvsq, zinv)
	yOut = new(big.Int).Mul(y, zinvsq)
	yOut.Mod(yOut, curve.P)
	return
}

func (curve sm2Curve) Add(x1, y1, x2, y2 *big.Int) (*big.Int, *big.Int) {
	z1 := zForAffine(x1, y1)
	z2 := zForAffine(x2, y2)
	return curve.AffineFromJacobian(curve.AddJacobian(x1, y1, z1, x2, y2, z2))
}

// addJacobian takes two points in Jacobian coordinates, (x1, y1, z1) and
// (x2, y2, z2) and returns their sum, also in Jacobian form.
func (curve sm2Curve) AddJacobian(x1, y1, z1, x2, y2, z2 *big.Int) (*big.Int, *big.Int, *big.Int) {
	// See https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html#addition-add-2007-bl
	x3, y3, z3 := new(big.Int), new(big.Int), new(big.Int)
	if z1.Sign() == 0 {
		x3.Set(x2)
		y3.Set(y2)
		z3.Set(z2)
		return x3, y3, z3
	}
	if z2.Sign() == 0 {
		x3.Set(x1)
		y3.Set(y1)
		z3.Set(z1)
		return x3, y3, z3
	}

	z1z1 := new(big.Int).Mul(z1, z1)
	z1z1.Mod(z1z1, curve.P)
	z2z2 := new(big.Int).Mul(z2, z2)
	z2z2.Mod(z2z2, curve.P)

	u1 := new(big.Int).Mul(x1, z2z2)
	u1.Mod(u1, curve.P)
	u2 := new(big.Int).Mul(x2, z1z1)
	u2.Mod(u2, curve.P)
	h := new(big.Int).Sub(u2, u1)
	xEqual := h.Sign() == 0
	if h.Sign() == -1 {
		h.Add(h, curve.P)
	}
	i := new(big.Int).Lsh(h, 1)
	i.Mul(i, i)
	j := new(big.Int).Mul(h, i)

	s1 := new(big.Int).Mul(y1, z2)
	s1.Mul(s1, z2z2)
	s1.Mod(s1, curve.P)
	s2 := new(big.Int).Mul(y2, z1)
	s2.Mul(s2, z1z1)
	s2.Mod(s2, curve.P)
	r := new(big.Int).Sub(s2, s1)
	if r.Sign() == -1 {
		r.Add(r, curve.P)
	}
	yEqual := r.Sign() == 0
	if xEqual && yEqual {
		return curve.DoubleJacobian(x1, y1, z1)
	}
	r.Lsh(r, 1)
	v := new(big.Int).Mul(u1, i)

	x3.Set(r)
	x3.Mul(x3, x3)
	x3.Sub(x3, j)
	x3.Sub(x3, v)
	x3.Sub(x3, v)
	x3.Mod(x3, curve.P)

	y3.Set(r)
	v.Sub(v, x3)
	y3.Mul(y3, v)
	s1.Mul(s1, j)
	s1.Lsh(s1, 1)
	y3.Sub(y3, s1)
	y3.Mod(y3, curve.P)

	z3.Add(z1, z2)
	z3.Mul(z3, z3)
	z3.Sub(z3, z1z1)
	z3.Sub(z3, z2z2)
	z3.Mul(z3, h)
	z3.Mod(z3, curve.P)

	return x3, y3, z3
}

func (curve sm2Curve) Double(x1, y1 *big.Int) (*big.Int, *big.Int) {
	z1 := zForAffine(x1, y1)
	return curve.AffineFromJacobian(curve.DoubleJacobian(x1, y1, z1))
}

// doubleJacobian takes a point in Jacobian coordinates, (x, y, z), and
// returns its double, also in Jacobian form.
func (curve sm2Curve) DoubleJacobian(x, y, z *big.Int) (*big.Int, *big.Int, *big.Int) {
	// See https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html#doubling-dbl-2001-b
	if z.Sign() == 0 {
		return x, y, z
	}
	delta := new(big.Int).Mul(z, z)
	delta.Mod(delta, curve.P)
	gamma := new(big.Int).Mul(y, y)
	gamma.Mod(gamma, curve.P)
	alpha := new(big.Int).Sub(x, delta)
	if alpha.Sign() == -1 {
		alpha.Add(alpha, curve.P)
	}
	alpha2 := new(big.Int).Add(x, delta)
	alpha.Mul(alpha, alpha2)
	alpha2.Set(alpha)
	alpha.Lsh(alpha, 1)
	alpha.Add(alpha, alpha2)

	beta := alpha2.Mul(x, gamma)

	x3 := new(big.Int).Mul(alpha, alpha)
	beta8 := new(big.Int).Lsh(beta, 3)
	beta8.Mod(beta8, curve.P)
	x3.Sub(x3, beta8)
	if x3.Sign() == -1 {
		x3.Add(x3, curve.P)
	}
	x3.Mod(x3, curve.P)

	z3 := new(big.Int).Add(y, z)
	z3.Mul(z3, z3)
	z3.Sub(z3, gamma)
	if z3.Sign() == -1 {
		z3.Add(z3, curve.P)
	}
	z3.Sub(z3, delta)
	if z3.Sign() == -1 {
		z3.Add(z3, curve.P)
	}
	z3.Mod(z3, curve.P)

	beta.Lsh(beta, 2)
	beta.Sub(beta, x3)
	if beta.Sign() == -1 {
		beta.Add(beta, curve.P)
	}
	y3 := alpha.Mul(alpha, beta)

	gamma.Mul(gamma, gamma)
	gamma.Lsh(gamma, 3)
	gamma.Mod(gamma, curve.P)

	y3.Sub(y3, gamma)
	if y3.Sign() == -1 {
		y3.Add(y3, curve.P)
	}
	y3.Mod(y3, curve.P)

	return x3, y3, z3
}

func (curve sm2Curve) ScalarMult(Bx, By *big.Int, k []byte) (*big.Int, *big.Int) {
	Bz := new(big.Int).SetInt64(1)
	x, y, z := new(big.Int), new(big.Int), new(big.Int)

	for _, byte := range k {
		for bitNum := 0; bitNum < 8; bitNum++ {
			x, y, z = curve.DoubleJacobian(x, y, z)
			if byte&0x80 == 0x80 {
				x, y, z = curve.AddJacobian(Bx, By, Bz, x, y, z)
			}
			byte <<= 1
		}
	}

	return curve.AffineFromJacobian(x, y, z)
}

func (curve sm2Curve) ScalarBaseMult(k []byte) (*big.Int, *big.Int) {
	return curve.ScalarMult(curve.Gx, curve.Gy, k)
}
