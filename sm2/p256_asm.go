// Copyright 2015 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// This file contains the Go wrapper for the constant-time, 64-bit assembly
// implementation of P256. The optimizations performed here are described in
// detail in:
// S.Gueron and V.Krasnov, "Fast prime field elliptic-curve cryptography with
//                          256-bit primes"
// https://link.springer.com/article/10.1007%2Fs13389-014-0090-x
// https://eprint.iacr.org/2013/816.pdf

// +build amd64 arm64

package sm2

import (
	"math/big"
	"sync"
)

type (
	p256Curve struct {
		*CurveParams
	}

	p256Point struct {
		xyz [12]uint64
	}
)

var (
	pSM2            p256Curve
	p256Precomputed *[43][32 * 8]uint64
	precomputeOnce  sync.Once
	p256MontOne     = []uint64{1, 0xffffffff, 0, 0x100000000}
	nRR             = []uint64{0x901192af7c114f20, 0x3464504ade6fa2fa, 0x620fc84c3affe0d4, 0x1eb5e412a22b3d3b}
)

// p256Mul operates in a Montgomery domain with R = 2^256 mod p, where p is the
// underlying field of the curve. (See initP256 for the value.) Thus rr here is
// R×R mod p. See comment in Inverse about how this is used.
var rr = []uint64{0x0000000200000003, 0x00000002ffffffff, 0x0000000100000001, 0x0000000400000002}

func initSM2() {
	// See FIPS 186-3, section D.2.3
	pSM2.CurveParams = sm2Params
}

func (curve p256Curve) Params() *CurveParams {
	return curve.CurveParams
}

// Functions implemented in p256_asm_*64.s
// Montgomery multiplication modulo P256
//go:noescape
func p256Mul(res, in1, in2 []uint64)

// Montgomery square modulo P256, repeated n times (n >= 1)
//go:noescape
func p256Sqr(res, in []uint64, n int)

// Montgomery multiplication by 1
//go:noescape
func p256FromMont(res, in []uint64)

// iff cond == 1  val <- -val
//go:noescape
func p256NegCond(val []uint64, cond int)

// if cond == 0 res <- b; else res <- a
//go:noescape
func p256MovCond(res, a, b []uint64, cond int)

// Endianness swap
//go:noescape
func p256BigToLittle(res []uint64, in []byte)

//go:noescape
func p256LittleToBig(res []byte, in []uint64)

// Constant time table access
//go:noescape
func p256Select(point, table []uint64, idx int)

//go:noescape
func p256SelectBase(point, table []uint64, idx int)

// Montgomery multiplication modulo Ord(G)
//go:noescape
func p256OrdMul(res, in1, in2 []uint64)

// Montgomery square modulo Ord(G), repeated n times
//go:noescape
func p256OrdSqr(res, in []uint64, n int)

// Point add with in2 being affine point
// If sign == 1 -> in2 = -in2
// If sel == 0 -> res = in1
// if zero == 0 -> res = in2
//go:noescape
func p256PointAddAffineAsm(res, in1, in2 []uint64, sign, sel, zero int)

// Point add. Returns one if the two input points were equal and zero
// otherwise. (Note that, due to the way that the equations work out, some
// representations of ∞ are considered equal to everything by this function.)
//go:noescape
func p256PointAddAsm(res, in1, in2 []uint64) int

// Point double
//go:noescape
func p256PointDoubleAsm(res, in []uint64)

func (curve p256Curve) Inverse(k *big.Int) *big.Int {
	if k.Sign() < 0 {
		// This should never happen.
		k = new(big.Int).Neg(k)
	}

	if k.Cmp(pSM2.N) >= 0 {
		// This should never happen.
		k = new(big.Int).Mod(k, pSM2.N)
	}
	return new(big.Int).ModInverse(k, curve.CurveParams.N)
	/*
		// table will store precomputed powers of x.
		var table [4 * 9]uint64
		var (
			_1      = table[4*0 : 4*1]
			_11     = table[4*1 : 4*2]
			_101    = table[4*2 : 4*3]
			_111    = table[4*3 : 4*4]
			_1111   = table[4*4 : 4*5]
			_10101  = table[4*5 : 4*6]
			_101111 = table[4*6 : 4*7]
			x       = table[4*7 : 4*8]
			t       = table[4*8 : 4*9]
		)

		fromBig(x[:], k)
		// This code operates in the Montgomery domain where R = 2^256 mod n
		// and n is the order of the scalar field. (See initP256 for the
		// value.) Elements in the Montgomery domain take the form a×R and
		// multiplication of x and y in the calculates (x × y × R^-1) mod n. RR
		// is R×R mod n thus the Montgomery multiplication x and RR gives x×R,
		// i.e. converts x into the Montgomery domain.
		// Window values borrowed from https://briansmith.org/ecc-inversion-addition-chains-01#p256_scalar_inversion
		p256OrdMul(_1, x, RR)      // _1
		p256OrdSqr(x, _1, 1)       // _10
		p256OrdMul(_11, x, _1)     // _11
		p256OrdMul(_101, x, _11)   // _101
		p256OrdMul(_111, x, _101)  // _111
		p256OrdSqr(x, _101, 1)     // _1010
		p256OrdMul(_1111, _101, x) // _1111

		p256OrdSqr(t, x, 1)          // _10100
		p256OrdMul(_10101, t, _1)    // _10101
		p256OrdSqr(x, _10101, 1)     // _101010
		p256OrdMul(_101111, _101, x) // _101111
		p256OrdMul(x, _10101, x)     // _111111 = x6
		p256OrdSqr(t, x, 2)          // _11111100
		p256OrdMul(t, t, _11)        // _11111111 = x8
		p256OrdSqr(x, t, 8)          // _ff00
		p256OrdMul(x, x, t)          // _ffff = x16
		p256OrdSqr(t, x, 16)         // _ffff0000
		p256OrdMul(t, t, x)          // _ffffffff = x32

		p256OrdSqr(x, t, 64)
		p256OrdMul(x, x, t)
		p256OrdSqr(x, x, 32)
		p256OrdMul(x, x, t)

		sqrs := []uint8{
			6, 5, 4, 5, 5,
			4, 3, 3, 5, 9,
			6, 2, 5, 6, 5,
			4, 5, 5, 3, 10,
			2, 5, 5, 3, 7, 6}
		muls := [][]uint64{
			_101111, _111, _11, _1111, _10101,
			_101, _101, _101, _111, _101111,
			_1111, _1, _1, _1111, _111,
			_111, _111, _101, _11, _101111,
			_11, _11, _11, _1, _10101, _1111}

		for i, s := range sqrs {
			p256OrdSqr(x, x, int(s))
			p256OrdMul(x, x, muls[i])
		}

		// Multiplying by one in the Montgomery domain converts a Montgomery
		// value out of the domain.
		one := []uint64{1, 0, 0, 0}
		p256OrdMul(x, x, one)

		xOut := make([]byte, 32)
		p256LittleToBig(xOut, x)
		return new(big.Int).SetBytes(xOut)
	*/
}

// fromBig converts a *big.Int into a format used by this code.
func fromBig(out []uint64, big *big.Int) {
	for i := range out {
		out[i] = 0
	}

	for i, v := range big.Bits() {
		out[i] = uint64(v)
	}
}

// toBig convert u64 slice format to *big.Int
func toBig(x []uint64) *big.Int {
	var res [4]big.Word
	res[0] = big.Word(x[0])
	res[1] = big.Word(x[1])
	res[2] = big.Word(x[2])
	res[3] = big.Word(x[3])
	return new(big.Int).SetBits(res[:])
}

// p256GetScalar endian-swaps the big-endian scalar value from in and writes it
// to out. If the scalar is equal or greater than the order of the group, it's
// reduced modulo that order.
func p256GetScalar(out []uint64, in []byte) {
	n := new(big.Int).SetBytes(in)

	if n.Cmp(pSM2.N) >= 0 {
		n.Mod(n, pSM2.N)
	}
	fromBig(out, n)
}

func maybeReduceModP(in *big.Int) *big.Int {
	if in.Cmp(pSM2.P) < 0 {
		return in
	}
	return new(big.Int).Mod(in, pSM2.P)
}

func (curve p256Curve) CombinedMult(bigX, bigY *big.Int, baseScalar, scalar []byte) (x, y *big.Int) {
	scalarReversed := make([]uint64, 4)
	var r1, r2 p256Point
	p256GetScalar(scalarReversed, baseScalar)
	r1IsInfinity := scalarIsZero(scalarReversed)
	r1.p256BaseMult(scalarReversed)

	p256GetScalar(scalarReversed, scalar)
	r2IsInfinity := scalarIsZero(scalarReversed)
	fromBig(r2.xyz[0:4], maybeReduceModP(bigX))
	fromBig(r2.xyz[4:8], maybeReduceModP(bigY))
	p256Mul(r2.xyz[0:4], r2.xyz[0:4], rr[:])
	p256Mul(r2.xyz[4:8], r2.xyz[4:8], rr[:])

	// This sets r2's Z value to 1, in the Montgomery domain.
	copy(r2.xyz[8:], p256MontOne)
	//r2.xyz[8] = 0x0000000000000001
	//r2.xyz[9] = 0xffffffff
	//r2.xyz[10] = 0
	//r2.xyz[11] = 0x100000000

	r2.p256ScalarMult(scalarReversed)

	var sum, double p256Point
	pointsEqual := p256PointAddAsm(sum.xyz[:], r1.xyz[:], r2.xyz[:])
	p256PointDoubleAsm(double.xyz[:], r1.xyz[:])
	sum.CopyConditional(&double, pointsEqual)
	sum.CopyConditional(&r1, r2IsInfinity)
	sum.CopyConditional(&r2, r1IsInfinity)

	return sum.p256PointToAffine()
}

func (curve p256Curve) ScalarBaseMult(scalar []byte) (x, y *big.Int) {
	scalarReversed := make([]uint64, 4)
	p256GetScalar(scalarReversed, scalar)

	var r p256Point
	r.p256BaseMult(scalarReversed)
	return r.p256PointToAffine()
}

func (curve p256Curve) ScalarMult(bigX, bigY *big.Int, scalar []byte) (x, y *big.Int) {
	scalarReversed := make([]uint64, 4)
	p256GetScalar(scalarReversed, scalar)

	var r p256Point
	fromBig(r.xyz[0:4], maybeReduceModP(bigX))
	fromBig(r.xyz[4:8], maybeReduceModP(bigY))
	p256Mul(r.xyz[0:4], r.xyz[0:4], rr[:])
	p256Mul(r.xyz[4:8], r.xyz[4:8], rr[:])
	// This sets r2's Z value to 1, in the Montgomery domain.
	copy(r.xyz[8:], p256MontOne)
	//r.xyz[8] = 0x0000000000000001
	//r.xyz[9] = 0xffffffff
	//r.xyz[10] = 0
	//r.xyz[11] = 0x100000000

	r.p256ScalarMult(scalarReversed)
	return r.p256PointToAffine()
}

// uint64IsZero returns 1 if x is zero and zero otherwise.
func uint64IsZero(x uint64) int {
	x = ^x
	x &= x >> 32
	x &= x >> 16
	x &= x >> 8
	x &= x >> 4
	x &= x >> 2
	x &= x >> 1
	return int(x & 1)
}

// scalarIsZero returns 1 if scalar represents the zero value, and zero
// otherwise.
func scalarIsZero(scalar []uint64) int {
	return uint64IsZero(scalar[0] | scalar[1] | scalar[2] | scalar[3])
}

func (p *p256Point) p256PointToAffine() (x, y *big.Int) {
	zInv := make([]uint64, 4)
	zInvSq := make([]uint64, 4)
	p256Inverse(zInv, p.xyz[8:12])
	p256Sqr(zInvSq, zInv, 1)
	p256Mul(zInv, zInv, zInvSq)

	p256Mul(zInvSq, p.xyz[0:4], zInvSq)
	p256Mul(zInv, p.xyz[4:8], zInv)

	p256FromMont(zInvSq, zInvSq)
	p256FromMont(zInv, zInv)

	xOut := make([]byte, 32)
	yOut := make([]byte, 32)
	p256LittleToBig(xOut, zInvSq)
	p256LittleToBig(yOut, zInv)

	return new(big.Int).SetBytes(xOut), new(big.Int).SetBytes(yOut)
}

// CopyConditional copies overwrites p with src if v == 1, and leaves p
// unchanged if v == 0.
func (p *p256Point) CopyConditional(src *p256Point, v int) {
	pMask := uint64(v) - 1
	srcMask := ^pMask

	for i, n := range p.xyz {
		p.xyz[i] = (n & pMask) | (src.xyz[i] & srcMask)
	}
}

// p256Inverse sets out to in^-1 mod p.
func p256Inverse(out, in []uint64) {
	var stack [6 * 4]uint64
	p2 := stack[4*0 : 4*0+4]
	p4 := stack[4*1 : 4*1+4]
	p8 := stack[4*2 : 4*2+4]
	p16 := stack[4*3 : 4*3+4]
	p32 := stack[4*4 : 4*4+4]
	p3 := stack[4*5 : 4*5+4]

	p256Sqr(out, in, 1)
	p256Mul(p2, out, in) // 3*p

	p256Sqr(out, p2, 1)
	p256Mul(p3, out, in) // 7*p

	p256Sqr(out, p2, 2)
	p256Mul(p4, out, p2) // f*p

	p256Sqr(out, p4, 4)
	p256Mul(p8, out, p4) // ff*p

	p256Sqr(out, p8, 8)
	p256Mul(p16, out, p8) // ffff*p

	p256Sqr(out, p16, 16)
	p256Mul(p32, out, p16) // ffffffff*p

	p256Sqr(out, p16, 8)
	p256Mul(out, out, p8) // ffffff*p
	p256Sqr(out, out, 4)
	p256Mul(out, out, p4) // fffffff*p
	p256Sqr(out, out, 3)
	p256Mul(out, out, p3)
	p256Sqr(out, out, 1) // fffffffe*p

	p256Sqr(out, out, 32)
	p256Mul(out, out, p32) // fffffffeffffffff

	p256Sqr(out, out, 32)
	p256Mul(out, out, p32)
	p256Sqr(out, out, 32)
	p256Mul(out, out, p32) // ... ffffffffffffffff
	p256Sqr(out, out, 32)
	p256Mul(out, out, p32)
	p256Sqr(out, out, 32) // ... ffffffff00000000

	p256Sqr(out, out, 32)
	p256Mul(out, out, p32)

	p256Sqr(out, out, 16)
	p256Mul(out, out, p16)

	p256Sqr(out, out, 8)
	p256Mul(out, out, p8)

	p256Sqr(out, out, 4)
	p256Mul(out, out, p4)

	p256Sqr(out, out, 2)
	p256Mul(out, out, p2)

	p256Sqr(out, out, 2)
	p256Mul(out, out, in)
}

func (p *p256Point) p256StorePoint(r *[16 * 4 * 3]uint64, index int) {
	copy(r[index*12:], p.xyz[:])
}

func boothW5(in uint) (int, int) {
	var s uint = ^((in >> 5) - 1)
	var d uint = (1 << 6) - in - 1
	d = (d & s) | (in & (^s))
	d = (d >> 1) + (d & 1)
	return int(d), int(s & 1)
}

func boothW6(in uint) (int, int) {
	var s uint = ^((in >> 6) - 1)
	var d uint = (1 << 7) - in - 1
	d = (d & s) | (in & (^s))
	d = (d >> 1) + (d & 1)
	return int(d), int(s & 1)
}

func initTable() {
	p256Precomputed = new([43][32 * 8]uint64)

	// basePoint in montgomery form
	basePoint := []uint64{
		0x61328990f418029e, 0x3e7981eddca6c050, 0xd6a1ed99ac24c3c3, 0x91167a5ee1c13b05,
		0xc1354e593c2d0ddd, 0xc1f5e5788d3295fa, 0x8d4cfb066e2a48f8, 0x63cd65d481d735bd,
		0x0000000000000001, 0xffffffff, 0, 0x100000000,
	}
	// convert Gx, Gy to montgomery form
	fromBig(basePoint[:4], sm2Params.Gx)
	fromBig(basePoint[4:8], sm2Params.Gy)
	p256Mul(basePoint[:4], basePoint[:4], rr)
	p256Mul(basePoint[4:8], basePoint[4:8], rr)
	copy(basePoint[8:], p256MontOne)
	t1 := make([]uint64, 12)
	t2 := make([]uint64, 12)
	copy(t2, basePoint)

	zInv := make([]uint64, 4)
	zInvSq := make([]uint64, 4)
	for j := 0; j < 32; j++ {
		copy(t1, t2)
		for i := 0; i < 43; i++ {
			// The window size is 6 so we need to double 6 times.
			if i != 0 {
				for k := 0; k < 6; k++ {
					p256PointDoubleAsm(t1, t1)
				}
			}
			// Convert the point to affine form. (Its values are
			// still in Montgomery form however.)
			p256Inverse(zInv, t1[8:12])
			p256Sqr(zInvSq, zInv, 1)
			p256Mul(zInv, zInv, zInvSq)

			p256Mul(t1[:4], t1[:4], zInvSq)
			p256Mul(t1[4:8], t1[4:8], zInv)

			copy(t1[8:12], basePoint[8:12])
			// Update the table entry
			copy(p256Precomputed[i][j*8:], t1[:8])
		}
		if j == 0 {
			p256PointDoubleAsm(t2, basePoint)
		} else {
			p256PointAddAsm(t2, t2, basePoint)
		}
	}
}

func (p *p256Point) p256BaseMult(scalar []uint64) {
	precomputeOnce.Do(initTable)

	wvalue := (scalar[0] << 1) & 0x7f
	sel, sign := boothW6(uint(wvalue))
	p256SelectBase(p.xyz[0:8], p256Precomputed[0][0:], sel)
	p256NegCond(p.xyz[4:8], sign)

	// (This is one, in the Montgomery domain.)
	copy(p.xyz[8:], p256MontOne)
	//p.xyz[8] = 0x0000000000000001
	//p.xyz[9] = 0xffffffff
	//p.xyz[10] = 0
	//p.xyz[11] = 0x100000000

	var t0 p256Point
	// (This is one, in the Montgomery domain.)
	copy(t0.xyz[8:], p256MontOne)
	//t0.xyz[8] = 0x0000000000000001
	//t0.xyz[9] = 0xffffffff
	//t0.xyz[10] = 0
	//t0.xyz[11] = 0x100000000

	index := uint(5)
	zero := sel

	for i := 1; i < 43; i++ {
		if index < 192 {
			wvalue = ((scalar[index/64] >> (index % 64)) + (scalar[index/64+1] << (64 - (index % 64)))) & 0x7f
		} else {
			wvalue = (scalar[index/64] >> (index % 64)) & 0x7f
		}
		index += 6
		sel, sign = boothW6(uint(wvalue))
		p256SelectBase(t0.xyz[0:8], p256Precomputed[i][0:], sel)
		p256PointAddAffineAsm(p.xyz[0:12], p.xyz[0:12], t0.xyz[0:8], sign, sel, zero)
		zero |= sel
	}
}

func (p *p256Point) p256ScalarMult(scalar []uint64) {
	// precomp is a table of precomputed points that stores powers of p
	// from p^1 to p^16.
	var precomp [16 * 4 * 3]uint64
	var t0, t1, t2, t3 p256Point

	// Prepare the table
	p.p256StorePoint(&precomp, 0) // 1

	p256PointDoubleAsm(t0.xyz[:], p.xyz[:])
	p256PointDoubleAsm(t1.xyz[:], t0.xyz[:])
	p256PointDoubleAsm(t2.xyz[:], t1.xyz[:])
	p256PointDoubleAsm(t3.xyz[:], t2.xyz[:])
	t0.p256StorePoint(&precomp, 1)  // 2
	t1.p256StorePoint(&precomp, 3)  // 4
	t2.p256StorePoint(&precomp, 7)  // 8
	t3.p256StorePoint(&precomp, 15) // 16

	p256PointAddAsm(t0.xyz[:], t0.xyz[:], p.xyz[:])
	p256PointAddAsm(t1.xyz[:], t1.xyz[:], p.xyz[:])
	p256PointAddAsm(t2.xyz[:], t2.xyz[:], p.xyz[:])
	t0.p256StorePoint(&precomp, 2) // 3
	t1.p256StorePoint(&precomp, 4) // 5
	t2.p256StorePoint(&precomp, 8) // 9

	p256PointDoubleAsm(t0.xyz[:], t0.xyz[:])
	p256PointDoubleAsm(t1.xyz[:], t1.xyz[:])
	t0.p256StorePoint(&precomp, 5) // 6
	t1.p256StorePoint(&precomp, 9) // 10

	p256PointAddAsm(t2.xyz[:], t0.xyz[:], p.xyz[:])
	p256PointAddAsm(t1.xyz[:], t1.xyz[:], p.xyz[:])
	t2.p256StorePoint(&precomp, 6)  // 7
	t1.p256StorePoint(&precomp, 10) // 11

	p256PointDoubleAsm(t0.xyz[:], t0.xyz[:])
	p256PointDoubleAsm(t2.xyz[:], t2.xyz[:])
	t0.p256StorePoint(&precomp, 11) // 12
	t2.p256StorePoint(&precomp, 13) // 14

	p256PointAddAsm(t0.xyz[:], t0.xyz[:], p.xyz[:])
	p256PointAddAsm(t2.xyz[:], t2.xyz[:], p.xyz[:])
	t0.p256StorePoint(&precomp, 12) // 13
	t2.p256StorePoint(&precomp, 14) // 15

	// Start scanning the window from top bit
	index := uint(254)
	var sel, sign int

	wvalue := (scalar[index/64] >> (index % 64)) & 0x3f
	sel, _ = boothW5(uint(wvalue))

	p256Select(p.xyz[0:12], precomp[0:], sel)
	zero := sel

	for index > 4 {
		index -= 5
		p256PointDoubleAsm(p.xyz[:], p.xyz[:])
		p256PointDoubleAsm(p.xyz[:], p.xyz[:])
		p256PointDoubleAsm(p.xyz[:], p.xyz[:])
		p256PointDoubleAsm(p.xyz[:], p.xyz[:])
		p256PointDoubleAsm(p.xyz[:], p.xyz[:])

		if index < 192 {
			wvalue = ((scalar[index/64] >> (index % 64)) + (scalar[index/64+1] << (64 - (index % 64)))) & 0x3f
		} else {
			wvalue = (scalar[index/64] >> (index % 64)) & 0x3f
		}

		sel, sign = boothW5(uint(wvalue))

		p256Select(t0.xyz[0:], precomp[0:], sel)
		p256NegCond(t0.xyz[4:8], sign)
		p256PointAddAsm(t1.xyz[:], p.xyz[:], t0.xyz[:])
		p256MovCond(t1.xyz[0:12], t1.xyz[0:12], p.xyz[0:12], sel)
		p256MovCond(p.xyz[0:12], t1.xyz[0:12], t0.xyz[0:12], zero)
		zero |= sel
	}

	p256PointDoubleAsm(p.xyz[:], p.xyz[:])
	p256PointDoubleAsm(p.xyz[:], p.xyz[:])
	p256PointDoubleAsm(p.xyz[:], p.xyz[:])
	p256PointDoubleAsm(p.xyz[:], p.xyz[:])
	p256PointDoubleAsm(p.xyz[:], p.xyz[:])

	wvalue = (scalar[0] << 1) & 0x3f
	sel, sign = boothW5(uint(wvalue))

	p256Select(t0.xyz[0:], precomp[0:], sel)
	p256NegCond(t0.xyz[4:8], sign)
	p256PointAddAsm(t1.xyz[:], p.xyz[:], t0.xyz[:])
	p256MovCond(t1.xyz[0:12], t1.xyz[0:12], p.xyz[0:12], sel)
	p256MovCond(p.xyz[0:12], t1.xyz[0:12], t0.xyz[0:12], zero)
}
