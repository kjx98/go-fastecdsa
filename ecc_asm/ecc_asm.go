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

package ecc

import (
	"math/big"
	"sync"
)

// Functions implemented in ecc_asm_*64.s
// multiplication modulo p
//go:noescape
func vli_mod_mult_fast(res, in1, in2, p []uint64, ndigits uint)

// Functions implemented in ecc_asm_*64.s
// multiplication
//go:noescape
func vli_mult(res, in1, in2 []uint64)

// Functions implemented in ecc_asm_*64.s
// in inverse mod prime p
//go:noescape
func vli_mod_inv(res, in, p []uint64)

// Function mod prime with barrett reduction
// mod MUST be prime following with mu
//go:noescape
func vli_mmod_barrett(res, prod, mod []uint64)

// Function calc div quo with barrett reduction
//go:noescape
func vli_div_barrett(res, prod, mu []uint64)

// Functions implemented in ecc_asm_*64.s
// multiplication modulo p
//go:noescape
func mont_MulMod(res, in1, in2, p, rr []uint64, k0 uint64)

// Function implemented in ecc_asm_*86.s
// Exp modulo prime p
//go:noescape
func mont_ExpMod(res, in1, in2, p, rr []uint64, k0 uint64)

/*
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
*/
