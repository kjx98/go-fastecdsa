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
	"unsafe"
)

// Function check CPU support FMA
func _vli_asm_acc() (res uint)

// Function get ecc_curve point
func _get_curve(curve_id uint) (curveH unsafe.Pointer)

func _get_curve_params(p, n, b, gx, gy unsafe.Pointer, curveH unsafe.Pointer)

func _point_double_jacobian(ptR, p unsafe.Pointer, curveH unsafe.Pointer)
func _point_double(ptR, p unsafe.Pointer, curveH unsafe.Pointer)
func _point_add_jacobian(ptR, p, q unsafe.Pointer, curveH unsafe.Pointer)
func _point_add(ptR, p, q unsafe.Pointer, curveH unsafe.Pointer)
func _affine_from_jacobian(x, y, pt unsafe.Pointer, curveH unsafe.Pointer)

func _point_mult(res, p, scalar, curveH unsafe.Pointer)

// Functions implemented in ecc_asm_*64.s
// multiplication modulo p
//go:noescape
func _vli_mod_mult_fast(res, in1, in2, p unsafe.Pointer, ndigits uint)

// Functions implemented in ecc_asm_*64.s
// multiplication
//go:noescape
func _vli_mult(res, in1, in2 unsafe.Pointer)

// Functions implemented in ecc_asm_*64.s
// in inverse mod prime p
//go:noescape
func _vli_mod_inv(res, in, p, buff unsafe.Pointer)

// Function mod prime with barrett reduction
// mod MUST be prime following with mu
//go:noescape
func _vli_mmod_barrett(res, prod, mod, buff unsafe.Pointer)

// Function calc div quo with barrett reduction
//go:noescape
func _vli_div_barrett(res, prod, mu unsafe.Pointer)

// Functions implemented in ecc_asm_*64.s
// multiplication modulo p
//go:noescape
func _mont_mod_mult(res, in1, in2, p, rr unsafe.Pointer, k0 uint64)

// Functions implemented in ecc_asm_*64.s
// sqr modulo p, dup n times
//go:noescape
func _mont_mod_sqr(res, in1, p, rr unsafe.Pointer, k0 uint64, n uint)

// Function implemented in ecc_asm_*86.s
// Exp modulo prime p
//go:noescape
func _mont_mod_exp(res, in1, in2, p, rr unsafe.Pointer, k0 uint64)

// Function implemented in ecc_asm_*86.s
// u64 multiply prime P of SM2
//go:noescape
func _vli_sm2_mult_p(res *[]big.Word, u uint64)

// Function implemented in ecc_asm_*86.s
// u64 multiply prime N of SM2
//go:noescape
func _vli_sm2_mult_n(res unsafe.Pointer, u uint64)

// Functions implemented in ecc_asm_*64.s
// multiplication modulo p for SM2
//go:noescape
func _mont_sm2_mod_mult_p(res, in1, in2, p, rr unsafe.Pointer)

// Functions implemented in ecc_asm_*64.s
// multiplication modulo p for SM2
//go:noescape
func _mont_sm2_mod_mult_n(res, in1, in2, p, rr unsafe.Pointer)

// Functions implemented in ecc_asm_*64.s
// Montgomery inverse modulo prime mod
func vliModInv(in, mod []big.Word) (result []big.Word) {
	var res [32]big.Word
	var x [4]big.Word
	var buff [32]big.Word
	copy(x[:], in)
	_vli_mod_inv(unsafe.Pointer(&res[0]), unsafe.Pointer(&x[0]),
		unsafe.Pointer(&mod[0]), unsafe.Pointer(&buff[0]))
	result = res[:4]
	return
}

// Montgomery multiplication modulo sm2
func vliModMult(left, right, mdU []big.Word) (result *big.Int) {
	var res [32]big.Word
	var prod [8]big.Word
	var lf, rt [4]big.Word
	var mod [9]big.Word // should be 9, 4 word for mod, 5 word for mu
	var buff [20]big.Word
	copy(lf[:], left)
	copy(rt[:], right)
	copy(mod[:], mdU)
	_vli_mult(unsafe.Pointer(&prod[0]), unsafe.Pointer(&lf[0]),
		unsafe.Pointer(&rt[0]))
	_vli_mmod_barrett(unsafe.Pointer(&res[0]), unsafe.Pointer(&prod[0]),
		unsafe.Pointer(&mod[0]), unsafe.Pointer(&buff[0]))
	result = new(big.Int).SetBits(res[:4])
	return
}

func vliModMultBarrett(left, right *big.Int, mdU []big.Word) *big.Int {
	var res [4]big.Word
	prod := new(big.Int).Mul(left, right)
	var prd [8]big.Word
	var mod [9]big.Word // should be 9, 4 word for mod, 5 word for mu
	var buff [20]big.Word
	copy(prd[:], prod.Bits())
	copy(mod[:], mdU)
	_vli_mmod_barrett(unsafe.Pointer(&res[0]), unsafe.Pointer(&prd[0]),
		unsafe.Pointer(&mod[0]), unsafe.Pointer(&buff[0]))
	return new(big.Int).SetBits(res[:4])
}

/*
func vliBarrettDiv(prod *big.Int, muB []big.Word) (result *big.Int) {
	var res [8]big.Word
	var prd [8]big.Word
	var mu [5]big.Word // should be 9, 4 word for mod, 5 word for mu
	copy(prd[:], prod.Bits())
	copy(mu[:], muB)
	_vli_div_barrett(unsafe.Pointer(&res[0]), unsafe.Pointer(&prd[0]),
		unsafe.Pointer(&mu[0]))
	result = new(big.Int).SetBits(res[:4])
	return
}
*/

func vliModMultMont(x, y, mod []big.Word, rr []uint64, k0 uint64) *big.Int {
	var r [32]big.Word
	_mont_mod_mult(unsafe.Pointer(&r[0]), unsafe.Pointer(&x[0]),
		unsafe.Pointer(&y[0]), unsafe.Pointer(&mod[0]), unsafe.Pointer(&rr[0]),
		k0)
	return new(big.Int).SetBits(r[:4])
}

func vliExpModMont(x, y, mod []big.Word, rr []uint64, k0 uint64) (res *big.Int) {
	var r [32]big.Word
	_mont_mod_exp(unsafe.Pointer(&r[0]), unsafe.Pointer(&x[0]),
		unsafe.Pointer(&y[0]), unsafe.Pointer(&mod[0]),
		unsafe.Pointer(&rr[0]), k0)
	res = new(big.Int).SetBits(r[:4])
	return
}

func vliModMultMontP(x, y, mod []big.Word, rr []uint64) *big.Int {
	var r [32]big.Word
	_mont_sm2_mod_mult_p(unsafe.Pointer(&r[0]), unsafe.Pointer(&x[0]),
		unsafe.Pointer(&y[0]), unsafe.Pointer(&mod[0]), unsafe.Pointer(&rr[0]))
	return new(big.Int).SetBits(r[:4])
}

func vliMultP(u uint64) *big.Int {
	res := make([]big.Word, 8)
	_vli_sm2_mult_p(&res, u)
	return new(big.Int).SetBits(res[:5])
}
