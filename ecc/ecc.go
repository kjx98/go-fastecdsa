// +build amd64 arm64

package ecc

// #cgo CXXFLAGS: -O2 -Wpedantic -Wall -std=gnu++11
// #include "ecc.h"
import "C"

import (
	"math/big"
	"unsafe"
)

// Functions implemented in ecc_asm_*64.s
// Montgomery inverse modulo P256
func vliModInv(input, modB []byte) (result []big.Word) {
	var res [4]big.Word
	in := vliFromBE64(input)
	mod := vliFromBE64(modB)
	C.vli_mod_inv((*C.u64)(unsafe.Pointer(&res[0])), &in[0], &mod[0], 4)
	result = res[:]
	return
}

// Montgomery multiplication modulo P256
func vliModMult(left, right, mdU []big.Word) (result *big.Int) {
	var res [4]big.Word
	var prod [8]big.Word
	var lf, rt [4]big.Word
	var mod [9]big.Word // should be 9, 4 word for mod, 5 word for mu
	copy(lf[:], left)
	copy(rt[:], right)
	copy(mod[:], mdU)
	C.vli_mult((*C.u64)(unsafe.Pointer(&prod[0])),
		(*C.u64)(unsafe.Pointer(&lf[0])),
		(*C.u64)(unsafe.Pointer(&rt[0])), 4)
	C.vli_mmod_barrett((*C.u64)((unsafe.Pointer)(&res[0])),
		(*C.u64)(unsafe.Pointer(&prod[0])),
		(*C.u64)(unsafe.Pointer(&mod[0])), 4)
	result = new(big.Int).SetBits(res[:4])
	return
}

func vliModMultBarrett(left, right *big.Int, mdU []big.Word) (result *big.Int) {
	var res [4]big.Word
	prod := new(big.Int).Mul(left, right)
	var prd [8]big.Word
	var mod [9]big.Word // should be 9, 4 word for mod, 5 word for mu
	copy(prd[:], prod.Bits())
	copy(mod[:], mdU)
	C.vli_mmod_barrett((*C.u64)((unsafe.Pointer)(&res[0])),
		(*C.u64)(unsafe.Pointer(&prd[0])),
		(*C.u64)(unsafe.Pointer(&mod[0])), 4)
	result = new(big.Int).SetBits(res[:4])
	return
}

func vliBarrettDiv(prod *big.Int, muB []big.Word) (result *big.Int) {
	var res [8]big.Word
	var prd [8]big.Word
	var mu [5]big.Word // should be 9, 4 word for mod, 5 word for mu
	copy(prd[:], prod.Bits())
	copy(mu[:], muB)
	C.vli_div_barrett((*C.u64)(unsafe.Pointer(&res[0])),
		(*C.u64)(unsafe.Pointer(&prd[0])),
		(*C.u64)(unsafe.Pointer(&mu[0])), 4)
	result = new(big.Int).SetBits(res[:4])
	return
}

func vliModMultMont(x, y, mod []big.Word, rr []uint64, k0 uint64) (res *big.Int) {
	var r [4]big.Word
	C.mont_MulMod((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])),
		(*C.u64)(unsafe.Pointer(&mod[0])), (*C.u64)(unsafe.Pointer(&rr[0])),
		C.u64(k0))
	res = new(big.Int).SetBits(r[:4])
	return
}

func vliFromBE64(src []byte) (dest []C.u64) {
	var res [4]C.u64
	var ss [32]byte
	copy(ss[:], src) // process 32 bytes
	C.vli_from_be64(&res[0], unsafe.Pointer(&ss[0]), 4)
	dest = res[:]
	return
}
