// +build amd64 arm64

package ecc

// #cgo CXXFLAGS: -O3 -Wpedantic -Wall -Wno-maybe-uninitialized -std=gnu++11
// #include "ecc.h"
import "C"

// for gcc 4.8.5, __builtin_add/sub_overflow improve 10% -- 20%
// -fvar-tracking-assignments
// cgo CXXFLAGS: -O2 -Wpedantic -Wall -std=gnu++11 -DNO_BUILTIN_OVERFLOW
// WITH_SM2_MULTP using shift instead of multiply, no improvement @x86_64/arm64
// cgo CXXFLAGS: -O2 -Wpedantic -Wall -std=gnu++11 -DWITH_SM2_MULTP

import (
	"math/big"
	"unsafe"
)

func vliTestFMA() bool {
	return C.vli_asm_acc() != 0
}

func toWordSlice(x C.fElem) []big.Word {
	pt := (*[4]big.Word)(unsafe.Pointer(&x))
	return pt[:]
}

func fromWordSlice(bits []big.Word) *C.fElem {
	var bb [4]big.Word
	copy(bb[:], bits)
	return (*C.fElem)(unsafe.Pointer(&bb[0]))
}

func fillMontParams(mod []big.Word, rr []uint64, k0 uint64) *C.montParams {
	var pa C.montParams
	for i := 0; i < 4; i++ {
		if i < len(mod) {
			pa.p[i] = C.u64(mod[i])
		}
		if i < len(rr) {
			pa.rr[i] = C.u64(rr[i])
		}
	}
	pa.k0 = C.u64(k0)
	return &pa
}

// Functions implemented in ecc_asm_*64.s
// Montgomery inverse modulo P256
func vliModInv(in, mod []big.Word) (result []big.Word) {
	var res C.fElem
	C.vli_mod_inv((*C.u64)(unsafe.Pointer(&res)),
		(*C.u64)(unsafe.Pointer(fromWordSlice(in))),
		(*C.u64)(unsafe.Pointer(fromWordSlice(mod))))
	result = toWordSlice(res)
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
		(*C.u64)(unsafe.Pointer(&rt[0])))
	C.vli_mmod_barrett((*C.u64)((unsafe.Pointer)(&res[0])),
		(*C.u64)(unsafe.Pointer(&prod[0])),
		(*C.u64)(unsafe.Pointer(&mod[0])))
	result = new(big.Int).SetBits(res[:4])
	return
}

func vliModMultBarrett(left, right *big.Int, mdU []big.Word) *big.Int {
	var res [4]big.Word
	prod := new(big.Int).Mul(left, right)
	var prd [8]big.Word
	var mod [9]big.Word // should be 9, 4 word for mod, 5 word for mu
	copy(prd[:], prod.Bits())
	copy(mod[:], mdU)
	C.vli_mmod_barrett((*C.u64)((unsafe.Pointer)(&res[0])),
		(*C.u64)(unsafe.Pointer(&prd[0])),
		(*C.u64)(unsafe.Pointer(&mod[0])))
	return new(big.Int).SetBits(res[:4])
}

func vliBarrettDiv(prod *big.Int, muB []big.Word) (result *big.Int) {
	var res [8]big.Word
	var prd [8]big.Word
	var mu [5]big.Word // should be 9, 4 word for mod, 5 word for mu
	copy(prd[:], prod.Bits())
	copy(mu[:], muB)
	C.vli_div_barrett((*C.u64)(unsafe.Pointer(&res[0])),
		(*C.u64)(unsafe.Pointer(&prd[0])),
		(*C.u64)(unsafe.Pointer(&mu[0])))
	result = new(big.Int).SetBits(res[:4])
	return
}

func vliModMultMont(x, y, mod []big.Word, rr []uint64, k0 uint64) *big.Int {
	var r [4]big.Word
	pa := fillMontParams(mod, rr, k0)
	C.mont_mod_mult((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])), pa)
	return new(big.Int).SetBits(r[:4])
}

func vliExpModMont(x, y, mod []big.Word, rr []uint64, k0 uint64) *big.Int {
	var r [4]big.Word
	pa := fillMontParams(mod, rr, k0)
	C.mont_mod_exp((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])), pa)
	return new(big.Int).SetBits(r[:4])
}

func vliSM2MultP(u uint64) *big.Int {
	r := make([]big.Word, 6)
	//can't convert []big.Word to slice_t
	//C.vli_sm2_mult_p(r[:], C.u64(u))
	C.vli_sm2_mult_p((*C.u64)(unsafe.Pointer(&r[0])), C.u64(len(r)), C.u64(u))
	if r[4] == 0 {
		return new(big.Int).SetBits(r[:4])
	}
	return new(big.Int).SetBits(r[:5])
}

/*
func vliFromBE64(src []byte) (dest []C.u64) {
	var res [4]C.u64
	var ss [32]byte
	copy(ss[:], src) // process 32 bytes
	C.vli_from_be64(&res[0], unsafe.Pointer(&ss[0]), 4)
	dest = res[:]
	return
}
*/
