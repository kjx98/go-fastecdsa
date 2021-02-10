// +build amd64 arm64

package ecc

// #cgo CXXFLAGS: -O3 -Wall -I../include -Wuninitialized -std=c++11
// #cgo CFLAGS: -O3 -Wpedantic -I../include -Wno-maybe-uninitialized -std=c11
// #include "ecc.h"
import "C"

// for gcc 4.8.5, __builtin_add/sub_overflow improve 10% -- 20%
// -fvar-tracking-assignments
// WITH_SM2_MULTP using shift instead of multiply, no improvement @x86_64/arm64
// cgo CXXFLAGS: -O2 -Wpedantic -Wall -std=gnu++11 -DWITH_SM2_MULTP

import (
	"math/big"
	"unsafe"
)

func vliTestFMA() bool {
	return C.vli_asm_acc() != 0
}

func toWordSlice(x C.bn_words_t) []big.Word {
	pt := *(*[4]big.Word)(unsafe.Pointer(&x))
	return pt[:]
}

func fromWordSlice(bits []big.Word) *C.bn_words_t {
	var bb [4]big.Word
	copy(bb[:], bits)
	return (*C.bn_words_t)(unsafe.Pointer(&bb[0]))
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
	var res C.bn_words_t
	C.vli_mod_inv((*C.u64)(unsafe.Pointer(&res)),
		(*C.u64)(unsafe.Pointer(fromWordSlice(in))),
		(*C.u64)(unsafe.Pointer(fromWordSlice(mod))))
	result = toWordSlice(res)
	return
}

func vliToMont(x, mod []big.Word, rr []uint64, k0 uint64) *big.Int {
	var r [4]big.Word
	pa := fillMontParams(mod, rr, k0)
	C.to_montgomery((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), pa)
	return new(big.Int).SetBits(r[:4])
}

func vliFromMont(y, mod []big.Word, rr []uint64, k0 uint64) *big.Int {
	var r [4]big.Word
	pa := fillMontParams(mod, rr, k0)
	C.from_montgomery((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&y[0])), pa)
	return new(big.Int).SetBits(r[:4])
}

func vliModMultMont(x, y, mod []big.Word, rr []uint64, k0 uint64) *big.Int {
	var r [4]big.Word
	pa := fillMontParams(mod, rr, k0)
	C.mont_mod_mult((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])), pa)
	return new(big.Int).SetBits(r[:4])
}

func vliModSqrMont(x, mod []big.Word, rr []uint64, k0 uint64) *big.Int {
	var r [4]big.Word
	pa := fillMontParams(mod, rr, k0)
	C.mont_mod_sqr((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), pa, 1)
	return new(big.Int).SetBits(r[:4])
}

func vliExpModMont(x, y, mod []big.Word, rr []uint64, k0 uint64) *big.Int {
	var r [4]big.Word
	pa := fillMontParams(mod, rr, k0)
	C.mont_mod_exp((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])), pa)
	return new(big.Int).SetBits(r[:4])
}

func vliSM2ModMultP(x, y []big.Word) *big.Int {
	var r [4]big.Word
	C.mont_sm2_mod_mult_p((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])))
	return new(big.Int).SetBits(r[:4])
}

func vliSM2ModMultN(x, y []big.Word) *big.Int {
	var r [4]big.Word
	C.mont_sm2_mod_mult_n((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])))
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
