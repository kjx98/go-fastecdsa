// +build amd64 arm64

package ecc

// #cgo CFLAGS: -O2 -Wpedantic -Wall -std=gnu11
// #include <sys/types.h>
// #include "ecc.h"
/*
static struct ecc_point *setEC_point(struct ecc_point *pt, u_int64_t *x, u_int64_t *y) {
	if (pt == NULL) return NULL;
	pt->x[0] = x[0]; pt->x[1] = x[1]; pt->x[2] = x[2]; pt->x[3] = x[3];
	pt->y[0] = y[0]; pt->y[1] = y[1]; pt->y[2] = y[2]; pt->y[3] = y[3];
	return pt;
}
*/
import "C"

import (
	"math/big"
	"unsafe"
)

// Functions implemented in ecc_asm_*64.s
// Montgomery addition modulo P256
func vliSub(left, right []uint64) (result []uint64, borrow bool) {
	var res [4]uint64
	ret := C.vli_sub((*C.u64)(&res[0]), (*C.u64)(&left[0]), (*C.u64)(&right[0]), 4)
	result = res[:]
	if ret != 0 {
		borrow = true
	}
	return
}

// Montgomery inverse modulo P256
func vliModInv(input, mod []byte) (result []uint64) {
	var res [4]uint64
	in := vliFromBE64(input)
	modP := vliFromBE64(mod)
	C.vli_mod_inv((*C.u64)(&res[0]), (*C.u64)(&in[0]), (*C.u64)(&modP[0]), 4)
	result = res[:]
	return
}

// Montgomery multiplication modulo P256
func vliModMult(left, right, mod []byte) (result *big.Int) {
	var res [4]big.Word
	lf := vliFromBE64(left)
	rt := vliFromBE64(right)
	modP := vliFromBE64(mod)
	C.vli_mod_mult_slow((*C.u64)(unsafe.Pointer(&res[0])), (*C.u64)(&lf[0]),
		(*C.u64)(&rt[0]), (*C.u64)(&modP[0]), 4)
	result = new(big.Int).SetBits(res[:])
	return
}

func vliModMultBarrett(left, right []byte, mod []big.Word) (result *big.Int) {
	var res [8]big.Word
	lf := vliFromBE64(left)
	rt := vliFromBE64(right)
	//modP := vliFromBE64(mod)
	C.vli_mod_mult_fast((*C.u64)((unsafe.Pointer)(&res[0])), (*C.u64)(&lf[0]),
		(*C.u64)(&rt[0]), (*C.u64)((unsafe.Pointer)(&mod[0])), 4)
	result = new(big.Int).SetBits(res[:4])
	return
}

func vliBarrettDiv(prod *big.Int, mu []big.Word) (result *big.Int) {
	var res [8]big.Word
	prd := vliFromBE64(prod.Bytes())
	//modP := vliFromBE64(mod)
	C.vli_div_barrett((*C.u64)((unsafe.Pointer)(&res[0])), (*C.u64)(&prd[0]),
		(*C.u64)((unsafe.Pointer)(&mu[0])), 4)
	result = new(big.Int).SetBits(res[:4])
	return
}

func vliFromBE64(src []byte) (dest []uint64) {
	var res [4]uint64
	C.vli_from_be64((*C.u64)(&res[0]), unsafe.Pointer(&src[0]), 4)
	dest = res[:]
	return
}
