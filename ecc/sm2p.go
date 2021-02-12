// +build sm2p

package ecc

// #cgo CXXFLAGS: -O3 -Wpedantic -I../include -Wno-uninitialized -std=c++11
// #cgo CFLAGS: -O3 -Wpedantic -I../include -Wno-uninitialized -std=c11
// #include "sm2p.h"
import "C"

// for gcc 4.8.5, __builtin_add/sub_overflow improve 10% -- 20%
// -fvar-tracking-assignments
// WITH_SM2_MULTP using shift instead of multiply, no improvement @x86_64/arm64
// cgo CXXFLAGS: -O2 -Wpedantic -Wall -std=gnu++11 -DWITH_SM2_MULTP

import (
	//"crypto/elliptic"
	"math/big"
	"unsafe"
)

/*
func sm2ModMult(x, y []big.Word) *big.Int {
	var r [4]big.Word
	xb := fromWordSlice(x)
	yb := fromWordSlice(y)
	C.sm2_mod_mul((*C.bn_words_t)(unsafe.Pointer(&r[0])),
		(*[4]C.u64)(unsafe.Pointer(xb)), (*[4]C.u64)(unsafe.Pointer(yb)))
	return new(big.Int).SetBits(r[:4])
}
*/

func newPoint(x, y, z *big.Int) *C.Point {
	var pt C.Point
	if z == nil {
		z = bigOne
	}
	pt.x = *fromWordSlice(x.Bits())
	pt.y = *fromWordSlice(y.Bits())
	pt.z = *fromWordSlice(z.Bits())
	return &pt
}

func sm2ModInv(x []big.Word) *big.Int {
	var r [4]big.Word
	in := fromWordSlice(x)
	C.sm2_mod_inv((*C.u64)(unsafe.Pointer(&r[0])), (*C.u64)(unsafe.Pointer(in)))
	return new(big.Int).SetBits(r[:4])
}

func sm2Add(x1, y1, x2, y2 *big.Int) (rx, ry *big.Int) {
	pt1 := newPoint(x1, y1, nil)
	pt2 := newPoint(x2, y2, nil)
	var pt C.Point
	C.sm2_point_add(&pt, pt1, pt2)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

func sm2AddJacobian(x1, y1, z1, x2, y2, z2 *big.Int) (rx, ry, rz *big.Int) {
	pt1 := newPoint(x1, y1, z1)
	pt2 := newPoint(x2, y2, z2)
	var pt C.Point
	C.sm2_point_add_jacobian(&pt, pt1, pt2)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	rz = new(big.Int).SetBits(toWordSlice(pt.z))
	return
}

func sm2Double(x, y *big.Int) (rx, ry *big.Int) {
	pt1 := newPoint(x, y, nil)
	var pt C.Point
	C.sm2_point_double(&pt, pt1)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

func sm2DoubleJacobian(x, y, z *big.Int) (rx, ry, rz *big.Int) {
	pt1 := newPoint(x, y, z)
	var pt C.Point
	C.sm2_point_double_jacobian(&pt, pt1)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	rz = new(big.Int).SetBits(toWordSlice(pt.z))
	return
}

func sm2ScalarBaseMult(k []byte) (rx, ry *big.Int) {
	var pt C.Point
	var bb [32]byte
	copy(bb[:], k)
	C.sm2_scalar_base_mult(&pt, (*C.u8)(unsafe.Pointer(&bb[0])))
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

func sm2ScalarMult(x, y *big.Int, k []byte) (rx, ry *big.Int) {
	pt1 := newPoint(x, y, nil)
	var pt C.Point
	scal := new(big.Int).SetBytes(k)
	var ss [4]big.Word
	copy(ss[:], scal.Bits())
	C.sm2_scalar_mult(&pt, pt1, (*C.u64)(unsafe.Pointer(&ss[0])))
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

/*
func (c eccCurve) DoubleJacobian(x, y, z *big.Int) (rx, ry, rz *big.Int) {
	pt1 := c.newPoint(x, y, z)
	var pt C.Point
	C.point_double_jacobian(&pt, pt1, c.hnd)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	rz = new(big.Int).SetBits(toWordSlice(pt.z))
	return
}


func (c eccCurve) AffineFromJacobian(x, y, z *big.Int) (xOut, yOut *big.Int) {
	var xb, yb [4]big.Word
	pt := c.newPoint(x, y, z)
	C.affine_from_jacobian((*C.u64)(unsafe.Pointer(&xb[0])),
		(*C.u64)(unsafe.Pointer(&yb[0])), pt, c.hnd)
	xOut = new(big.Int).SetBits(xb[:])
	yOut = new(big.Int).SetBits(yb[:])
	return
}
*/
