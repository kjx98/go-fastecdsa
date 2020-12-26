// +build sm2p256

package ecc

// #cgo CXXFLAGS: -O3 -Wpedantic -Wall -Wno-maybe-uninitialized -std=gnu++11
// #include "sm2p.h"
import "C"

// for gcc 4.8.5, __builtin_add/sub_overflow improve 10% -- 20%
// -fvar-tracking-assignments
// cgo CXXFLAGS: -O2 -Wpedantic -Wall -std=gnu++11 -DNO_BUILTIN_OVERFLOW
// WITH_SM2_MULTP using shift instead of multiply, no improvement @x86_64/arm64
// cgo CXXFLAGS: -O2 -Wpedantic -Wall -std=gnu++11 -DWITH_SM2_MULTP

import (
//"crypto/elliptic"
//"math/big"
//"unsafe"
)

/*
func vliModMultMontP(x, y []big.Word) *big.Int {
	var r [4]big.Word
	C.mont_sm2_mod_mult_p((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])))
	return new(big.Int).SetBits(r[:4])
}

func vliModMultMontN(x, y []big.Word) *big.Int {
	var r [4]big.Word
	C.mont_sm2_mod_mult_n((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])))
	return new(big.Int).SetBits(r[:4])
}

func getCurveParams(curveId uint) *CurveParams {
	var cHnd C.CURVE_HND
	if curveId == 0 {
		cHnd = sm2c.hnd
	} else {
		cHnd = C.get_curve(C.uint(curveId))
	}
	if cHnd == C.CURVE_HND(uintptr(0)) {
		return nil
	}
	var p, n, b, gx, gy [4]big.Word
	C.get_curve_params((*C.u64)(unsafe.Pointer(&p[0])),
		(*C.u64)(unsafe.Pointer(&n[0])), (*C.u64)(unsafe.Pointer(&b[0])),
		(*C.u64)(unsafe.Pointer(&gx[0])), (*C.u64)(unsafe.Pointer(&gy[0])),
		cHnd)
	sm2Params := &CurveParams{Name: "SM2-c"}
	sm2Params.P = new(big.Int).SetBits(p[:])
	sm2Params.N = new(big.Int).SetBits(n[:])
	sm2Params.B = new(big.Int).SetBits(b[:])
	sm2Params.Gx = new(big.Int).SetBits(gx[:])
	sm2Params.Gy = new(big.Int).SetBits(gy[:])
	sm2Params.BitSize = 256
	return sm2Params
}

func (c eccCurve) newPoint(x, y, z *big.Int) *C.Point {
	var pt C.Point
	if z == nil {
		z = bigOne
	}
	pt.x = *fromWordSlice(x.Bits())
	pt.y = *fromWordSlice(y.Bits())
	pt.z = *fromWordSlice(z.Bits())
	return &pt
}

func (c eccCurve) Add(x1, y1, x2, y2 *big.Int) (rx, ry *big.Int) {
	pt1 := c.newPoint(x1, y1, nil)
	pt2 := c.newPoint(x2, y2, nil)
	var pt C.Point
	C.point_add(&pt, pt1, pt2, c.hnd)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

func (c eccCurve) AddJacobian(x1, y1, z1, x2, y2, z2 *big.Int) (rx, ry, rz *big.Int) {
	pt1 := c.newPoint(x1, y1, z1)
	pt2 := c.newPoint(x2, y2, z2)
	var pt C.Point
	C.point_add_jacobian(&pt, pt1, pt2, c.hnd)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	rz = new(big.Int).SetBits(toWordSlice(pt.z))
	return
}

func (c eccCurve) DoubleJacobian(x, y, z *big.Int) (rx, ry, rz *big.Int) {
	pt1 := c.newPoint(x, y, z)
	var pt C.Point
	C.point_double_jacobian(&pt, pt1, c.hnd)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	rz = new(big.Int).SetBits(toWordSlice(pt.z))
	return
}

func (c eccCurve) Double(x, y *big.Int) (rx, ry *big.Int) {
	pt1 := c.newPoint(x, y, nil)
	var pt C.Point
	C.point_double(&pt, pt1, c.hnd)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
	return
}

func (c eccCurve) ScalarMult(x, y *big.Int, k []byte) (rx, ry *big.Int) {
	pt1 := c.newPoint(x, y, nil)
	var pt C.Point
	scal := new(big.Int).SetBytes(k)
	var ss [4]big.Word
	copy(ss[:], scal.Bits())
	C.point_mult(&pt, pt1, (*C.u64)(unsafe.Pointer(&ss[0])), c.hnd)
	rx = new(big.Int).SetBits(toWordSlice(pt.x))
	ry = new(big.Int).SetBits(toWordSlice(pt.y))
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
