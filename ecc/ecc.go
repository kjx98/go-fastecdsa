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
	"crypto/elliptic"
	"math/big"
	"unsafe"
)

// A Curve represents a short-form Weierstrass curve with a=-3.
// See https://www.hyperelliptic.org/EFD/g1p/auto-shortw.html
type Curve = elliptic.Curve

// CurveParams contains the parameters of an elliptic curve and also provides
// a generic, non-constant time implementation of Curve.
type CurveParams = elliptic.CurveParams

// eccCurve contains the parameters of an elliptic curve and also provides
// a generic, non-constant time implementation of Curve.
type eccCurve struct {
	*CurveParams
	hnd    C.CURVE_HND
	inited bool
}

var sm2c eccCurve
var bigOne *big.Int

func init() {
	if sm2c.inited {
		return
	}
	sm2c.inited = true
	bigOne = big.NewInt(1)
	cHnd := C.get_curve(C.ECC_CURVE_SM2)
	if cHnd == C.CURVE_HND(uintptr(0)) {
		return
	}
	sm2c.hnd = cHnd
	sm2c.CurveParams = getCurveParams(0)
}

func SM2C() *eccCurve {
	if sm2c.CurveParams == nil {
		return nil
	}
	return &sm2c
}

func vliTestFMA() bool {
	return C.vli_asm_acc() != 0
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

func toWordSlice(x C.felem) []big.Word {
	pt := (*[4]big.Word)(unsafe.Pointer(&x))
	return pt[:]
}

func fromWordSlice(bits []big.Word) *C.felem {
	var bb [4]big.Word
	copy(bb[:], bits)
	return (*C.felem)(unsafe.Pointer(&bb[0]))
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

/*
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
*/

func (c eccCurve) AffineFromJacobian(x, y, z *big.Int) (xOut, yOut *big.Int) {
	var xb, yb [4]big.Word
	pt := c.newPoint(x, y, z)
	C.affine_from_jacobian((*C.u64)(unsafe.Pointer(&xb[0])),
		(*C.u64)(unsafe.Pointer(&yb[0])), pt, c.hnd)
	xOut = new(big.Int).SetBits(xb[:])
	yOut = new(big.Int).SetBits(yb[:])
	return
}

// Functions implemented in ecc_asm_*64.s
// Montgomery inverse modulo P256
func vliModInv(in, mod []big.Word) (result []big.Word) {
	var res C.felem
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
	C.mont_mod_mult((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])),
		(*C.u64)(unsafe.Pointer(&mod[0])), (*C.u64)(unsafe.Pointer(&rr[0])),
		C.u64(k0))
	return new(big.Int).SetBits(r[:4])
}

func vliExpModMont(x, y, mod []big.Word, rr []uint64, k0 uint64) *big.Int {
	var r [4]big.Word
	C.mont_mod_exp((*C.u64)(unsafe.Pointer(&r[0])),
		(*C.u64)(unsafe.Pointer(&x[0])), (*C.u64)(unsafe.Pointer(&y[0])),
		(*C.u64)(unsafe.Pointer(&mod[0])), (*C.u64)(unsafe.Pointer(&rr[0])),
		C.u64(k0))
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
