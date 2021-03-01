// Copyright 2010 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package elliptic implements several standard elliptic curves over prime
// fields.
package sm2

// This package operates, internally, on Jacobian coordinates. For a given
// (x, y) position on the curve, the Jacobian coordinates are (x1, y1, z1)
// where x = x1/z1² and y = y1/z1³. The greatest speedups come when the whole
// calculation can be performed within the transform (as in ScalarMult and
// ScalarBaseMult). But even for Add and Double, it's faster to apply and
// reverse the transform than to operate in affine coordinates.

import (
	"crypto/elliptic"
	"math/big"
	"sync"
)

// A Curve represents a short-form Weierstrass curve with a=-3.
// See https://www.hyperelliptic.org/EFD/g1p/auto-shortw.html
type Curve = elliptic.Curve

// CurveParams contains the parameters of an elliptic curve and also provides
// a generic, non-constant time implementation of Curve.
type CurveParams = elliptic.CurveParams

// zForAffine returns a Jacobian Z value for the affine point (x, y). If x and
// y are zero, it assumes that they represent the point at infinity because (0,
// 0) is not on the any of the curves handled here.
func zForAffine(x, y *big.Int) *big.Int {
	z := new(big.Int)
	if x.Sign() != 0 || y.Sign() != 0 {
		z.SetInt64(1)
	}
	return z
}

// Marshal converts a point into the uncompressed form specified in section 4.3.6 of ANSI X9.62.
func Marshal(curve Curve, x, y *big.Int) []byte {
	byteLen := (curve.Params().BitSize + 7) >> 3

	ret := make([]byte, 1+2*byteLen)
	ret[0] = 4 // uncompressed point

	xBytes := x.Bytes()
	copy(ret[1+byteLen-len(xBytes):], xBytes)
	yBytes := y.Bytes()
	copy(ret[1+2*byteLen-len(yBytes):], yBytes)
	return ret
}

// Unmarshal converts a point, serialized by Marshal, into an x, y pair.
// It is an error if the point is not in uncompressed form or is not on the curve.
// On error, x = nil.
func Unmarshal(curve Curve, data []byte) (x, y *big.Int) {
	byteLen := (curve.Params().BitSize + 7) >> 3
	if len(data) == 1+byteLen {
		if data[0] != 0x2 && data[0] != 0x3 {
			return
		}
		p := curve.Params().P
		x = new(big.Int).SetBytes(data[1:])
		if x.Cmp(p) >= 0 {
			return nil, nil
		}
		x3 := new(big.Int).Mul(x, x)
		x3.Mul(x3, x)
		threeX := new(big.Int).Lsh(x, 1)
		threeX.Add(threeX, x)
		x3.Sub(x3, threeX)
		x3.Add(x3, curve.Params().B)
		x3.Mod(x3, p)
		if y = new(big.Int).ModSqrt(x3, p); y == nil {
			return nil, nil
		}
		switch data[0] {
		case 0x2: // should be even
			if y.Bit(0) != 0 {
				y = y.Sub(p, y)
			}
		case 0x3: // should be odd
			if y.Bit(0) == 0 {
				y = y.Sub(p, y)
			}
		}
		if !curve.IsOnCurve(x, y) {
			return nil, nil
		}
		return
	}
	if len(data) != 1+2*byteLen {
		return
	}
	if data[0] != 4 { // uncompressed form
		return
	}
	p := curve.Params().P
	x = new(big.Int).SetBytes(data[1 : 1+byteLen])
	y = new(big.Int).SetBytes(data[1+byteLen:])
	if x.Cmp(p) >= 0 || y.Cmp(p) >= 0 {
		return nil, nil
	}
	if !curve.IsOnCurve(x, y) {
		return nil, nil
	}
	return
}

var initonce sync.Once
var sm2Params = &CurveParams{Name: "SM2"}

var secp256k1Params = &CurveParams{Name: "BC"}

func init() {
	sm2Params.P, _ = new(big.Int).SetString("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF", 16)
	sm2Params.N, _ = new(big.Int).SetString("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123", 16)
	sm2Params.B, _ = new(big.Int).SetString("28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93", 16)
	sm2Params.Gx, _ = new(big.Int).SetString("32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7", 16)
	sm2Params.Gy, _ = new(big.Int).SetString("BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0", 16)
	sm2Params.BitSize = 256
	secp256k1Params.P, _ = new(big.Int).SetString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16)
	secp256k1Params.N, _ = new(big.Int).SetString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16)
	secp256k1Params.B, _ = new(big.Int).SetString("07", 16)
	secp256k1Params.Gx, _ = new(big.Int).SetString("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16)
	secp256k1Params.Gy, _ = new(big.Int).SetString("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16)
	secp256k1Params.BitSize = 256
}

func initAll() {
	initSM2()
	initSM2go()
}

// SM2 returns a Curve which implements sm2
//
// The cryptographic operations are implemented using constant-time algorithms.
func SM2() Curve {
	initonce.Do(initAll)
	return pSM2
}

func SM2asm() p256Curve {
	initonce.Do(initAll)
	return pSM2
}

func SM2go() Curve {
	initonce.Do(initAll)
	return sm2g
}

// P256 returns a Curve which implements sm2 via elliptic package
//
// The cryptographic operations do not use constant-time algorithms.
func P256() Curve {
	initonce.Do(initAll)
	return sm2Params
}

func BTC() Curve {
	initonce.Do(initAll)
	return secp256k1Params
}
