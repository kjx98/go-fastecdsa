// Copyright 2010 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package elliptic implements several standard elliptic curves over prime
// fields.
package btc

// This package operates, internally, on Jacobian coordinates. For a given
// (x, y) position on the curve, the Jacobian coordinates are (x1, y1, z1)
// where x = x1/z1² and y = y1/z1³. The greatest speedups come when the whole
// calculation can be performed within the transform (as in ScalarMult and
// ScalarBaseMult). But even for Add and Double, it's faster to apply and
// reverse the transform than to operate in affine coordinates.

import (
	"errors"
	"math/big"
)

const (
	maxBits  = 256
	maxWords = 4
)

var (
	errTooLarge = errors.New("Diff too large, > 64Bits")
)

// CalcMu   mu = b^2k / p
func CalcMu(p *big.Int) *big.Int {
	n2k := new(big.Int).SetUint64(1)
	n2k.Lsh(n2k, maxBits*2)
	res := new(big.Int).Div(n2k, p)
	return res
}

func fastDivBaseExp(x *big.Int, nExp int) *big.Int {
	var res big.Int
	qb := x.Bits()
	if len(qb) < nExp {
		return x
	}
	res.SetBits(qb[nExp:])
	return &res
}

func BarrettDiv(prod, mu *big.Int) (r *big.Int) {
	var qq big.Int
	q1 := fastDivBaseExp(prod, maxWords-1)
	qq.Mul(q1, mu)
	r = fastDivBaseExp(&qq, maxWords+1)
	return
}

func DiffInt(x, y *big.Int) (delta int64, err error) {
	var dd big.Int
	dd.Sub(x, y)
	if len(dd.Bits()) > 1 {
		return 0, errTooLarge
	}
	delta = dd.Int64()
	return
}
