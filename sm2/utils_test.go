package sm2

import (
	"math/big"
	"testing"
)

func TestCalcMu(t *testing.T) {
	n256 := new(big.Int).SetUint64(1)
	n256.Lsh(n256, 256)
	n512 := new(big.Int).Mul(n256, n256)

	cParams := P256().Params()
	p := cParams.P

	mu := CalcMu(p)
	pmu := new(big.Int).Mul(mu, p)
	if pmu.Cmp(n512) >= 0 {
		t.Log("p * mu > 2^512")
		t.Fail()
	}
	pmu.Add(pmu, p)
	if pmu.Cmp(n512) < 0 {
		t.Log("p * mu + p < 2^512")
		t.Fail()
	}
	ww := mu.Bits()
	t.Logf("mu: %x %x %x %x %x", ww[0], ww[1], ww[2], ww[3], ww[4])
	mu1 := new(big.Int).Div(n512, p)
	if mu.Cmp(mu1) != 0 {
		t.Log("mu diff")
		t.Fail()
	}
}

func TestBarrettDiv(t *testing.T) {
	cParams := P256().Params()
	goGx := cParams.Gx
	goGy := cParams.Gy
	prod := new(big.Int).Mul(goGx, d1)
	n64 := new(big.Int).SetUint64(1)
	n64.Lsh(n64, 64)
	tt1 := fastDivBaseExp(prod, 1)
	tt2 := new(big.Int).Div(prod, n64)
	if tt1.Cmp(tt2) != 0 {
		t.Log("fastDivBaseExp diff Div")
		t.Fail()
	}
	p := cParams.P
	mu := CalcMu(p)
	ww := mu.Bits()
	t.Logf("mu: %x %x %x %x %x", ww[0], ww[1], ww[2], ww[3], ww[4])
	q := BarrettDiv(prod, mu)
	qe := new(big.Int).Div(prod, p)
	if dd, err := DiffInt(qe, q); err != nil {
		t.Log("Diff error", err)
		t.Fail()
	} else if dd < 0 {
		t.Log("delta prod/p, q < 0, delta: ", dd)
		t.Fail()
	} else if dd > 2 {
		t.Log("delta > 2, delta: ", dd)
		t.Fail()
	}
	prod = prod.Mul(goGy, d1)
	q = BarrettDiv(prod, mu)
	qe = qe.Div(prod, p)
	if dd, err := DiffInt(qe, q); err != nil {
		t.Log("Diff error", err)
		t.Fail()
	} else if dd < 0 || dd > 2 {
		t.Log("delta prod/p, q < 0 or > 2, delta: ", dd)
		t.Fail()
	}
}
