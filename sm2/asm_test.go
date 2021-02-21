// +build amd64 arm64

package sm2

import (
	"crypto/rand"
	"gitee.com/jkuang/go-fastecdsa"
	"math/big"
	"testing"
)

//func p256Sqr(res, in []uint64, n int)
// p256Inverse sets out to in^-1 mod p.
//func p256Inverse(out, in []uint64)
func init() {
	SM2go()
}

func toBig(x []uint64) *big.Int {
	var		res	[4]big.Word
	res[0] = big.Word(x[0])
	res[1] = big.Word(x[1])
	res[2] = big.Word(x[2])
	res[3] = big.Word(x[3])
	return new(big.Int).SetBits(res[:])
}

func asmMontMult(x, y *big.Int) *big.Int {
	var xp, yp		[4]uint64
	var	res			[4]uint64
	fromBig(xp[:], x)
	fromBig(yp[:], y)
	p256Mul(res[:], xp[:], yp[:])
	return toBig(res[:])
}

func asmMontSqr(x *big.Int) *big.Int {
	var xp,rp	[4]uint64
	var	res		[4]uint64
	fromBig(xp[:], x)
	fromBig(rp[:], sm2g.rr)
	p256Mul(res[:], xp[:], rp[:])
	p256Sqr(res[:], res[:], 1)
	p256FromMont(res[:], res[:])
	return toBig(res[:])
}

func asmMontRed(y *big.Int) *big.Int {
	var yp		[4]uint64
	var	res		[4]uint64
	fromBig(yp[:], y)
	p256FromMont(res[:], yp[:])
	return toBig(res[:])
}

func asmMontMul(x,y *big.Int) *big.Int {
	xp := asmMontMult(x, sm2g.rr)
	yp := asmMontMult(y, sm2g.rr)
	res := asmMontMult(xp, yp)
	return  asmMontRed(res)
}

func TestAsmMontMul(t *testing.T) {
	prod := new(big.Int).Mul(x1, y1)
	m1 := new(big.Int).Mod(prod, sm2g.P)
	m2 := asmMontMul(x1, y1)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontMulMod step 1 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	/*
		ww := x1.Bits()
		t.Logf("x1: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
		ww = y1.Bits()
		t.Logf("y1: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
		ww = m1.Bits()
		t.Logf("xy1mod: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	*/
	prod = new(big.Int).Mul(x2, y2)
	m1 = new(big.Int).Mod(prod, sm2g.P)
	m2 = asmMontMul(x2, y2)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontMulMod step2 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	/*
		ww = x2.Bits()
		t.Logf("x2: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
		ww = y2.Bits()
		t.Logf("y2: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
		ww = m1.Bits()
		t.Logf("xy2mod: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	*/
}

func TestAsmMontSqr(t *testing.T) {
	prod := new(big.Int).Mul(x1, x1)
	m1 := new(big.Int).Mod(prod, sm2g.P)
	m2 := asmMontSqr(x1)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontSqrMod step 1 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	prod = new(big.Int).Mul(x2, x2)
	m1 = new(big.Int).Mod(prod, sm2g.P)
	m2 = asmMontSqr(x2)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontSqrMod step2 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
}

func TestAsmInverse(t *testing.T) {
	p := sm2g.P
	RR := new(big.Int).SetUint64(1)
	RR.Lsh(RR, 257)
	RR.Mod(RR, p)
	Rinv := new(big.Int).ModInverse(RR, p)
	var		res, yy	[4]uint64
	xp := asmMontMult(RR, sm2g.rr)
	fromBig(yy[:], xp)
	p256Inverse(res[:], yy[:])
	p256FromMont(res[:], res[:])
	RinvA := toBig(res[:])
	if Rinv.Cmp(RinvA) != 0 {
		t.Logf("p256Inverse diff:\n%s vs\n%s", Rinv.Text(16), RinvA.Text(16))
		t.Fail()
	}
}

func BenchmarkAsmInverse(b *testing.B) {
	b.ResetTimer()
	priv, _ := fastecdsa.GenerateKey(P256(), rand.Reader)
	var		res, yy		[4]uint64
	fromBig(yy[:], priv.PublicKey.X)

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256Inverse(res[:], yy[:])
		}
	})
}

func BenchmarkAsmMontModMul(b *testing.B) {
	b.ResetTimer()
	c := sm2g

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = asmMontMult(x1, c.rr)
	}
}

func BenchmarkAsmECMULT(b *testing.B) {
	b.ResetTimer()
	Curve := SM2()
	goGx := Curve.Params().Gx
	goGy := Curve.Params().Gy

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = Curve.ScalarMult(goGx, goGy, d1.Bytes())
		}
	})
}
