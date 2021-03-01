// +build amd64 arm64

package btc

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
	BTCgo()
}

func asmMontMult(x, y *big.Int) *big.Int {
	var xp, yp [4]uint64
	var res [4]uint64
	fromBig(xp[:], x)
	fromBig(yp[:], y)
	p256Mul(res[:], xp[:], yp[:])
	return toBig(res[:])
}

func asmMontSqr(x *big.Int) *big.Int {
	var xp [4]uint64
	var res [4]uint64
	fromBig(xp[:], x)
	// toMontgomery form
	p256Mul(xp[:], xp[:], rr)
	p256Sqr(res[:], xp[:], 1)
	p256FromMont(res[:], res[:])
	return toBig(res[:])
}

func asmOrdMult(x, y *big.Int) *big.Int {
	var xp, yp [4]uint64
	var res [4]uint64
	fromBig(xp[:], x)
	fromBig(yp[:], y)
	p256OrdMul(xp[:], xp[:], nRR)
	p256OrdMul(yp[:], yp[:], nRR)
	p256OrdMul(res[:], xp[:], yp[:])
	p256OrdMul(res[:], res[:], one)
	return toBig(res[:])
}

func asmOrdSqr(x *big.Int) *big.Int {
	var xp [4]uint64
	var res [4]uint64
	fromBig(xp[:], x)
	p256OrdMul(xp[:], xp[:], nRR)
	p256OrdSqr(res[:], xp[:], 1)
	p256OrdMul(res[:], res[:], one)
	return toBig(res[:])
}

func asmMontRed(y *big.Int) *big.Int {
	var yp [4]uint64
	var res [4]uint64
	fromBig(yp[:], y)
	p256FromMont(res[:], yp[:])
	return toBig(res[:])
}

func asmMontMul(x, y *big.Int) *big.Int {
	xp := asmMontMult(x, btcg.rr)
	yp := asmMontMult(y, btcg.rr)
	res := asmMontMult(xp, yp)
	return asmMontRed(res)
}

func asmSqrt(y *big.Int) (rt *big.Int) {
	var yp [4]uint64
	var res [4]uint64
	fromBig(yp[:], y)
	if p256Sqrt(res[:], yp[:]) {
		rt = toBig(res[:])
	}
	return
}

func TestAsmMontMul(t *testing.T) {
	prod := new(big.Int).Mul(x1, y1)
	m1 := new(big.Int).Mod(prod, btcg.P)
	m2 := asmMontMul(x1, y1)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontMulMod step 1 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	prod = new(big.Int).Mul(x2, y2)
	m1 = new(big.Int).Mod(prod, btcg.P)
	m2 = asmMontMul(x2, y2)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontMulMod step2 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
}

func TestAsmMontSqr(t *testing.T) {
	prod := new(big.Int).Mul(x1, x1)
	m1 := new(big.Int).Mod(prod, btcg.P)
	m2 := asmMontSqr(x1)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontSqrMod step 1 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	if m3 := asmSqrt(m2); m3 == nil {
		t.Log("Can't get Sqrt")
		t.Fail()
	} else if m3.Cmp(x1) != 0 {
		t.Logf("ModSqrt step 1 diff:\n%s vs\n%s", x1.Text(16), m3.Text(16))
		t.Fail()
	}
	prod = new(big.Int).Mul(x2, x2)
	m1 = new(big.Int).Mod(prod, btcg.P)
	m2 = asmMontSqr(x2)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontSqrMod step2 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	if m3 := asmSqrt(m2); m3 == nil {
		t.Log("Can't get Sqrt")
		t.Fail()
	} else if m3.Cmp(x2) != 0 {
		t.Logf("ModSqrt step 2 diff:\n%s vs\n%s", x1.Text(16), m3.Text(16))
		t.Fail()
	}
}

func TestAsmOrdMul(t *testing.T) {
	prod := new(big.Int).Mul(x1, y1)
	m1 := new(big.Int).Mod(prod, btcg.N)
	m2 := asmOrdMult(x1, y1)
	if m1.Cmp(m2) != 0 {
		t.Logf("OrdMulMod step 1 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	prod = new(big.Int).Mul(x2, y2)
	m1 = new(big.Int).Mod(prod, btcg.N)
	m2 = asmOrdMult(x2, y2)
	if m1.Cmp(m2) != 0 {
		t.Logf("OrdMulMod step2 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
}

func TestAsmOrdSqr(t *testing.T) {
	prod := new(big.Int).Mul(x1, x1)
	m1 := new(big.Int).Mod(prod, btcg.N)
	//m2 := asmOrdMult(x1, x1)
	m2 := asmOrdSqr(x1)
	if m1.Cmp(m2) != 0 {
		t.Logf("OrdSqrMod step 1 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	prod = new(big.Int).Mul(x2, x2)
	m1 = new(big.Int).Mod(prod, btcg.N)
	//m2 = asmOrdMult(x2, x2)
	m2 = asmOrdSqr(x2)
	if m1.Cmp(m2) != 0 {
		t.Logf("OrdSqrMod step2 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
}

func TestAsmInverse(t *testing.T) {
	p := btcg.P
	RR := new(big.Int).SetUint64(1)
	RR.Lsh(RR, 257)
	RR.Mod(RR, p)
	Rinv := new(big.Int).ModInverse(RR, p)
	var res, yy [4]uint64
	xp := asmMontMult(RR, btcg.rr)
	fromBig(yy[:], xp)
	p256Inverse(res[:], yy[:])
	p256FromMont(res[:], res[:])
	RinvA := toBig(res[:])
	if Rinv.Cmp(RinvA) != 0 {
		t.Logf("p256Inverse diff:\n%s vs\n%s", Rinv.Text(16), RinvA.Text(16))
		t.Fail()
	}
}

func TestPointRecover(t *testing.T) {
	c := pBTC
	px, py := c.ScalarBaseMult(d1.Bytes())
	v := py.Bit(0)
	if py2, err := RecoverPoint(px, v); err != nil {
		t.Log("Can't recover pointY, error:", err)
		t.Fail()
	} else if py2.Cmp(py) != 0 {
		t.Logf("RecoverPoint diff:\n%s vs\n%s", py.Text(16), py2.Text(16))
		t.Fail()
	}
	px, py = c.ScalarBaseMult(d2.Bytes())
	v = py.Bit(0)
	if py2, err := RecoverPoint(px, v); err != nil {
		t.Log("Can't recover step2 pointY, error:", err)
		t.Fail()
	} else if py2.Cmp(py) != 0 {
		t.Logf("RecoverPoint step2 diff:\n%s vs\n%s", py.Text(16), py2.Text(16))
		t.Fail()
	}
}

func BenchmarkAsmInverse(b *testing.B) {
	priv, _ := fastecdsa.GenerateKey(P256(), rand.Reader)
	var res, yy [4]uint64
	fromBig(yy[:], priv.PublicKey.X)

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256Inverse(res[:], yy[:])
		}
	})
}

func BenchmarkAsmMontModMul(b *testing.B) {
	//c := btcg
	xp := x1 //toMont(x1)

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_ = asmMontMult(xp, xp)
		}
	})
}

func BenchmarkAsmMontSqr(b *testing.B) {
	var res, xp [4]uint64
	fromBig(xp[:], x1)
	p256Mul(res[:], xp[:], rr)
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256Sqr(res[:], res[:], 1)
		}
	})
}

func BenchmarkAsmOrdMul(b *testing.B) {
	var res, xp [4]uint64
	fromBig(xp[:], x1)
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256OrdMul(res[:], xp[:], nRR)
		}
	})
}

func BenchmarkAsmOrdSqr(b *testing.B) {
	var res, xp [4]uint64
	fromBig(xp[:], x1)
	p256OrdMul(res[:], xp[:], nRR)
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256OrdSqr(res[:], res[:], 1)
		}
	})
}

func BenchmarkPointRecover(b *testing.B) {
	c := pBTC
	px, py := c.ScalarBaseMult(d1.Bytes())
	v := py.Bit(0)
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			RecoverPoint(px, v)
		}
	})
}

func BenchmarkAsmECMULT(b *testing.B) {
	Curve := BTC()
	goGx := Curve.Params().Gx
	goGy := Curve.Params().Gy

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = Curve.ScalarMult(goGx, goGy, d1.Bytes())
		}
	})
}

func BenchmarkAsmECGMULT(b *testing.B) {
	Curve := BTC()

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = Curve.ScalarBaseMult(d1.Bytes())
		}
	})
}
