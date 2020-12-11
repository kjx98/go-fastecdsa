package sm2

import (
	"crypto/elliptic"
	"crypto/rand"
	"gitee.com/jkuang/go-fastecdsa"
	"math/big"
	"testing"
)

var (
	x1, y1 *big.Int
	x2, y2 *big.Int
	d1, d2 *big.Int
)

func init() {
	x1, _ = new(big.Int).SetString("12504736780261570232047283225790923734936736733904083456330261180763076903173", 16)
	y1, _ = new(big.Int).SetString("81307578426096630754000956949482484966368962694792140728234975018259774469569", 16)
	x2, _ = new(big.Int).SetString("48121564271922987841895377752074498583012355812029682461364979458450873405695", 16)
	y2, _ = new(big.Int).SetString("37446874645719659508108418738243030372422533994743756508347851178597805285404", 16)
	d1, _ = new(big.Int).SetString("44284103725881469857708921760250942555493926042739521955667532690943419533349", 16)
	d2, _ = new(big.Int).SetString("30959114451459308172620308816717028819652930302076989559150120900216006593354", 16)
}

func TestRRbyP256(t *testing.T) {
	n256 := new(big.Int).SetUint64(1)
	n256.Lsh(n256, 256)
	cParams := elliptic.P256().Params()
	n := cParams.N
	R := new(big.Int).Mod(n256, n)
	RR := new(big.Int).Mul(R, R)
	RR.Mod(RR, n)
	ww := RR.Bits()
	t.Logf("RR is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	ww = n.Bits()
	t.Logf("N(order) is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])

	p := cParams.P
	r := new(big.Int).Mod(n256, p)
	rr := new(big.Int).Mul(r, r)
	rr.Mod(rr, p)
	ww = rr.Bits()
	t.Logf("rr is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	Rinv := new(big.Int).SetUint64(1)
	Rinv.Lsh(Rinv, 257)
	Rinv.Mod(Rinv, p)
	Rinv.ModInverse(Rinv, p)
	t.Log("RInverse is ", Rinv.Text(16))
	K0 := new(big.Int).SetUint64(0xccd1c8aaee00bc4f)
	N0 := new(big.Int).Mul(K0, n)
	ww = N0.Bits()
	t.Logf("K0: %s, n*K0: %x %x %x %x", K0.Text(16), ww[0], ww[1], ww[2], ww[3])
	K0.SetUint64(1)
	K0.Lsh(K0, 64)
	N0 = N0.ModInverse(n, K0)
	if N0 == nil {
		t.Log("Ca'nt calc N0")
	} else {
		if N0.Cmp(K0) >= 0 {
			t.Log("SHOULD NEVER OCCUR")
			N0 = N0.Mod(K0, N0)
		} else {
			N0 = N0.Sub(K0, N0)
		}
		t.Logf("N0: %x", N0.Bits()[0])
	}
	ww = p.Bits()
	t.Logf("P: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	ww = cParams.B.Bits()
	t.Logf("B: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	ww = cParams.Gx.Bits()
	t.Logf("Gx: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	ww = cParams.Gy.Bits()
	t.Logf("Gy: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
}

func TestRRbySM2(t *testing.T) {
	n256 := new(big.Int).SetUint64(1)
	n256.Lsh(n256, 256)
	n512 := new(big.Int).Mul(n256, n256)
	cParams := P256().Params()
	n := cParams.N
	R := new(big.Int).Mod(n256, n)
	RR := new(big.Int).Mul(R, R)
	RR.Mod(RR, n)
	ww := RR.Bits()
	t.Logf("RR of sm2 is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	ww = n.Bits()
	t.Logf("N(order) is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])

	p := cParams.P
	mu := new(big.Int).Div(n512, p)
	ww = mu.Bits()
	t.Logf("mu: %x %x %x %x %x", ww[0], ww[1], ww[2], ww[3], ww[4])
	r := new(big.Int).Mod(n256, p)
	rr := new(big.Int).Mul(r, r)
	rr.Mod(rr, p)
	ww = rr.Bits()
	t.Logf("rr is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	Rinv := new(big.Int).SetUint64(1)
	Rinv.Lsh(Rinv, 257)
	Rinv.Mod(Rinv, p)
	Rinv.ModInverse(Rinv, p)
	t.Log("RInverse of sm2 is ", Rinv.Text(16))
	K0 := new(big.Int).SetUint64(0x327f9e8872350975)
	N0 := new(big.Int).Mul(K0, n)
	ww = N0.Bits()
	t.Logf("K0: %s, n*K0: %x %x %x %x", K0.Text(16), ww[0], ww[1], ww[2], ww[3])
	K0.SetUint64(1)
	K0.Lsh(K0, 64)
	N0 = N0.ModInverse(n, K0)
	if N0 == nil {
		t.Log("Ca'nt calc N0")
	} else {
		if N0.Cmp(K0) >= 0 {
			t.Log("SHOULD NEVER OCCUR")
			N0 = N0.Mod(K0, N0)
		} else {
			N0 = N0.Sub(K0, N0)
		}
		t.Logf("N0: %x", N0.Bits()[0])
	}
	ww = p.Bits()
	t.Logf("P: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	ww = cParams.B.Bits()
	t.Logf("B: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	ww = cParams.Gx.Bits()
	t.Logf("Gx: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	ww = cParams.Gy.Bits()
	t.Logf("Gy: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
}

func TestSM2AsmGo(t *testing.T) {
	goCurve := SM2go()
	asmCurve := P256()
	goGx := goCurve.Params().Gx
	goGy := goCurve.Params().Gy
	gx, gy := goCurve.ScalarMult(goGx, goGy, d1.Bytes())
	aGx := asmCurve.Params().Gx
	aGy := asmCurve.Params().Gy
	ax, ay := asmCurve.ScalarMult(aGx, aGy, d1.Bytes())
	if ax.Cmp(gx) != 0 {
		t.Log("mult X diff")
		t.Fail()
	}
	if ay.Cmp(gy) != 0 {
		t.Log("mult Y diff")
		t.Fail()
	}
}

func BenchmarkInverse(b *testing.B) {
	b.ResetTimer()
	priv, _ := fastecdsa.GenerateKey(P256(), rand.Reader)
	p := P256().Params().P
	res := new(big.Int)

	b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_ = res.ModInverse(priv.PublicKey.X, p)
		}
	})
}

func BenchmarkMul(b *testing.B) {
	b.ResetTimer()
	//p := P256().Params().P
	res := new(big.Int)

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = res.Mul(x1, x2)
		//_ = res.Mod(res, p)
	}
}

func BenchmarkECADD(b *testing.B) {
	b.ResetTimer()
	curve := P256()

	b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.Add(x1, y1, x2, y2)
		}
	})
}

func BenchmarkECDBL(b *testing.B) {
	b.ResetTimer()
	curve := P256()

	b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.Double(x1, y1)
		}
	})
}
