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
	one    *big.Int
	n256   *big.Int
)

func init() {
	x1, _ = new(big.Int).SetString("afe524a88091d6daf7f9188477f6086a3c3c6c2f18c0e55b68e8ded5dbec39ca", 16)
	y1, _ = new(big.Int).SetString("8d837895a7dc179ef176066831aad5c2af60d71184a4ca536f18b74046a55994", 16)
	x2, _ = new(big.Int).SetString("684446fbee5167ef8baaa8adfa83e4606c0a05bb2f2c9125ca9ad478f1770dad", 16)
	y2, _ = new(big.Int).SetString("c7337843bdb886bff9965b00c6d87aff04f8d6bfa6a6846c0f28e513642bf309", 16)
	d1, _ = new(big.Int).SetString("44960d13c3ae7889e7fdfc0c48f4ac1da4e68fd3a5be28ad3f53eddad6d9c892", 16)
	d2, _ = new(big.Int).SetString("b68c5c25852521c647d7d0eddd09494949602ebaa885202a5573bb6ec8c5d96f", 16)
	one = new(big.Int).SetUint64(1)
	n256 = new(big.Int).Lsh(one, 256)
}

func calcK0(p *big.Int) uint64 {
	t := uint64(1)
	N := uint64(p.Bits()[0])
	for i := 1; i < 64; i++ {
		t = t * t * N // mod 2^64
	}
	return -t
}

func calcRR(p *big.Int) *big.Int {
	t := new(big.Int).Sub(n256, p)
	for i := 256; i < 512; i++ {
		t.Add(t, t)
		if t.Cmp(p) >= 0 {
			t.Sub(t, p)
		}
	}
	return t
}

func TestRRbyP256(t *testing.T) {
	cParams := elliptic.P256().Params()
	n := cParams.N
	R := new(big.Int).Mod(n256, n)
	RR := new(big.Int).Mul(R, R)
	RR.Mod(RR, n)
	ww := RR.Bits()
	t.Logf("NIST P256 RR is %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	cRR := calcRR(n)
	if cRR.Cmp(RR) != 0 {
		t.Logf("calcRR diff, %s", cRR.Text(16))
	} else {
		t.Log("calcRR NIST P256 n works")
	}
	ww = n.Bits()
	t.Logf("N(order) is %X %X %X %X", ww[0], ww[1], ww[2], ww[3])

	p := cParams.P
	// verify Prime poly
	n96 := new(big.Int).Lsh(one, 96)
	n224 := new(big.Int).Lsh(one, 224)
	n192 := new(big.Int).Lsh(one, 192)
	polyP := new(big.Int).Sub(n256, n224)
	polyP.Add(polyP, n192)
	polyP.Add(polyP, n96)
	polyP.Sub(polyP, one)
	if polyP.Cmp(p) != 0 {
		ww = polyP.Bits()
		t.Logf("P256 polyP diff P: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
		t.Fail()
	}
	r := new(big.Int).Mod(n256, p)
	rr := new(big.Int).Mul(r, r)
	rr.Mod(rr, p)
	ww = rr.Bits()
	t.Logf("rr is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	crr := calcRR(p)
	if crr.Cmp(rr) != 0 {
		t.Logf("calcRR diff, %s", crr.Text(16))
	} else {
		t.Log("calcRR NIST P256 p works")
	}
	Rinv := new(big.Int).SetUint64(1)
	Rinv.Lsh(Rinv, 257)
	Rinv.Mod(Rinv, p)
	Rinv.ModInverse(Rinv, p)
	t.Log("RInverse is ", Rinv.Text(16))
	K0 := new(big.Int).SetUint64(0xccd1c8aaee00bc4f)
	N0 := new(big.Int).Mul(K0, n)
	ww = N0.Bits()
	ck := calcK0(n)
	t.Logf("K0: %s (%x), n*K0: %x %x %x %x", K0.Text(16), ck, ww[0], ww[1], ww[2], ww[3])
	K0.SetUint64(1)
	K0.Lsh(K0, 64)
	N0 = N0.ModInverse(n, K0)
	if N0 == nil {
		t.Log("Can't calc N0")
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

func TestRRbyBTC(t *testing.T) {
	cParams := BTC().Params()
	n := cParams.N
	R := new(big.Int).Mod(n256, n)
	RR := new(big.Int).Mul(R, R)
	RR.Mod(RR, n)
	ww := RR.Bits()
	t.Logf("BTC RR is %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	cRR := calcRR(n)
	if cRR.Cmp(RR) != 0 {
		t.Logf("calcRR diff, %s", cRR.Text(16))
	} else {
		t.Log("calcRR secp256k1 n works")
	}
	ww = n.Bits()
	t.Logf("N(order) is %X %X %X %X", ww[0], ww[1], ww[2], ww[3])

	p := cParams.P
	// verify Prime poly
	btcPC := big.NewInt(0x1000003d1)
	polyP := new(big.Int).Sub(n256, btcPC)
	if polyP.Cmp(p) != 0 {
		ww = polyP.Bits()
		t.Logf("BTC polyP diff P: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
		t.Fail()
	}
	r := new(big.Int).Mod(n256, p)
	rr := new(big.Int).Mul(r, r)
	rr.Mod(rr, p)
	ww = rr.Bits()
	//t.Logf("rr is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	t.Logf("rr is %d words, %x %x", len(ww), ww[0], ww[1])
	crr := calcRR(p)
	if crr.Cmp(rr) != 0 {
		t.Logf("calcRR diff, %s", crr.Text(16))
	} else {
		t.Log("calcRR secp256k1 p works")
	}
	Rinv := new(big.Int).SetUint64(1)
	Rinv.Lsh(Rinv, 257)
	Rinv.Mod(Rinv, p)
	Rinv.ModInverse(Rinv, p)
	t.Log("RInverse is ", Rinv.Text(16))
	K0 := new(big.Int).SetUint64(0x4b0dff665588b13f)
	N0 := new(big.Int).Mul(K0, n)
	ww = N0.Bits()
	ck := calcK0(n)
	t.Logf("K0: %s (%x), n*K0: %x %x %x %x", K0.Text(16), ck, ww[0], ww[1], ww[2], ww[3])
	K0.SetUint64(1)
	K0.Lsh(K0, 64)
	N0 = N0.ModInverse(n, K0)
	if N0 == nil {
		t.Log("Can't calc N0")
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
	ww = cParams.Gx.Bits()
	t.Logf("Gx: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	ww = cParams.Gy.Bits()
	t.Logf("Gy: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
}

func TestRRbySM2(t *testing.T) {
	n96 := new(big.Int).Lsh(one, 96)
	n64 := new(big.Int).Lsh(one, 64)
	//n128 := new(big.Int).Lsh(one, 128)
	n224 := new(big.Int).Lsh(one, 224)
	smPP := new(big.Int).Sub(n256, n224)
	//smPP.Sub(smPP, n128)
	smPP.Sub(smPP, n96)
	smPP.Add(smPP, n64)
	smPP.Sub(smPP, one)
	n512 := new(big.Int).Mul(n256, n256)
	cParams := sm2g.Params()
	n := cParams.N
	R := new(big.Int).Mod(n256, n)
	RR := new(big.Int).Mul(R, R)
	RR.Mod(RR, n)
	ww := RR.Bits()
	t.Logf("RR mod N of sm2 is %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	cRR := calcRR(n)
	if cRR.Cmp(RR) != 0 {
		t.Logf("calcRR diff, %s", cRR.Text(16))
	} else {
		t.Log("calcRR sm2 n works")
	}
	ww = n.Bits()
	t.Logf("N(order) is %X %X %X %X", ww[0], ww[1], ww[2], ww[3])

	p := cParams.P
	mu := new(big.Int).Div(n512, p)
	ww = mu.Bits()
	t.Logf("mu: %X %X %X %X %X", ww[0], ww[1], ww[2], ww[3], ww[4])
	muN := new(big.Int).Div(n512, n)
	ww = muN.Bits()
	t.Logf("muN: %X %X %X %X %X", ww[0], ww[1], ww[2], ww[3], ww[4])
	r := new(big.Int).Mod(n256, p)
	rr := new(big.Int).Mul(r, r)
	rr.Mod(rr, p)
	ww = rr.Bits()
	t.Logf("rr mod P is %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	crr := calcRR(p)
	if crr.Cmp(rr) != 0 {
		t.Logf("calcRR diff, %s", crr.Text(16))
	} else {
		t.Log("calcRR sm2 p works")
	}
	Rinv := new(big.Int).SetUint64(1)
	Rinv.Lsh(Rinv, 257)
	Rinv.Mod(Rinv, p)
	Rinv.ModInverse(Rinv, p)
	t.Log("RInverse of sm2 is ", Rinv.Text(16))
	K0 := new(big.Int).SetUint64(0x327f9e8872350975)
	N0 := new(big.Int).Mul(K0, n)
	ww = N0.Bits()
	ck := calcK0(n)
	t.Logf("K0: %s (%x), n*K0: %X %X %X %X", K0.Text(16), ck, ww[0], ww[1], ww[2], ww[3])
	K0.SetUint64(1)
	K0.Lsh(K0, 64)
	N0 = N0.ModInverse(n, K0)
	if N0 == nil {
		t.Log("Can't calc N0")
	} else {
		if N0.Cmp(K0) >= 0 {
			t.Log("SHOULD NEVER OCCUR")
			N0 = N0.Mod(K0, N0)
		} else {
			N0 = N0.Sub(K0, N0)
		}
		t.Logf("N0: %x", N0.Bits()[0])
	}
	K0.Lsh(K0, 192)
	N0 = N0.ModInverse(p, K0)
	ck = calcK0(p)
	if N0 == nil {
		t.Log("Can't calc N0")
	} else {
		if N0.Cmp(K0) >= 0 {
			t.Log("SHOULD NEVER OCCUR")
			N0 = N0.Mod(K0, N0)
		} else {
			N0 = N0.Sub(K0, N0)
		}
		t.Logf("new N0: %x, mont K0: %x (%x)", N0.Bits()[0], sm2g.montK0(), ck)
	}
	polyP := sm2g.multP(1)
	if polyP.Cmp(p) != 0 {
		ww = polyP.Bits()
		t.Logf("multP polyP diff P: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
		t.Fail()
	} else {
		t.Log("polynomial Prime OK")
	}
	if smPP.Cmp(p) != 0 {
		ww = smPP.Bits()
		t.Logf("sm2 polyP diff P: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
		t.Fail()
	} else {
		t.Log("Lsh simu polynomial Prime OK")
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
	//asm version SM2 not work yet
	//asmCurve := SM2()
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

func TestMontMulMod(t *testing.T) {
	SM2go()
	c := sm2g
	prod := new(big.Int).Mul(x1, y1)
	m1 := new(big.Int).Mod(prod, c.P)
	m2 := c.montModMul(x1, y1)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontMulMod step 1 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
	prod = new(big.Int).Mul(x2, y2)
	m1 = new(big.Int).Mod(prod, c.P)
	m2 = c.montModMul(x2, y2)
	if m1.Cmp(m2) != 0 {
		t.Logf("MontMulMod step2 diff:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	}
}

func TestBarrettMod(t *testing.T) {
	SM2go()
	c := sm2g
	mu := c.GetMu()
	if mu == nil {
		t.Log("Can't get mu")
		t.Fail()
	} else {
		ww := mu.Bits()
		t.Logf("SM2 mu: %x %x %x %x %x", ww[0], ww[1], ww[2], ww[3], ww[4])
	}
	prod := new(big.Int).Mul(x1, y1)
	m1 := new(big.Int).Mod(prod, c.P)
	m2 := c.BarrettMod(prod)
	if m1.Cmp(m2) != 0 {
		t.Logf("step1 m1 diff m2:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	} else {
		t.Log("BarrettMod step1 ok")
	}
	prod = new(big.Int).Mul(x2, y2)
	m1 = new(big.Int).Mod(prod, c.P)
	m2 = c.BarrettMod(prod)
	if m1.Cmp(m2) != 0 {
		t.Logf("step2 m1 diff m2:\n%s vs\n%s", m1.Text(16), m2.Text(16))
		t.Fail()
	} else {
		t.Log("BarrettMod step2 ok")
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

func BenchmarkModMul(b *testing.B) {
	b.ResetTimer()
	p := P256().Params().P
	res := new(big.Int)

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = res.Mul(x1, y1)
		_ = res.Mod(res, p)
	}
}

func BenchmarkBarrettModMul(b *testing.B) {
	b.ResetTimer()
	c := sm2g
	res := new(big.Int)

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = res.Mul(x1, y1)
		_ = c.BarrettMod(res)
	}
}

func BenchmarkMontModMul(b *testing.B) {
	b.ResetTimer()
	c := sm2g

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = c.montModMul(x1, y1)
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

func BenchmarkECMULT(b *testing.B) {
	b.ResetTimer()
	Curve := SM2go()
	goGx := Curve.Params().Gx
	goGy := Curve.Params().Gy

	b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = Curve.ScalarMult(goGx, goGy, d1.Bytes())
		}
	})
}
