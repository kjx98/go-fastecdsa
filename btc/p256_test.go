package btc

import (
	//"crypto/elliptic"
	"math/big"
	"testing"
)

var (
	x1, y1 *big.Int
	x2, y2 *big.Int
	d1, d2 *big.Int
	n256   *big.Int
)

func init() {
	x1, _ = new(big.Int).SetString("afe524a88091d6daf7f9188477f6086a3c3c6c2f18c0e55b68e8ded5dbec39ca", 16)
	y1, _ = new(big.Int).SetString("8d837895a7dc179ef176066831aad5c2af60d71184a4ca536f18b74046a55994", 16)
	x2, _ = new(big.Int).SetString("684446fbee5167ef8baaa8adfa83e4606c0a05bb2f2c9125ca9ad478f1770dad", 16)
	y2, _ = new(big.Int).SetString("c7337843bdb886bff9965b00c6d87aff04f8d6bfa6a6846c0f28e513642bf309", 16)
	d1, _ = new(big.Int).SetString("44960d13c3ae7889e7fdfc0c48f4ac1da4e68fd3a5be28ad3f53eddad6d9c892", 16)
	d2, _ = new(big.Int).SetString("b68c5c25852521c647d7d0eddd09494949602ebaa885202a5573bb6ec8c5d96f", 16)
	n256 = new(big.Int).Lsh(bigOne, 256)
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

func calcK0a(p *big.Int) uint64 {
	t := uint64(1)
	N := uint64(p.Bits()[0])
	for i := 1; i < 52; i++ {
		t = t * t * N // mod 2^52
	}
	return -t
}

func calcRRa(p *big.Int) *big.Int {
	n260 := new(big.Int).Lsh(bigOne, 260)
	t := new(big.Int).Sub(n260, p)
	for i := 260; i < 520; i++ {
		t.Add(t, t)
		if t.Cmp(p) >= 0 {
			t.Sub(t, p)
		}
	}
	return t
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
	if crr.Cmp(btcg.rr) != 0 {
		t.Logf("calcRR diff btcg.rr, %s", crr.Text(16))
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
	K0 = new(big.Int).SetUint64(0xd838091dd2253531)
	N0 = new(big.Int).Mul(K0, p)
	ww = N0.Bits()
	ck = calcK0(p)
	t.Logf("K0: %s (%x), p*K0: %x %x %x %x", K0.Text(16), ck, ww[0], ww[1], ww[2], ww[3])
	ww = p.Bits()
	t.Logf("P: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	ww = cParams.Gx.Bits()
	t.Logf("Gx: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
	ww = cParams.Gy.Bits()
	t.Logf("Gy: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
}

func TestBTCAsmGo(t *testing.T) {
	goCurve := BTCgo()
	//asm version BTC works
	asmCurve := BTC()
	//asmCurve := P256()
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

func BenchmarkModMul(b *testing.B) {
	p := P256().Params().P
	res := new(big.Int)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = res.Mul(x1, y1)
		_ = res.Mod(res, p)
	}
}

func BenchmarkECADD(b *testing.B) {
	curve := P256()

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.Add(x1, y1, x2, y2)
		}
	})
}

func BenchmarkECDBL(b *testing.B) {
	curve := P256()

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.Double(x1, y1)
		}
	})
}

func BenchmarkECMULT(b *testing.B) {
	Curve := BTCgo()
	goGx := Curve.Params().Gx
	goGy := Curve.Params().Gy

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = Curve.ScalarMult(goGx, goGy, d1.Bytes())
		}
	})
}
