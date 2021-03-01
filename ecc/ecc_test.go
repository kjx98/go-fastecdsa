package ecc

import (
	"gitee.com/jkuang/go-fastecdsa/sm2"
	"math/big"
	"testing"
)

var (
	x1, y1  *big.Int
	x2, y2  *big.Int
	d1, d2  *big.Int
	bTwo    *big.Int
	rr, rrN []uint64
	k0N     uint64
)

func init() {
	x1, _ = new(big.Int).SetString("12504736780261570232047283225790923734936736733904083456330261180763076903173", 10)
	y1, _ = new(big.Int).SetString("81307578426096630754000956949482484966368962694792140728234975018259774469569", 10)
	x2, _ = new(big.Int).SetString("48121564271922987841895377752074498583012355812029682461364979458450873405695", 10)
	y2, _ = new(big.Int).SetString("37446874645719659508108418738243030372422533994743756508347851178597805285404", 10)
	d1, _ = new(big.Int).SetString("44960d13c3ae7889e7fdfc0c48f4ac1da4e68fd3a5be28ad3f53eddad6d9c892", 16)
	d2, _ = new(big.Int).SetString("b68c5c25852521c647d7d0eddd09494949602ebaa885202a5573bb6ec8c5d96f", 16)
	rr = []uint64{0x200000003, 0x2ffffffff, 0x100000001, 0x400000002}
	rrN = []uint64{0x901192af7c114f20, 0x3464504ade6fa2fa,
		0x620fc84c3affe0d4, 0x1eb5e412a22b3d3b}
	k0N = 0x327f9e8872350975
	bTwo = big.NewInt(2)
}

func dumpBits(v *big.Int, t *testing.T) {
	bLen := v.BitLen()
	if bLen == 0 {
		t.Log("big.Int v is zero")
		return
	}
	cOne := true
	bitsC := 0
	for i := bLen - 1; i >= 0; i-- {
		if cOne {
			if v.Bit(i) == 0 {
				t.Logf("%d one ", bitsC)
				bitsC = 1
				cOne = false
			} else {
				bitsC++
			}
		} else {
			if v.Bit(i) != 0 {
				t.Logf("%d zero ", bitsC)
				bitsC = 1
				cOne = true
			} else {
				bitsC++
			}
		}
	}
	if bitsC > 0 {
		ss := "one"
		if !cOne {
			ss = "zero"
		}
		t.Logf("%d %s\n", bitsC, ss)
	}
}

func TestEccMMod(t *testing.T) {
	bFMA := vliTestFMA()
	t.Log("CPU support FMA: ", bFMA)
	t.Log("CPU support CMOV: ", vliTestCMOV())
	p := sm2.P256().Params().P
	//xy := new(big.Int).Mul(x1, y1)
	//xyMod := new(big.Int).Mod(xy, p)
	pb := make([]big.Word, 9)
	copy(pb, p.Bits())
	mu := new(big.Int).SetUint64(1)
	mu.Lsh(mu, 512)
	mu.Div(mu, p)
	copy(pb[4:], mu.Bits())
	t.Logf("P0: %x, P3: %x, Mu0: %x, Mu3: %x", pb[0], pb[3], pb[4], pb[7])
}

func TestSM2MultP(t *testing.T) {
	p := sm2.P256().Params().P
	polyP := vliSM2MultP(1)
	if polyP.Cmp(p) != 0 {
		ww := polyP.Bits()
		t.Logf("sm2 polyP diff P: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
		ww = p.Bits()
		t.Logf("sm2 P: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
		t.Fail()
	} else {
		t.Log("polynomial Prime OK")
	}
	mp := new(big.Int).Mul(p, big.NewInt(0x33445566))
	polyP = vliSM2MultP(0x33445566)
	if polyP.Cmp(mp) != 0 {
		ww := polyP.Bits()
		t.Logf("sm2 polyP diff P: %X %X %X %X", ww[0], ww[1], ww[2], ww[3])
		ww = mp.Bits()
		t.Logf("sm2 P: %X %X %X %X %X", ww[0], ww[1], ww[2], ww[3], ww[4])
		t.Fail()
	} else {
		t.Log("polynomial Prime OK")
		/*
			polyP.Sub(p, bTwo)
			t.Log("Dump P - 2:")
			dumpBits(polyP, t)
			polyP.Sub(sm2.P256().Params().N, bTwo)
			t.Log("Dump N - 2:")
			dumpBits(polyP, t)
		*/
	}
}

func TestMontMultMod(t *testing.T) {
	p := sm2.P256().Params().P
	xy := new(big.Int).Mul(x1, y1)
	xyMod := new(big.Int).Mod(xy, p)
	xp := vliToMont(x1.Bits(), p.Bits(), rr, 1)
	yp := vliToMont(y1.Bits(), p.Bits(), rr, 1)
	bMod1 := vliModMultMont(xp.Bits(), yp.Bits(), p.Bits(), rr, 1)
	bMod := vliFromMont(bMod1.Bits(), p.Bits(), rr, 1)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step1 big.mulmod diff ModMultMont:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	xy = new(big.Int).Mul(x2, y2)
	xyMod = new(big.Int).Mod(xy, p)
	xp = vliToMont(x2.Bits(), p.Bits(), rr, 1)
	yp = vliToMont(y2.Bits(), p.Bits(), rr, 1)
	bMod1 = vliModMultMont(xp.Bits(), yp.Bits(), p.Bits(), rr, 1)
	bMod = vliFromMont(bMod1.Bits(), p.Bits(), rr, 1)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step2 big.mulmod diff ModMultMont:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	n := sm2.P256().Params().N
	xy = new(big.Int).Mul(x1, y1)
	xyMod = new(big.Int).Mod(xy, n)
	xp = vliToMont(x1.Bits(), n.Bits(), rrN, k0N)
	yp = vliToMont(y1.Bits(), n.Bits(), rrN, k0N)
	bMod1 = vliModMultMont(xp.Bits(), yp.Bits(), n.Bits(), rrN, k0N)
	bMod = vliFromMont(bMod1.Bits(), n.Bits(), rrN, k0N)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step3 big.mulmod diff ModMultMont:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	xy = new(big.Int).Mul(x2, y2)
	xyMod = new(big.Int).Mod(xy, n)
	xp = vliToMont(x2.Bits(), n.Bits(), rrN, k0N)
	yp = vliToMont(y2.Bits(), n.Bits(), rrN, k0N)
	bMod1 = vliModMultMont(xp.Bits(), yp.Bits(), n.Bits(), rrN, k0N)
	bMod = vliFromMont(bMod1.Bits(), n.Bits(), rrN, k0N)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step4 big.mulmod diff ModMultMont:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
}

func TestMontSqrMod(t *testing.T) {
	p := sm2.P256().Params().P
	xy := new(big.Int).Mul(x1, x1)
	xyMod := new(big.Int).Mod(xy, p)
	xp := vliToMont(x1.Bits(), p.Bits(), rr, 1)
	bMod1 := vliModSqrMont(xp.Bits(), p.Bits(), rr, 1)
	bMod := vliFromMont(bMod1.Bits(), p.Bits(), rr, 1)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step1 big.mulsqr diff ModSqrMont:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	xy = new(big.Int).Mul(x2, x2)
	xyMod = new(big.Int).Mod(xy, p)
	xp = vliToMont(x2.Bits(), p.Bits(), rr, 1)
	bMod1 = vliModSqrMont(xp.Bits(), p.Bits(), rr, 1)
	bMod = vliFromMont(bMod1.Bits(), p.Bits(), rr, 1)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step2 big.mulsqr diff ModSqrMont:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	n := sm2.P256().Params().N
	xy = new(big.Int).Mul(x1, x1)
	xyMod = new(big.Int).Mod(xy, n)
	xp = vliToMont(x1.Bits(), n.Bits(), rrN, k0N)
	bMod1 = vliModSqrMont(xp.Bits(), n.Bits(), rrN, k0N)
	bMod = vliFromMont(bMod1.Bits(), n.Bits(), rrN, k0N)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step3 big.mulsqr diff ModSqrMont:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	xy = new(big.Int).Mul(x2, x2)
	xyMod = new(big.Int).Mod(xy, n)
	xp = vliToMont(x2.Bits(), n.Bits(), rrN, k0N)
	bMod1 = vliModSqrMont(xp.Bits(), n.Bits(), rrN, k0N)
	bMod = vliFromMont(bMod1.Bits(), n.Bits(), rrN, k0N)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step4 big.mulsqr diff ModSqrMont:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
}

func TestMontExpMod(t *testing.T) {
	p := sm2.P256().Params().P
	xyMod := new(big.Int).Exp(x1, y1, p)
	bMod := vliExpModMont(x1.Bits(), y1.Bits(), p.Bits(), rr, 1)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step1 big.mulmod diff ModMultMont:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	xyMod = new(big.Int).Exp(x2, y2, p)
	bMod = vliExpModMont(x2.Bits(), y2.Bits(), p.Bits(), rr, 1)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step2 big.mulmod diff ModMultMont:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
}

func TestSqrtMod(t *testing.T) {
	p := sm2.P256().Params().P
	yy := new(big.Int).Mul(y1, y1)
	yy.Mod(yy, p)
	ySqrt := new(big.Int).ModSqrt(yy, p)
	if ySqrt == nil {
		t.Log("can't find sqrt of yy")
		t.Fail()
	} else {
		ww := yy.Bits()
		t.Logf("sqr y1: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	}
	t0 := new(big.Int).Rsh(p, 2)
	a1 := new(big.Int).Exp(yy, t0, p)
	a0 := new(big.Int).Mul(a1, yy)
	a0.Mod(a0, p)
	tt := new(big.Int).Mul(a0, a1)
	tt.Mod(tt, p)
	bOne := big.NewInt(1)
	if tt.Cmp(bOne) != 0 {
		t.Log("no sqrt for yy via yy^(p/4)")
		t.Fail()
	} else {
		ww := t0.Bits()
		t.Logf("quadP : %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
		dumpBits(t0, t)
	}
	if a0.Cmp(y1) == 0 {
		t.Log("sqrt = a * a ^ (p/4) ")
		ww := a1.Bits()
		t.Logf("y1 ^ quadP: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	}
	yy = yy.Exp(y2, big.NewInt(2), p)
	ySqrt = ySqrt.ModSqrt(yy, p)
	if ySqrt == nil {
		t.Log("can't find sqrt of yy")
		t.Fail()
	} else {
		ww := yy.Bits()
		t.Logf("sqr y2: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	}
	a1 = a1.Exp(yy, t0, p)
	a0 = a0.Mul(a1, yy)
	a0.Mod(a0, p)
	tt = tt.Mul(a0, a1)
	tt.Mod(tt, p)
	if tt.Cmp(bOne) != 0 {
		t.Log("no sqrt for yy via yy^(p/4)")
		t.Fail()
	}
	if a0.Cmp(y2) != 0 {
		// y2 should be even
		if y2.Bit(0) != 0 {
			t.Log("error, y2 not even")
			t.Fail()
		}
		a0.Sub(p, a0)
	}
	if a0.Cmp(y2) == 0 {
		t.Log("sqrt = a * a ^ (p/4) ")
		ww := a1.Bits()
		t.Logf("y2 ^ quadP: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	}
}

func TestSM2MontMultMod(t *testing.T) {
	p := sm2.P256().Params().P
	xy := new(big.Int).Mul(x1, y1)
	xyMod := new(big.Int).Mod(xy, p)
	bMod := vliSM2ModMultP(x1.Bits(), y1.Bits())
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step1 big.mulmod diff SM2ModMultP:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	xy = new(big.Int).Mul(x2, y2)
	xyMod = new(big.Int).Mod(xy, p)
	bMod = vliSM2ModMultP(x2.Bits(), y2.Bits())
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step2 big.mulmod diff SM2ModMultP:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	n := sm2.P256().Params().N
	xy = new(big.Int).Mul(x1, y1)
	xyMod = new(big.Int).Mod(xy, n)
	bMod = vliSM2ModMultN(x1.Bits(), y1.Bits())
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step3 big.mulmod diff SM2ModMultN:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	xy = new(big.Int).Mul(x2, y2)
	xyMod = new(big.Int).Mod(xy, n)
	bMod = vliSM2ModMultN(x2.Bits(), y2.Bits())
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step4 big.mulmod diff SM2ModMultN:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
}

func TestEccInverse(t *testing.T) {
	p := sm2.P256().Params().P
	two := big.NewInt(2)
	nMinus2 := new(big.Int).Sub(p, two)
	eInv := new(big.Int).Exp(x1, nMinus2, p)
	mInv := vliExpModMont(x1.Bits(), nMinus2.Bits(), p.Bits(), rr, 1)
	inv := new(big.Int).ModInverse(x1, p)
	if inv.Cmp(eInv) != 0 {
		t.Logf("vliModInv diff: ModInverse vs fermatInverse\n%s\n%s",
			inv.Text(16), eInv.Text(16))
		t.Fail()
	}
	if inv.Cmp(mInv) != 0 {
		t.Logf("vliMontInv diff: ModInverse vs fermatInverse\n%s\n%s",
			inv.Text(16), mInv.Text(16))
		t.Fail()
	}
	ww := vliModInv(x1.Bits(), p.Bits())
	cInv := new(big.Int).SetBits(ww)
	if cInv.Cmp(inv) != 0 {
		t.Logf("vliModInv diff: ModInverse vs vliModInv\n%s\n%s",
			inv.Text(16), cInv.Text(16))
		t.Fail()
	}
	eInv = new(big.Int).Exp(x2, nMinus2, p)
	inv = new(big.Int).ModInverse(x2, p)
	if inv.Cmp(eInv) != 0 {
		t.Logf("vliModInv diff: ModInverse vs fermatInverse\n%s\n%s",
			inv.Text(16), eInv.Text(16))
		t.Fail()
	}
	ww = vliModInv(x2.Bits(), p.Bits())
	cInv = new(big.Int).SetBits(ww)
	if cInv.Cmp(inv) != 0 {
		t.Logf("vliModInv diff: ModInverse vs vliModInv\n%s\n%s",
			inv.Text(16), cInv.Text(16))
		t.Fail()
	}
}

func BenchmarkInverse(b *testing.B) {
	p := sm2.P256().Params().P.Bits()

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_ = vliModInv(x1.Bits(), p)
		}
	})
}

func BenchmarkInverseGo(b *testing.B) {
	p := sm2.P256().Params().P

	inv := new(big.Int)
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			inv.ModInverse(x1, p)
		}
	})
}

func BenchmarkInverseFermat(b *testing.B) {
	p := sm2.P256().Params().P
	two := big.NewInt(2)
	nMinus2 := new(big.Int).Sub(p, two)

	inv := new(big.Int)
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			inv.Exp(x1, nMinus2, p)
			//inv.ModInverse(x1, p)
		}
	})
}

func BenchmarkMontModMul(b *testing.B) {
	p := sm2.P256().Params().P

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_ = vliModMultMont(x1.Bits(), y1.Bits(), p.Bits(), rr, 1)
		}
	})
}

func BenchmarkMontModSqr(b *testing.B) {
	p := sm2.P256().Params().P

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_ = vliModSqrMont(x1.Bits(), p.Bits(), rr, 1)
		}
	})
}

func BenchmarkSM2ModMulP(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = vliSM2ModMultP(x1.Bits(), y1.Bits())
	}
}

func BenchmarkMontExpMod(b *testing.B) {
	p := sm2.P256().Params().P

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = vliExpModMont(x1.Bits(), y1.Bits(), p.Bits(), rr, 1)
	}
}

func BenchmarkExpMod(b *testing.B) {
	p := sm2.P256().Params().P
	xymod := new(big.Int)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		xymod.Exp(x1, y1, p)
	}
}

func BenchmarkSM2MultP(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = vliSM2MultP(0x55557777)
	}
}
