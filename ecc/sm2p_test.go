// +build sm2p

package ecc

import (
	"gitee.com/jkuang/go-fastecdsa/sm2"
	"math/big"
	"testing"
)

/*
func TestMontMultModP(t *testing.T) {
	p := sm2.P256().Params().P
	xy := new(big.Int).Mul(x1, y1)
	xyMod := new(big.Int).Mod(xy, p)
	bMod := vliModMultMontP(x1.Bits(), y1.Bits())
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step1 big.mulmod diff ModMultMontP:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
	xy = new(big.Int).Mul(x2, y2)
	xyMod = new(big.Int).Mod(xy, p)
	bMod = vliModMultMontP(x2.Bits(), y2.Bits())
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("step2 big.mulmod diff ModMultMontP:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
}
*/

func TestSM2Inverse(t *testing.T) {
	p := sm2.P256().Params().P
	inv := new(big.Int).ModInverse(x1, p)
	cInv := sm2ModInv(x1.Bits())
	if cInv.Cmp(inv) != 0 {
		t.Logf("vliModInv diff: ModInverse vs vliModInv\n%s\n%s",
			inv.Text(16), cInv.Text(16))
		t.Fail()
	}
	inv = new(big.Int).ModInverse(x2, p)
	cInv = sm2ModInv(x2.Bits())
	if cInv.Cmp(inv) != 0 {
		t.Logf("vliModInv diff: ModInverse vs vliModInv\n%s\n%s",
			inv.Text(16), cInv.Text(16))
		t.Fail()
	}
}

func TestSM2Add(t *testing.T) {
	c := sm2.SM2go()
	x3, y3 := c.Add(x1, y1, x2, y2)
	x3a, y3a := sm2Add(x1, y1, x2, y2)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.Add step1 diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.AddJ step1 diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
}

func TestSM2AddJacobian(t *testing.T) {
	c := sm2.SM2go()
	var jc jacobianIntf
	if opt, ok := c.(jacobianIntf); ok {
		jc = opt
	} else {
		t.Log("not jacobian interface")
		t.Fail()
	}
	x3, y3, z3 := jc.AddJacobian(x1, y1, bigOne, x2, y2, bigOne)
	x3a, y3a, z3a := sm2AddJacobian(x1, y1, bigOne, x2, y2, bigOne)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.AddJ step1 diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.AddJ step1 diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
	if z3.Cmp(z3a) != 0 {
		t.Logf("sm2c.AddJ step1 diff sm2.Double z:\n%s\n%s",
			z3.Text(16), z3a.Text(16))
		t.Fail()
	}
}

func TestSM2Double(t *testing.T) {
	c := sm2.SM2go()
	x3, y3 := c.Double(x1, y1)
	x3a, y3a := sm2Double(x1, y1)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.DoubleJ step1 diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.DoubleJ step1 diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
	x3, y3 = c.Double(x2, y2)
	x3a, y3a = sm2Double(x2, y2)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.DoubleJ step2 diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.DoubleJ step2 diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
}

func TestSM2DoubleJacobian(t *testing.T) {
	c := sm2.SM2go()
	var jc jacobianIntf
	if opt, ok := c.(jacobianIntf); ok {
		jc = opt
	} else {
		t.Log("not jacobian interface")
		t.Fail()
	}
	x3, y3, z3 := jc.DoubleJacobian(x1, y1, bigOne)
	x3a, y3a, z3a := sm2DoubleJacobian(x1, y1, bigOne)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.DoubleJ step1 diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.DoubleJ step1 diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
	if z3.Cmp(z3a) != 0 {
		t.Logf("sm2c.DoubleJ step1 diff sm2.Double z:\n%s\n%s",
			z3.Text(16), z3a.Text(16))
		t.Fail()
	}
	x3, y3, z3 = jc.DoubleJacobian(x2, y2, bigOne)
	x3a, y3a, z3a = sm2DoubleJacobian(x2, y2, bigOne)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.DoubleJ step2 diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.DoubleJ step2 diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
	if z3.Cmp(z3a) != 0 {
		t.Logf("sm2c.DoubleJ step2 diff sm2.Double z:\n%s\n%s",
			z3.Text(16), z3a.Text(16))
		t.Fail()
	}
}

/*
func TestCurveMult(t *testing.T) {
	goCurve := sm2.SM2go()
	cCurve := sm2c
	goGx := goCurve.Params().Gx
	goGy := goCurve.Params().Gy
	gx, gy := goCurve.ScalarMult(goGx, goGy, d1.Bytes())
	aGx := cCurve.Params().Gx
	aGy := cCurve.Params().Gy
	ax, ay := cCurve.ScalarMult(aGx, aGy, d1.Bytes())
	if ax.Cmp(gx) != 0 {
		t.Log("mult X diff")
		t.Fail()
	}
	if ay.Cmp(gy) != 0 {
		t.Log("mult Y diff")
		t.Fail()
	}
}

func BenchmarkMontMultModP(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = vliModMultMontP(x1.Bits(), y1.Bits())
	}
}
*/

func BenchmarkSM2ModInv(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = sm2ModInv(x1.Bits())
	}
}

func BenchmarkSM2ECADD(b *testing.B) {
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = sm2Add(x1, y1, x2, y2)
		}
	})
}

func BenchmarkSM2ECADDJac(b *testing.B) {
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _, _ = sm2AddJacobian(x1, y1, bigOne, x2, y2, bigOne)
		}
	})
}

func BenchmarkSM2ECDBLJac(b *testing.B) {
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _, _ = sm2DoubleJacobian(x1, y1, bigOne)
		}
	})
}

/*
func BenchmarkECMULT(b *testing.B) {
	b.ResetTimer()
	curve := sm2c
	aGx := curve.Params().Gx
	aGy := curve.Params().Gy

	b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.ScalarMult(aGx, aGy, d1.Bytes())
		}
	})
}
*/
