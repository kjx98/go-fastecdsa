package ecc

import (
	"gitee.com/jkuang/go-fastecdsa/sm2"
	"math/big"
	"testing"
)

type jacobianIntf interface {
	// affineFromJacobian reverses the Jacobian transform
	AffineFromJacobian(x, y, z *big.Int) (xOut, yOut *big.Int)
	// addJacobian takes two points in Jacobian coordinates
	AddJacobian(x1, y1, z1, x2, y2, z2 *big.Int) (*big.Int, *big.Int, *big.Int)
	// doubleJacobian takes a point in Jacobian coordinates
	DoubleJacobian(x, y, z *big.Int) (*big.Int, *big.Int, *big.Int)
}

func TestEccCurveParam(t *testing.T) {
	c := getCurveParams(0)
	if c == nil {
		t.Log("Can't get ECC Curve Params")
		t.Fail()
	}
	t.Logf("ECC Curve Params: %s\nP: %s\nN: %s\nB: %s\nGx: %s\nGy: %s\n",
		c.Name, c.P.Text(16), c.N.Text(16), c.B.Text(16),
		c.Gx.Text(16), c.Gy.Text(16))
	ww := c.B.Bits()
	t.Logf("Ecc B: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
}

func TestMontMultModP(t *testing.T) {
	//p := sm2.P256().Params().P
	p := sm2.SM2().Params().P
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

func TestCurveAdd(t *testing.T) {
	//c := sm2.SM2go()
	c := sm2.SM2()
	x3, y3 := c.Add(x1, y1, x2, y2)
	x3a, y3a := sm2c.Add(x1, y1, x2, y2)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.Add diff sm2.Add x:\n%s\n%s", x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.Add diff sm2.Add y:\n%s\n%s", y3.Text(16), y3a.Text(16))
		t.Fail()
	}
}

func TestCurveAddJacobian(t *testing.T) {
	c := sm2.SM2go()
	var jc jacobianIntf
	if opt, ok := c.(jacobianIntf); ok {
		jc = opt
	} else {
		t.Log("not jacobian interface")
		t.Fail()
	}
	x3, y3, z3 := jc.AddJacobian(x1, y1, bigOne, x2, y2, bigOne)
	x3a, y3a, z3a := sm2c.addJacobian(x1, y1, bigOne, x2, y2, bigOne)
	xx3, yy3 := jc.AffineFromJacobian(x3, y3, z3)
	xx3a, yy3a := sm2c.affineFromJacobian(x3a, y3a, z3a)
	if xx3.Cmp(xx3a) != 0 {
		t.Logf("sm2c.AddJ step1 diff sm2.Add x:\n%s\n%s", xx3.Text(16), xx3a.Text(16))
		t.Fail()
	}
	if yy3.Cmp(yy3a) != 0 {
		t.Logf("sm2c.AddJ step1 diff sm2.Add y:\n%s\n%s", yy3.Text(16), yy3a.Text(16))
		t.Fail()
	}
	if z3.Cmp(z3a) != 0 {
		t.Logf("sm2c.AddJ step1 diff sm2.Add z:\n%s\n%s", z3.Text(16), z3a.Text(16))
		t.Logf("sm2c.AddJ sm2c x3a: %s\ny3a: %s", x3a.Text(16), y3a.Text(16))
		ww := x3a.Bits()
		t.Logf("x3a: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
		ww = y3a.Bits()
		t.Logf("y3a: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
		ww = z3a.Bits()
		t.Logf("z3a: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	} else {
		t.Log("sm2c.AddJ step1 z is same")
	}
	x4, y4, z4 := jc.AddJacobian(x1, y1, bigOne, x3, y3, z3)
	x4a, y4a, z4a := sm2c.addJacobian(x1, y1, bigOne, x3a, y3a, z3a)
	xx3, yy3 = jc.AffineFromJacobian(x4, y4, z4)
	xx3a, yy3a = sm2c.affineFromJacobian(x4a, y4a, z4a)
	if xx3.Cmp(xx3a) != 0 {
		t.Logf("sm2c.AddJ diff step1b sm2.Add x:\n%s\n%s", xx3.Text(16), xx3a.Text(16))
		t.Fail()
	}
	if yy3.Cmp(yy3a) != 0 {
		t.Logf("sm2c.AddJ diff step1b sm2.Add y:\n%s\n%s", yy3.Text(16), yy3a.Text(16))
		t.Fail()
	}
	if z4.Cmp(z4a) != 0 {
		t.Logf("sm2c.AddJ step1b diff sm2.Add z:\n%s\n%s", z4.Text(16), z4a.Text(16))
	} else {
		t.Log("sm2c.AddJ step1b z is same")
	}
	x5, y5, z5 := jc.AddJacobian(x4, y4, z4, x3, y3, z3)
	x5a, y5a, z5a := sm2c.addJacobian(x4a, y4a, z4a, x3a, y3a, z3a)
	xx3, yy3 = jc.AffineFromJacobian(x5, y5, z5)
	xx3a, yy3a = sm2c.affineFromJacobian(x5a, y5a, z5a)
	if xx3.Cmp(xx3a) != 0 {
		t.Logf("sm2c.AddJ diff step2 sm2.Add x:\n%s\n%s", xx3.Text(16), xx3a.Text(16))
		t.Fail()
	}
	if yy3.Cmp(yy3a) != 0 {
		t.Logf("sm2c.AddJ diff step2 sm2.Add y:\n%s\n%s", yy3.Text(16), yy3a.Text(16))
		t.Fail()
	}
	if z5.Cmp(z5a) != 0 {
		t.Logf("sm2c.AddJ step2 diff sm2.Add z:\n%s\n%s", z4.Text(16), z4a.Text(16))
	} else {
		t.Log("sm2c.AddJ step2 z is same")
	}
}

func TestCurveDouble(t *testing.T) {
	//c := sm2.SM2go()
	c := sm2.SM2()
	x3, y3 := c.Double(x1, y1)
	x3a, y3a := sm2c.Double(x1, y1)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.Double step1 diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.Double step1 diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
	x3, y3 = c.Double(x3, y3)
	x3a, y3a = sm2c.Double(x3a, y3a)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.Double step1b diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.Double step1b diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
	x3, y3 = c.Double(x2, y2)
	x3a, y3a = sm2c.Double(x2, y2)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.Double step2 diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.Double step2 diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
	x3, y3 = c.Double(x3, y3)
	x3a, y3a = sm2c.Double(x3a, y3a)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.Double step2b diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.Double step2b diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
}

func TestCurveDoubleJacobian(t *testing.T) {
	c := sm2.SM2go()
	var jc jacobianIntf
	if opt, ok := c.(jacobianIntf); ok {
		jc = opt
	} else {
		t.Log("not jacobian interface")
		t.Fail()
	}
	x3, y3, z3 := jc.DoubleJacobian(x1, y1, bigOne)
	x3a, y3a, z3a := sm2c.doubleJacobian(x1, y1, bigOne)
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
	} else {
		t.Logf("sm2c.DoubleJ step1 z: %s", z3.Text(16))
	}
	ww := x3a.Bits()
	t.Logf("x3a: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	ww = y3a.Bits()
	t.Logf("y3a: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	ww = z3a.Bits()
	t.Logf("z3a: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	x3, y3, z3 = jc.DoubleJacobian(x3, y3, z3)
	x3a, y3a, z3a = sm2c.doubleJacobian(x3a, y3a, z3a)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.DoubleJ step1b diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.DoubleJ step1b diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
	if z3.Cmp(z3a) != 0 {
		t.Logf("sm2c.DoubleJ step1b diff sm2.Double z:\n%s\n%s",
			z3.Text(16), z3a.Text(16))
		t.Fail()
	}
	x3, y3, z3 = jc.DoubleJacobian(x2, y2, bigOne)
	x3a, y3a, z3a = sm2c.doubleJacobian(x2, y2, bigOne)
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
	x3, y3, z3 = jc.DoubleJacobian(x3, y3, z3)
	x3a, y3a, z3a = sm2c.doubleJacobian(x3a, y3a, z3a)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.DoubleJ step2b diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.DoubleJ step2b diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
	if z3.Cmp(z3a) != 0 {
		t.Logf("sm2c.DoubleJ step2b diff sm2.Double z:\n%s\n%s",
			z3.Text(16), z3a.Text(16))
		t.Fail()
	}
}

func TestCurveMult(t *testing.T) {
	//goCurve := sm2.SM2go()
	goCurve := sm2.SM2()
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

func TestCurveBaseMult(t *testing.T) {
	//goCurve := sm2.SM2go()
	goCurve := sm2.SM2()
	cCurve := sm2c
	gx, gy := goCurve.ScalarBaseMult(d1.Bytes())
	ax, ay := cCurve.ScalarBaseMult(d1.Bytes())
	if ax.Cmp(gx) != 0 {
		t.Log("d1 basemult X diff")
		t.Fail()
	}
	if ay.Cmp(gy) != 0 {
		t.Log("d1 basemult Y diff")
		t.Fail()
	}
	gx, gy = goCurve.ScalarBaseMult(d2.Bytes())
	ax, ay = cCurve.ScalarBaseMult(d2.Bytes())
	if ax.Cmp(gx) != 0 {
		t.Log("d2 basemult X diff")
		t.Fail()
	}
	if ay.Cmp(gy) != 0 {
		t.Log("d2 basemult Y diff")
		t.Fail()
	}
}

func BenchmarkMontMultModP(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_ = vliModMultMontP(x1.Bits(), y1.Bits())
	}
}

func BenchmarkECADD(b *testing.B) {
	curve := sm2c

	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.Add(x1, y1, x2, y2)
		}
	})
}

func BenchmarkECADDJac(b *testing.B) {
	curve := sm2c

	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _, _ = curve.addJacobian(x1, y1, bigOne, x2, y2, bigOne)
		}
	})
}

func BenchmarkECDBL(b *testing.B) {
	curve := sm2c

	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.Double(x1, y1)
		}
	})
}

func BenchmarkECDBLJac(b *testing.B) {
	curve := sm2c

	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _, _ = curve.doubleJacobian(x1, y1, bigOne)
		}
	})
}

func BenchmarkECMULT(b *testing.B) {
	curve := sm2c
	aGx := curve.Params().Gx
	aGy := curve.Params().Gy

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.ScalarMult(aGx, aGy, d1.Bytes())
		}
	})
}

func BenchmarkECCombinedMULT(b *testing.B) {
	curve := sm2c
	aGx := curve.Params().Gx
	aGy := curve.Params().Gy
	sBuff := make([]byte, 2048)

	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.cMult(aGx, aGy, d1.Bytes(), d2.Bytes(), sBuff)
		}
	})
}

func BenchmarkECGMULT(b *testing.B) {
	curve := sm2c

	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = curve.ScalarBaseMult(d1.Bytes())
		}
	})
}
