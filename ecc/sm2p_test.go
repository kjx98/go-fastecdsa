// +build sm2p

package ecc

import (
	"github.com/kjx98/go-fastecdsa/sm2"
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
	x3, y3 = c.Add(x1, y1, x3, y3)
	x3a, y3a = sm2Add(x1, y1, x3a, y3a)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.Add step1b diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.AddJ step1b diff sm2.Double y:\n%s\n%s",
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
	x3, y3, z3 = jc.AddJacobian(x1, y1, bigOne, x3, y3, z3)
	x3a, y3a, z3a = sm2AddJacobian(x1, y1, bigOne, x3a, y3a, z3a)
	if x3.Cmp(x3a) != 0 {
		t.Logf("sm2c.AddJ step1b diff sm2.Double x:\n%s\n%s",
			x3.Text(16), x3a.Text(16))
		t.Fail()
	}
	if y3.Cmp(y3a) != 0 {
		t.Logf("sm2c.AddJ step1b diff sm2.Double y:\n%s\n%s",
			y3.Text(16), y3a.Text(16))
		t.Fail()
	}
	if z3.Cmp(z3a) != 0 {
		t.Logf("sm2c.AddJ step1b diff sm2.Double z:\n%s\n%s",
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
	x3, y3, z3 = jc.DoubleJacobian(x3, y3, z3)
	x3a, y3a, z3a = sm2DoubleJacobian(x3a, y3a, z3a)
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
	x3, y3, z3 = jc.DoubleJacobian(x3, y3, z3)
	x3a, y3a, z3a = sm2DoubleJacobian(x3a, y3a, z3a)
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

func TestSM2CurveBaseMult(t *testing.T) {
	goCurve := sm2.SM2go()
	gx, gy := goCurve.ScalarBaseMult(d1.Bytes())
	ax, ay := sm2ScalarBaseMult(d1.Bytes())
	if ax.Cmp(gx) != 0 {
		t.Logf("BaseMult step1 X diff:\n%s\n%s", gx.Text(16), ax.Text(16))
		t.Fail()
	}
	if ay.Cmp(gy) != 0 {
		t.Logf("BaseMult step1 Y diff:\n%s\n%s", gy.Text(16), ay.Text(16))
		t.Fail()
	}
	ww := gx.Bits()
	t.Logf("d1G x: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	ww = gy.Bits()
	t.Logf("d1G y: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	gx, gy = goCurve.ScalarBaseMult(d2.Bytes())
	ax, ay = sm2ScalarBaseMult(d2.Bytes())
	if ax.Cmp(gx) != 0 {
		t.Logf("BaseMult step2 X diff:\n%s\n%s", gx.Text(16), ax.Text(16))
		t.Fail()
	}
	if ay.Cmp(gy) != 0 {
		t.Logf("BaseMult step2 Y diff:\n%s\n%s", gy.Text(16), ay.Text(16))
		t.Fail()
	}
	ww = gx.Bits()
	t.Logf("d1G x: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	ww = gy.Bits()
	t.Logf("d1G y: %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
}

func TestSM2CurveMult(t *testing.T) {
	goCurve := sm2.SM2go()
	gx, gy := goCurve.ScalarMult(x1, y1, d1.Bytes())
	ax, ay := sm2ScalarMult(x1, y1, d1.Bytes())
	if ax.Cmp(gx) != 0 {
		t.Logf("BaseMult step1 X diff:\n%s\n%s", gx.Text(16), ax.Text(16))
		t.Fail()
	}
	if ay.Cmp(gy) != 0 {
		t.Logf("BaseMult step1 Y diff:\n%s\n%s", gy.Text(16), ay.Text(16))
		t.Fail()
	}
	gx, gy = goCurve.ScalarMult(x2, y2, d2.Bytes())
	ax, ay = sm2ScalarMult(x2, y2, d2.Bytes())
	if ax.Cmp(gx) != 0 {
		t.Logf("BaseMult step2 X diff:\n%s\n%s", gx.Text(16), ax.Text(16))
		t.Fail()
	}
	if ay.Cmp(gy) != 0 {
		t.Logf("BaseMult step2 Y diff:\n%s\n%s", gy.Text(16), ay.Text(16))
		t.Fail()
	}
}

/*
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

func BenchmarkSM2ECMULT(b *testing.B) {
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = sm2ScalarMult(x1, y1, d1.Bytes())
		}
	})
}

func BenchmarkSM2ECGMULT(b *testing.B) {
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _ = sm2ScalarBaseMult(d1.Bytes())
		}
	})
}
