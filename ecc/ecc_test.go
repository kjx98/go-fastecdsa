package ecc

import (
	"gitee.com/jkuang/go-fastecdsa/sm2"
	"math/big"
	"testing"
)

var (
	x1, y1 *big.Int
	x2, y2 *big.Int
)

func init() {
	x1, _ = new(big.Int).SetString("12504736780261570232047283225790923734936736733904083456330261180763076903173", 16)
	y1, _ = new(big.Int).SetString("81307578426096630754000956949482484966368962694792140728234975018259774469569", 16)
	x2, _ = new(big.Int).SetString("48121564271922987841895377752074498583012355812029682461364979458450873405695", 16)
	y2, _ = new(big.Int).SetString("37446874645719659508108418738243030372422533994743756508347851178597805285404", 16)
}

func TestEccMMod(t *testing.T) {
	p := sm2.P256().Params().P
	xy := new(big.Int).Mul(x1, y1)
	xyMod := new(big.Int).Mod(xy, p)
	bMod := vliModMult(x1.Bytes(), y1.Bytes(), p.Bytes())
	if bMod.Cmp(xyMod) != 0 {
		t.Log("big.MulMod diff vliModMult")
		t.Fail()
	}
	pb := make([]big.Word, 9)
	copy(pb, p.Bits())
	mu := new(big.Int).SetUint64(1)
	mu.Lsh(mu, 512)
	mu.Div(mu, p)
	copy(pb[4:], mu.Bits())
	t.Logf("P0: %x, P3: %x, Mu0: %x, Mu3: %x", pb[0], pb[3], pb[4], pb[7])
	bMod = vliModMultBarrett(x1.Bytes(), y1.Bytes(), pb)
	if bMod.Cmp(xyMod) != 0 {
		t.Logf("big.mulmod diff barrett mulmod:\n%s vs\n%s\n",
			xyMod.Text(16), bMod.Text(16))
		t.Fail()
	}
}

func BenchmarkInverse(b *testing.B) {
	b.ResetTimer()
	p := sm2.P256().Params().P.Bytes()

	b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_ = vliModInv(x1.Bytes(), p)
		}
	})
}

func BenchmarkMul(b *testing.B) {
	b.ResetTimer()
	p := sm2.P256().Params().P.Bytes()

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = vliModMult(x1.Bytes(), x2.Bytes(), p)
	}
}

/*
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
*/
