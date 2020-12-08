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
)

func init() {
	x1, _ = new(big.Int).SetString("12504736780261570232047283225790923734936736733904083456330261180763076903173", 16)
	y1, _ = new(big.Int).SetString("81307578426096630754000956949482484966368962694792140728234975018259774469569", 16)
	x2, _ = new(big.Int).SetString("48121564271922987841895377752074498583012355812029682461364979458450873405695", 16)
	y2, _ = new(big.Int).SetString("37446874645719659508108418738243030372422533994743756508347851178597805285404", 16)
}

func TestRRbyP256(t *testing.T) {
	n256 := new(big.Int).SetUint64(1)
	n256.Lsh(n256, 256)
	n := elliptic.P256().Params().N
	R := new(big.Int).Mod(n256, n)
	RR := new(big.Int).Mul(R, R)
	RR.Mod(RR, n)
	ww := RR.Bits()
	t.Logf("RR is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	ww = n.Bits()
	t.Logf("N(order) is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])

	p := elliptic.P256().Params().P
	Rinv := new(big.Int).SetUint64(1)
	Rinv.Lsh(Rinv, 257)
	Rinv.Mod(Rinv, p)
	Rinv.ModInverse(Rinv, p)
	t.Log("RInverse is ", Rinv.Text(16))
	t.Logf("x1: %s, y1: %s", x1.Text(16), y1.Text(16))
}

func TestRRbySM2(t *testing.T) {
	n256 := new(big.Int).SetUint64(1)
	n256.Lsh(n256, 256)
	n := P256().Params().N
	R := new(big.Int).Mod(n256, n)
	RR := new(big.Int).Mul(R, R)
	RR.Mod(RR, n)
	ww := RR.Bits()
	t.Logf("RR of sm2 is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])
	ww = n.Bits()
	t.Logf("N(order) is %x %x %x %x", ww[0], ww[1], ww[2], ww[3])

	p := P256().Params().P
	Rinv := new(big.Int).SetUint64(1)
	Rinv.Lsh(Rinv, 257)
	Rinv.Mod(Rinv, p)
	Rinv.ModInverse(Rinv, p)
	t.Log("RInverse of sm2 is ", Rinv.Text(16))
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
