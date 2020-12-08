package sm2

import (
	"crypto/elliptic"
	"crypto/rand"
	"gitee.com/jkuang/go-fastecdsa"
	"math/big"
	"testing"
)

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
	p256 := elliptic.P256()
	priv, _ := fastecdsa.GenerateKey(p256, rand.Reader)
	p := p256.Params().P
	res := new(big.Int)

	b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_ = res.Mul(priv.PublicKey.X, priv.PublicKey.Y)
			_ = res.Mod(priv.PublicKey.X, p)
		}
	})
}
