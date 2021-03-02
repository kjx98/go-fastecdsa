// Copyright 2011 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fastecdsa

import (
	"bufio"
	"compress/bzip2"
	"crypto/elliptic"
	"crypto/rand"
	"crypto/sha1"
	"crypto/sha256"
	"crypto/sha512"
	"encoding/hex"
	"github.com/kjx98/go-fastecdsa/btc"
	"github.com/kjx98/go-fastecdsa/ecc"
	"github.com/kjx98/go-fastecdsa/sm2"
	"hash"
	"io"
	"math/big"
	"os"
	"strings"
	"testing"
)

func testKeyGeneration(t *testing.T, c elliptic.Curve, tag string) {
	priv, err := GenerateKey(c, rand.Reader)
	if err != nil {
		t.Errorf("%s: error: %s", tag, err)
		return
	}
	if !c.IsOnCurve(priv.PublicKey.X, priv.PublicKey.Y) {
		t.Errorf("%s: public key invalid: %v", tag, err)
	}
	t.Log(tag, " KeyGen test PASS✅")
	/*
		t.Logf("%s priv D: %s", tag, priv.D.Text(16))
		t.Logf("%s curve X: %s, Y: %s", tag, priv.PublicKey.X.Text(16),
			priv.PublicKey.Y.Text(16))
	*/
}

func TestKeyGeneration(t *testing.T) {
	testKeyGeneration(t, elliptic.P224(), "p224")
	if testing.Short() {
		return
	}
	testKeyGeneration(t, elliptic.P256(), "p256")
	testKeyGeneration(t, elliptic.P384(), "p384")
	testKeyGeneration(t, sm2.P256(), "sm2")
	testKeyGeneration(t, sm2.SM2go(), "sm2go")
	testKeyGeneration(t, ecc.SM2C(), "sm2c")
	testKeyGeneration(t, sm2.SM2(), "sm2Asm")
	testKeyGeneration(t, btc.BTC(), "btcAsm")
}

func BenchmarkSignP256(b *testing.B) {
	//b.ResetTimer()
	p256 := elliptic.P256()
	hashed := []byte("testing")
	priv, _ := GenerateKey(p256, rand.Reader)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _, _ = Sign(rand.Reader, priv, hashed)
		}
	})
}

func BenchmarkSignSM2(b *testing.B) {
	//b.ResetTimer()
	//p256 := sm2.P256()
	p256 := sm2.SM2()
	hashed := []byte("testing")
	priv, _ := GenerateKey(p256, rand.Reader)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _, _ = Sign(rand.Reader, priv, hashed)
		}
	})
}

func BenchmarkSignSM2C(b *testing.B) {
	//b.ResetTimer()
	p256 := ecc.SM2C()
	hashed := []byte("testing")
	priv, _ := GenerateKey(p256, rand.Reader)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _, _ = Sign(rand.Reader, priv, hashed)
		}
	})
}

func BenchmarkSignBTC(b *testing.B) {
	//b.ResetTimer()
	p256 := btc.BTC()
	hashed := []byte("testing")
	priv, _ := GenerateKey(p256, rand.Reader)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			_, _, _ = Sign(rand.Reader, priv, hashed)
		}
	})
}

func BenchmarkVerifyP256(b *testing.B) {
	//b.ResetTimer()
	p256 := elliptic.P256()
	hashed := []byte("testing")
	priv, _ := GenerateKey(p256, rand.Reader)
	r, s, _ := Sign(rand.Reader, priv, hashed)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			Verify(&priv.PublicKey, hashed, r, s)
		}
	})
}

func BenchmarkVerifySM2(b *testing.B) {
	//b.ResetTimer()
	//p256 := sm2.P256()
	p256 := sm2.SM2()
	hashed := []byte("testing")
	priv, _ := GenerateKey(p256, rand.Reader)
	r, s, _ := Sign(rand.Reader, priv, hashed)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			Verify(&priv.PublicKey, hashed, r, s)
		}
	})
}

func BenchmarkVerifySM2C(b *testing.B) {
	//b.ResetTimer()
	p256 := ecc.SM2C()
	hashed := []byte("testing")
	priv, _ := GenerateKey(p256, rand.Reader)
	r, s, _ := Sign(rand.Reader, priv, hashed)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			Verify(&priv.PublicKey, hashed, r, s)
		}
	})
}

func BenchmarkVerifyBTC(b *testing.B) {
	//b.ResetTimer()
	p256 := btc.BTC()
	hashed := []byte("testing")
	priv, _ := GenerateKey(p256, rand.Reader)
	r, s, _ := Sign(rand.Reader, priv, hashed)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			Verify(&priv.PublicKey, hashed, r, s)
		}
	})
}

func BenchmarkRecoverSM2(b *testing.B) {
	//b.ResetTimer()
	//p256 := sm2.P256()
	p256 := sm2.SM2asm()
	hashed := []byte("testing")
	e := hashToInt(hashed, p256)
	priv, _ := GenerateKey(p256, rand.Reader)
	r, s, v, _ := p256.Sign(rand.Reader, e, priv.D)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256.Recover(r, s, e, v)
		}
	})
	/*
		for i := 0; i < b.N; i++ {
			p256.Recover(r, s, e, v)
		}
	*/
}

func BenchmarkRecoverSM2C(b *testing.B) {
	//b.ResetTimer()
	p256 := ecc.SM2C()
	hashed := []byte("testing")
	e := hashToInt(hashed, p256)
	priv, _ := GenerateKey(p256, rand.Reader)
	r, s, v, _ := p256.Sign(rand.Reader, e, priv.D)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256.Recover(r, s, e, v)
		}
	})
}

func BenchmarkRecoverBTC(b *testing.B) {
	//b.ResetTimer()
	p256 := btc.BTCasm()
	hashed := []byte("testing")
	e := hashToInt(hashed, p256)
	priv, _ := GenerateKey(p256, rand.Reader)
	r, s, v, _ := p256.Sign(rand.Reader, e, priv.D)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			p256.Recover(r, s, e, v)
		}
	})
}

func BenchmarkVerifySM2go(b *testing.B) {
	//b.ResetTimer()
	p256 := sm2.SM2go()
	hashed := []byte("testing")
	priv, _ := GenerateKey(p256, rand.Reader)
	r, s, _ := Sign(rand.Reader, priv, hashed)

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			Verify(&priv.PublicKey, hashed, r, s)
		}
	})
}

func BenchmarkKeyGeneration(b *testing.B) {
	//b.ResetTimer()
	p256 := elliptic.P256()

	//b.ReportAllocs()
	b.ResetTimer()
	b.RunParallel(func(pb *testing.PB) {
		for pb.Next() {
			GenerateKey(p256, rand.Reader)
		}
	})
}

func testSignAndVerify(t *testing.T, c elliptic.Curve, tag string) {
	priv, _ := GenerateKey(c, rand.Reader)

	hashed := []byte("testing")
	r, s, err := Sign(rand.Reader, priv, hashed)
	if err != nil {
		t.Errorf("%s: error signing: %s", tag, err)
		return
	}

	if !Verify(&priv.PublicKey, hashed, r, s) {
		t.Errorf("%s: Verify failed", tag)
		t.Logf("Pubkey X: %s\nY: %s\nr: %s", priv.PublicKey.X.Text(16),
			priv.PublicKey.Y.Text(16), r.Text(16))
	}

	hashed[0] ^= 0xff
	if Verify(&priv.PublicKey, hashed, r, s) {
		t.Errorf("%s: Verify always works!", tag)
	}
	t.Log(tag, " SignAndVerify test ✔")
}

func testSignAndVerify2(t *testing.T, c, c1 elliptic.Curve, tag string) {
	priv, _ := GenerateKey(c, rand.Reader)

	// Sign via Curve c
	hashed := []byte("testing")
	r, s, err := Sign(rand.Reader, priv, hashed)
	if err != nil {
		t.Errorf("%s: error signing: %s", tag, err)
		return
	}

	// Verify via Curve c1
	priv.PublicKey.Curve = c1
	if !Verify(&priv.PublicKey, hashed, r, s) {
		t.Errorf("%s: Verify failed", tag)
		t.Logf("Pubkey X: %s\nY: %s\nr: %s", priv.PublicKey.X.Text(16),
			priv.PublicKey.Y.Text(16), r.Text(16))
		return
	}

	hashed[0] ^= 0xff
	if Verify(&priv.PublicKey, hashed, r, s) {
		t.Errorf("%s: Verify always works!", tag)
	}
	t.Log(tag, " SignAndVerify2 test PASS ✅")
}

func testSignAndRecover(t *testing.T, c elliptic.Curve, tag string) {
	var opt signIntf
	if in, ok := c.(signIntf); !ok {
		t.Errorf("%s: not support sign/recover interface", tag)
	} else {
		opt = in
	}
	priv, _ := GenerateKey(c, rand.Reader)

	// Sign via Curve c
	hashed := []byte("testing")
	e := hashToInt(hashed, c)
	for i := 0; i < 10; i++ {
		r, s, v, err := opt.Sign(rand.Reader, e, priv.D)
		if err != nil {
			t.Errorf("%s: error signing: %s", tag, err)
			return
		}

		// recover
		if px, py, err := opt.Recover(r, s, e, v); err != nil {
			t.Errorf("%s: Recover failed", tag)
			t.Logf("Pubkey X: %s\nY: %s\nr: %s", priv.PublicKey.X.Text(16),
				priv.PublicKey.Y.Text(16), r.Text(16))
		} else if px.Cmp(priv.PublicKey.X) != 0 || py.Cmp(priv.PublicKey.Y) != 0 {
			t.Logf("%s: Recover point diff:\nX1: %s\nX2: %s\nY1: %s\nY2: %s", tag,
				priv.PublicKey.X.Text(16), px.Text(16),
				priv.PublicKey.Y.Text(16), py.Text(16))
		}
	}
	t.Log(tag, " SignAndRecover test PASS ✅")
}

func testSignAndRecover2(t *testing.T, c, c1 elliptic.Curve, tag string) {
	var opt, opt1 signIntf
	if in, ok := c.(signIntf); !ok {
		t.Errorf("%s: curve c not support sign/recover interface", tag)
		t.Fail()
	} else {
		opt = in
	}
	if in, ok := c.(signIntf); !ok {
		t.Errorf("%s: curve c not support sign/recover interface", tag)
		t.Fail()
	} else {
		opt1 = in
	}
	priv, _ := GenerateKey(c, rand.Reader)

	// Sign via Curve c
	hashed := []byte("testing")
	e := hashToInt(hashed, c)
	r, s, v, err := opt.Sign(rand.Reader, e, priv.D)
	if err != nil {
		t.Errorf("%s: error signing: %s", tag, err)
		return
	}

	// recover
	if px, py, err := opt1.Recover(r, s, e, v); err != nil {
		t.Errorf("%s: Recover failed", tag)
		t.Logf("Pubkey X: %s\nY: %s\nr: %s", priv.PublicKey.X.Text(16),
			priv.PublicKey.Y.Text(16), r.Text(16))
	} else if px.Cmp(priv.PublicKey.X) != 0 || py.Cmp(priv.PublicKey.Y) != 0 {
		t.Logf("%s: Recover point diff:\nX1: %s\nX2: %s\nY1: %s\nY2: %s", tag,
			priv.PublicKey.X.Text(16), px.Text(16),
			priv.PublicKey.Y.Text(16), py.Text(16))
	}
	t.Log(tag, " SignAndRecover2 test PASS ✅")
}

func TestSignAndVerify(t *testing.T) {
	testSignAndVerify(t, elliptic.P224(), "p224")
	if testing.Short() {
		return
	}
	testSignAndVerify(t, elliptic.P256(), "p256")
	testSignAndVerify(t, elliptic.P384(), "p384")
	testSignAndVerify(t, elliptic.P521(), "p521")
	testSignAndVerify(t, sm2.SM2(), "SM2asm")
	testSignAndVerify(t, sm2.P256(), "SM2go")
	testSignAndVerify(t, ecc.SM2C(), "SM2C")
	testSignAndVerify(t, btc.BTC(), "BTCasm")
	testSignAndVerify2(t, sm2.SM2(), ecc.SM2C(), "SM2asm vs SM2C")
	testSignAndVerify2(t, ecc.SM2C(), sm2.SM2(), "SM2C vs SM2asm")
	// sm2.P256() using ecdsa while sm2.SM2(), ecc.SM2C() using sm2 sign/verify
	//testSignAndVerify2(t, sm2.SM2(), sm2.P256(), "SM2asm vs SM2go")
	//testSignAndVerify2(t, sm2.P256(), sm2.SM2(), "SM2go vs SM2asm")
	testSignAndVerify2(t, sm2.SM2go(), sm2.P256(), "SM2p vs SM2go")
	testSignAndVerify2(t, sm2.P256(), sm2.SM2go(), "SM2go vs SM2p")
	testSignAndVerify2(t, btc.BTCgo(), btc.BTC(), "BTCgo vs BTCasm")
	testSignAndVerify2(t, btc.BTC(), btc.BTCgo(), "BTCasm vs BTCgo")
}

func TestSignAndRecover(t *testing.T) {
	testSignAndRecover(t, sm2.SM2(), "SM2asm")
	testSignAndRecover(t, ecc.SM2C(), "SM2C")
	testSignAndRecover(t, btc.BTC(), "BTCasm")
	testSignAndRecover2(t, sm2.SM2(), ecc.SM2C(), "SM2asm vs SM2C")
	testSignAndRecover2(t, ecc.SM2C(), sm2.SM2(), "SM2C vs SM2asm")
}

func testNonceSafety(t *testing.T, c elliptic.Curve, tag string) {
	priv, _ := GenerateKey(c, rand.Reader)

	hashed := []byte("testing")
	r0, s0, err := Sign(zeroReader, priv, hashed)
	if err != nil {
		t.Errorf("%s: error signing: %s", tag, err)
		return
	}

	hashed = []byte("testing...")
	r1, s1, err := Sign(zeroReader, priv, hashed)
	if err != nil {
		t.Errorf("%s: error signing: %s", tag, err)
		return
	}

	if s0.Cmp(s1) == 0 {
		// This should never happen.
		t.Errorf("%s: the signatures on two different messages were the same", tag)
	}

	if r0.Cmp(r1) == 0 {
		t.Errorf("%s: the nonce used for two different messages was the same", tag)
	}
}

func TestNonceSafety(t *testing.T) {
	testNonceSafety(t, elliptic.P224(), "p224")
	if testing.Short() {
		return
	}
	testNonceSafety(t, elliptic.P256(), "p256")
	testNonceSafety(t, elliptic.P384(), "p384")
	testNonceSafety(t, elliptic.P521(), "p521")
}

func testINDCCA(t *testing.T, c elliptic.Curve, tag string) {
	priv, _ := GenerateKey(c, rand.Reader)

	hashed := []byte("testing")
	r0, s0, err := Sign(rand.Reader, priv, hashed)
	if err != nil {
		t.Errorf("%s: error signing: %s", tag, err)
		return
	}

	r1, s1, err := Sign(rand.Reader, priv, hashed)
	if err != nil {
		t.Errorf("%s: error signing: %s", tag, err)
		return
	}

	if s0.Cmp(s1) == 0 {
		t.Errorf("%s: two signatures of the same message produced the same result", tag)
	}

	if r0.Cmp(r1) == 0 {
		t.Errorf("%s: two signatures of the same message produced the same nonce", tag)
	}
}

func TestINDCCA(t *testing.T) {
	testINDCCA(t, elliptic.P224(), "p224")
	if testing.Short() {
		return
	}
	testINDCCA(t, elliptic.P256(), "p256")
	testINDCCA(t, elliptic.P384(), "p384")
	testINDCCA(t, elliptic.P521(), "p521")
}

func fromHex(s string) *big.Int {
	r, ok := new(big.Int).SetString(s, 16)
	if !ok {
		panic("bad hex")
	}
	return r
}

func TestVectors(t *testing.T) {
	// This test runs the full set of NIST test vectors from
	// https://csrc.nist.gov/groups/STM/cavp/documents/dss/186-3ecdsatestvectors.zip
	//
	// The SigVer.rsp file has been edited to remove test vectors for
	// unsupported algorithms and has been compressed.

	if testing.Short() {
		return
	}

	f, err := os.Open("testdata/SigVer.rsp.bz2")
	if err != nil {
		t.Fatal(err)
	}

	buf := bufio.NewReader(bzip2.NewReader(f))

	lineNo := 1
	var h hash.Hash
	var msg []byte
	var hashed []byte
	var r, s *big.Int
	pub := new(PublicKey)

	for {
		line, err := buf.ReadString('\n')
		if len(line) == 0 {
			if err == io.EOF {
				break
			}
			t.Fatalf("error reading from input: %s", err)
		}
		lineNo++
		// Need to remove \r\n from the end of the line.
		if !strings.HasSuffix(line, "\r\n") {
			t.Fatalf("bad line ending (expected \\r\\n) on line %d", lineNo)
		}
		line = line[:len(line)-2]

		if len(line) == 0 || line[0] == '#' {
			continue
		}

		if line[0] == '[' {
			line = line[1 : len(line)-1]
			parts := strings.SplitN(line, ",", 2)

			switch parts[0] {
			case "P-224":
				pub.Curve = elliptic.P224()
			case "P-256":
				pub.Curve = elliptic.P256()
			case "P-384":
				pub.Curve = elliptic.P384()
			case "P-521":
				pub.Curve = elliptic.P521()
			default:
				pub.Curve = nil
			}

			switch parts[1] {
			case "SHA-1":
				h = sha1.New()
			case "SHA-224":
				h = sha256.New224()
			case "SHA-256":
				h = sha256.New()
			case "SHA-384":
				h = sha512.New384()
			case "SHA-512":
				h = sha512.New()
			default:
				h = nil
			}

			continue
		}

		if h == nil || pub.Curve == nil {
			continue
		}

		switch {
		case strings.HasPrefix(line, "Msg = "):
			if msg, err = hex.DecodeString(line[6:]); err != nil {
				t.Fatalf("failed to decode message on line %d: %s", lineNo, err)
			}
		case strings.HasPrefix(line, "Qx = "):
			pub.X = fromHex(line[5:])
		case strings.HasPrefix(line, "Qy = "):
			pub.Y = fromHex(line[5:])
		case strings.HasPrefix(line, "R = "):
			r = fromHex(line[4:])
		case strings.HasPrefix(line, "S = "):
			s = fromHex(line[4:])
		case strings.HasPrefix(line, "Result = "):
			expected := line[9] == 'P'
			h.Reset()
			h.Write(msg)
			hashed := h.Sum(hashed[:0])
			if Verify(pub, hashed, r, s) != expected {
				t.Fatalf("incorrect result on line %d", lineNo)
			}
		default:
			t.Fatalf("unknown variable on line %d: %s", lineNo, line)
		}
	}
}

func testNegativeInputs(t *testing.T, curve elliptic.Curve, tag string) {
	key, err := GenerateKey(curve, rand.Reader)
	if err != nil {
		t.Errorf("failed to generate key for %q", tag)
	}

	var hash [32]byte
	r := new(big.Int).SetInt64(1)
	r.Lsh(r, 550 /* larger than any supported curve */)
	r.Neg(r)

	if Verify(&key.PublicKey, hash[:], r, r) {
		t.Errorf("bogus signature accepted for %q", tag)
	}
}

func TestNegativeInputs(t *testing.T) {
	testNegativeInputs(t, elliptic.P224(), "p224")
	testNegativeInputs(t, elliptic.P256(), "p256")
	testNegativeInputs(t, elliptic.P384(), "p384")
	testNegativeInputs(t, elliptic.P521(), "p521")
}

func TestZeroHashSignature(t *testing.T) {
	zeroHash := make([]byte, 64)

	for _, curve := range []elliptic.Curve{elliptic.P224(), elliptic.P256(), elliptic.P384(), elliptic.P521()} {
		privKey, err := GenerateKey(curve, rand.Reader)
		if err != nil {
			panic(err)
		}

		// Sign a hash consisting of all zeros.
		r, s, err := Sign(rand.Reader, privKey, zeroHash)
		if err != nil {
			panic(err)
		}

		// Confirm that it can be verified.
		if !Verify(&privKey.PublicKey, zeroHash, r, s) {
			t.Errorf("zero hash signature verify failed for %T", curve)
		}
	}
}
