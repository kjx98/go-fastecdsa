// Copyright 2018 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fastecdsa

import (
	"crypto/elliptic"
	"crypto/rand"
	"crypto/sha256"
	"fmt"
)

func Example() {
	privateKey, err := GenerateKey(elliptic.P256(), rand.Reader)
	if err != nil {
		panic(err)
	}

	msg := "hello, world"
	hash := sha256.Sum256([]byte(msg))

	r, s, err := Sign(rand.Reader, privateKey, hash[:])
	if err != nil {
		panic(err)
	}
	fmt.Printf("signature: (0x%x, 0x%x)\n", r, s)

	valid := Verify(&privateKey.PublicKey, hash[:], r, s)
	fmt.Println("signature verified:", valid)
}
