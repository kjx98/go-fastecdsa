/*
 * Copyright(c) 2020 Jesse Kuang  <jkuang@21cn.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#pragma once
#ifndef __VLI_HPP__
#define __VLI_HPP__

#include <stdint.h>
#include <stdbool.h>
#include <endian.h>
#include <cassert>
#include <type_traits>
#include "cdefs.h"

#if	__cplusplus < 201103L
# error "C++ std MUST at least c++11"
#endif

#ifndef	ARRAY_SIZE
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))
#endif

template<typename T> forceinline
T max(const T a, const T b) noexcept {
	static_assert(std::is_integral<T>::value, "Integral required.");
	return (a>b)?a:b;
}


template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static void vli_clear(u64 *vli)
{
	uint i;
	for (i = 0; i < ndigits; i++)
	vli[i] = 0;
}

/**
 * vli_is_zero() - Determine is vli is zero
 *
 * @vli:		vli to check.
 * @ndigits:		length of the @vli
 */
/* Returns true if vli == 0, false otherwise. */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static bool vli_is_zero(const u64 *vli)
{
	uint i;
//#pragma unroll 4
	for (i = 0; i < ndigits; i++) {
		if (vli[i]) return false;
	}
	return true;
}

/* Sets dest = src. */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static void vli_set(u64 *dest, const u64 *src)
{
	uint i;

	for (i = 0; i < ndigits; i++)
		dest[i] = src[i];
}

/**
 * vli_cmp() - compare left and right vlis
 *
 * @left:		vli
 * @right:		vli
 * @ndigits:		length of both vlis
 *
 * Returns sign of @left - @right, i.e. -1 if @left < @right,
 * 0 if @left == @right, 1 if @left > @right.
 */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
int vli_cmp(const u64 *left, const u64 *right)
{
	int i;
//#pragma unroll 4
	for (i = ndigits - 1; i >= 0; i--) {
		if (left[i] > right[i]) return 1;
		else if (left[i] < right[i]) return -1;
	}
	return 0;
}

/* Computes result = in << c, returning carry. Can modify in place
 * (if result == in). 0 < shift < 64.
 */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static u64 vli_lshift(u64 *result, const u64 *in, unsigned int shift)
{
	u64 carry = 0;
	int i;

	for (i = 0; i < ndigits; i++) {
		u64 temp = in[i];

		result[i] = (temp << shift) | carry;
		carry = temp >> (64 - shift);
	}

	return carry;
}

/* Computes vli = vli >> 1. */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static void vli_rshift1(u64 *vli)
{
	u64 *end = vli;
	u64 carry = 0;

	vli += ndigits;

	while (vli-- > end) {
		u64 temp = *vli;
		*vli = (temp >> 1) | carry;
		carry = temp << 63;
	}
}

/* Computes result = left + right, returning carry. Can modify in place. */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static u64 vli_add(u64 *result, const u64 *left, const u64 *right)
{
	u64 carry = 0;
	int i;

	for (i = 0; i < ndigits; i++) {
		u64 sum;

		sum = left[i] + right[i] + carry;
		if (sum != left[i])
			carry = (sum < left[i]);

		result[i] = sum;
	}

	return carry;
}

/* Computes result = left + right, returning carry. Can modify in place. */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static u64 vli_uadd(u64 *result, const u64 *left, u64 right)
{
	u64 carry = right;
	int i;

	for (i = 0; i < ndigits; i++) {
		u64 sum;

		sum = left[i] + carry;
		if (sum != left[i])
			carry = (sum < left[i]);
		else
			carry = !!carry;

		result[i] = sum;
	}

	return carry;
}

/**
 * vli_sub() - Subtracts right from left
 *
 * @result:		where to write result
 * @left:		vli
 * @right		vli
 * @ndigits:		length of all vlis
 *
 * Note: can modify in-place.
 *
 * Return: carry bit.
 */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static u64 vli_sub(u64 *result, const u64 *left, const u64 *right)
{
	u64 borrow = 0;
	int i;

//#pragma unroll 4
	for (i = 0; i < ndigits; i++) {
		u64 diff;
		diff = left[i] - right[i] - borrow;
		if (diff != left[i]) borrow = (diff > left[i]);
		result[i] = diff;
	}

	return borrow;
}

/* Computes result = left - right, returning borrow. Can modify in place. */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static u64 vli_usub(u64 *result, const u64 *left, u64 right)
{
	u64 borrow = right;
	int i;

	for (i = 0; i < ndigits; i++) {
		u64 diff;

		diff = left[i] - borrow;
		if (diff != left[i])
			borrow = (diff > left[i]);

		result[i] = diff;
	}

	return borrow;
}

/* Counts the number of 64-bit "digits" in vli. */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static uint vli_num_digits(const u64 *vli)
{
	int i;

	/* Search from the end until we find a non-zero digit.
	 * We do it in reverse because we expect that most digits will
	 * be nonzero.
	 */
	for (i = ndigits - 1; i >= 0 && vli[i] == 0; i--);

	return (i + 1);
}

/* Counts the number of bits required for vli. */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
static uint vli_num_bits(const u64 *vli)
{
	uint i, num_digits;
	u64 digit;

	num_digits = vli_num_digits<ndigits>(vli);
	if (num_digits == 0) return 0;

	digit = vli[num_digits - 1];
	i = 64 - __builtin_clzl(digit);

	return ((num_digits - 1) * 64 + i);
}

/**
 * vli_from_be64() - Load vli from big-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 BE values
 * @ndigits:		length of both vli and array
 */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
void vli_from_be64(u64 *dest, const void *src)
{
	int i;
	const u64 *from = src;
//#pragma unroll 4
	for (i = 0; i < ndigits; i++)
		dest[i] = be64toh(from[ndigits - 1 - i]);
}

/**
 * vli_from_le64() - Load vli from little-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 LE values
 * @ndigits:		length of both vli and array
 */
template<uint ndigits> forceinline
__attribute__((optimize("unroll-loops")))
void vli_from_le64(u64 *dest, const void *src)
{
	int i;
	const u64 *from = src;
//#pragma unroll 4
	for (i = 0; i < ndigits; i++)
		dest[i] = le64toh(from[ndigits - 1 - i]);
}

#endif	//	__VLI_HPP__
