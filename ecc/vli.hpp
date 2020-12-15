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
#ifndef __VLI_H__
#define __VLI_H__

#include <stdint.h>
#include <stdbool.h>
#include <cassert>
#include <type_traits>
#include <endian.h>
typedef	__uint32_t u32;
typedef	__uint64_t u64;
typedef __uint64_t be64;
typedef __uint8_t u8;
typedef	unsigned int	uint;

#if	__cplusplus < 201103L
# error "C++ std MUST at least c++11"
#endif

#ifdef	__GNUC__
# define unlikely(cond)	__builtin_expect ((cond), 0)
# define likely(cond)	__builtin_expect (!!(cond), 1)
#define forceinline __inline__ __attribute__((always_inline))
#else
# error "MUST compiled by gcc/clang"
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
static void vli_clear(u64 *vli)
{
	uint i;
#pragma unroll 4
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
static bool vli_is_zero(const u64 *vli)
{
	uint i;
#pragma unroll 4
	for (i = 0; i < ndigits; i++) {
		if (vli[i]) return false;
	}
	return true;
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
int vli_cmp(const u64 *left, const u64 *right)
{
	int i;
#pragma unroll 4
	for (i = ndigits - 1; i >= 0; i--) {
		if (left[i] > right[i]) return 1;
		else if (left[i] < right[i]) return -1;
	}
	return 0;
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
u64 vli_sub(u64 *result, const u64 *left, const u64 *right)
{
	u64 borrow = 0;
	int i;

#pragma unroll 4
	for (i = 0; i < ndigits; i++) {
		u64 diff;
		diff = left[i] - right[i] - borrow;
		if (diff != left[i]) borrow = (diff > left[i]);
		result[i] = diff;
	}

	return borrow;
}

/**
 * vli_from_be64() - Load vli from big-endian u64 array
 *
 * @dest:		destination vli
 * @src:		source array of u64 BE values
 * @ndigits:		length of both vli and array
 */
template<uint ndigits> forceinline
void vli_from_be64(u64 *dest, const void *src)
{
	int i;
	const u64 *from = src;
#pragma unroll 4
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
void vli_from_le64(u64 *dest, const void *src)
{
	int i;
	const u64 *from = src;
#pragma unroll 4
	for (i = 0; i < ndigits; i++)
		dest[i] = le64toh(from[ndigits - 1 - i]);
}

#endif
