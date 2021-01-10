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
#ifndef __CDEFS_H__
#define __CDEFS_H__

#include <stdint.h>
#include <stdbool.h>
#include <sys/types.h>

typedef	__uint32_t u32;
typedef	__uint64_t u64;
typedef __uint64_t be64;
typedef __uint8_t u8;
typedef __int64_t s64;
typedef	unsigned int	uint;

// gcc before 6 w/out builtin_usubl_overflow...
#if	__GNUC__ > 5 || __clang_major__ > 6
// has buildin_usubl_overflow...
#elif defined(__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
  /* even with gcc, the typedef won't work for 32-bit platforms */
# define	NO_BUILTIN_OVERFLOW
#else
#  error "Need GCC 3.1 or later to define type uint128_t"
#endif

# define unlikely(cond)	__builtin_expect ((cond), 0)
# define likely(cond)	__builtin_expect (!!(cond), 1)
# define forceinline __inline__ __attribute__((always_inline))

typedef u64	bn_words_t[4];

typedef	struct {
	bn_words_t	x;
	bn_words_t	y;
	bn_words_t	z;
}	Point;

#endif	//	__CDEFS_H__
