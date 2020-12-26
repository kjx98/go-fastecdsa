// +build sm2p256

#pragma once
#ifndef	__SM2P_H__
#define	__SM2P_H__

# include <stdint.h>
# include <stdbool.h>
# include <sys/types.h>
# include <string.h>

# if defined(__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
  /* even with gcc, the typedef won't work for 32-bit platforms */
typedef __uint128_t uint128_t;  /* nonstandard; implemented by gcc on 64-bit
                                 * platforms */
typedef __int128_t int128_t;
# else
#  error "Need GCC 3.1 or later to define type uint128_t"
# endif

typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t s64;

void sm2_mod_inv(u64 *result, const u64 *in, const u64 *mod);
void sm2_mod_mul(u64 *result, const u64 *x, const u64 *y);
void sm2_point_add(u64 *result, const u64 *p, const u64 *q);

#endif	//	__SM2P_H__
