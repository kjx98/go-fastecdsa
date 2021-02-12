#pragma once
#ifndef	__SM2P_H__
#define	__SM2P_H__

# include <string.h>
# include "cdefs.h"

# if defined(__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
  /* even with gcc, the typedef won't work for 32-bit platforms */
typedef __uint128_t uint128_t;  /* nonstandard; implemented by gcc on 64-bit
                                 * platforms */
typedef __int128_t int128_t;
# else
#  error "Need GCC 3.1 or later to define type uint128_t"
# endif


void sm2_mod_inv(u64 *result, const u64 *in);
//void sm2_mod_mul(u64 *result, const u64 *x, const u64 *y);
void sm2_point_add(Point *result, const Point *p, const Point *q);
void sm2_point_add_jacobian(Point *result, const Point *p, const Point *q);
void sm2_point_double(Point *result, const Point *p);
void sm2_point_double_jacobian(Point *result, const Point *p);
void sm2_scalar_base_mult(Point *result, const u8 *scalar);
void sm2_scalar_mult(Point *result,const Point *p,  const bn_words_t scalar);

#endif	//	__SM2P_H__
