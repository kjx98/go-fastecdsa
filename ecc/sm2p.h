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


void sm2_mod_inv(bn_words *result, const bn_words *in);
void sm2_mod_mul(bn_words *result, const bn_words *x, const bn_words *y);
void sm2_point_add(Point *result, const Point *p, const Point *q);
void sm2_point_add_jacobian(Point *result, const Point *p, const Point *q);
void sm2_point_double(Point *result, const Point *p);
void sm2_point_double_jacobian(Point *result, const Point *p);

#endif	//	__SM2P_H__
