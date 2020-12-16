/*
 * Copyright (c) 2013, Kenneth MacKay
 * All rights reserved.
 * Copyright(c) 2020 Jesse Kuang  <jkuang@21cn.com>
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
#ifndef __ECC_IMPL_H__
#define __ECC_IMPL_H__

#include "cdefs.h"

/* One digit is u64 qword. */
#define ECC_MAX_DIGITS             (512 / 64)

/**
 * struct ecc_curve - definition of elliptic curve
 *
 * @name:	Short name of the curve.
 * @g:		Generator point of the curve.
 * @p:		Prime number, if Barrett's reduction is used for this curve
 *		pre-calculated value 'mu' is appended to the @p after ndigits.
 *		Use of Barrett's reduction is heuristically determined in
 *		vli_mmod_fast().
 * @n:		Order of the curve group.
 * @a:		Curve parameter a.
 * @b:		Curve parameter b.
 */
struct ecc_curve {
	ecc_curve(const char *_name, const u64 *_gx, const u64 *_gy, const u64 *_p,
				const u64 *_n, const u64 *_a, const u64 *_b, const uint ndig=4,
				const bool bBat=false) : name(_name), gx(_gx), gy(_gy),
			p(_p), n(_n), a(_a), b(_b),
			ndigits(ndig), use_barrett(bBat) {}
	const char *name;
	const u64 *gx;
	const u64 *gy;
	const u64 *p;
	const u64 *n;
	const u64 *a;
	const u64 *b;
	const uint ndigits = 4;
	const bool use_barrett = false;
};

#endif	//	__ECC_IMPL_H__
