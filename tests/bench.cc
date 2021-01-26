#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <benchmark/benchmark.h>
#include "vli_bn.hpp"
#include "curve_const.hpp"
#include "curve_defs.hpp"
#include "mont.hpp"

using namespace vli;
using namespace ecc;

#include "testData.hpp"
//#define	sm2_p256	(*sm2_p256p)

static void test_montMult(benchmark::State &state)
{
	bignum<4>	xp;
	bignum<4>	bx1(dx1);
	for (auto _ : state) {
		xp.mont_mult(bx1, rr, prime, sm2_p_k0);
	}
}
BENCHMARK(test_montMult);

static void test_montMultP(benchmark::State &state)
{
	bignum<4>	xp;
	bignum<4>	bx1(dx1);
	for (auto _ : state) {
		mont_mult<sm2_p_k0>(xp, bx1, rr, prime);
	}
}
BENCHMARK(test_montMultP);

static void test_montSqr(benchmark::State &state)
{
	bignum<4>	xp;
	bignum<4>	bx1(dx1);
	for (auto _ : state) {
		xp.mont_sqr(bx1, prime, sm2_p_k0);
	}
}
BENCHMARK(test_montSqr);

static void test_inverse(benchmark::State &state)
{
	u64	res[4];
	for (auto _ : state) {
		vli_mod_inv<4>(res, dx1, sm2_p);
	}
}
BENCHMARK(test_inverse);

static void test_inverseNew(benchmark::State &state)
{
	u64	res[4];
	for (auto _ : state) {
		vli_mod_inv_new<4>(res, dx1, sm2_p);
	}
}
BENCHMARK(test_inverseNew);

static void test_mult2(benchmark::State &state)
{
	bignum<4>	x3(dx3), res;
	for (auto _ : state) {
		sm2_p256.mont_mult2(res, x3);
	}
}
BENCHMARK(test_mult2);

static void test_mult8(benchmark::State &state)
{
	bignum<4>	x3(dx3), res;
	for (auto _ : state) {
		res = x3;
		sm2_p256.mont_mult8(res);
	}
}
BENCHMARK(test_mult8);

static void test_mult2k(benchmark::State &state)
{
	bignum<4>	x3(dx3), res;
	for (auto _ : state) {
		sm2_k256.mont_mult2(res, x3);
	}
}
BENCHMARK(test_mult2k);

static void test_mult8k(benchmark::State &state)
{
	bignum<4>	x3(dx3), res;
	for (auto _ : state) {
		res = x3;
		sm2_k256.mont_mult8(res);
	}
}
BENCHMARK(test_mult8k);

static void test_div2(benchmark::State &state)
{
	bignum<4>	x3(dx3), res;
	for (auto _ : state) {
		res = x3;
		sm2_p256.mont_div2(res);
	}
}
BENCHMARK(test_div2);

static void test_sqrt(benchmark::State &state)
{
	bignum<4>	yy1(dy1y1), res;
	for (auto _ : state) {
		sm2_p256.mod_sqrt(res, yy1);
	}
}
BENCHMARK(test_sqrt);

static void test_ksqrt(benchmark::State &state)
{
	bignum<4>	yy1(dy1y1), res;
	for (auto _ : state) {
		sm2_k256.mod_sqrt(res, yy1);
	}
}
BENCHMARK(test_ksqrt);

static void test_ptRecovery(benchmark::State &state)
{
	bignum<4>	x1(d1Gx);
	bignum<4>	y1(d1Gy), res;
	for (auto _ : state) {
		point_recovery(sm2_p256, res, x1, y1.is_even());
	}
}
BENCHMARK(test_ptRecovery);

static void test_ptRecoveryK(benchmark::State &state)
{
	bignum<4>	x1(d1Gx);
	bignum<4>	y1(d1Gy), res;
	for (auto _ : state) {
		point_recovery(sm2_k256, res, x1, y1.is_even());
	}
}
BENCHMARK(test_ptRecoveryK);

static void test_ECADDJac(benchmark::State &state)
{
	u64		xx3[4], yy3[4], zz3[4];
	for (auto _ : state) {
		sm2_p256.point_add_jacobian(xx3, yy3, zz3, dx1, dy1,
						bigOne.data(), dx2, dy2);
	}
}
BENCHMARK(test_ECADDJac);

static void test_ECDBLJac(benchmark::State &state)
{
	bignum<4>	bigOne(1);
	u64		xx3[4], yy3[4], zz3[4];
	for (auto _ : state) {
		//sm2_p256.point_double_jacobian(xx3, yy3, zz3, dx3, dy3, dz3);
		sm2_p256.point_double_jacobian(xx3, yy3, zz3, dx1, dy1);
	}
}
BENCHMARK(test_ECDBLJac);

static void test_ECADDJacK(benchmark::State &state)
{
	u64		xx3[4], yy3[4], zz3[4];
	for (auto _ : state) {
		sm2_k256.point_add_jacobian(xx3, yy3, zz3, dx1, dy1,
						bigOne.data(), dx2, dy2);
	}
}
BENCHMARK(test_ECADDJacK);

static void test_ECDBLJacK(benchmark::State &state)
{
	bignum<4>	bigOne(1);
	u64		xx3[4], yy3[4], zz3[4];
	for (auto _ : state) {
		//sm2_k256.point_double_jacobian(xx3, yy3, zz3, dx3, dy3, dz3);
		sm2_k256.point_double_jacobian(xx3, yy3, zz3, dx1, dy1);
	}
}
BENCHMARK(test_ECDBLJacK);

static void test_ECScalarMultBase(benchmark::State &state)
{
	bignum<4>	d1(d1d);
	point_t<4>	res;
	for (auto _ : state) {
		sm2_k256.scalar_mult_base(res, d1);
	}
}
BENCHMARK(test_ECScalarMultBase);

static void test_ECScalarMultBaseN(benchmark::State &state)
{
	bignum<4>	d1(d1d);
	point_t<4>	res;
	for (auto _ : state) {
		sm2_p256.scalar_mult_base(res, d1);
	}
}
BENCHMARK(test_ECScalarMultBaseN);

static void test_ECScalarMult(benchmark::State &state)
{
	bignum<4>	d1(d1d);
	point_t<4>	res;
	point_t<4>	gg(sm2_gx, sm2_gy);
	for (auto _ : state) {
		sm2_p256.scalar_mult(res, gg, d1);
	}
}
BENCHMARK(test_ECScalarMult);

static void test_ECScalarMultNAF2(benchmark::State &state)
{
	bignum<4>	d1(d1d);
	point_t<4>	res;
	point_t<4>	gg(sm2_gx, sm2_gy);
	for (auto _ : state) {
		sm2_p256.scalar_multNAF2(res, gg, d1);
	}
}
BENCHMARK(test_ECScalarMultNAF2);

static void test_ECScalarMult256(benchmark::State &state)
{
	bignum<4>	d1(d1d);
	point_t<4>	res;
	point_t<4>	gg(sm2_gx, sm2_gy);
	bignum<4>	bn_zero(0ul);
	for (auto _ : state) {
		sm2_k256.combined_mult(res, gg, d1, bn_zero);
	}
}
BENCHMARK(test_ECScalarMult256);

static void test_ECScalarCMult(benchmark::State &state)
{
	point_t<4>	pt1(dx1, dy1);
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	for (auto _ : state) {
		sm2_k256.combined_mult(res, pt1, d1, d2);
	}
}
BENCHMARK(test_ECScalarCMult);

static void test_ECScalarCMultNAF(benchmark::State &state)
{
	point_t<4>	pt1(dx1, dy1);
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	for (auto _ : state) {
		sm2_p256.combined_mult(res, pt1, d1, d2);
	}
}
BENCHMARK(test_ECScalarCMultNAF);

int main(int argc, char ** argv) {
	//sm2_p256.init();
	//sm2_k256.init();
	benchmark::Initialize(&argc, argv);
	benchmark::RunSpecifiedBenchmarks();
}
