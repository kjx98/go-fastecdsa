#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <benchmark/benchmark.h>
#include "vli_bn.hpp"
#include "curve_const.hpp"
#include "curve_defs.hpp"
#include "mont.hpp"

using namespace vli;

#include "testData.hpp"

static void test_montMult(benchmark::State &state)
{
	bignum<4>	xp;
	bignum<4>	bx1(dx1);
	for (auto _ : state) {
		xp.mont_mult(bx1, rr, prime, sm2_p_k0);
	}
}
BENCHMARK(test_montMult);

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


int main(int argc, char ** argv) {
	sm2_p256.init();
	benchmark::Initialize(&argc, argv);
	benchmark::RunSpecifiedBenchmarks();
}
