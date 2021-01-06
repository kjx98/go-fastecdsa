#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <benchmark/benchmark.h>
#include "vli_bn.hpp"
#include "curve_const.hpp"
#include "curve_defs.hpp"
#include "mont.hpp"

using namespace vli;

#include "testData.h"

static void test_montMult(benchmark::State &state)
{
	bignum<4>	xp;
	for (auto _ : state) {
		xp.mont_mult(bx1, rr, prime, sm2_p_k0);
	}
}
BENCHMARK(test_montMult);

static void test_montSqr(benchmark::State &state)
{
	bignum<4>	xp;
	for (auto _ : state) {
		xp.mont_sqr(bx1, prime, sm2_p_k0);
	}
}
BENCHMARK(test_montSqr);

static void test_ECADDJac(benchmark::State &state)
{
	u64		x3[4], y3[4], z3[4];
	for (auto _ : state) {
		sm2_p256.point_add_jacobian(x3, y3, z3, bx1.data(), by1.data(),
						bigOne.data(), bx2.data(), by2.data());
	}
}
BENCHMARK(test_ECADDJac);

static void test_ECDBLJac(benchmark::State &state)
{
	bignum<4>	bigOne(1);
	u64		x3[4], y3[4], z3[4];
	for (auto _ : state) {
		sm2_p256.point_double_jacobian(x3, y3, z3, bx1.data(), by1.data());
	}
}
BENCHMARK(test_ECDBLJac);



int main(int argc, char ** argv) {
	initData();
	sm2_p256.init();
	benchmark::Initialize(&argc, argv);
	benchmark::RunSpecifiedBenchmarks();
}
