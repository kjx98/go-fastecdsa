#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <benchmark/benchmark.h>
#include "vli_bn.hpp"
#include "curve_const.hpp"
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


int main(int argc, char ** argv) {
	initData();
	benchmark::Initialize(&argc, argv);
	benchmark::RunSpecifiedBenchmarks();
}
