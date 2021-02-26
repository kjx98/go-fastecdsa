#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <benchmark/benchmark.h>
#include "vli_bn.hpp"
#include "curve_const.hpp"
#include "curve_defs.hpp"
#include "ecc_key.hpp"
#include "mont.hpp"

using namespace vli;
using namespace ecc;

#include "testData.hpp"
//#define	sm2_p256	(*sm2_p256p)
bignum<4>	tt;

static void test_sm2MultP(benchmark::State &state)
{
	u64		r[6];
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		vli_sm2_multP(r, 0x33445566);
		tt = bignum<4>(r);
	}
}
BENCHMARK(test_sm2MultP);

static void test_montRed(benchmark::State &state)
{
	bignum<4>	xp;
	bignum<4>	bx1(dx1);
	xp.mont_mult(bx1, rr, prime, sm2_p_k0);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		bx1.mont_reduction(xp, prime, sm2_p_k0);
	}
	tt = bx1;
}
BENCHMARK(test_montRed);

static void test_montRedP(benchmark::State &state)
{
	bignum<4>	xp;
	bignum<4>	bx1(dx1);
	xp.mont_mult(bx1, rr, prime, sm2_p_k0);
	u64   *resp = reinterpret_cast<u64 *>(&bx1);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		mont_reduction<4,sm2_p_k0>(resp, xp.data(), prime.data());
	}
	tt = bx1;
}
BENCHMARK(test_montRedP);

static void test_montSM2Red(benchmark::State &state)
{
	bignum<4>	xp;
	bignum<4>	bx1(dx1);
	u64		res[4];
	xp.mont_mult(bx1, rr, prime, sm2_p_k0);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		sm2p_reduction(res, xp.data());
	}
	tt = bignum<4>(res);
}
BENCHMARK(test_montSM2Red);

static void test_montSM2RedN(benchmark::State &state)
{
	bignum<4>	xp;
	bignum<4>	bx1(dx1);
	u64		res[4];
	xp.mont_mult(bx1, rr, prime, sm2_p_k0);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		sm2p_reductionN(res, xp.data());
		tt = bignum<4>(res);
	}
}
BENCHMARK(test_montSM2RedN);

static void test_montMult(benchmark::State &state)
{
	bignum<4>	xp;
	bignum<4>	bx1(dx1);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		xp.mont_mult(bx1, rr, prime, sm2_p_k0);
	}
	tt.mont_reduction(xp, prime, sm2_p_k0);
}
BENCHMARK(test_montMult);

static void test_montMultP(benchmark::State &state)
{
	bignum<4>	xp;
	//bignum<4>	bx1(dx1);
	u64   *resp = reinterpret_cast<u64 *>(&xp);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		mont_mult<4,sm2_p_k0>(resp, dx1, rr.data(), prime.data());
	}
	tt.mont_reduction(xp, prime, sm2_p_k0);
}
BENCHMARK(test_montMultP);

static void test_montMultK01(benchmark::State &state)
{
	u64	xp[4];
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		sm2p_mult(xp, dx1, rr.data());
	}
	tt.mont_reduction(xp, prime, 1);
}
BENCHMARK(test_montMultK01);

static void test_montMultK01N(benchmark::State &state)
{
	u64	xp[4];
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		sm2p_multN(xp, dx1, rr.data());
	}
	tt.mont_reduction(xp, prime, 1);
}
BENCHMARK(test_montMultK01N);

static void test_montSqr(benchmark::State &state)
{
	bignum<4>	xp, bp;
	bignum<4>	bx1(dx1);
	bp.mont_mult(bx1, rr, prime, sm2_p_k0);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		xp.mont_sqr(bp, prime, sm2_p_k0);
	}
	tt.mont_reduction(xp, prime, sm2_p_k0);
}
BENCHMARK(test_montSqr);

static void test_montSqrP(benchmark::State &state)
{
	bignum<4>	xp, bp;
	bignum<4>	bx1(dx1);
	u64   *resp = reinterpret_cast<u64 *>(&xp);
	bp.mont_mult(bx1, rr, prime, sm2_p_k0);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		mont_sqr<4,sm2_p_k0>(resp, bp.data(), prime.data());
	}
	tt.mont_reduction(xp, prime, sm2_p_k0);
}
BENCHMARK(test_montSqrP);

static void test_montSqrN(benchmark::State &state)
{
	bignum<4>	bp;
	bignum<4>	bx1(dx1);
	u64			res[4];
	bp.mont_mult(bx1, rr, prime, sm2_p_k0);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
			mont_sqrN<4,sm2_p_k0>(res, bp.data(), prime.data());
	}
	bignum<4>	xp(res);
	tt.mont_reduction(xp, prime, sm2_p_k0);
}
BENCHMARK(test_montSqrN);

static void test_sm2montSqrN(benchmark::State &state)
{
	bignum<4>	bp;
	bignum<4>	bx1(dx1);
	u64			res[4];
	bp.mont_mult(bx1, rr, prime, sm2_p_k0);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
			sm2p_sqrN(res, bp.data());
	}
	bignum<4>	xp(res);
	tt.mont_reduction(xp, prime, sm2_p_k0);
}
BENCHMARK(test_sm2montSqrN);

static void test_montSqrK01(benchmark::State &state)
{
	u64		res[4];
	u64		bx1[4];
	sm2p_mult(bx1, dx1, rr.data());
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		sm2p_mult(res, bx1, bx1);
	}
	tt.mont_reduction(res, prime, 1);
}
BENCHMARK(test_montSqrK01);

static void test_montBtcRed(benchmark::State &state)
{
	bignum<4>	xp;
	bignum<4>	bx1(dx1);
	u64		res[4];
	xp.mont_mult(bx1, btc_rr, btc_prime, secp256k1_p_k0);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
		btc_reduction(res, xp.data());
		tt = bignum<4>(res);
	}
}
BENCHMARK(test_montBtcRed);

static void test_BtcMontSqrN(benchmark::State &state)
{
	bignum<4>	bp;
	bignum<4>	bx1(dx1);
	u64			res[4];
	bp.mont_mult(bx1, btc_rr, btc_prime, secp256k1_p_k0);
	for (auto _ : state) {
		for (int i=0; i<1000; ++i)
			btc_sqrN(res, bp.data());
	}
	bignum<4>	xp(res);
	tt.mont_reduction(xp, btc_prime, secp256k1_p_k0);
}
BENCHMARK(test_BtcMontSqrN);


static void test_bnMult(benchmark::State &state)
{
	bignum<4>	bx1(dx1);
	bignum<4>	bx2(dx2);
	bn_prod<4>	xp;
	for (auto _ : state) {
		for (int i=0; i<1000; ++i) xp.mult(bx1, bx1);
	}
	tt = xp.bn256();
}
BENCHMARK(test_bnMult);

static void test_vliSqr(benchmark::State &state)
{
	u64		res[8];
	for (auto _ : state) {
		for (int i=0; i<1000; ++i) vli_squareOld<4>(res, dx1);
	}
	tt = bignum<4>(res);
}
BENCHMARK(test_vliSqr);

static void test_vliSqrN(benchmark::State &state)
{
	u64		res[8];
	for (auto _ : state) {
		for (int i=0; i<1000; ++i) vli_square<4>(res, dx1);
	}
	tt = bignum<4>(res);
}
BENCHMARK(test_vliSqrN);

static void test_bnSqr(benchmark::State &state)
{
	bignum<4>	bx1(dx1);
	bn_prod<4>	xp;
	for (auto _ : state) {
		for (int i=0; i<1000; ++i) xp.squareOld(bx1);
	}
	tt = xp.bn256();
}
BENCHMARK(test_bnSqr);

static void test_bnSqrN(benchmark::State &state)
{
	bignum<4>	bx1(dx1);
	bn_prod<4>	xp;
	for (auto _ : state) {
		for (int i=0; i<1000; ++i) xp.square(bx1);
	}
	tt = xp.bn256();
}
BENCHMARK(test_bnSqrN);

static void test_inverse(benchmark::State &state)
{
	u64	res[4];
	for (auto _ : state) {
		vli_mod_inv<4>(res, dx1, sm2_p);
		tt = bignum<4>(res);
	}
}
BENCHMARK(test_inverse);

static void test_inverseNew(benchmark::State &state)
{
	u64	res[4];
	for (auto _ : state) {
		vli_mod_inv_new<4>(res, dx1, sm2_p);
		tt = bignum<4>(res);
	}
}
BENCHMARK(test_inverseNew);

static void test_mult2(benchmark::State &state)
{
	bignum<4>	x3(dx3), res;
	for (auto _ : state) {
		for (int i=0; i<1000; ++i) {
			sm2_p256.mont_mult2(res, x3);
			tt = res;
		}
	}
}
BENCHMARK(test_mult2);

static void test_mult8(benchmark::State &state)
{
	bignum<4>	x3(dx3), res;
	for (auto _ : state) {
		for (int i=0; i<1000; ++i) {
			res = x3;
			sm2_p256.mont_mult8(res);
			tt = res;
		}
	}
}
BENCHMARK(test_mult8);

static void test_mult2k(benchmark::State &state)
{
	bignum<4>	x3(dx3), res;
	for (auto _ : state) {
		for (int i=0; i<1000; ++i) {
			sm2_k256.mont_mult2(res, x3);
			tt = res;
		}
	}
}
BENCHMARK(test_mult2k);

static void test_mult8k(benchmark::State &state)
{
	bignum<4>	x3(dx3), res;
	for (auto _ : state) {
		for (int i=0; i<1000; ++i) {
			res = x3;
			sm2_k256.mont_mult8(res);
			tt = res;
		}
	}
}
BENCHMARK(test_mult8k);

static void test_div2(benchmark::State &state)
{
	bignum<4>	x3(dx3), res;
	for (auto _ : state) {
		for (int i=0; i<1000; ++i) {
			res = x3;
			sm2_p256.mont_div2(res);
			tt = res;
		}
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

static void test_ptRecover(benchmark::State &state)
{
	bignum<4>	x1(d1Gx);
	bignum<4>	y1(d1Gy), res;
	for (auto _ : state) {
		pointY_recover(sm2_p256, res, x1, y1.is_odd());
	}
}
BENCHMARK(test_ptRecover);

static void test_ptRecoverK(benchmark::State &state)
{
	bignum<4>	x1(d1Gx);
	bignum<4>	y1(d1Gy), res;
	for (auto _ : state) {
		pointY_recover(sm2_k256, res, x1, y1.is_odd());
	}
}
BENCHMARK(test_ptRecoverK);

static void test_ECADDJac(benchmark::State &state)
{
	point_t<4>	pt1, pt2, pt3;
	sm2_p256.to_affined(pt1, dx1, dy1);
	sm2_p256.to_affined(pt2, dx2, dy2);
	for (auto _ : state) {
		sm2_p256.point_add(pt3, pt1, pt2.x, pt2.y);
	}
}
BENCHMARK(test_ECADDJac);

static void test_ECDBLJac(benchmark::State &state)
{
	point_t<4>	pt1, pt3;
	sm2_p256.to_affined(pt1, dx1, dy1);
	for (auto _ : state) {
		sm2_p256.point_double(pt3, pt1.x, pt1.y);
	}
}
BENCHMARK(test_ECDBLJac);

static void test_ECADDJacK(benchmark::State &state)
{
	point_t<4>	pt1, pt2, pt3;
	sm2_p256.to_affined(pt1, dx1, dy1);
	sm2_p256.to_affined(pt2, dx2, dy2);
	for (auto _ : state) {
		sm2_k256.point_add(pt3, pt1, pt2.x, pt2.y);
	}
}
BENCHMARK(test_ECADDJacK);

static void test_ECDBLJacK(benchmark::State &state)
{
	point_t<4>	pt1, pt3;
	sm2_p256.to_affined(pt1, dx1, dy1);
	for (auto _ : state) {
		sm2_k256.point_double(pt3, pt1.x, pt1.y);
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

static void test_ECScalarCMultK(benchmark::State &state)
{
	point_t<4>	pt1(dx1, dy1);
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	for (auto _ : state) {
		sm2_k256.combined_mult(res, pt1, d1, d2);
	}
}
BENCHMARK(test_ECScalarCMultK);

static void test_ECScalarCMultP(benchmark::State &state)
{
	point_t<4>	pt1(dx1, dy1);
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	for (auto _ : state) {
		sm2_p256.combined_mult(res, pt1, d1, d2);
	}
}
BENCHMARK(test_ECScalarCMultP);

static void test_ECGenKey(benchmark::State &state)
{
	for (auto _ : state) {
		private_key<4>	priv(sm2_k256);
	}
}
BENCHMARK(test_ECGenKey);

static void test_ECGenKeyK(benchmark::State &state)
{
	bignum<4>	px, py, dd;
	for (auto _ : state) {
		gen_keypair(sm2_k256, dd, px, py);
	}
}
BENCHMARK(test_ECGenKeyK);

static void test_ECSign(benchmark::State &state)
{
	private_key<4>	priv(sm2_k256);
	bignum<4>	msg = bn_random<4>::Instance().get_random();
	bignum<4>	r, s;
	int		ecInd __attribute__ ((unused));
	for (auto _ : state) {
		ecInd = ec_sign(sm2_k256, r, s, priv, msg);
	}
}
BENCHMARK(test_ECSign);

static void test_ECVerify(benchmark::State &state)
{
	private_key<4>	priv(sm2_k256);
	bignum<4>	msg = bn_random<4>::Instance().get_random();
	bignum<4>	r, s;
	int		ecInd __attribute__ ((unused));
	ecInd = ec_sign(sm2_k256, r, s, priv, msg);
	for (auto _ : state) {
		ec_verify(sm2_k256, r, s, priv.PubKey(), msg);
	}
}
BENCHMARK(test_ECVerify);

static void test_ECRecover(benchmark::State &state)
{
	private_key<4>	priv(sm2_k256);
	bignum<4>	msg = bn_random<4>::Instance().get_random();
	bignum<4>	r, s;
	int		ecInd __attribute__ ((unused));
	ecInd = ec_sign(sm2_k256, r, s, priv, msg);
	spoint_t<4>	Ps;
	for (auto _ : state) {
		ec_recover(sm2_k256, Ps, r, s, ecInd, msg);
	}
}
BENCHMARK(test_ECRecover);

int main(int argc, char ** argv) {
	benchmark::Initialize(&argc, argv);
	benchmark::RunSpecifiedBenchmarks();
}
