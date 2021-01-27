#include <stdlib.h>
#include <string>
#include <chrono>
#include <gtest/gtest.h>
#include "vli_bn.hpp"
#include "curve_const.hpp"
#include "curve_defs.hpp"
#include "mont.hpp"


using namespace vli;
using namespace ecc;
using namespace std::chrono;

#include "testData.hpp"
//#define	sm2_p256	(*sm2_p256p)

TEST(testEcc, TestCalcK0RR)
{
	bignum<4>	trr;
	ASSERT_EQ(calcK0<4>(prime), sm2_p_k0);
	ASSERT_TRUE(rr.cmp(calcRR<4>(trr, prime)) == 0);
	bignum<4>	primeN(sm2_n);
	bignum<4>	rrN(sm2_n_rr);
	ASSERT_EQ(calcK0<4>(primeN), sm2_n_k0);
	ASSERT_TRUE(rrN.cmp(calcRR<4>(trr, primeN)) == 0);
}

static void mont_mul(bignum<4>& res, const bignum<4>& x, const bignum<4>& y)
{
	bignum<4>	xp, yp;
	xp.mont_mult(x, rr, prime, sm2_p_k0);
	yp.mont_mult(y, rr, prime, sm2_p_k0);
	res.mont_mult(xp, yp, prime, sm2_p_k0);
	res.mont_reduction(res, prime, sm2_p_k0);
}

static void mont_mulp(bignum<4>& res, const bignum<4>& x, const bignum<4>& y)
{
	bignum<4>	xp, yp;
	mont_mult<sm2_p_k0>(xp, x, rr, prime);
	mont_mult<sm2_p_k0>(yp, y, rr, prime);
	mont_mult<sm2_p_k0>(res, xp, yp, prime);
	mont_reduction<sm2_p_k0>(res, res, prime);
}

TEST(testEcc, TestMontMult)
{
	bignum<4>	res;
	bignum<4>	bx1(dx1), by1(dy1);
	mont_mul(res, bx1, by1);
	EXPECT_TRUE(res.cmp(xy1mod) == 0);
	bignum<4>	bx2(dx2), by2(dy2);
	mont_mul(res, bx2, by2);
	EXPECT_TRUE(res.cmp(xy2mod) == 0);
}

TEST(testEcc, TestMontMultP)
{
	bignum<4>	res;
	bignum<4>	bx1(dx1), by1(dy1);
	mont_mulp(res, bx1, by1);
	EXPECT_TRUE(res.cmp(xy1mod) == 0);
	bignum<4>	bx2(dx2), by2(dy2);
	mont_mulp(res, bx2, by2);
	EXPECT_TRUE(res.cmp(xy2mod) == 0);
}

TEST(TestEcc, TestSM2MultR)
{
	u64		ur[4];
	vli_sm2_multR(ur, 1);
	bignum<4>	res;
	res.clear();
	res.sub_from(prime);
	EXPECT_EQ(res.cmp(ur), 0);
	vli_sm2_multR(ur, 8);
	EXPECT_EQ(res.lshift(3), 0);
	EXPECT_EQ(res.cmp(ur), 0);
	vli_sm2_multR(ur, 64);
	EXPECT_EQ(res.lshift(3), 0);
	EXPECT_EQ(res.cmp(ur), 0);
	u64		u=dx1[0] & 0xffffff;
	vli_sm2_multR(ur, u);
	res.clear();
	res.sub_from(prime);
	bn_prod<4>	prod;
	prod.umult(res, u);
	EXPECT_TRUE(vli_is_zero<4>(prod.data()+4));
	EXPECT_EQ(vli_cmp<4>(prod.data(), ur), 0);
}


TEST(testVli, TestBignumz)
{
	bignumz<4>	p2(2), p3(3);
	bignumz<4>	n2(-2), n3(-3);
	bignumz<4>	bn_z(0L);
	ASSERT_TRUE(p2 < p3);
	ASSERT_TRUE(n3 < n2);
	ASSERT_TRUE(n2 < p2);
	ASSERT_FALSE(p2 < bn_z);
	ASSERT_TRUE(n2 < bn_z);
	bignumz<4>	res;
	res.sub(p2, p3);
	EXPECT_TRUE(res < bn_z);
	res.sub(n2, n3);
	EXPECT_TRUE(bn_z < res);
	res.add(n2, p3);
	EXPECT_TRUE(bn_z < res);
	res.add(n2, p2);
	EXPECT_TRUE(res == bn_z);
	res.sub(n2, p2);
	ASSERT_TRUE(res < n3);
	res.rshift1();
	ASSERT_TRUE(res == n2);
}

TEST(testvli, TestIsZero)
{
	ASSERT_EQ(u64IsZero(0), 1);
	ASSERT_EQ(u64IsZero(1), 0);
	EXPECT_EQ(u64IsZero(176453), 0);
}

TEST(testvli, TestIsOne)
{
	ASSERT_EQ(u64IsOne(0), 0);
	ASSERT_EQ(u64IsOne(1), 1);
	EXPECT_EQ(u64IsOne(176453), 0);
	EXPECT_EQ(u64IsOne(-1), 0);
}

TEST(testvli, TestNumBits)
{
	bignum<4>	p8(8);
	ASSERT_EQ(p8.num_bits(), 4);
	ASSERT_EQ(prime.num_bits(), 256);
	EXPECT_EQ(bigOne.num_bits(), 1);
	EXPECT_EQ(rr.num_bits(), 227);
	bignum<4> d1(d1d);
	EXPECT_EQ(d1.num_bits(), 255);
}

TEST(testvli, TestGetBits)
{
	uint	res;
	res = vli_get_bits<4,4>(d1d, 0);
	EXPECT_EQ(res, 2);
	res = vli_get_bits<4,4>(d1d, 20);
	EXPECT_EQ(res, 0xd);
	EXPECT_EQ((vli_get_bits<4, 2>(d1d, 20)), 1);
	EXPECT_EQ((vli_get_bits<4, 4>(d1d, 64)), 0xd);
	EXPECT_EQ((vli_get_bits<4, 4>(d1d, 124)), 0xa);
	EXPECT_EQ((vli_get_bits<4, 2>(d1d, 124)), 2);
	EXPECT_EQ((vli_get_bits<4, 4>(d1d, 188)), 0xe);
}

TEST(testVli, TestInverse)
{
	u64	res[4];
	vli_mod_inv<4>(res, dx1, sm2_p);
	EXPECT_EQ(vli_cmp<4>(res, x1_inv), 0);
	vli_mod_inv<4>(res, dx2, sm2_p);
	EXPECT_EQ(vli_cmp<4>(res, x2_inv), 0);
}

TEST(testVli, TestInverseNew)
{
	u64	res[4];
	vli_mod_inv_new<4>(res, dx1, sm2_p);
	EXPECT_EQ(vli_cmp<4>(res, x1_inv), 0);
	bignum<4>	xinv(res);
	std::cout << "inv_new x1_inv: " << xinv << std::endl;
	vli_mod_inv_new<4>(res, dx2, sm2_p);
	EXPECT_EQ(vli_cmp<4>(res, x2_inv), 0);
}


template<typename bnT, typename curveT>
forceinline static
void mont_sqr(const curveT& curve, bnT& res, const bnT &x1) noexcept
{
	bnT		tmp;
	curve.to_montgomery(tmp, x1);
	curve.mont_msqr(res, tmp);
	curve.from_montgomery(res, res);
}


TEST(testEcc, TestModSqrt)
{
	bignum<4>	res, sqrT;
	bignum<4>	tt;
	bignum<4>	by1(dy1), by2(dy2);
	ASSERT_EQ(sm2_p256.quadP(), quadPrime);
	mont_sqr(sm2_p256, res, by1);
	ASSERT_EQ(res.cmp(dy1y1), 0);
	sm2_p256.mod_exp(tt, res, quadPrime);
	ASSERT_EQ(tt.cmp(dy1Quad), 0);
	EXPECT_TRUE(sm2_p256.mod_sqrt(sqrT, res));
	EXPECT_TRUE(sqrT.cmp(by1) == 0);
	mont_sqr(sm2_p256, res, by2);
	sm2_p256.mod_exp(tt, res, quadPrime);
	ASSERT_EQ(tt.cmp(dy2Quad), 0);
	ASSERT_EQ(res.cmp(dy2y2), 0);
	EXPECT_TRUE(sm2_p256.mod_sqrt(sqrT, res));
	if (by2.is_even()) sqrT.sub(prime, sqrT);
	EXPECT_TRUE(sqrT.cmp(by2) == 0);
	ASSERT_EQ(sm2_p256.quadP(), quadPrime);
	mont_sqr(sm2_p256, res, by1);
	ASSERT_EQ(res.cmp(dy1y1), 0);
	sm2_k256.mod_exp(tt, res, quadPrime);
	ASSERT_EQ(tt.cmp(dy1Quad), 0);
	EXPECT_TRUE(sm2_k256.mod_sqrt(sqrT, res));
	EXPECT_TRUE(sqrT.cmp(by1) == 0);
	mont_sqr(sm2_p256, res, by2);
	sm2_k256.mod_exp(tt, res, quadPrime);
	ASSERT_EQ(tt.cmp(dy2Quad), 0);
	ASSERT_EQ(res.cmp(dy2y2), 0);
	EXPECT_TRUE(sm2_k256.mod_sqrt(sqrT, res));
	if (by2.is_even()) sqrT.sub(prime, sqrT);
	EXPECT_TRUE(sqrT.cmp(by2) == 0);
}

TEST(testEcc, TestMult248)
{
	bignum<4>	x3(dx3), y3(dy3), res1, res2, hh;
	sm2_p256.mont_mult2(res1, x3);
	hh = res1;
	sm2_p256.mont_div2(hh);
	EXPECT_EQ(hh, x3);
	sm2_p256.mont_mult2(res1, res1);
	res2 = x3;
	sm2_p256.mont_mult4(res2);
	EXPECT_EQ(res1.cmp(res2), 0);
	sm2_p256.mont_mult2(res1, res1);
	hh = res1;
	sm2_p256.mont_div2(hh);
	EXPECT_EQ(hh, res2);
	res2 = x3;
	sm2_p256.mont_mult8(res2);
	EXPECT_EQ(res1.cmp(res2), 0);
	sm2_p256.mont_mult2(res1, y3);
	sm2_p256.mont_mult2(res1, res1);
	res2 = y3;
	sm2_p256.mont_mult4(res2);
	EXPECT_EQ(res1.cmp(res2), 0);
	sm2_p256.mont_mult2(res1, res1);
	res2 = y3;
	sm2_p256.mont_mult8(res2);
	EXPECT_EQ(res1.cmp(res2), 0);
}

TEST(testEcc, TestOnCurve)
{
	bignum<4>	x1(d1Gx), x2(d2Gx);
	bignum<4>	y1(d1Gy), y2(d2Gy);
	bignum<4>	res;
	res.sub(sm2_p256.paramP(), sm2_p256.paramA());
	ASSERT_TRUE(res.is_u64());
	EXPECT_EQ(res.data()[0], 3);
	sm2_p256.mont_mult2(res, sm2_p256.mont_one());
	sm2_p256.mod_add_to(res, sm2_p256.mont_one());
	ASSERT_EQ(res, sm2_p256.mont_three());
	std::cerr << "mont_one: " << sm2_p256.mont_one() << std::endl;
	std::cerr << "mont_three: " << sm2_p256.mont_three() << std::endl;
	// res is 3 in montgomery form
	sm2_p256.mod_add_to(res, sm2_p256.montParamA());
	if (res.cmp(sm2_p256.paramP()) >= 0) res.sub_from(sm2_p256.paramP());
	EXPECT_TRUE(res.is_zero());
	std::cerr << "mont_A mont(-3): " << sm2_p256.montParamA() << std::endl;
	ASSERT_TRUE(sm2_p256.is_on_curve(x1, y1));
	ASSERT_TRUE(sm2_p256.is_on_curve(x2, y2));
}

TEST(testEcc, TestPointRecovery)
{
	bignum<4>	x1(d1Gx), x2(d2Gx);
	bignum<4>	y1(d1Gy), y2(d2Gy);
	bignum<4>	res;
	EXPECT_TRUE(pointY_recovery(sm2_p256, res, x1, !vli_is_even(d1Gy)));
	EXPECT_EQ(res, y1);
	EXPECT_TRUE(pointY_recovery(sm2_p256, res, x2, y2.is_odd()));
	EXPECT_EQ(res.cmp(d2Gy), 0);
}

TEST(testEcc, TestPointNeg)
{
	point_t<4>	pt1(dx1, dy1);
	point_t<4>	pt2(dx2, dy2);
	point_t<4>	res, q;
	sm2_p256.point_neg(q, pt1);
	sm2_p256.point_add(res, pt1, q);
	EXPECT_TRUE(res.is_zero());
	sm2_p256.point_neg(q, pt2);
	sm2_p256.point_add(res, pt2, q);
	EXPECT_TRUE(res.is_zero());
}

TEST(testEcc, TestPreCompute)
{
	point_t<4>	res, res1;
	point_t<4>	gg;
	sm2_p256.to_montgomery(gg.x, sm2_gx);
	sm2_p256.to_montgomery(gg.y, sm2_gy);
	gg.z = sm2_p256.mont_one();
	point_t<4>	pre_comps[wSize];
	steady_clock::time_point t1 = steady_clock::now();
	pre_compute<4>(sm2_p256, pre_comps, gg);
	steady_clock::time_point t2 = steady_clock::now();
	std::chrono::duration<double> time_span1;
	time_span1 = (t2 - t1);
	std::cerr << "pre_compute<4>() cost " << time_span1.count() * 1e6
			<< " us" << std::endl;
	res = gg;
	EXPECT_TRUE(pre_comps[0] == res);
	sm2_p256.point_double(res, res);
	EXPECT_TRUE(pre_comps[1] == res);
	sm2_p256.point_add(res1, res, gg.x, gg.y);
	EXPECT_TRUE(sm2_p256.point_eq(pre_comps[2],res1));
	sm2_p256.point_double(res, res);
	EXPECT_TRUE(sm2_p256.point_eq(pre_comps[3],res));
	sm2_p256.point_double(res, res);
	EXPECT_TRUE(pre_comps[7] == res);
	if (wSize < 16) return;
	sm2_p256.point_add(res1, res, gg.x, gg.y);
	EXPECT_TRUE(pre_comps[8] == res1);
	sm2_p256.point_double(res, res);
	EXPECT_TRUE(pre_comps[15] == res);
}

TEST(testEcc, TestECADD)
{
	u64		xx3[4], yy3[4], zz3[4];
	bignum<4>	x3(dx3), y3(dy3), z3(dz3);
	sm2_p256.point_add_jacobian(xx3, yy3, zz3, dx1, dy1,
						bigOne.data(), dx2, dy2);
	//vli_mod_inv<4>(z, z3, sm2_p256.paramP().data());
	//sm2_p256.apply_z(xx3, yy3, z);
	EXPECT_EQ(x3.cmp(xx3), 0);
	EXPECT_EQ(y3.cmp(yy3), 0);
	EXPECT_EQ(z3.cmp(zz3), 0);
}

TEST(testEcc, TestECDBL)
{
	u64		xx3[4], yy3[4], zz3[4];
	bignum<4>	x3(x1x1), y3(y1y1), z3(z1z1);
	sm2_p256.point_double_jacobian(xx3, yy3, zz3, dx1, dy1);
	//vli_mod_inv<4>(z, z3, sm2_p256.paramP().data());
	//sm2_p256.apply_z(xx3, yy3, z);
	EXPECT_EQ(x3.cmp(xx3), 0);
	EXPECT_EQ(y3.cmp(yy3), 0);
	EXPECT_EQ(z3.cmp(zz3), 0);
}

TEST(testEcc, TestECADDk256)
{
	u64		xx3[4], yy3[4], zz3[4];
	bignum<4>	x3(dx3), y3(dy3), z3(dz3);
	sm2_k256.point_add_jacobian(xx3, yy3, zz3, dx1, dy1,
						bigOne.data(), dx2, dy2);
	//vli_mod_inv<4>(z, z3, sm2_p256.paramP().data());
	//sm2_p256.apply_z(xx3, yy3, z);
	EXPECT_EQ(x3.cmp(xx3), 0);
	EXPECT_EQ(y3.cmp(yy3), 0);
	EXPECT_EQ(z3.cmp(zz3), 0);
}

TEST(testEcc, TestECDBLk256)
{
	u64		xx3[4], yy3[4], zz3[4];
	bignum<4>	x3(x1x1), y3(y1y1), z3(z1z1);
	sm2_k256.point_double_jacobian(xx3, yy3, zz3, dx1, dy1);
	//vli_mod_inv<4>(z, z3, sm2_p256.paramP().data());
	//sm2_p256.apply_z(xx3, yy3, z);
	EXPECT_EQ(x3.cmp(xx3), 0);
	EXPECT_EQ(y3.cmp(yy3), 0);
	EXPECT_EQ(z3.cmp(zz3), 0);
}

TEST(testEcc, TestScalarMult)
{
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	point_t<4>	gg(sm2_gx, sm2_gy);
	sm2_p256.scalar_mult(res, gg, d1);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d1Gx), 0);
	EXPECT_EQ(res.y.cmp(d1Gy), 0);
#ifdef	ommit
	std::cout << "res x1: " << res.x << std::endl;
	std::cout << "res y1: " << res.y << std::endl;
#endif
	sm2_p256.scalar_mult(res, gg, d2);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d2Gx), 0);
	EXPECT_EQ(res.y.cmp(d2Gy), 0);
#ifdef	ommit
	std::cout << "res x2: " << res.x << std::endl;
	std::cout << "res y2: " << res.y << std::endl;
#endif
}

TEST(testEcc, TestBaseNAF)
{
	ASSERT_EQ(BaseW, 6);
	EXPECT_EQ(nwBaseNAF<4>(), 43);
}

TEST(testEcc, TestPreComputeBaseN)
{
	point_t<4>	res, res1;
	point_t<4>	gg;
	sm2_p256.to_montgomery(gg.x, sm2_gx);
	sm2_p256.to_montgomery(gg.y, sm2_gy);
	gg.z = sm2_p256.mont_one();
	auto& pre_comps = sm2_p256.select_base_NAF();
	res = gg;
	EXPECT_TRUE(pre_comps[0] == res);
	sm2_p256.point_double(res, res);
	sm2_p256.apply_z_mont(res);
	EXPECT_TRUE(pre_comps[1] == res);
	sm2_p256.point_add(res1, res, gg.x, gg.y);
	sm2_p256.apply_z_mont(res1);
	EXPECT_TRUE(pre_comps[2] == res1);
	sm2_p256.point_double(res, res);
	sm2_p256.apply_z_mont(res);
	EXPECT_TRUE(pre_comps[3] == res);
	sm2_p256.point_double(res, res);
	sm2_p256.apply_z_mont(res);
	EXPECT_TRUE(pre_comps[7] == res);
	sm2_p256.point_add(res1, res, gg.x, gg.y);
	sm2_p256.apply_z_mont(res1);
	EXPECT_TRUE(pre_comps[8] == res1);
	sm2_p256.point_double(res, res);
	sm2_p256.apply_z_mont(res);
	EXPECT_TRUE(pre_comps[15] == res);
	sm2_p256.point_add(res1, res, res1);
	sm2_p256.apply_z_mont(res1);
	EXPECT_TRUE(sm2_p256.point_eq(pre_comps[24],res1));
	auto& pre_comps1 = sm2_p256.select_base_NAF(1);
	point_t<4>	resN, resN1;
	sm2_p256.point_double(resN, res.x, res.y);
	sm2_p256.point_double(resN1, res1.x, res1.y);
	for (int i=1; i < 6; ++i) {
		sm2_p256.point_double(resN, resN);
		sm2_p256.point_double(resN1, resN1);
	}
	sm2_p256.apply_z_mont(resN);
	EXPECT_TRUE(pre_comps1[15] == resN);
	sm2_p256.apply_z_mont(resN1);
	EXPECT_TRUE(sm2_p256.point_eq(pre_comps1[24],resN1));
	spoint_t<4>	rr;
	ASSERT_TRUE(sm2_p256.select_base_point(rr, 24, 1));
	EXPECT_TRUE(rr == resN1);
}

TEST(testEcc, TestScalarMultBaseNAF)
{
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	// initTable moved int protected
#ifdef	ommit
	steady_clock::time_point t1 = steady_clock::now();
	bool iniRes = sm2_p256.initTable();
	steady_clock::time_point t2 = steady_clock::now();
	std::chrono::duration<double> time_span1;
	time_span1 = (t2 - t1);
	std::cerr << "initTable() cost " << time_span1.count() * 1e6
			<< " us" << std::endl;
	ASSERT_TRUE(iniRes);
#endif
	sm2_p256.scalar_mult_base(res, d1);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d1Gx), 0);
	EXPECT_EQ(res.y.cmp(d1Gy), 0);
#ifdef	ommit
	std::cout << "res x1: " << res.x << std::endl;
	std::cout << "res y1: " << res.y << std::endl;
#endif
	sm2_p256.scalar_mult_base(res, d2);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d2Gx), 0);
	EXPECT_EQ(res.y.cmp(d2Gy), 0);
#ifdef	ommit
	std::cout << "res x2: " << res.x << std::endl;
	std::cout << "res y2: " << res.y << std::endl;
#endif
}

TEST(testEcc, TestScalar256Mult)
{
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	bignum<4>	bn_zero(0ul);
	point_t<4>	gg(sm2_gx, sm2_gy);
	sm2_k256.combined_mult(res, gg, d1, bn_zero);
	EXPECT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d1Gx), 0);
	EXPECT_EQ(res.y.cmp(d1Gy), 0);
#ifdef	ommit
	std::cout << "res x1: " << res.x << std::endl;
	std::cout << "res y1: " << res.y << std::endl;
#endif
	sm2_k256.combined_mult(res, gg, d2, bn_zero);
	EXPECT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d2Gx), 0);
	EXPECT_EQ(res.y.cmp(d2Gy), 0);
#ifdef	ommit
	std::cout << "res x2: " << res.x << std::endl;
	std::cout << "res y2: " << res.y << std::endl;
#endif
	sm2_k256.combined_mult(res, gg, bn_zero, d1);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d1Gx), 0);
	EXPECT_EQ(res.y.cmp(d1Gy), 0);
#ifdef	ommit
	std::cout << "res x1: " << res.x << std::endl;
	std::cout << "res y1: " << res.y << std::endl;
#endif
	sm2_k256.combined_mult(res, gg, bn_zero, d2);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d2Gx), 0);
	EXPECT_EQ(res.y.cmp(d2Gy), 0);
#ifdef	ommit
	std::cout << "res x2: " << res.x << std::endl;
	std::cout << "res y2: " << res.y << std::endl;
#endif
}

TEST(testEcc, TestScalarMultNAF2)
{
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	point_t<4>	gg(sm2_gx, sm2_gy);
	point_t<4>	pre_comps[wSize];
	sm2_p256.to_montgomery(res.x, gg.x);
	sm2_p256.to_montgomery(res.y, gg.y);
	sm2_p256.to_montgomery(res.z, gg.z);
	pre_compute<4>(sm2_p256, pre_comps, res);
	for (int i=0; i<wSize; ++i) {
		sm2_p256.from_montgomery(pre_comps[i].x, pre_comps[i].x);
		sm2_p256.from_montgomery(pre_comps[i].y, pre_comps[i].y);
		sm2_p256.from_montgomery(pre_comps[i].z, pre_comps[i].z);
	}
	{
		bignum<4>	ss(3);
		sm2_p256.scalar_multNAF2(res, gg, ss);
		EXPECT_TRUE(sm2_p256.point_eq(res, pre_comps[2]));
	}
	{
		bignum<4>	ss(4);
		sm2_p256.scalar_multNAF2(res, gg, ss);
		EXPECT_TRUE(sm2_p256.point_eq(res, pre_comps[3]));
	}
	{
		bignum<4>	ss(8);
		sm2_p256.scalar_multNAF2(res, gg, ss);
		EXPECT_TRUE(sm2_p256.point_eq(res, pre_comps[7]));
	}
	for(uint i=0; i<wSize; ++i) {
		bignum<4>	ss(i+1);
		sm2_p256.scalar_multNAF2(res, gg, ss);
		auto bCheck = sm2_p256.point_eq(res, pre_comps[i]);
		if (!bCheck) std::cerr << "diff multNAF2 index: " << i << std::endl;
		EXPECT_TRUE(sm2_p256.point_eq(res, pre_comps[i]));
	}
	sm2_p256.scalar_multNAF2(res, gg, d1);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d1Gx), 0);
	EXPECT_EQ(res.y.cmp(d1Gy), 0);
#ifdef	ommit
	std::cout << "res x1: " << res.x << std::endl;
	std::cout << "res y1: " << res.y << std::endl;
#endif
	sm2_p256.scalar_multNAF2(res, gg, d2);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d2Gx), 0);
	EXPECT_EQ(res.y.cmp(d2Gy), 0);
#ifdef	ommit
	std::cout << "res x2: " << res.x << std::endl;
	std::cout << "res y2: " << res.y << std::endl;
#endif
}

TEST(testEcc, TestScalarMultBase)
{
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	sm2_k256.scalar_mult_base(res, d1);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d1Gx), 0);
	EXPECT_EQ(res.y.cmp(d1Gy), 0);
#ifdef	ommit
	std::cout << "res x1: " << res.x << std::endl;
	std::cout << "res y1: " << res.y << std::endl;
#endif
	sm2_k256.scalar_mult_base(res, d2);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(d2Gx), 0);
	EXPECT_EQ(res.y.cmp(d2Gy), 0);
#ifdef	ommit
	std::cout << "res x1: " << res.x << std::endl;
	std::cout << "res y1: " << res.y << std::endl;
#endif
}

TEST(testEcc, TestScalarCombinedMult)
{
	point_t<4>	pt1(dx1, dy1);
	bignum<4>	d1(d1d), d2(d2d);
	u64		x3[4], y3[4], z3[4];
	point_t<4>	res, res1, res2;
	sm2_k256.combined_mult(res, pt1, d1, d2);
	ASSERT_TRUE(res.z.is_one());

	sm2_p256.scalar_mult(res1, pt1, d1);
	ASSERT_TRUE(res1.z.is_one());
	sm2_k256.scalar_mult_base(res2, d2);
	ASSERT_TRUE(res2.z.is_one());
	sm2_p256.point_add_jacobian(x3, y3, z3, res1.x.data(), res1.y.data(),
						bigOne.data(), res2.x.data(), res2.y.data());
	point_t<4>	pt3(x3, y3, z3);
	sm2_p256.apply_z(pt3);
	EXPECT_EQ(res, pt3);
	std::cout << "res x3: " << res.x << std::endl;
	std::cout << "res y3: " << res.y << std::endl;
	bignum<4>	xx3(x3);
	bignum<4>	yy3(y3);
	std::cout << "res xx3: " << xx3 << std::endl;
	std::cout << "res yy3: " << yy3 << std::endl;
}

TEST(testEcc, TestScalarCombinedMultN)
{
	point_t<4>	pt1(dx1, dy1);
	bignum<4>	d1(d1d), d2(d2d);
	u64		x3[4], y3[4], z3[4];
	point_t<4>	res, res1, res2;
	sm2_p256.combined_mult(res, pt1, d1, d2);
	ASSERT_TRUE(res.z.is_one());

	sm2_p256.scalar_mult(res1, pt1, d1);
	ASSERT_TRUE(res1.z.is_one());
	sm2_p256.scalar_mult_base(res2, d2);
	ASSERT_TRUE(res2.z.is_one());
	sm2_p256.point_add_jacobian(x3, y3, z3, res1.x.data(), res1.y.data(),
						res1.z.data(), res2.x.data(), res2.y.data());
	point_t<4>	pt3(x3, y3, z3);
	sm2_p256.apply_z(pt3);
	EXPECT_EQ(res, pt3);
	std::cout << "res x3: " << res.x << std::endl;
	std::cout << "res y3: " << res.y << std::endl;
	bignum<4>	xx3(x3);
	bignum<4>	yy3(y3);
	std::cout << "res xx3: " << xx3 << std::endl;
	std::cout << "res yy3: " << yy3 << std::endl;
}

int main(int argc, char *argv[])
{
	//sm2_p256.init();
    testing::InitGoogleTest(&argc, argv);//将命令行参数传递给gtest
    return RUN_ALL_TESTS();   //RUN_ALL_TESTS()运行所有测试案例
}
