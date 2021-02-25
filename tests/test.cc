#include <stdlib.h>
#include <string>
#include <chrono>
#include <gtest/gtest.h>
#include "vli_bn.hpp"
#include "curve_const.hpp"
#include "curve_defs.hpp"
#include "ecc_key.hpp"
#include "ecdsa_key.hpp"
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
	EXPECT_TRUE(rr.cmp(calcRR<4>(trr, prime)) == 0);
	bignum<4>	primeN(sm2_n);
	bignum<4>	rrN(sm2_n_rr);
	ASSERT_EQ(calcK0<4>(primeN), sm2_n_k0);
	EXPECT_TRUE(rrN.cmp(calcRR<4>(trr, primeN)) == 0);
}

TEST(testEcc, TestSM2umultP)
{
	u64		r[6];
	vli_sm2_multP(r, 1);
	ASSERT_EQ(prime.cmp(r), 0);
	u64	v3456[]={0xFFFFFFFFCCBBAA9Aull, 0xCCBBAA9A33445565ull,
		   			0xFFFFFFFFFFFFFFFFull, 0xCCBBAA99FFFFFFFFull, 0x33445565};
	bignum<5>	bv(v3456);
	vli_sm2_multP(r, 0x33445566);
	EXPECT_EQ(bv.cmp(r), 0);
}

TEST(testEcc, TestBTCumultP)
{
	u64		r[6];
	vli_btc_multP(r, 1);
	ASSERT_EQ(btc_prime.cmp(r), 0);
	u64	v3456[]={0xccbba9d6583615baull, 0xffffffffffffffffull,
		   			0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull, 0x33445565};
	bignum<5>	bv(v3456);
	vli_btc_multP(r, 0x33445566);
	bignum<5>	rv(r);
	//std::cerr << "secp256k1p * 33445566: " << rv << std::endl;
	EXPECT_EQ(bv.cmp(r), 0);
}

TEST(TestVli, TestUmult)
{
	u64	x[4]={0xffffffff00112233ull, 0xffffffffffffffffull, 0xaabbccdd88776655ull,
			0x1122334455667788ull};
	u64	u = -1;
	x[0] = -1;
	x[1] = -1;
	x[3] = -1;
	u64	r0[8], r1[6];
	vli_umult<4>(r0, x, u);
	vli_umult2<4>(r1, x, u);
	bignum<5>	re0(r0), re1(r1);
	ASSERT_EQ(re0, re1);
}

TEST(testVli, TestU64Addc)
{
	u64		x=12, y= 10, carry=1;
	ASSERT_EQ(u64_addc(x, y, carry), 23);
	ASSERT_EQ(carry, 0);
	y = -10;
	ASSERT_EQ(u64_addc(x, y, carry), 2);
	ASSERT_EQ(carry, 1);
	ASSERT_EQ(u64_addc(x, y, carry), 3);
	ASSERT_EQ(carry, 1);
	ASSERT_EQ(u64_addcz(x, carry), 13);
	ASSERT_EQ(carry, 0);
	ASSERT_EQ(u64_addcz(x, carry), 12);
	ASSERT_EQ(carry, 0);
	ASSERT_EQ(u64_addcz(y, carry), -10);
	ASSERT_EQ(carry, 0);
	carry = 1;
	y = -1;
	ASSERT_EQ(u64_addcz(y, carry), 0);
	ASSERT_EQ(carry, 1);
}

TEST(testVli, TestU64Subc)
{
	u64		x=12, y= 10, carry=1;
	ASSERT_EQ(u64_subc(x, y, carry), 1);
	ASSERT_EQ(carry, 0);
	y = -14;
	ASSERT_EQ(u64_subc(x, y, carry), 26);
	ASSERT_EQ(carry, 1);
	ASSERT_EQ(u64_subc(x, y, carry), 25);
	ASSERT_EQ(carry, 1);
	ASSERT_EQ(u64_subcz(x, carry), 11);
	ASSERT_EQ(carry, 0);
	ASSERT_EQ(u64_subcz(x, carry), 12);
	ASSERT_EQ(carry, 0);
	ASSERT_EQ(u64_subcz(y, carry), -14);
	ASSERT_EQ(carry, 0);
	carry = 1;
	y = 1;
	ASSERT_EQ(u64_subcz(y, carry), 0);
	ASSERT_EQ(carry, 0);
}


TEST(testVli, TestSquare)
{
	bignum<4>	bx1(dx1), by1(dy1);
	bignum<4>	bx2(dx2);
	bn_prod<4>	res, res1;
	res.squareOld(bx1);
	res1.square(bx1);
	EXPECT_EQ(vli_cmp<8>(res.data(), res1.data()), 0);
	ASSERT_EQ(res, res1);
	res.squareOld(by1);
	res1.square(by1);
	EXPECT_EQ(vli_cmp<8>(res.data(), res1.data()), 0);
	ASSERT_EQ(res, res1);
	res.squareOld(bx2);
	res1.square(bx2);
	EXPECT_EQ(res.cmp(res1), 0);
	ASSERT_EQ(res, res1);
	for(int i=0; i<10; ++i) {
		auto&  rd = bn_random<4>::Instance();
		bignum<4>	tmp = rd.get_random();
		res.squareOld(tmp);
		res1.square(tmp);
		ASSERT_EQ(res, res1);
	}
}

TEST(testVli, TestSquareN)
{
	bignum<4>	bx1(dx1), by1(dy1);
	bignum<4>	bx2(dx2);
	bn_prod<4>	res;
	u64			res1[8];
	res.squareOld(bx1);
	vli_square<4>(res1, dx1);
	EXPECT_EQ(vli_cmp<8>(res.data(), res1), 0);
	ASSERT_TRUE(res == res1);
	res.squareOld(by1);
	vli_square<4>(res1, dy1);
	EXPECT_EQ(vli_cmp<8>(res.data(), res1), 0);
	ASSERT_TRUE(res == res1);
	res.squareOld(bx2);
	vli_square<4>(res1, dx2);
	EXPECT_EQ(res.cmp(res1), 0);
	ASSERT_TRUE(res == res1);
}

static void mont_mul(bignum<4>& res, const bignum<4>& x, const bignum<4>& y)
{
	bignum<4>	xp, yp;
	xp.mont_mult(x, rr, prime, sm2_p_k0);
	yp.mont_mult(y, rr, prime, sm2_p_k0);
	res.mont_mult(xp, yp, prime, sm2_p_k0);
	res.mont_reduction(res, prime, sm2_p_k0);
}

static void mont_mulPr(bignum<4>& res, const bignum<4>& x, const bignum<4>& y)
{
	bignum<4>	xp, yp;
	xp.mont_mult(x, rr, prime, sm2_p_k0);
	yp.mont_mult(y, rr, prime, sm2_p_k0);
	bn_prod<4>	pd;
	pd.mult(xp, yp);
	res.mont_reduction(pd.m_low(), prime, 1);
	res.mod_add_to(pd.m_high(), prime);
	res.mont_reduction(res, prime, sm2_p_k0);
}

static void mont_mulp(bignum<4>& res, const bignum<4>& x, const bignum<4>& y)
{
	u64	xp[4], yp[4];
	u64   *resp = reinterpret_cast<u64 *>(&res);
	mont_mult<4,sm2_p_k0>(xp, x.data(), rr.data(), prime.data());
	mont_mult<4,sm2_p_k0>(yp, y.data(), rr.data(), prime.data());
	mont_mult<4,sm2_p_k0>(resp, xp, yp, prime.data());
	mont_reduction<4,sm2_p_k0>(resp, resp, prime.data());
}

static void mont_mulK01(bignum<4>& res, const bignum<4>& x, const bignum<4>& y)
{
	u64			xp[4], yp[4];
	sm2p_mult(xp, x.data(), rr.data());
	sm2p_mult(yp, y.data(), rr.data());
	sm2p_mult(xp, xp, yp);
	res.mont_reduction(xp, prime, 1);
}

static void mont_mulK01N(bignum<4>& res, const bignum<4>& x, const bignum<4>& y)
{
	u64			xp[4], yp[4];
	sm2p_multN(xp, x.data(), rr.data());
	sm2p_multN(yp, y.data(), rr.data());
	sm2p_multN(xp, xp, yp);
	res.mont_reduction(xp, prime, 1);
}

static void mont_sqrN(u64 *res, const bignum<4>& x)
{
	bignum<4>	xp;
	xp.mont_mult(x, rr, prime, sm2_p_k0);
	mont_sqrN<4,sm2_p_k0>(res, xp.data(), prime.data());
	mont_reduction<4,sm2_p_k0>(res, res, prime.data());
}

static void sm2p_mont_sqrN(u64 *res, const bignum<4>& x)
{
	bignum<4>	xp;
	xp.mont_mult(x, rr, prime, sm2_p_k0);
	sm2p_sqrN(res, xp.data());
	mont_reduction<4,sm2_p_k0>(res, res, prime.data());
}

TEST(TestEcc, TestMontMultBase)
{
	u64	x[4]={0xffffffff00112233ull, 0xffffffffffffffffull, 0xaabbccdd88776655ull,
			0x1122334455667788ull};
	u64		re0[4];
	bignum<4>	res;
	sm2p_mult(re0, x, x);
	bignum<4>	xx(x);
	res.mont_mult(xx, xx, prime, sm2_p_k0);
	EXPECT_TRUE(res == re0);
	sm2p_multN(re0, x, x);
	EXPECT_TRUE(res == re0);
	for(int i=0; i<10; ++i) {
		auto&  rd = bn_random<4>::Instance();
		bignum<4>	tmp = rd.get_random();
		bignum<4>	tmp1 = rd.get_random();
		sm2p_mult(re0, tmp.data(), tmp1.data());
		res.mont_mult(tmp, tmp1, prime, sm2_p_k0);
		bignum<4>	res0(re0);
		ASSERT_EQ(res, res0);
		sm2p_multN(re0, tmp.data(), tmp1.data());
		bignum<4>	res1(re0);
		EXPECT_EQ(res0, res1);
	}
}

TEST(testEcc, TestMontRedK01)
{
	auto&  rd = bn_random<4>::Instance();
	bignum<4>	xp, res;
	for (int i=0; i<10; ++i) {
		bignum<4>	tmp = rd.get_random();
		u64   *resp = reinterpret_cast<u64 *>(&xp);
		mont_mult<4,sm2_p_k0>(resp, tmp.data(), rr.data(), prime.data());
		res.mont_reduction(xp, prime, 1);
		EXPECT_EQ(res, tmp);
	}
}

TEST(testEcc, TestSM2pRed)
{
	auto&  rd = bn_random<4>::Instance();
	bignum<4>	xp;
	u64		res[4];
	for (int i=0; i<10; ++i) {
		bignum<4>	tmp = rd.get_random();
		u64   *resp = reinterpret_cast<u64 *>(&xp);
		mont_mult<4, sm2_p_k0>(resp, tmp.data(), rr.data(), prime.data());
		sm2p_reduction(res, xp.data());
		EXPECT_TRUE(tmp == res);
	}
}

TEST(testEcc, TestSM2pRedN)
{
	auto&  rd = bn_random<4>::Instance();
	bignum<4>	xp;
	u64		res[4];
	for (int i=0; i<10; ++i) {
		bignum<4>	tmp = rd.get_random();
		u64   *resp = reinterpret_cast<u64 *>(&xp);
		mont_mult<4, sm2_p_k0>(resp, tmp.data(), rr.data(), prime.data());
		sm2p_reductionN(res, xp.data());
		EXPECT_TRUE(tmp == res);
	}
}

TEST(testEcc, TestMontMult)
{
	bignum<4>	res, res2;
	bignum<4>	bx1(dx1), by1(dy1);
	mont_mul(res, bx1, by1);
	EXPECT_TRUE(res.cmp(xy1mod) == 0);
	mont_mulPr(res2, bx1, by1);
	EXPECT_EQ(res2, res);
	bignum<4>	bx2(dx2), by2(dy2);
	mont_mul(res, bx2, by2);
	EXPECT_TRUE(res.cmp(xy2mod) == 0);
	mont_mulPr(res2, bx2, by2);
	EXPECT_EQ(res2, res);
}

TEST(testEcc, TestMontSqr)
{
	bignum<4>	res;
	u64			res2[4];
	bignum<4>	bx1(dx1), by1(dy1);
	mont_mul(res, bx1, bx1);
	mont_sqrN(res2, bx1);
	EXPECT_TRUE(res.cmp(res2) == 0);
	sm2p_mont_sqrN(res2, bx1);
	EXPECT_TRUE(res.cmp(res2) == 0);
	mont_mul(res, by1, by1);
	mont_sqrN(res2, by1);
	EXPECT_TRUE(res.cmp(res2) == 0);
	sm2p_mont_sqrN(res2, by1);
	EXPECT_TRUE(res.cmp(res2) == 0);
}

TEST(testEcc, TestMontMultK01)
{
	bignum<4>	res;
	bignum<4>	bx1(dx1), by1(dy1);
	mont_mulK01(res, bx1, by1);
	EXPECT_TRUE(res.cmp(xy1mod) == 0);
	bignum<4>	bx2(dx2), by2(dy2);
	mont_mulK01(res, bx2, by2);
	EXPECT_TRUE(res.cmp(xy2mod) == 0);
}

TEST(testEcc, TestSM2pMultN)
{
	bignum<4>	res;
	bignum<4>	bx1(dx1), by1(dy1);
	mont_mulK01(res, bx1, by1);
	EXPECT_TRUE(res.cmp(xy1mod) == 0);
	bignum<4>	bx2(dx2), by2(dy2);
	mont_mulK01N(res, bx2, by2);
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

TEST(testEcc, TestBigNumRandom)
{
	auto&  rd = bn_random<4>::Instance();
	bignum<4>	tmp = rd.get_random();
	EXPECT_FALSE(tmp.is_zero());
	std::cerr << "bn_random: " << tmp << std::endl;
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

TEST(testEcc, TestPointRecover)
{
	bignum<4>	x1(d1Gx), x2(d2Gx);
	bignum<4>	y1(d1Gy), y2(d2Gy);
	bignum<4>	res;
	EXPECT_TRUE(pointY_recover(sm2_p256, res, x1, !vli_is_even(d1Gy)));
	EXPECT_EQ(res, y1);
	EXPECT_TRUE(pointY_recover(sm2_p256, res, x2, y2.is_odd()));
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
	bignum<4>	x3(dx3), y3(dy3), z3(dz3);
	point_t<4>	pt1, pt2, pt3;
	sm2_p256.to_affined(pt1, dx1, dy1);
	sm2_p256.to_affined(pt2, dx2, dy2);
	sm2_p256.point_add(pt3, pt1, pt2.x, pt2.y);
	sm2_p256.from_montgomery(pt3.x, pt3.x);
	sm2_p256.from_montgomery(pt3.y, pt3.y);
	sm2_p256.from_montgomery(pt3.z, pt3.z);
	EXPECT_EQ(x3.cmp(pt3.x), 0);
	EXPECT_EQ(y3.cmp(pt3.y), 0);
	EXPECT_EQ(z3.cmp(pt3.z), 0);
}

TEST(testEcc, TestECDBL)
{
	bignum<4>	x3(x1x1), y3(y1y1), z3(z1z1);
	point_t<4>	pt1, pt3;
	sm2_p256.to_affined(pt1, dx1, dy1);
	sm2_p256.point_double(pt3, pt1.x, pt1.y);
	sm2_p256.from_montgomery(pt3.x, pt3.x);
	sm2_p256.from_montgomery(pt3.y, pt3.y);
	sm2_p256.from_montgomery(pt3.z, pt3.z);
	EXPECT_EQ(x3.cmp(pt3.x), 0);
	EXPECT_EQ(y3.cmp(pt3.y), 0);
	EXPECT_EQ(z3.cmp(pt3.z), 0);
}

TEST(testEcc, TestECADDk256)
{
	bignum<4>	x3(dx3), y3(dy3), z3(dz3);
	point_t<4>	pt1, pt2, pt3;
	sm2_k256.to_affined(pt1, dx1, dy1);
	sm2_k256.to_affined(pt2, dx2, dy2);
	sm2_k256.point_add(pt3, pt1, pt2.x, pt2.y);
	sm2_p256.from_montgomery(pt3.x, pt3.x);
	sm2_p256.from_montgomery(pt3.y, pt3.y);
	sm2_p256.from_montgomery(pt3.z, pt3.z);
	EXPECT_EQ(x3.cmp(pt3.x), 0);
	EXPECT_EQ(y3.cmp(pt3.y), 0);
	EXPECT_EQ(z3.cmp(pt3.z), 0);
}

TEST(testEcc, TestECDBLk256)
{
	bignum<4>	x3(x1x1), y3(y1y1), z3(z1z1);
	point_t<4>	pt1, pt3;
	sm2_k256.to_affined(pt1, dx1, dy1);
	sm2_k256.point_double(pt3, pt1.x, pt1.y);
	sm2_p256.from_montgomery(pt3.x, pt3.x);
	sm2_p256.from_montgomery(pt3.y, pt3.y);
	sm2_p256.from_montgomery(pt3.z, pt3.z);
	EXPECT_EQ(x3.cmp(pt3.x), 0);
	EXPECT_EQ(y3.cmp(pt3.y), 0);
	EXPECT_EQ(z3.cmp(pt3.z), 0);
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
	sm2_p256.scalar_mult(res, gg, sm2_p256.paramN());
	EXPECT_TRUE(res.is_zero());
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
	sm2_k256.combined_mult(res, gg, bn_zero, sm2_p256.paramN());
	EXPECT_TRUE(res.is_zero());
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
	sm2_k256.scalar_mult_base(res, sm2_p256.paramN());
	EXPECT_TRUE(res.is_zero());
}

TEST(testEcc, TestScalarCombinedMult)
{
	point_t<4>	pt1(dx1, dy1);
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res, res1, res2;
	sm2_k256.combined_mult(res, pt1, d1, d2);
	ASSERT_TRUE(res.z.is_one());

	sm2_p256.scalar_mult(res1, pt1, d1);
	ASSERT_TRUE(res1.z.is_one());
	sm2_p256.scalar_mult_base(res2, d2);
	ASSERT_TRUE(res2.z.is_one());
	sm2_p256.to_montgomery(res1.x, res1.x);
	sm2_p256.to_montgomery(res1.y, res1.y);
	sm2_p256.to_montgomery(res2.x, res2.x);
	sm2_p256.to_montgomery(res2.y, res2.y);
	res1.z = sm2_p256.mont_one();
	point_t<4>	pt3;
	bignum<4>	xx3, yy3;
	sm2_p256.point_add(pt3, res1, res2.x, res2.y);
	sm2_p256.apply_z(xx3, yy3, pt3);
	EXPECT_EQ(res.x, xx3);
	EXPECT_EQ(res.y, yy3);
	std::cout << "res x3: " << res.x << std::endl;
	std::cout << "res y3: " << res.y << std::endl;
	std::cout << "res xx3: " << xx3 << std::endl;
	std::cout << "res yy3: " << yy3 << std::endl;
}

TEST(testEcc, TestScalarCombinedMultN)
{
	point_t<4>	pt1(dx1, dy1);
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res, res1, res2;
	sm2_p256.combined_mult(res, pt1, d1, d2);
	ASSERT_TRUE(res.z.is_one());

	sm2_p256.scalar_mult(res1, pt1, d1);
	ASSERT_TRUE(res1.z.is_one());
	sm2_p256.scalar_mult_base(res2, d2);
	ASSERT_TRUE(res2.z.is_one());
	sm2_p256.to_montgomery(res1.x, res1.x);
	sm2_p256.to_montgomery(res1.y, res1.y);
	sm2_p256.to_montgomery(res2.x, res2.x);
	sm2_p256.to_montgomery(res2.y, res2.y);
	res1.z = sm2_p256.mont_one();
	point_t<4>	pt3;
	bignum<4>	xx3, yy3;
	sm2_p256.point_add(pt3, res1, res2.x, res2.y);
	sm2_p256.apply_z(xx3, yy3, pt3);
	EXPECT_EQ(res.x, xx3);
	EXPECT_EQ(res.y, yy3);
	std::cout << "res x3: " << res.x << std::endl;
	std::cout << "res y3: " << res.y << std::endl;
	std::cout << "res xx3: " << xx3 << std::endl;
	std::cout << "res yy3: " << yy3 << std::endl;
}

TEST(TestECDSA, TestPrivateKey)
{
	private_key<4>	priv(sm2_p256);
	// operator bool() not support by gtest-1.6, need gtest-1.8 or above
	//ASSERT_TRUE(priv);
	if (! priv ) std::cerr << "private_key MUST be true" << std::endl;
	auto&	pk = priv.PubKey();
	ASSERT_TRUE(sm2_p256.is_on_curve(pk.x, pk.y));
	std::cerr << "PrivateKey: " << priv.D() << std::endl;
	std::cerr << "PubX: " << pk.x << std::endl;
	std::cerr << "PubY: " << pk.y << std::endl;
	bignum<4>	res;
	EXPECT_TRUE(pointY_recover(sm2_k256, res, pk.x, pk.y.is_odd()));
	EXPECT_EQ(res, pk.y);
	bignum<4>	px, py, dd;
	for (int i=0; i<10; ++i) {
		gen_keypair(sm2_k256, dd, px, py);
		ASSERT_TRUE(sm2_k256.is_on_curve(px, py));
		EXPECT_TRUE(pointY_recover(sm2_k256, res, px, py.is_odd()));
		EXPECT_EQ(res, py);
		std::cerr << "GenPrivKeyPair test round: " << i << "\tPointY odd: "
				<< (bool)py.is_odd() << std::endl;
	}
}

TEST(TestCurve, TestPresBuff)
{
	ASSERT_EQ(sizeof(point_t<4>), 4*8*3);
	EXPECT_EQ(sm2_k256.presBuffSize(), sizeof(point_t<4>) * ecc::wSize);
	std::cerr << "Sizeof wNAF precompute buff: "
			<< sizeof(point_t<4>) * ecc::wSize
			<< " bytes" << std::endl;
}

TEST(TestECDSA, TestSign)
{
	private_key<4>	priv(sm2_p256);
	bignum<4>	msg = bn_random<4>::Instance().get_random();
	bignum<4>	r, s;
	int		ecInd;
	for (int i=0; i<10; i++) {
		ecInd = ec_sign(sm2_p256, r, s, priv, msg);
		ASSERT_FALSE(r.is_zero());
		ASSERT_FALSE(s.is_zero());
		std::cerr << "sign test round: " << i << "\tsign ecInd: "
				<< ecInd << std::endl;
		bignum<4>	t;
		t.mod_add(r, s, sm2_p256.paramN());
		ASSERT_FALSE(t.is_zero());
		point_t<4>	q, p(priv.PubKey().x, priv.PubKey().y);
		sm2_p256.combined_mult(q, p, t, s);
		bignum<4>	tmp;
		tmp.mod_add(q.x, msg, sm2_p256.paramN());
		ASSERT_EQ(tmp, r);
	}
	steady_clock::time_point t1 = steady_clock::now();
	for (int i=0; i<1000; ++i)
		ecInd = ec_sign(sm2_k256, r, s, priv, msg);
	steady_clock::time_point t2 = steady_clock::now();
	std::chrono::duration<double> time_span1;
	time_span1 = (t2 - t1);
	std::cerr << "1000 sign cost " << time_span1.count() * 1e3
			<< " ms" << std::endl;
	std::cerr << "sign " << (int)(1000/time_span1.count()) << " msg/sec"
			<< std::endl;
}

TEST(TestECDSA, TestVerify)
{
	private_key<4>	priv(sm2_p256);
	bignum<4>	msg = bn_random<4>::Instance().get_random();
	bignum<4>	r, s;
	int		ecInd;
	ecInd = ec_sign(sm2_p256, r, s, priv, msg);
	ASSERT_FALSE(r.is_zero());
	ASSERT_FALSE(s.is_zero());
	std::cerr << "signed ret: " << ecInd << std::endl;
	bignum<4>	t;
	t.mod_add(r, s, sm2_p256.paramN());
	ASSERT_FALSE(t.is_zero());
	ASSERT_TRUE(ec_verify(sm2_p256, r, s, priv.PubKey(), msg));
}

TEST(TestECDSA, TestRecover)
{
	private_key<4>	priv(sm2_p256);
	bignum<4>	msg = bn_random<4>::Instance().get_random();
	bignum<4>	r, s;
	int		ecInd;
	ecInd = ec_sign(sm2_p256, r, s, priv, msg);
	ASSERT_FALSE(r.is_zero());
	ASSERT_FALSE(s.is_zero());
	std::cerr << "signed ret: " << ecInd << std::endl;
	bignum<4>	t;
	t.mod_add(r, s, sm2_p256.paramN());
	ASSERT_FALSE(t.is_zero());
	spoint_t<4>	Ps;
	ASSERT_TRUE(ec_recover(sm2_p256, Ps, r, s, ecInd, msg));
	EXPECT_EQ(Ps.x, priv.PubKey().x);
	EXPECT_EQ(Ps.y, priv.PubKey().y);
	steady_clock::time_point t1 = steady_clock::now();
	for (int i=0; i<1000; ++i)
		ec_recover(sm2_k256, Ps, r, s, ecInd, msg);
	steady_clock::time_point t2 = steady_clock::now();
	std::chrono::duration<double> time_span1;
	time_span1 = (t2 - t1);
	std::cerr << "1000 recover pubkey from signed msg cost "
			<< time_span1.count() * 1e3 << " ms" << std::endl;
	std::cerr << "recover " << (int)(1000/time_span1.count()) << " msg/sec"
			<< std::endl;
}

#ifdef	WITH_ECDSA
TEST(TestECDSA, TestEcdsaSign)
{
	private_key<4>	priv(sm2_p256);
	bignum<4>	msg = bn_random<4>::Instance().get_random();
	bignum<4>	r, s;
	int		ecInd;
	for (int i=0; i<10; i++) {
		ecInd = ecdsa_sign(sm2_k256, r, s, priv, msg);
		ASSERT_FALSE(r.is_zero());
		ASSERT_FALSE(s.is_zero());
		std::cerr << "ECDSA sign test round: " << i << "\tsign ecInd: "
				<< ecInd << std::endl;
		bignum<4>	u,v, sInv;
		point_t<4>	q, p(priv.PubKey().x, priv.PubKey().y);
		mod_inv<4>(sInv, s, sm2_k256.paramN());
		sm2_k256.to_montgomeryN(sInv, sInv);
		sm2_k256.to_montgomeryN(u, msg);
		sm2_k256.to_montgomeryN(v, r);
		sm2_k256.mont_nmult(u, u, sInv);
		sm2_k256.mont_nmult(v, v, sInv);
		sm2_k256.from_montgomeryN(u, u);
		sm2_k256.from_montgomeryN(v, v);
		EXPECT_FALSE(u.is_zero());
		EXPECT_FALSE(v.is_zero());
		sm2_k256.combined_mult(q, p, v, u);
		ASSERT_EQ(q.x, r);
	}
	steady_clock::time_point t1 = steady_clock::now();
	for (int i=0; i<1000; ++i)
		ecInd = ecdsa_sign(sm2_k256, r, s, priv, msg);
	steady_clock::time_point t2 = steady_clock::now();
	std::chrono::duration<double> time_span1;
	time_span1 = (t2 - t1);
	std::cerr << "1000 ECDSA sign cost " << time_span1.count() * 1e3
			<< " ms" << std::endl;
	std::cerr << "ECDSA sign " << (int)(1000/time_span1.count()) << " msg/sec"
			<< std::endl;
}

TEST(TestECDSA, TestEcdsaVerify)
{
	private_key<4>	priv(sm2_p256);
	bignum<4>	msg = bn_random<4>::Instance().get_random();
	bignum<4>	r, s;
	int		ecInd;
	ecInd = ecdsa_sign(sm2_p256, r, s, priv, msg);
	ASSERT_FALSE(r.is_zero());
	ASSERT_FALSE(s.is_zero());
	std::cerr << "ECDSA sign ret: " << ecInd << std::endl;
	ASSERT_TRUE(ecdsa_verify(sm2_p256, r, s, priv.PubKey(), msg));
}

TEST(TestECDSA, TestEcdsaRecover)
{
	private_key<4>	priv(sm2_p256);
	bignum<4>	msg = bn_random<4>::Instance().get_random();
	bignum<4>	r, s;
	int		ecInd;
	ecInd = ecdsa_sign(sm2_p256, r, s, priv, msg);
	ASSERT_FALSE(r.is_zero());
	ASSERT_FALSE(s.is_zero());
	std::cerr << "ECDSA sign ret: " << ecInd << std::endl;
	spoint_t<4>	Ps;
	ASSERT_TRUE(ecdsa_recover(sm2_p256, Ps, r, s, ecInd, msg));
#ifdef	WITH_HALF_N_ommit
	if (Ps.x != priv.PubKey().x) {
		s.sub(sm2_p256.paramN(), s);
		ASSERT_TRUE(ecdsa_recover(sm2_p256, Ps, r, s, ecInd, msg));
	}
#endif
	EXPECT_EQ(Ps.x, priv.PubKey().x);
	EXPECT_EQ(Ps.y, priv.PubKey().y);
	steady_clock::time_point t1 = steady_clock::now();
	for (int i=0; i<1000; ++i)
		ecdsa_recover(sm2_k256, Ps, r, s, ecInd, msg);
	steady_clock::time_point t2 = steady_clock::now();
	std::chrono::duration<double> time_span1;
	time_span1 = (t2 - t1);
	std::cerr << "1000 ECDSA recover pubkey from signed msg cost "
			<< time_span1.count() * 1e3 << " ms" << std::endl;
	std::cerr << "ECDSA recover " << (int)(1000/time_span1.count())
			<< " msg/sec" << std::endl;
}
#endif

int main(int argc, char *argv[])
{
	//sm2_p256.init();
    testing::InitGoogleTest(&argc, argv);//将命令行参数传递给gtest
#ifdef	__clang__
	std::cerr << "Built via clang " << __clang_major__ << "."
			<< __clang_minor__ << std::endl;
#elif	defined(__GNUC__)
	std::cerr << "Built via GCC " << __GNUC__ << "."
			<< __GNUC_MINOR__ << std::endl;
#endif
    return RUN_ALL_TESTS();   //RUN_ALL_TESTS()运行所有测试案例
}
