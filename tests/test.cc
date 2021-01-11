#include <stdlib.h>
#include <string>
#include <gtest/gtest.h>
#include "vli_bn.hpp"
#include "curve_const.hpp"
#include "curve_defs.hpp"
#include "mont.hpp"


using namespace vli;

#include "testData.hpp"

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

TEST(testEcc, TestMult248)
{
	bignum<4>	x3(dx3), y3(dy3), res1, res2;
	sm2_p256.mont_mult2(res1, x3);
	sm2_p256.mont_mult2(res1, res1);
	res2 = x3;
	sm2_p256.mont_mult4(res2);
	EXPECT_EQ(res1.cmp(res2), 0);
	sm2_p256.mont_mult2(res1, res1);
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


static bool select_gmul(point_t<4>& pt, const uint idx)
{
	if (idx >= 2 * 16) return false;
	pt.x = bignum<4>(gmul[idx>>4][idx&0xf][0]);
	pt.y = bignum<4>(gmul[idx>>4][idx&0xf][1]);
	pt.z = bignum<4>(gmul[idx>>4][idx&0xf][2]);
	return true;
}

TEST(testEcc, TestPreCompute)
{
	sm2_k256.pre_compute_base();
	point_t<4>	gp, gm;
	for(int i=0;i<32;i++) {
		EXPECT_TRUE(sm2_k256.select_base_point(gp, i));
		EXPECT_TRUE(select_gmul(gm, i));
		ASSERT_TRUE(gp == gm);
	}
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

TEST(testEcc, TestScalartBase)
{
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	point_t<4>	gg(sm2_gx, sm2_gy, bigOne.data());
	sm2_p256.scalar_mul(res, gg, d1);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(dx1), 0);
	EXPECT_EQ(res.y.cmp(dy1), 0);
	sm2_p256.scalar_mul(res, gg, d2);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(dx2), 0);
	EXPECT_EQ(res.y.cmp(dy2), 0);
}

// not work yet
#ifdef	ommit
TEST(testEcc, TestScalarMultBase)
{
	bignum<4>	d1(d1d), d2(d2d);
	point_t<4>	res;
	sm2_k256.pre_compute_base();
	sm2_k256.scalar_mult_base(res, d1);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(dx1), 0);
	EXPECT_EQ(res.y.cmp(dy1), 0);
	sm2_k256.scalar_mult_base(res, d2);
	ASSERT_TRUE(res.z.is_one());
	EXPECT_EQ(res.x.cmp(dx2), 0);
	EXPECT_EQ(res.y.cmp(dy2), 0);
}
#endif

int main(int argc, char *argv[])
{
	sm2_p256.init();
	sm2_k256.init();
    testing::InitGoogleTest(&argc, argv);//将命令行参数传递给gtest
    return RUN_ALL_TESTS();   //RUN_ALL_TESTS()运行所有测试案例
}
