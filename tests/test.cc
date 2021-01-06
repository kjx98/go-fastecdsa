#include <stdlib.h>
#include <string>
#include <gtest/gtest.h>
#include "vli_bn.hpp"
#include "curve_const.hpp"
#include "curve_defs.hpp"
#include "mont.hpp"


using namespace vli;

#include "testData.h"

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

TEST(testEcc, TestMontMult)
{
	bignum<4>	res;
	mont_mul(res, bx1, by1);
	EXPECT_TRUE(res.cmp(xy1mod) == 0);
	mont_mul(res, bx2, by2);
	EXPECT_TRUE(res.cmp(xy2mod) == 0);
}

TEST(testCurve, TestECADD)
{
	u64		xx3[4], yy3[4], zz3[4];
	ASSERT_TRUE(sm2_p256.init());
	sm2_p256.point_add_jacobian(xx3, yy3, zz3, bx1.data(), by1.data(),
						bigOne.data(), bx2.data(), by2.data(), bigOne.data());
	//vli_mod_inv<4>(z, z3, sm2_p256.paramP().data());
	//sm2_p256.apply_z(xx3, yy3, z);
	EXPECT_EQ(x3.cmp(xx3), 0);
	EXPECT_EQ(y3.cmp(yy3), 0);
	EXPECT_EQ(z3.cmp(zz3), 0);
}

int main(int argc, char *argv[])
{
	initData();
	sm2_p256.init();
    testing::InitGoogleTest(&argc, argv);//将命令行参数传递给gtest
    return RUN_ALL_TESTS();   //RUN_ALL_TESTS()运行所有测试案例
}
