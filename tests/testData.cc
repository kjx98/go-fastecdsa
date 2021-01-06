#include "vli_bn.hpp"
#include "curve_const.hpp"
#include "mont.hpp"


using namespace vli;

bignum<4>	bx1({0x68e8ded5dbec39caull, 0x3c3c6c2f18c0e55bull,
		   		0xf7f9188477f6086aull, 0xafe524a88091d6daull});
bignum<4>	by1({0x6f18b74046a55994ull, 0xaf60d71184a4ca53ull,
				0xf176066831aad5c2ull, 0x8d837895a7dc179eull});
bignum<4>	bx2({0xca9ad478f1770dadull, 0x6c0a05bb2f2c9125ull,
				0x8baaa8adfa83e460ull, 0x684446fbee5167efull});
bignum<4>	by2({0xf28e513642bf309ull, 0x4f8d6bfa6a6846cull,
				0xf9965b00c6d87affull, 0xc7337843bdb886bfull});
bignum<4>	x3({0xd7731e82c7bad690ull, 0x02df40e54cf5b276ull,
				0x5ff5c92c30e3444cull, 0xf201d844ba3e9888ull});
bignum<4>	y3({0x8599d7a8a10da37bull, 0x785a72ba5b4dd000ull,
				0x0461357c00065ca1ull, 0xca0f618bf18fd27bull});
bignum<4>	z3({0x85007478729b0ffaull, 0xdb74944f5ba5a8f9ull,
				0xaf8c3aa0924f1dadull, 0x4ebe68f4001c6cc0ull});
bignum<4>	xy1mod({0x2c81b3f085aebe77ull, 0x3de394579a4dd4fbull,
				0xbceab95bcfab140dull, 0x6a23e4075ba931c2ull});
bignum<4>	xy2mod({0x8990cba9fcf0e413ull, 0x7ee32096736ae93eull,
				0xd274495d2719cfceull, 0x3467d7471311224bull});
bignum<4>	prime(sm2_p);
bignum<4>	rr(sm2_p_rr);
bignum<4>	bigOne(1);

void initData()
{
#ifdef	ommit
	const u64	*p = sm2_p;
	prime = bignum<4>(p);
	p = sm2_p_rr;
	rr = bignum<4>(p);
#endif
}
