#include "vli_bn.hpp"
#include "curve_const.hpp"

constexpr u64 dx1[] = {0x456490a734454105ull, 0x434764fd47d6fddull,
				0x1a3fc4b41b76d16eull, 0x1ba56e21c7de3566ull};
constexpr u64 dy1[] = {0xeec161d863c875c1ull, 0x94455772c85a4af3ull,
				0x145072ff099a874eull, 0xb3c27563f4adad8eull};
constexpr u64 dx2[] = {0xca65051fa6e050ffull, 0xdfa90a9f302318d6ull,
				0xc9cbff54adc5ef1bull, 0x6a63d715c7faa226ull};
constexpr u64 dy2[] = {0x158ea8d9c0a5601cull, 0x706e149f68929555ull,
				0xc595111bdd1514dull, 0x52ca2e119967eeefull};
constexpr u64 dx3[] = {0xd7731e82c7bad690ull, 0x02df40e54cf5b276ull,
				0x5ff5c92c30e3444cull, 0xf201d844ba3e9888ull};
constexpr u64 dy3[] = {0x8599d7a8a10da37bull, 0x785a72ba5b4dd000ull,
				0x0461357c00065ca1ull, 0xca0f618bf18fd27bull};
constexpr u64 dz3[] ={0x85007478729b0ffaull, 0xdb74944f5ba5a8f9ull,
				0xaf8c3aa0924f1dadull, 0x4ebe68f4001c6cc0ull};
constexpr u64 x1x1[] = {0xc5e0b83c849c9973ull, 0xe5069289576049efull,
				0xd5ebe7339cb40052ull, 0xe15dcc87234dc74dull};
constexpr u64 y1y1[] = {0x544dcd8f0cdf4ec8ull, 0xc4870641ced55afcull,
				0xf6a42b3e5ca6afc8ull, 0x575fea61de0e329full};
constexpr u64 z1z1[] = {0xdd82c3b0c790eb83ull, 0x288aaee690b495e6ull,
				0x28a0e5fe13350e9dull, 0x6784eac8e95b5b1cull};
constexpr u64 xy1mod[] = {0x8818da276776c620ull, 0xb20a4220e5254038ull,
				0x8e5a907ec15b6583ull, 0x4b23c3d16b70bd7ull};
constexpr u64 xy2mod[] = {0x5f2ac54e9268d236ull, 0xef1910202dbf4528ull,
				0x772534c53fbe9f4bull, 0x75654a24b92dc172ull};
constexpr u64 x1_inv[] = {0xc8e6f24a12b6be53ull, 0xd38d81d6e013b90aull,
				0x9c1cb188f2a91d75ull, 0xb8936229cd327d3aull};
constexpr u64 x2_inv[] = {0xe258fda37352c178ull, 0x3e87218e3e2cfd07ull,
				0xfa88a9e1fabcd344ull, 0x9bab218bea302819ull};
bignum<4>	prime(sm2_p);
bignum<4>	rr(sm2_p_rr);
bignum<4>	bigOne(1);
