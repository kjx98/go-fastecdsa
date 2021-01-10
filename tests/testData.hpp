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

/*-
 * Base point pre computation
 * --------------------------
 *
 * Two different sorts of precomputed tables are used in the following code.
 * Each contain various points on the curve, where each point is three field
 * elements (x, y, z).
 *
 * For the base point table, z is usually 1 (0 for the point at infinity).
 * This table has 2 * 16 elements, starting with the following:
 * index | bits    | point
 * ------+---------+------------------------------
 *     0 | 0 0 0 0 | 0G
 *     1 | 0 0 0 1 | 1G
 *     2 | 0 0 1 0 | 2^64G
 *     3 | 0 0 1 1 | (2^64 + 1)G
 *     4 | 0 1 0 0 | 2^128G
 *     5 | 0 1 0 1 | (2^128 + 1)G
 *     6 | 0 1 1 0 | (2^128 + 2^64)G
 *     7 | 0 1 1 1 | (2^128 + 2^64 + 1)G
 *     8 | 1 0 0 0 | 2^192G
 *     9 | 1 0 0 1 | (2^192 + 1)G
 *    10 | 1 0 1 0 | (2^192 + 2^64)G
 *    11 | 1 0 1 1 | (2^192 + 2^64 + 1)G
 *    12 | 1 1 0 0 | (2^192 + 2^128)G
 *    13 | 1 1 0 1 | (2^192 + 2^128 + 1)G
 *    14 | 1 1 1 0 | (2^192 + 2^128 + 2^64)G
 *    15 | 1 1 1 1 | (2^192 + 2^128 + 2^64 + 1)G
 * followed by a copy of this with each element multiplied by 2^32.
 *
 * The reason for this is so that we can clock bits into four different
 * locations when doing simple scalar multiplies against the base point,
 * and then another four locations using the second 16 elements.
 *
 * Tables for other points have table[i] = iG for i in 0 .. 16. */

/* gmul is the table of precomputed base points */
constexpr bn_words_t gmul[2][16][3] = {
    {{{0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}},
     {{0x715a4589334c74c7, 0x8fe30bbff2660be1, 0x5f9904466a39c994,
       0x32c4ae2c1f198119},
      {0x2df32e52139f0a0, 0xd0a9877cc62a4740, 0x59bdcee36b692153,
       0xbc3736a2f4f6779c},
      {1, 0, 0, 0}},
     {{0xe18bd546b5824517, 0x673891d791caa486, 0xba220b99df9f9a14,
       0x95afbd1155c1da54},
      {0x8e4450eb334acdcb, 0xc3c7d1898a53f20d, 0x2eee750f4053017c,
       0xe8a6d82c517388c2},
      {1, 0, 0, 0}},
     {{0xf81c8da9b99fba55, 0x137f6c6149feef6e, 0xcb129aa494da9ad4,
       0x82a0f5407d123db6},
      {0xfdeca00772c4dbc9, 0xa961b58f0cf58373, 0xecacab94e973f9c3,
       0xf12fa4696a22ca3f},
      {1, 0, 0, 0}},
     {{0xeae3d9a9d13a42ed, 0x2b2308f6484e1b38, 0x3db7b24888c21f3a,
       0xb692e5b574d55da9},
      {0xd186469de295e5ab, 0xdb61ac1773438e6d, 0x5a924f85544926f9,
       0xa175051b0f3fb613},
      {1, 0, 0, 0}},
     {{0xa72d084f62c8d58b, 0xe3d6467deaf48fd7, 0x8fe75e5a128a56a7,
       0xc0023fe7ff2b68bd},
      {0x64f67782316815f9, 0xb52b6d9b19a69cd2, 0x5d1ed6fa89cbbade,
       0x796c910ee7f4ccdb},
      {1, 0, 0, 0}},
     {{0x1b2150c1c5f13015, 0xdaaba91b5d952c9b, 0xe8cc24c3f546142,
       0x75a34b243705f260},
      {0x77d195421cef1339, 0x636644aa0c3a0623, 0x4683df176eeb2444,
       0x642ce3bd3535e74d},
      {1, 0, 0, 0}},
     {{0x4a59ac2c6e7ecc08, 0xaf2b71164f191d63, 0x3622a87fb284554f,
       0xd9eb397b441e9cd0},
      {0xa66b8a4893b6a54d, 0x26fb89a40b4a663a, 0xafa87501eedfc9f4,
       0xf3f000bc66f98108},
      {1, 0, 0, 0}},
     {{0xad8bc68ce031d616, 0x16888d8ee4003187, 0x44c0757f3bb8b600,
       0x793fae7af0164245},
      {0x210cd042973f333b, 0x8666ff52dbd25f9, 0x65c5b129f5f7ad5d,
       0xe03d7a8d19b3219a},
      {1, 0, 0, 0}},
     {{0xd68bfbace0e00392, 0x261014f7d3445dc7, 0xd9f46b2714a071ee,
       0x1b200af30810b682},
      {0xd91d8b12ae69bcd, 0x74a08f17bf8cd981, 0xd822913cf0d2b82d,
       0x248b7af0b05bfad2},
      {1, 0, 0, 0}},
     {{0xba119a049e62f2e2, 0xf278e8a34df05ae5, 0xd269f3564eb5d180,
       0x8e74ad0f4f957cb1},
      {0x112ff4dabd76e2dd, 0x91373f20630fdb7f, 0xf43eab474992904c,
       0x55a5ccc7af3b6db4},
      {1, 0, 0, 0}},
     {{0x5ad104a8bdd23de9, 0xf5a9e515eb71c2c1, 0x390542a0ba95c174,
       0x4c55fb20426491bf},
      {0x91525735ef626289, 0xd2ed977f88f09635, 0xfd48731b7a8a8521,
       0x8f89a03b8fdebea},
      {1, 0, 0, 0}},
     {{0x7e8e61ea35eb8e2e, 0x1bb2700db98a762c, 0xd81ea23b7738c17c,
       0xf9def2a46dba26a3},
      {0x183a7912d05e329f, 0x34664a0896ccde0e, 0x56c22652614283bb,
       0x91692899d5ff0513},
      {1, 0, 0, 0}},
     {{0x449d48d8f3bdbe19, 0xab95de03cc8510cb, 0xaef159463f8bfb25,
       0xda72c379dae3ca8b},
      {0xcba9315ce82cc3ea, 0x4e524bac38a58020, 0x36ba2752538e348c,
       0xb170d0da75ed450f},
      {1, 0, 0, 0}},
     {{0x947af0f52b4f8da6, 0x7eda17d917827976, 0x5ba79a0c705853a0,
       0xa5d9873b3fb2ddc7},
      {0xc2a48162a5fd9ce9, 0x80ee8ae526f25f02, 0xf60c8ef6633be6a9,
       0xe2e23f0229a84a35},
      {1, 0, 0, 0}},
     {{0xbc4945bd86bb6afb, 0x237eb711eba46fee, 0x7c1db58b7b86eb33,
       0xd94eb728273b3ac7},
      {0xbe1717e59568d0a4, 0x4a6067cc45f70212, 0x19b32eb5afc2fb17,
       0xbe3c1e7ac3ac9d3c},
      {1, 0, 0, 0}}},
    {{{0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}},
     {{0x68a88405ae53c1e9, 0x51e46707fd558656, 0x71e834cf86896c10,
       0x3d251b54e10d581f},
      {0x1884d5b0eeb19032, 0xeeaf729853e526fe, 0x5931f6831a8d8c11,
       0x87891d33fb98b4d8},
      {1, 0, 0, 0}},
     {{0x9047673fcac14893, 0xf5df5d83bfb58659, 0xa6230c81642e71a,
       0xef14b33800777791},
      {0xcf1e99afa3386fca, 0x7ace937791313d53, 0x36fe159b6dcd01bb,
       0xc9bc50d02e2b960a},
      {1, 0, 0, 0}},
     {{0x716e5a7ee12e162d, 0xbbf9bb2c62dd5a00, 0xca235ccb4144dd05,
       0xbcb7de0f8f70520e},
      {0x981e8964947cb8eb, 0x53c7102ea04de08d, 0xe9076332afc6a10d,
       0x93d90f776b58c35d},
      {1, 0, 0, 0}},
     {{0x834dbff6678337ee, 0xc607e811fef0785a, 0xaaefc62be30a298b,
       0xeb5ca335326afad3},
      {0x9774fe1384af54a8, 0xca4b6ef5785388b4, 0x1346c82d66f6c642,
       0xedcc0c2aaa2d53ce},
      {1, 0, 0, 0}},
     {{0xb896b3f764b9e6f4, 0x47e4018c736fb3d0, 0xfc2fc86707413920,
       0x1a8526428e1aeae7},
      {0x1386802650e2ae60, 0x7474dedc995384d0, 0x2c4cc396dd43b011,
       0x63b0e9c7141de1b0},
      {1, 0, 0, 0}},
     {{0xeb5fb3b369d17771, 0x1fe07b18933ed257, 0xdfc4c81ce3673912,
       0x913614c66a91a647},
      {0x18aee853c0ba877f, 0x3109c2deceff091, 0x8532307e7e4ee08c,
       0xcef0791a6e6ce0bb},
      {1, 0, 0, 0}},
     {{0xf0e9f5d8057a4a0f, 0xbbf7f8b49f125aa9, 0x51e8fdd6283187c2,
       0xe0997d4759d36298},
      {0x67ec3c5c6f4221c3, 0x3ea275dbc860722f, 0x152d01e23859f5e2,
       0xfb57404312680f44},
      {1, 0, 0, 0}},
     {{0x21ac3df849be2a1f, 0x11006e9fc51d112f, 0x9151aa584775c857,
       0x5159d218ba04a8d9},
      {0x98b7d1a925fd1866, 0x8f4753cafc2ad9d8, 0x8eb91ec1569c05a9,
       0x4abbd1ae27e13f11},
      {1, 0, 0, 0}},
     {{0x616f6644b2c11f4c, 0x251cd7140e540758, 0xf927a40110f02017,
       0x92ff3cc3c1c941b6},
      {0x3249906213f565fe, 0x4633e3ddeb9dbd4e, 0xea9a9d1ec402e6c2,
       0xdc84ce34b14bb7cf},
      {1, 0, 0, 0}},
     {{0xa93e23e5436ff69a, 0x52dcb0a79b63efce, 0x34f6538a9e90cb41,
       0x9cac08f200234bc0},
      {0x6661825b5174a02d, 0x7d4d06de036be57, 0x589d74610ae6bd27,
       0xa296f5577fc91a93},
      {1, 0, 0, 0}},
     {{0x10acefa9d29721d0, 0x8b0f6b8bb5bcd340, 0x921d318c3d86785c,
       0xd6916f3bc16aa378},
      {0x2a0d646a7ad84a0e, 0x7b93256c2fe7e97a, 0x5765e27626479e41,
       0xae9da2272daaced3},
      {1, 0, 0, 0}},
     {{0x56fdc215f7f34ac5, 0xebcb4ff2da3877d3, 0x1eb96792aba6b832,
       0x807ce6bea24741aa},
      {0xff1c10109c721fb4, 0xd187d4bc796353a7, 0x7639ae749af2d303,
       0xaff6d783d56c9286},
      {1, 0, 0, 0}},
     {{0x6002d51b6290dd01, 0xcba3ab0099a836a5, 0x71776611e00d2528,
       0xfaf2cb8c87fce119},
      {0xd445228bdf6882ae, 0xcbbfade17cbce919, 0x837b6335a2eb2453,
       0x11ad7c4b8597f6b6},
      {1, 0, 0, 0}},
     {{0x48de8f368cf2e399, 0x7ae3d25630a74277, 0xdef1a9a6c505323f,
       0xe55f203b4b8d9672},
      {0xc58d8f0d9a1e6e97, 0xe160e6d4b2737a76, 0xd60bd087d47cbdd8,
       0x687d41364d5fef53},
      {1, 0, 0, 0}},
     {{0x83f21bbe056bbf9b, 0x4c2a9d120b4ba5ab, 0xff383d1845b64e4f,
       0x8f13cc8d06dd7867},
      {0xf3a292d8424f0995, 0xfd2546eae7cbe44b, 0x67d14dee6c1e75a3,
       0x53b49e6cc93fb5a8},
      {1, 0, 0, 0}}}
};

