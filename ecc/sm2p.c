// +build sm2p

#include "sm2p.h"

#include "sm2p256.inc"

static void bn_to_felem(felem out, const bn_words in)
{
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
    out[3] = in[3];
}

void sm2_mod_inv(bn_words *result, const bn_words *in)
{
    felem tmp;
	smallfelem *out=(smallfelem*)result;

    bn_to_felem(tmp, *in);
    felem_inv(tmp, tmp);
    felem_contract(*out, tmp);
}

void sm2_mod_mul(bn_words *result, const bn_words *x, const bn_words *y)
{
}

void sm2_point_add(Point *result, const Point *p, const Point *q)
{
    felem felem_x3, felem_y3, felem_z3;
    felem felem_x1, felem_y1, felem_z1;
    smallfelem_expand(felem_x1, p->x);
    smallfelem_expand(felem_y1, p->y);
    smallfelem_expand(felem_z1, p->z);
    point_add(felem_x3, felem_y3, felem_z3, felem_x1, felem_y1, felem_z1, 0,
              q->x, q->y, q->z);
    felem_shrink(result->x, felem_x3);
    felem_shrink(result->y, felem_y3);
    felem_shrink(result->z, felem_z3);
}

void sm2_point_double(Point *result, const Point *p)
{
    felem felem_x_out, felem_y_out, felem_z_out;
    felem felem_x_in, felem_y_in, felem_z_in;

    smallfelem_expand(felem_x_in, p->x);
    smallfelem_expand(felem_y_in, p->y);
    smallfelem_expand(felem_z_in, p->z);
    point_double(felem_x_out, felem_y_out, felem_z_out,
                 felem_x_in, felem_y_in, felem_z_in);
    felem_shrink(result->x, felem_x_out);
    felem_shrink(result->y, felem_y_out);
    felem_shrink(result->z, felem_z_out);
}
