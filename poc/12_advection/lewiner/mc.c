#include <float.h>
#include <math.h>
#include <memory.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "mc.h"
#include "table.h"

#define ALLOC_SIZE (10)
#define SIZE(a) (int)(sizeof(a)/sizeof(*(a)))

enum { X, Y, Z };
enum { OOO,
    IOO,
    IIO,
    OIO,
    OOI,
    IOI,
    III,
    OII
};
static int CuDir[] = { IOO, OIO, OOI };

static int CuOfset[][3] = {
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 1, 1, 0 },
    { 0, 1, 0 },
    { 0, 0, 1 },
    { 1, 0, 1 },
    { 1, 1, 1 },
    { 0, 1, 1 },
};

static int COff[][4] = {
    { 0, 0, 0, X },
    { 1, 0, 0, Y },
    { 0, 1, 0, X },
    { 0, 0, 0, Y },
    { 0, 0, 1, X },
    { 1, 0, 1, Y },
    { 0, 1, 1, X },
    { 0, 0, 1, Y },
    { 0, 0, 0, Z },
    { 1, 0, 0, Z },
    { 1, 1, 0, Z },
    { 0, 1, 0, Z },
};

static int add_c_vertex(void);
static void add_tri2(const int *trig, int n);
static void add_tri3(const int *trig, int n, int v12);
static int apply(void);
static double get_data(int i, int j, int k);
static int get_vert(int, int, int, int);
static void intersection(void);
static void process_cu(void);
static void realloc_tri(void);
static void realloc_ver(void);
static void set_vert(int, int, int, int);
static int test_face(int face);
static int test_interior(int s);
/* weighed average */
static double wavg(double a, double b, double t);
/* relative position */
static double rpos(double a, double b);

struct MarchLewiner {
    double cu[8];
    double *data;
    int Case;
    int config;
    int I;
    int J;
    int K;
    int N;
    int nt;
    int Nt;
    int nv;
    int Nv;
    int subconfig;
    int *verts;
    int x;
    int y;
    int z;
    int *tri;
    double *ver;
} *q;

struct MarchLewiner *
march_lewiner_ini(int x0, int y0, int z0)
{
    int i;
    int N;
    struct MarchLewiner *q;

    q = malloc(sizeof(*q));
    q->x = x0;
    q->y = y0;
    q->z = z0;
    N = x0 * y0 * z0;
    q->verts = malloc(3 * N * sizeof(*q->verts));
    for (i = 0; i < 3 * N; i++)
	q->verts[i] = -1;
    q->nv = q->nt = 0;
    q->Nv = q->Nt = ALLOC_SIZE;
    q->ver = malloc(3 * q->Nv * sizeof(*q->ver));
    q->tri = malloc(3 * q->Nt * sizeof(*q->tri));
    return q;
}

int
march_lewiner_fin(struct MarchLewiner *q)
{
    free(q->verts);
    free(q->ver);
    free(q->tri);
    free(q);
    return 0;
}

int
march_lewiner_apply(struct MarchLewiner *q0, double *data, int *nv,
		    double **ver, int *nt, int **tri)
{
    int status;

    q = q0;
    q->data = data;
    status = apply();

    *nv = q->nv;
    *ver = q->ver;
    *nt = q->nt;
    *tri = q->tri;
    return status;
}

static int
apply(void)
{
    int lut_entry, p, *o, i, j, k;
    int x, y, z;
    double *cu;

    x = q->x;
    y = q->y;
    z = q->z;
    cu = q->cu;

    intersection();
    for (k = 0; k < z - 1; k++)
	for (j = 0; j < y - 1; j++)
	    for (i = 0; i < x - 1; i++) {
		lut_entry = 0;
		for (p = 0; p < SIZE(CuOfset); ++p) {
		    o = CuOfset[p];
		    cu[p] = get_data(i + o[X], j + o[Y], k + o[Z]);
		    if (cu[p] > 0)
			lut_entry += 1 << p;
		}
		q->I = i;
		q->J = j;
		q->K = k;
		q->Case = cases[lut_entry][0];
		q->config = cases[lut_entry][1];
		process_cu();
	    }
    return 0;
}

static double
wavg(double a, double b, double t)
{
    return a + (b - a) * t;
}

static double
rpos(double a, double b)
{
    return a / (a - b);
}

static double
get_data(int i, int j, int k)
{
    double ans;
    int y, z;

    y = q->y;
    z = q->z;
    ans = q->data[i * y * z + j * y + k];
    return fabs(ans) < FLT_EPSILON ? FLT_EPSILON : ans;
}

static int
get_vert(int D, int i, int j, int k)
{
    int x, y;

    x = q->x;
    y = q->y;
    return q->verts[3 * (i + j * x + k * x * y) + D];
}

static void
set_vert(int D, int i, int j, int k)
{
    double u;
    int d;
    int x, y;
    double *cu;

    x = q->x;
    y = q->y;

    cu = q->cu;
    d = CuDir[D];
    u = rpos(cu[OOO], cu[d]);
    realloc_ver();
    q->ver[3 * q->nv + X] = i;
    q->ver[3 * q->nv + Y] = j;
    q->ver[3 * q->nv + Z] = k;
    q->ver[3 * q->nv + D] += u;
    q->verts[3 * (i + j * x + k * x * y) + D] = q->nv++;
}

static void
intersection(void)
{
    int i, j, k;
    int x, y, z;
    double *cu;

    x = q->x;
    y = q->y;
    z = q->z;
    cu = q->cu;

    for (k = 0; k < z; k++)
	for (j = 0; j < y; j++)
	    for (i = 0; i < x; i++) {
		cu[OOO] = get_data(i, j, k);
		if (i < x - 1)
		    cu[IOO] = get_data(i + 1, j, k);
		else
		    cu[IOO] = cu[OOO];
		if (j < y - 1)
		    cu[OIO] = get_data(i, j + 1, k);
		else
		    cu[OIO] = cu[OOO];
		if (k < z - 1)
		    cu[OOI] = get_data(i, j, k + 1);
		else
		    cu[OOI] = cu[OOO];
		if (cu[OOO] < 0) {
		    if (cu[IOO] > 0)
			set_vert(X, i, j, k);
		    if (cu[OIO] > 0)
			set_vert(Y, i, j, k);
		    if (cu[OOI] > 0)
			set_vert(Z, i, j, k);
		} else {
		    if (cu[IOO] < 0)
			set_vert(X, i, j, k);
		    if (cu[OIO] < 0)
			set_vert(Y, i, j, k);
		    if (cu[OOI] < 0)
			set_vert(Z, i, j, k);
		}
	    }
}

static int
test_face(int face)
{
    double A, B, C, D;
    double *cu;

    cu = q->cu;
    switch (face) {
    case -1:
    case 1:
	A = cu[OOO];
	B = cu[OOI];
	C = cu[IOI];
	D = cu[IOO];
	break;
    case -2:
    case 2:
	A = cu[IOO];
	B = cu[IOI];
	C = cu[III];
	D = cu[IIO];
	break;
    case -3:
    case 3:
	A = cu[IIO];
	B = cu[III];
	C = cu[OII];
	D = cu[OIO];
	break;
    case -4:
    case 4:
	A = cu[OIO];
	B = cu[OII];
	C = cu[OOI];
	D = cu[OOO];
	break;
    case -5:
    case 5:
	A = cu[OOO];
	B = cu[OIO];
	C = cu[IIO];
	D = cu[IOO];
	break;
    case -6:
    case 6:
	A = cu[OOI];
	B = cu[OII];
	C = cu[III];
	D = cu[IOI];
	break;
    default:
	fprintf(stderr, "Invalid face code %d\n", face);
	exit(2);
    };
    if (fabs(A * C - B * D) < FLT_EPSILON)
	return face >= 0;
    return face * A * (A * C - B * D) >= 0;
}

static int
test_interior(int s)
{
    double t, At = 0, Bt = 0, Ct = 0, Dt = 0, a, b;
    int test = 0;
    int edge = -1;
    double *cu;

    cu = q->cu;
    switch (q->Case) {
    case 4:
    case 10:
	a = (cu[OOI] - cu[OOO]) * (cu[III] - cu[IIO]) -
	    (cu[OII] - cu[OIO]) * (cu[IOI] - cu[IOO]);
	b = cu[IIO] * (cu[OOI] - cu[OOO]) + cu[OOO] * (cu[III] - cu[IIO])
	    - cu[IOO] * (cu[OII] - cu[OIO]) -
	    cu[OIO] * (cu[IOI] - cu[IOO]);
	t = -b / (2 * a);
	if (t < 0 || t > 1)
	    return s > 0;
	At = wavg(cu[OOO], cu[OOI], t);
	Bt = wavg(cu[OIO], cu[OII], t);
	Ct = wavg(cu[IIO], cu[III], t);
	Dt = wavg(cu[IOO], cu[IOI], t);
	break;
    case 6:
    case 7:
    case 12:
    case 13:
	switch (q->Case) {
	case 6:
	    edge = test6[q->config][2];
	    break;
	case 7:
	    edge = test7[q->config][4];
	    break;
	case 12:
	    edge = test12[q->config][3];
	    break;
	case 13:
	    edge = tiling13_5_1[q->config][q->subconfig][0];
	    break;
	}
	switch (edge) {
	case 0:
	    t = rpos(cu[OOO], cu[IOO]);
	    At = 0;
	    Bt = wavg(cu[OIO], cu[IIO], t);
	    Ct = wavg(cu[OII], cu[III], t);
	    Dt = wavg(cu[OOI], cu[IOI], t);
	    break;
	case 1:
	    t = rpos(cu[IOO], cu[IIO]);
	    At = 0;
	    Bt = wavg(cu[OOO], cu[OIO], t);
	    Ct = wavg(cu[OOI], cu[OII], t);
	    Dt = wavg(cu[IOI], cu[III], t);
	    break;
	case 2:
	    t = rpos(cu[IIO], cu[OIO]);
	    At = 0;
	    Bt = wavg(cu[IOO], cu[OOO], t);
	    Ct = wavg(cu[IOI], cu[OOI], t);
	    Dt = wavg(cu[III], cu[OII], t);
	    break;
	case 3:
	    t = rpos(cu[OIO], cu[OOO]);
	    At = 0;
	    Bt = wavg(cu[IIO], cu[IOO], t);
	    Ct = wavg(cu[III], cu[IOI], t);
	    Dt = wavg(cu[OII], cu[OOI], t);
	    break;
	case 4:
	    t = rpos(cu[OOI], cu[IOI]);
	    At = 0;
	    Bt = wavg(cu[OII], cu[III], t);
	    Ct = wavg(cu[OIO], cu[IIO], t);
	    Dt = wavg(cu[OOO], cu[IOO], t);
	    break;
	case 5:
	    t = rpos(cu[IOI], cu[III]);
	    At = 0;
	    Bt = wavg(cu[OOI], cu[OII], t);
	    Ct = wavg(cu[OOO], cu[OIO], t);
	    Dt = wavg(cu[IOO], cu[IIO], t);
	    break;
	case 6:
	    t = rpos(cu[III], cu[OII]);
	    At = 0;
	    Bt = wavg(cu[IOI], cu[OOI], t);
	    Ct = wavg(cu[IOO], cu[OOO], t);
	    Dt = wavg(cu[IIO], cu[OIO], t);
	    break;
	case 7:
	    t = rpos(cu[OII], cu[OOI]);
	    At = 0;
	    Bt = wavg(cu[III], cu[IOI], t);
	    Ct = wavg(cu[IIO], cu[IOO], t);
	    Dt = wavg(cu[OIO], cu[OOO], t);
	    break;
	case 8:
	    t = rpos(cu[OOO], cu[OOI]);
	    At = 0;
	    Bt = wavg(cu[OIO], cu[OII], t);
	    Ct = wavg(cu[IIO], cu[III], t);
	    Dt = wavg(cu[IOO], cu[IOI], t);
	    break;
	case 9:
	    t = rpos(cu[IOO], cu[IOI]);
	    At = 0;
	    Bt = wavg(cu[OOO], cu[OOI], t);
	    Ct = wavg(cu[OIO], cu[OII], t);
	    Dt = wavg(cu[IIO], cu[III], t);
	    break;
	case 10:
	    t = rpos(cu[IIO], cu[III]);
	    At = 0;
	    Bt = wavg(cu[IOO], cu[IOI], t);
	    Ct = wavg(cu[OOO], cu[OOI], t);
	    Dt = wavg(cu[OIO], cu[OII], t);
	    break;
	case 11:
	    t = rpos(cu[OIO], cu[OII]);
	    At = 0;
	    Bt = wavg(cu[IIO], cu[III], t);
	    Ct = wavg(cu[IOO], cu[IOI], t);
	    Dt = wavg(cu[OOO], cu[OOI], t);
	    break;
	default:
	    printf("Invalid edge %d\n", edge);
	    break;
	}
	break;
    default:
	printf("Invalid ambiguous case %d\n", q->Case);
	break;
    }
    if (At >= 0)
	test++;
    if (Bt >= 0)
	test += 2;
    if (Ct >= 0)
	test += 4;
    if (Dt >= 0)
	test += 8;
    switch (test) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 6:
    case 8:
    case 9:
    case 12:
	return s > 0;
    case 7:
    case 11:
    case 13:
    case 14:
    case 15:
	return s < 0;
    case 5:
	if (At * Ct - Bt * Dt < FLT_EPSILON)
	    return s > 0;
	break;
    case 10:
	if (At * Ct - Bt * Dt >= FLT_EPSILON)
	    return s > 0;
	break;
    }
    return s < 0;
}

static void
process_cu(void)
{
    int v12 = -1;

    q->subconfig = 0;
    switch (q->Case) {
    case 0:
	break;
    case 1:
	add_tri2(tiling1[q->config], 1);
	break;
    case 2:
	add_tri2(tiling2[q->config], 2);
	break;
    case 3:
	if (test_face(test3[q->config]))
	    add_tri2(tiling3_2[q->config], 4);
	else
	    add_tri2(tiling3_1[q->config], 2);
	break;
    case 4:
	if (test_interior(test4[q->config]))
	    add_tri2(tiling4_1[q->config], 2);
	else
	    add_tri2(tiling4_2[q->config], 6);
	break;
    case 5:
	add_tri2(tiling5[q->config], 3);
	break;
    case 6:
	if (test_face(test6[q->config][0]))
	    add_tri2(tiling6_2[q->config], 5);
	else {
	    if (test_interior(test6[q->config][1]))
		add_tri2(tiling6_1_1[q->config], 3);
	    else {
		v12 = add_c_vertex();
		add_tri3(tiling6_1_2[q->config], 9, v12);
	    }
	}
	break;
    case 7:
	if (test_face(test7[q->config][0]))
	    q->subconfig += 1;
	if (test_face(test7[q->config][1]))
	    q->subconfig += 2;
	if (test_face(test7[q->config][2]))
	    q->subconfig += 4;
	switch (q->subconfig) {
	case 0:
	    add_tri2(tiling7_1[q->config], 3);
	    break;
	case 1:
	    add_tri2(tiling7_2[q->config][0], 5);
	    break;
	case 2:
	    add_tri2(tiling7_2[q->config][1], 5);
	    break;
	case 3:
	    v12 = add_c_vertex();
	    add_tri3(tiling7_3[q->config][0], 9, v12);
	    break;
	case 4:
	    add_tri2(tiling7_2[q->config][2], 5);
	    break;
	case 5:
	    v12 = add_c_vertex();
	    add_tri3(tiling7_3[q->config][1], 9, v12);
	    break;
	case 6:
	    v12 = add_c_vertex();
	    add_tri3(tiling7_3[q->config][2], 9, v12);
	    break;
	case 7:
	    if (test_interior(test7[q->config][3]))
		add_tri2(tiling7_4_2[q->config], 9);
	    else
		add_tri2(tiling7_4_1[q->config], 5);
	    break;
	};
	break;
    case 8:
	add_tri2(tiling8[q->config], 2);
	break;
    case 9:
	add_tri2(tiling9[q->config], 4);
	break;
    case 10:
	if (test_face(test10[q->config][0])) {
	    if (test_face(test10[q->config][1]))
		add_tri2(tiling10_1_1_[q->config], 4);
	    else {
		v12 = add_c_vertex();
		add_tri3(tiling10_2[q->config], 8, v12);
	    }
	} else {
	    if (test_face(test10[q->config][1])) {
		v12 = add_c_vertex();
		add_tri3(tiling10_2_[q->config], 8, v12);
	    } else {
		if (test_interior(test10[q->config][2]))
		    add_tri2(tiling10_1_1[q->config], 4);
		else
		    add_tri2(tiling10_1_2[q->config], 8);
	    }
	}
	break;
    case 11:
	add_tri2(tiling11[q->config], 4);
	break;
    case 12:
	if (test_face(test12[q->config][0])) {
	    if (test_face(test12[q->config][1]))
		add_tri2(tiling12_1_1_[q->config], 4);
	    else {
		v12 = add_c_vertex();
		add_tri3(tiling12_2[q->config], 8, v12);
	    }
	} else {
	    if (test_face(test12[q->config][1])) {
		v12 = add_c_vertex();
		add_tri3(tiling12_2_[q->config], 8, v12);
	    } else {
		if (test_interior(test12[q->config][2]))
		    add_tri2(tiling12_1_1[q->config], 4);
		else
		    add_tri2(tiling12_1_2[q->config], 8);
	    }
	}
	break;
    case 13:
	if (test_face(test13[q->config][0]))
	    q->subconfig += 1;
	if (test_face(test13[q->config][1]))
	    q->subconfig += 2;
	if (test_face(test13[q->config][2]))
	    q->subconfig += 4;
	if (test_face(test13[q->config][3]))
	    q->subconfig += 8;
	if (test_face(test13[q->config][4]))
	    q->subconfig += 16;
	if (test_face(test13[q->config][5]))
	    q->subconfig += 32;
	switch (subconfig13[q->subconfig]) {
	case 0:
	    add_tri2(tiling13_1[q->config], 4);
	    break;
	case 1:
	    add_tri2(tiling13_2[q->config][0], 6);
	    break;
	case 2:
	    add_tri2(tiling13_2[q->config][1], 6);
	    break;
	case 3:
	    add_tri2(tiling13_2[q->config][2], 6);
	    break;
	case 4:
	    add_tri2(tiling13_2[q->config][3], 6);
	    break;
	case 5:
	    add_tri2(tiling13_2[q->config][4], 6);
	    break;
	case 6:
	    add_tri2(tiling13_2[q->config][5], 6);
	    break;
	case 7:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][0], 10, v12);
	    break;
	case 8:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][1], 10, v12);
	    break;
	case 9:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][2], 10, v12);
	    break;
	case 10:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][3], 10, v12);
	    break;
	case 11:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][4], 10, v12);
	    break;
	case 12:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][5], 10, v12);
	    break;
	case 13:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][6], 10, v12);
	    break;
	case 14:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][7], 10, v12);
	    break;
	case 15:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][8], 10, v12);
	    break;
	case 16:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][9], 10, v12);
	    break;
	case 17:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][10], 10, v12);
	    break;
	case 18:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[q->config][11], 10, v12);
	    break;
	case 19:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_4[q->config][0], 12, v12);
	    break;
	case 20:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_4[q->config][1], 12, v12);
	    break;
	case 21:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_4[q->config][2], 12, v12);
	    break;
	case 22:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_4[q->config][3], 12, v12);
	    break;
	case 23:
	    q->subconfig = 0;
	    if (test_interior(test13[q->config][6]))
		add_tri2(tiling13_5_1[q->config][0], 6);
	    else
		add_tri2(tiling13_5_2[q->config][0], 10);
	    break;
	case 24:
	    q->subconfig = 1;
	    if (test_interior(test13[q->config][6]))
		add_tri2(tiling13_5_1[q->config][1], 6);
	    else
		add_tri2(tiling13_5_2[q->config][1], 10);
	    break;
	case 25:
	    q->subconfig = 2;
	    if (test_interior(test13[q->config][6]))
		add_tri2(tiling13_5_1[q->config][2], 6);
	    else
		add_tri2(tiling13_5_2[q->config][2], 10);
	    break;
	case 26:
	    q->subconfig = 3;
	    if (test_interior(test13[q->config][6]))
		add_tri2(tiling13_5_1[q->config][3], 6);
	    else
		add_tri2(tiling13_5_2[q->config][3], 10);
	    break;
	case 27:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][0], 10, v12);
	    break;
	case 28:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][1], 10, v12);
	    break;
	case 29:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][2], 10, v12);
	    break;
	case 30:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][3], 10, v12);
	    break;
	case 31:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][4], 10, v12);
	    break;
	case 32:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][5], 10, v12);
	    break;
	case 33:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][6], 10, v12);
	    break;
	case 34:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][7], 10, v12);
	    break;
	case 35:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][8], 10, v12);
	    break;
	case 36:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][9], 10, v12);
	    break;
	case 37:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][10], 10, v12);
	    break;
	case 38:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[q->config][11], 10, v12);
	    break;
	case 39:
	    add_tri2(tiling13_2_[q->config][0], 6);
	    break;
	case 40:
	    add_tri2(tiling13_2_[q->config][1], 6);
	    break;
	case 41:
	    add_tri2(tiling13_2_[q->config][2], 6);
	    break;
	case 42:
	    add_tri2(tiling13_2_[q->config][3], 6);
	    break;
	case 43:
	    add_tri2(tiling13_2_[q->config][4], 6);
	    break;
	case 44:
	    add_tri2(tiling13_2_[q->config][5], 6);
	    break;
	case 45:
	    add_tri2(tiling13_1_[q->config], 4);
	    break;
	default:
	    printf("Impossible case 13?\n");
	}
	break;
    case 14:
	add_tri2(tiling14[q->config], 4);
	break;
    };
}

static void
add_tri2(const int *trig, int n)
{
    add_tri3(trig, n, -1);
}

static void
add_tri3(const int *trig, int n, int v12)
{
    int tv[3];
    int t;
    int g, *o;
    int I, J, K;

    I = q->I;
    J = q->J;
    K = q->K;

    for (t = 0; t < 3 * n; t++) {
	g = trig[t];
	if (g < 12) {
	    o = COff[g];
	    tv[t % 3] = get_vert(o[3], I + o[X], J + o[Y], K + o[Z]);
	} else
	    tv[t % 3] = v12;
	if (tv[t % 3] == -1) {
	    printf("Marching Cubes: invalid triangle %d\n", q->nt + 1);
	}
	if (t % 3 == 2) {
	    realloc_tri();
	    q->tri[3 * q->nt + X] = tv[0];
	    q->tri[3 * q->nt + Y] = tv[1];
	    q->tri[3 * q->nt + Z] = tv[2];
	    q->nt++;
	}
    }
}

static void
realloc_ver(void)
{
    if (q->nv >= q->Nv) {
	q->Nv *= 2;
	q->ver = realloc(q->ver, 3 * q->Nv * sizeof(*q->ver));
	fprintf(stderr, "%d allocated ver\n", q->Nv);
    }
}

static void
realloc_tri(void)
{
    if (q->nt >= q->Nt) {
	q->Nt *= 2;
	q->tri = realloc(q->tri, 3 * q->Nt * sizeof(*q->tri));
	fprintf(stderr, "%d allocated tri\n", q->Nt);
    }
}

static int
add_c_vertex(void)
{
    double u, rx, ry, rz;
    int g, m, *o;
    int I, J, K;

    I = q->I;
    J = q->J;
    K = q->K;

    u = 0;
    rx = ry = rz = 0;
    for (g = 0; g < SIZE(COff); g++) {
	o = COff[g];
	m = get_vert(o[3], I + o[X], J + o[Y], K + o[Z]);
	if (m != -1) {
	    u++;
	    rx += q->ver[3 * m + X];
	    ry += q->ver[3 * m + Y];
	    rz += q->ver[3 * m + Z];
	}
    }
    realloc_ver();
    q->ver[3 * q->nv + X] = rx / u;
    q->ver[3 * q->nv + Y] = ry / u;
    q->ver[3 * q->nv + Z] = rz / u;
    return q->nv++;
}
