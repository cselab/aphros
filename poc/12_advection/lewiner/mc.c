#include <float.h>
#include <math.h>
#include <memory.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "mc.h"
#include "table.h"

#define ALLOC_SIZE (99999)
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

static double get_data(int i, int j, int k);
static void process_cu(void);
static int test_face(int face);
static int test_interior(int s);
static void intersection(void);
static void add_tri3(const int *trig, int n, int v12);
static void add_tri2(const int *trig, int n);
static void test_vertex_addition(void);
static int add_c_vertex(void);
static int get_vert(int, int, int, int);
static void set_vert(int, int, int, int);

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

void
ini(int x0, int y0, int z0, double *data0)
{
    int i;
    int N;

    x = x0;
    y = y0;
    z = z0;
    data = data0;
    N = x * y * z;
    verts = malloc(3 * N * sizeof(*verts));
    for (i = 0; i < 3 * N; i++)
	verts[i] = -1;
    nv = nt = 0;
    Nv = Nt = ALLOC_SIZE;
    ver = malloc(3 * Nv * sizeof(*ver));
    tri = malloc(3 * Nt * sizeof(*tri));
}

void
run(void)
{
    int lut_entry, p, *o, i, j, k;

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
		I = i;
		J = j;
		K = k;
		Case = cases[lut_entry][0];
		config = cases[lut_entry][1];
		process_cu();
	    }
}

static double
get_data(int i, int j, int k)
{
    double ans;

    ans = data[i + j * x + k * x * y];
    return fabs(ans) < FLT_EPSILON ? FLT_EPSILON : ans;
}

static int
get_vert(int D, int i, int j, int k)
{
    return verts[3 * (i + j * x + k * x * y) + D];
}

static void
set_vert(int D, int i, int j, int k)
{
    double u;
    int d;

    test_vertex_addition();
    d = CuDir[D];
    u = (cu[OOO]) / (cu[OOO] - cu[d]);
    ver[3 * nv + X] = i;
    ver[3 * nv + Y] = j;
    ver[3 * nv + Z] = k;
    ver[3 * nv + D] += u;
    verts[3 * (i + j * x + k * x * y) + D] = nv++;
}

static void
intersection(void)
{
    int i, j, k;

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

    switch (Case) {
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
	At = cu[OOO] + (cu[OOI] - cu[OOO]) * t;
	Bt = cu[OIO] + (cu[OII] - cu[OIO]) * t;
	Ct = cu[IIO] + (cu[III] - cu[IIO]) * t;
	Dt = cu[IOO] + (cu[IOI] - cu[IOO]) * t;
	break;
    case 6:
    case 7:
    case 12:
    case 13:
	switch (Case) {
	case 6:
	    edge = test6[config][2];
	    break;
	case 7:
	    edge = test7[config][4];
	    break;
	case 12:
	    edge = test12[config][3];
	    break;
	case 13:
	    edge = tiling13_5_1[config][subconfig][0];
	    break;
	}
	switch (edge) {
	case 0:
	    t = cu[OOO] / (cu[OOO] - cu[IOO]);
	    At = 0;
	    Bt = cu[OIO] + (cu[IIO] - cu[OIO]) * t;
	    Ct = cu[OII] + (cu[III] - cu[OII]) * t;
	    Dt = cu[OOI] + (cu[IOI] - cu[OOI]) * t;
	    break;
	case 1:
	    t = cu[IOO] / (cu[IOO] - cu[IIO]);
	    At = 0;
	    Bt = cu[OOO] + (cu[OIO] - cu[OOO]) * t;
	    Ct = cu[OOI] + (cu[OII] - cu[OOI]) * t;
	    Dt = cu[IOI] + (cu[III] - cu[IOI]) * t;
	    break;
	case 2:
	    t = cu[IIO] / (cu[IIO] - cu[OIO]);
	    At = 0;
	    Bt = cu[IOO] + (cu[OOO] - cu[IOO]) * t;
	    Ct = cu[IOI] + (cu[OOI] - cu[IOI]) * t;
	    Dt = cu[III] + (cu[OII] - cu[III]) * t;
	    break;
	case 3:
	    t = cu[OIO] / (cu[OIO] - cu[OOO]);
	    At = 0;
	    Bt = cu[IIO] + (cu[IOO] - cu[IIO]) * t;
	    Ct = cu[III] + (cu[IOI] - cu[III]) * t;
	    Dt = cu[OII] + (cu[OOI] - cu[OII]) * t;
	    break;
	case 4:
	    t = cu[OOI] / (cu[OOI] - cu[IOI]);
	    At = 0;
	    Bt = cu[OII] + (cu[III] - cu[OII]) * t;
	    Ct = cu[OIO] + (cu[IIO] - cu[OIO]) * t;
	    Dt = cu[OOO] + (cu[IOO] - cu[OOO]) * t;
	    break;
	case 5:
	    t = cu[IOI] / (cu[IOI] - cu[III]);
	    At = 0;
	    Bt = cu[OOI] + (cu[OII] - cu[OOI]) * t;
	    Ct = cu[OOO] + (cu[OIO] - cu[OOO]) * t;
	    Dt = cu[IOO] + (cu[IIO] - cu[IOO]) * t;
	    break;
	case 6:
	    t = cu[III] / (cu[III] - cu[OII]);
	    At = 0;
	    Bt = cu[IOI] + (cu[OOI] - cu[IOI]) * t;
	    Ct = cu[IOO] + (cu[OOO] - cu[IOO]) * t;
	    Dt = cu[IIO] + (cu[OIO] - cu[IIO]) * t;
	    break;
	case 7:
	    t = cu[OII] / (cu[OII] - cu[OOI]);
	    At = 0;
	    Bt = cu[III] + (cu[IOI] - cu[III]) * t;
	    Ct = cu[IIO] + (cu[IOO] - cu[IIO]) * t;
	    Dt = cu[OIO] + (cu[OOO] - cu[OIO]) * t;
	    break;
	case 8:
	    t = cu[OOO] / (cu[OOO] - cu[OOI]);
	    At = 0;
	    Bt = cu[OIO] + (cu[OII] - cu[OIO]) * t;
	    Ct = cu[IIO] + (cu[III] - cu[IIO]) * t;
	    Dt = cu[IOO] + (cu[IOI] - cu[IOO]) * t;
	    break;
	case 9:
	    t = cu[IOO] / (cu[IOO] - cu[IOI]);
	    At = 0;
	    Bt = cu[OOO] + (cu[OOI] - cu[OOO]) * t;
	    Ct = cu[OIO] + (cu[OII] - cu[OIO]) * t;
	    Dt = cu[IIO] + (cu[III] - cu[IIO]) * t;
	    break;
	case 10:
	    t = cu[IIO] / (cu[IIO] - cu[III]);
	    At = 0;
	    Bt = cu[IOO] + (cu[IOI] - cu[IOO]) * t;
	    Ct = cu[OOO] + (cu[OOI] - cu[OOO]) * t;
	    Dt = cu[OIO] + (cu[OII] - cu[OIO]) * t;
	    break;
	case 11:
	    t = cu[OIO] / (cu[OIO] - cu[OII]);
	    At = 0;
	    Bt = cu[IIO] + (cu[III] - cu[IIO]) * t;
	    Ct = cu[IOO] + (cu[IOI] - cu[IOO]) * t;
	    Dt = cu[OOO] + (cu[OOI] - cu[OOO]) * t;
	    break;
	default:
	    printf("Invalid edge %d\n", edge);
	    break;
	}
	break;
    default:
	printf("Invalid ambiguous case %d\n", Case);
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

    subconfig = 0;
    switch (Case) {
    case 0:
	break;
    case 1:
	add_tri2(tiling1[config], 1);
	break;
    case 2:
	add_tri2(tiling2[config], 2);
	break;
    case 3:
	if (test_face(test3[config]))
	    add_tri2(tiling3_2[config], 4);
	else
	    add_tri2(tiling3_1[config], 2);
	break;
    case 4:
	if (test_interior(test4[config]))
	    add_tri2(tiling4_1[config], 2);
	else
	    add_tri2(tiling4_2[config], 6);
	break;
    case 5:
	add_tri2(tiling5[config], 3);
	break;
    case 6:
	if (test_face(test6[config][0]))
	    add_tri2(tiling6_2[config], 5);
	else {
	    if (test_interior(test6[config][1]))
		add_tri2(tiling6_1_1[config], 3);
	    else {
		v12 = add_c_vertex();
		add_tri3(tiling6_1_2[config], 9, v12);
	    }
	}
	break;
    case 7:
	if (test_face(test7[config][0]))
	    subconfig += 1;
	if (test_face(test7[config][1]))
	    subconfig += 2;
	if (test_face(test7[config][2]))
	    subconfig += 4;
	switch (subconfig) {
	case 0:
	    add_tri2(tiling7_1[config], 3);
	    break;
	case 1:
	    add_tri2(tiling7_2[config][0], 5);
	    break;
	case 2:
	    add_tri2(tiling7_2[config][1], 5);
	    break;
	case 3:
	    v12 = add_c_vertex();
	    add_tri3(tiling7_3[config][0], 9, v12);
	    break;
	case 4:
	    add_tri2(tiling7_2[config][2], 5);
	    break;
	case 5:
	    v12 = add_c_vertex();
	    add_tri3(tiling7_3[config][1], 9, v12);
	    break;
	case 6:
	    v12 = add_c_vertex();
	    add_tri3(tiling7_3[config][2], 9, v12);
	    break;
	case 7:
	    if (test_interior(test7[config][3]))
		add_tri2(tiling7_4_2[config], 9);
	    else
		add_tri2(tiling7_4_1[config], 5);
	    break;
	};
	break;
    case 8:
	add_tri2(tiling8[config], 2);
	break;
    case 9:
	add_tri2(tiling9[config], 4);
	break;
    case 10:
	if (test_face(test10[config][0])) {
	    if (test_face(test10[config][1]))
		add_tri2(tiling10_1_1_[config], 4);
	    else {
		v12 = add_c_vertex();
		add_tri3(tiling10_2[config], 8, v12);
	    }
	} else {
	    if (test_face(test10[config][1])) {
		v12 = add_c_vertex();
		add_tri3(tiling10_2_[config], 8, v12);
	    } else {
		if (test_interior(test10[config][2]))
		    add_tri2(tiling10_1_1[config], 4);
		else
		    add_tri2(tiling10_1_2[config], 8);
	    }
	}
	break;
    case 11:
	add_tri2(tiling11[config], 4);
	break;
    case 12:
	if (test_face(test12[config][0])) {
	    if (test_face(test12[config][1]))
		add_tri2(tiling12_1_1_[config], 4);
	    else {
		v12 = add_c_vertex();
		add_tri3(tiling12_2[config], 8, v12);
	    }
	} else {
	    if (test_face(test12[config][1])) {
		v12 = add_c_vertex();
		add_tri3(tiling12_2_[config], 8, v12);
	    } else {
		if (test_interior(test12[config][2]))
		    add_tri2(tiling12_1_1[config], 4);
		else
		    add_tri2(tiling12_1_2[config], 8);
	    }
	}
	break;
    case 13:
	if (test_face(test13[config][0]))
	    subconfig += 1;
	if (test_face(test13[config][1]))
	    subconfig += 2;
	if (test_face(test13[config][2]))
	    subconfig += 4;
	if (test_face(test13[config][3]))
	    subconfig += 8;
	if (test_face(test13[config][4]))
	    subconfig += 16;
	if (test_face(test13[config][5]))
	    subconfig += 32;
	switch (subconfig13[subconfig]) {
	case 0:
	    add_tri2(tiling13_1[config], 4);
	    break;
	case 1:
	    add_tri2(tiling13_2[config][0], 6);
	    break;
	case 2:
	    add_tri2(tiling13_2[config][1], 6);
	    break;
	case 3:
	    add_tri2(tiling13_2[config][2], 6);
	    break;
	case 4:
	    add_tri2(tiling13_2[config][3], 6);
	    break;
	case 5:
	    add_tri2(tiling13_2[config][4], 6);
	    break;
	case 6:
	    add_tri2(tiling13_2[config][5], 6);
	    break;
	case 7:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][0], 10, v12);
	    break;
	case 8:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][1], 10, v12);
	    break;
	case 9:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][2], 10, v12);
	    break;
	case 10:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][3], 10, v12);
	    break;
	case 11:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][4], 10, v12);
	    break;
	case 12:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][5], 10, v12);
	    break;
	case 13:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][6], 10, v12);
	    break;
	case 14:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][7], 10, v12);
	    break;
	case 15:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][8], 10, v12);
	    break;
	case 16:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][9], 10, v12);
	    break;
	case 17:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][10], 10, v12);
	    break;
	case 18:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3[config][11], 10, v12);
	    break;
	case 19:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_4[config][0], 12, v12);
	    break;
	case 20:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_4[config][1], 12, v12);
	    break;
	case 21:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_4[config][2], 12, v12);
	    break;
	case 22:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_4[config][3], 12, v12);
	    break;
	case 23:
	    subconfig = 0;
	    if (test_interior(test13[config][6]))
		add_tri2(tiling13_5_1[config][0], 6);
	    else
		add_tri2(tiling13_5_2[config][0], 10);
	    break;
	case 24:
	    subconfig = 1;
	    if (test_interior(test13[config][6]))
		add_tri2(tiling13_5_1[config][1], 6);
	    else
		add_tri2(tiling13_5_2[config][1], 10);
	    break;
	case 25:
	    subconfig = 2;
	    if (test_interior(test13[config][6]))
		add_tri2(tiling13_5_1[config][2], 6);
	    else
		add_tri2(tiling13_5_2[config][2], 10);
	    break;
	case 26:
	    subconfig = 3;
	    if (test_interior(test13[config][6]))
		add_tri2(tiling13_5_1[config][3], 6);
	    else
		add_tri2(tiling13_5_2[config][3], 10);
	    break;
	case 27:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][0], 10, v12);
	    break;
	case 28:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][1], 10, v12);
	    break;
	case 29:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][2], 10, v12);
	    break;
	case 30:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][3], 10, v12);
	    break;
	case 31:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][4], 10, v12);
	    break;
	case 32:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][5], 10, v12);
	    break;
	case 33:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][6], 10, v12);
	    break;
	case 34:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][7], 10, v12);
	    break;
	case 35:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][8], 10, v12);
	    break;
	case 36:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][9], 10, v12);
	    break;
	case 37:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][10], 10, v12);
	    break;
	case 38:
	    v12 = add_c_vertex();
	    add_tri3(tiling13_3_[config][11], 10, v12);
	    break;
	case 39:
	    add_tri2(tiling13_2_[config][0], 6);
	    break;
	case 40:
	    add_tri2(tiling13_2_[config][1], 6);
	    break;
	case 41:
	    add_tri2(tiling13_2_[config][2], 6);
	    break;
	case 42:
	    add_tri2(tiling13_2_[config][3], 6);
	    break;
	case 43:
	    add_tri2(tiling13_2_[config][4], 6);
	    break;
	case 44:
	    add_tri2(tiling13_2_[config][5], 6);
	    break;
	case 45:
	    add_tri2(tiling13_1_[config], 4);
	    break;
	default:
	    printf("Marching Cus: Impossible case 13?\n");
	}
	break;
    case 14:
	add_tri2(tiling14[config], 4);
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

    for (t = 0; t < 3 * n; t++) {
	g = trig[t];
	if (g < 12) {
	    o = COff[g];
	    tv[t % 3] = get_vert(o[3], I + o[X], J + o[Y], K + o[Z]);
	} else
	    tv[t % 3] = v12;
	if (tv[t % 3] == -1) {
	    printf("Marching Cubes: invalid triangle %d\n", nt + 1);
	}
	if (t % 3 == 2) {
	    if (nt >= Nt) {
		Nt *= 2;
		tri = realloc(tri, 3 * Nt * sizeof(*tri));
		fprintf(stderr, "%d allocated triangles\n", Nt);
	    }
	    tri[3 * nt + X] = tv[0];
	    tri[3 * nt + Y] = tv[1];
	    tri[3 * nt + Z] = tv[2];
	    nt++;
	}
    }
}

static void
test_vertex_addition()
{
    if (nv >= Nv) {
	Nv *= 2;
	ver = realloc(ver, 3 * Nv * sizeof(*ver));
	fprintf(stderr, "%d allocated ver\n", Nv);
    }
}

static int
add_c_vertex()
{
    double u, rx, ry, rz;
    int g, m, *o;

    u = 0;
    test_vertex_addition();
    rx = ry = rz = 0;
    for (g = 0; g < SIZE(COff); g++) {
	o = COff[g];
	m = get_vert(o[3], I + o[X], J + o[Y], K + o[Z]);
	if (m != -1) {
	    u++;
	    rx += ver[3 * m + X];
	    ry += ver[3 * m + Y];
	    rz += ver[3 * m + Z];
	}
    }
    ver[3 * nv + X] = rx / u;
    ver[3 * nv + Y] = ry / u;
    ver[3 * nv + Z] = rz / u;
    return nv++;
}

void
obj(void)
{
    int t;
    int v;
    int i;

    t = nt;
    v = nv;
    printf("# File type: ASCII OBJ\n");
    for (i = 0; i < v; i++)
	printf("v %.16g %.16g %.16g\n", ver[3 * i + X], ver[3 * i + Y],
	       ver[3 * i + Z]);
    for (i = 0; i < t; i++)
	printf("f %d %d %d\n", tri[3 * i + X] + 1, tri[3 * i + Y] + 1,
	       tri[3 * i + Z] + 1);
}
