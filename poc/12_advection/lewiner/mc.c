#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <float.h>
#include "mc.h"
#include "table.h"

#define ALLOC_SIZE 65536
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
static int CubeDir[] = { IOO, OIO, OOI };
static int CubeOff[][3] = {
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
    { 0, 1, 0, X },
    { 0, 1, 1, X },
    { 0, 0, 0, X },
    { 0, 0, 1, X },
    { 1, 0, 0, Y },
    { 1, 0, 1, Y },
    { 0, 0, 0, Y },
    { 0, 0, 1, Y },
    { 1, 1, 0, Z },
    { 1, 0, 0, Z },
    { 0, 1, 0, Z },
    { 0, 0, 0, Z }
};


static double get_data(int i, int j, int k);
static void process_cube();
static int test_face(int face);
static int test_interior(int s);
static void compute_intersection_points();
static void add_triangle3(const int *trig, int n, int v12);
static void add_triangle2(const int *trig, int n);
static void test_vertex_addition();
static int add_c_vertex();
static int get_vert(int, int, int, int);
static void set_vert(int, int, int, int);

double cube[8];
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
int *triangles;
double *vertices;

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
    return verts[3*(i + j * x + k * x * y) + D];
}

static void
set_vert(int D, int i, int j, int k)
{
    double u;
    int d;

    test_vertex_addition();
    d = CubeDir[D];
    u = (cube[OOO]) / (cube[OOO] - cube[d]);
    vertices[3 * nv + X] = i;
    vertices[3 * nv + Y] = j;
    vertices[3 * nv + Z] = k;
    vertices[3 * nv + D] += u;
    verts[3*(i + j * x + k * x * y) + D] = nv++;
}

void
MarchingCubes(int x0, int y0, int z0, double *data0)
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
    vertices = malloc(N * sizeof(*vertices));
    triangles = malloc(N * sizeof(*triangles));
}

void
run()
{
    int lut_entry, p, *o, i, j, k;

    compute_intersection_points();
    for (k = 0; k < z - 1; k++)
	for (j = 0; j < y - 1; j++)
	    for (i = 0; i < x - 1; i++) {
		lut_entry = 0;
		for (p = 0; p < SIZE(CubeOff); ++p) {
		    o = CubeOff[p];
		    cube[p] = get_data(i + o[X], j + o[Y], k + o[Z]);
		    if (cube[p] > 0)
			lut_entry += 1 << p;
		}
		I = i; J = j; K = k;
		Case = cases[lut_entry][0];
		config = cases[lut_entry][1];
		process_cube();
	    }
}

static void
compute_intersection_points()
{
    int i, j, k;

    for (k = 0; k < z; k++)
	for (j = 0; j < y; j++)
	    for (i = 0; i < x; i++) {
		cube[OOO] = get_data(i, j, k);
		if (i < x - 1)
		    cube[IOO] = get_data(i + 1, j, k);
		else
		    cube[IOO] = cube[OOO];
		if (j < y - 1)
		    cube[OIO] = get_data(i, j + 1, k);
		else
		    cube[OIO] = cube[OOO];
		if (k < z - 1)
		    cube[OOI] = get_data(i, j, k + 1);
		else
		    cube[OOI] = cube[OOO];
		if (cube[OOO] < 0) {
		    if (cube[IOO] > 0)
			set_vert(X, i, j, k);
		    if (cube[OIO] > 0)
			set_vert(Y, i, j, k);
		    if (cube[OOI] > 0)
			set_vert(Z, i, j, k);
		} else {
		    if (cube[IOO] < 0)
			set_vert(X, i, j, k);
		    if (cube[OIO] < 0)
			set_vert(Y, i, j, k);
		    if (cube[OOI] < 0)
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
	A = cube[OOO];
	B = cube[OOI];
	C = cube[IOI];
	D = cube[IOO];
	break;
    case -2:
    case 2:
	A = cube[IOO];
	B = cube[IOI];
	C = cube[III];
	D = cube[IIO];
	break;
    case -3:
    case 3:
	A = cube[IIO];
	B = cube[III];
	C = cube[OII];
	D = cube[OIO];
	break;
    case -4:
    case 4:
	A = cube[OIO];
	B = cube[OII];
	C = cube[OOI];
	D = cube[OOO];
	break;
    case -5:
    case 5:
	A = cube[OOO];
	B = cube[OIO];
	C = cube[IIO];
	D = cube[IOO];
	break;
    case -6:
    case 6:
	A = cube[OOI];
	B = cube[OII];
	C = cube[III];
	D = cube[IOI];
	break;
    default:
	printf("Invalid face code %d\n", face);
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
	a = (cube[OOI] - cube[OOO]) * (cube[III] - cube[IIO]) -
	    (cube[OII] - cube[OIO]) * (cube[IOI] - cube[IOO]);
	b = cube[IIO] * (cube[OOI] - cube[OOO]) + cube[OOO] * (cube[III] -
							       cube[IIO])
	    - cube[IOO] * (cube[OII] - cube[OIO]) -
	    cube[OIO] * (cube[IOI] - cube[IOO]);
	t = -b / (2 * a);
	if (t < 0 || t > 1)
	    return s > 0;
	At = cube[OOO] + (cube[OOI] - cube[OOO]) * t;
	Bt = cube[OIO] + (cube[OII] - cube[OIO]) * t;
	Ct = cube[IIO] + (cube[III] - cube[IIO]) * t;
	Dt = cube[IOO] + (cube[IOI] - cube[IOO]) * t;
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
	    t = cube[OOO] / (cube[OOO] - cube[IOO]);
	    At = 0;
	    Bt = cube[OIO] + (cube[IIO] - cube[OIO]) * t;
	    Ct = cube[OII] + (cube[III] - cube[OII]) * t;
	    Dt = cube[OOI] + (cube[IOI] - cube[OOI]) * t;
	    break;
	case 1:
	    t = cube[IOO] / (cube[IOO] - cube[IIO]);
	    At = 0;
	    Bt = cube[OOO] + (cube[OIO] - cube[OOO]) * t;
	    Ct = cube[OOI] + (cube[OII] - cube[OOI]) * t;
	    Dt = cube[IOI] + (cube[III] - cube[IOI]) * t;
	    break;
	case 2:
	    t = cube[IIO] / (cube[IIO] - cube[OIO]);
	    At = 0;
	    Bt = cube[IOO] + (cube[OOO] - cube[IOO]) * t;
	    Ct = cube[IOI] + (cube[OOI] - cube[IOI]) * t;
	    Dt = cube[III] + (cube[OII] - cube[III]) * t;
	    break;
	case 3:
	    t = cube[OIO] / (cube[OIO] - cube[OOO]);
	    At = 0;
	    Bt = cube[IIO] + (cube[IOO] - cube[IIO]) * t;
	    Ct = cube[III] + (cube[IOI] - cube[III]) * t;
	    Dt = cube[OII] + (cube[OOI] - cube[OII]) * t;
	    break;
	case 4:
	    t = cube[OOI] / (cube[OOI] - cube[IOI]);
	    At = 0;
	    Bt = cube[OII] + (cube[III] - cube[OII]) * t;
	    Ct = cube[OIO] + (cube[IIO] - cube[OIO]) * t;
	    Dt = cube[OOO] + (cube[IOO] - cube[OOO]) * t;
	    break;
	case 5:
	    t = cube[IOI] / (cube[IOI] - cube[III]);
	    At = 0;
	    Bt = cube[OOI] + (cube[OII] - cube[OOI]) * t;
	    Ct = cube[OOO] + (cube[OIO] - cube[OOO]) * t;
	    Dt = cube[IOO] + (cube[IIO] - cube[IOO]) * t;
	    break;
	case 6:
	    t = cube[III] / (cube[III] - cube[OII]);
	    At = 0;
	    Bt = cube[IOI] + (cube[OOI] - cube[IOI]) * t;
	    Ct = cube[IOO] + (cube[OOO] - cube[IOO]) * t;
	    Dt = cube[IIO] + (cube[OIO] - cube[IIO]) * t;
	    break;
	case 7:
	    t = cube[OII] / (cube[OII] - cube[OOI]);
	    At = 0;
	    Bt = cube[III] + (cube[IOI] - cube[III]) * t;
	    Ct = cube[IIO] + (cube[IOO] - cube[IIO]) * t;
	    Dt = cube[OIO] + (cube[OOO] - cube[OIO]) * t;
	    break;
	case 8:
	    t = cube[OOO] / (cube[OOO] - cube[OOI]);
	    At = 0;
	    Bt = cube[OIO] + (cube[OII] - cube[OIO]) * t;
	    Ct = cube[IIO] + (cube[III] - cube[IIO]) * t;
	    Dt = cube[IOO] + (cube[IOI] - cube[IOO]) * t;
	    break;
	case 9:
	    t = cube[IOO] / (cube[IOO] - cube[IOI]);
	    At = 0;
	    Bt = cube[OOO] + (cube[OOI] - cube[OOO]) * t;
	    Ct = cube[OIO] + (cube[OII] - cube[OIO]) * t;
	    Dt = cube[IIO] + (cube[III] - cube[IIO]) * t;
	    break;
	case 10:
	    t = cube[IIO] / (cube[IIO] - cube[III]);
	    At = 0;
	    Bt = cube[IOO] + (cube[IOI] - cube[IOO]) * t;
	    Ct = cube[OOO] + (cube[OOI] - cube[OOO]) * t;
	    Dt = cube[OIO] + (cube[OII] - cube[OIO]) * t;
	    break;
	case 11:
	    t = cube[OIO] / (cube[OIO] - cube[OII]);
	    At = 0;
	    Bt = cube[IIO] + (cube[III] - cube[IIO]) * t;
	    Ct = cube[IOO] + (cube[IOI] - cube[IOO]) * t;
	    Dt = cube[OOO] + (cube[OOI] - cube[OOO]) * t;
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
process_cube()
{
    int v12 = -1;

    subconfig = 0;
    switch (Case) {
    case 0:
	break;
    case 1:
	add_triangle2(tiling1[config], 1);
	break;
    case 2:
	add_triangle2(tiling2[config], 2);
	break;
    case 3:
	if (test_face(test3[config]))
	    add_triangle2(tiling3_2[config], 4);
	else
	    add_triangle2(tiling3_1[config], 2);
	break;
    case 4:
	if (test_interior(test4[config]))
	    add_triangle2(tiling4_1[config], 2);
	else
	    add_triangle2(tiling4_2[config], 6);
	break;
    case 5:
	add_triangle2(tiling5[config], 3);
	break;
    case 6:
	if (test_face(test6[config][0]))
	    add_triangle2(tiling6_2[config], 5);
	else {
	    if (test_interior(test6[config][1]))
		add_triangle2(tiling6_1_1[config], 3);
	    else {
		v12 = add_c_vertex();
		add_triangle3(tiling6_1_2[config], 9, v12);
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
	    add_triangle2(tiling7_1[config], 3);
	    break;
	case 1:
	    add_triangle2(tiling7_2[config][0], 5);
	    break;
	case 2:
	    add_triangle2(tiling7_2[config][1], 5);
	    break;
	case 3:
	    v12 = add_c_vertex();
	    add_triangle3(tiling7_3[config][0], 9, v12);
	    break;
	case 4:
	    add_triangle2(tiling7_2[config][2], 5);
	    break;
	case 5:
	    v12 = add_c_vertex();
	    add_triangle3(tiling7_3[config][1], 9, v12);
	    break;
	case 6:
	    v12 = add_c_vertex();
	    add_triangle3(tiling7_3[config][2], 9, v12);
	    break;
	case 7:
	    if (test_interior(test7[config][3]))
		add_triangle2(tiling7_4_2[config], 9);
	    else
		add_triangle2(tiling7_4_1[config], 5);
	    break;
	};
	break;
    case 8:
	add_triangle2(tiling8[config], 2);
	break;
    case 9:
	add_triangle2(tiling9[config], 4);
	break;
    case 10:
	if (test_face(test10[config][0])) {
	    if (test_face(test10[config][1]))
		add_triangle2(tiling10_1_1_[config], 4);
	    else {
		v12 = add_c_vertex();
		add_triangle3(tiling10_2[config], 8, v12);
	    }
	} else {
	    if (test_face(test10[config][1])) {
		v12 = add_c_vertex();
		add_triangle3(tiling10_2_[config], 8, v12);
	    } else {
		if (test_interior(test10[config][2]))
		    add_triangle2(tiling10_1_1[config], 4);
		else
		    add_triangle2(tiling10_1_2[config], 8);
	    }
	}
	break;
    case 11:
	add_triangle2(tiling11[config], 4);
	break;
    case 12:
	if (test_face(test12[config][0])) {
	    if (test_face(test12[config][1]))
		add_triangle2(tiling12_1_1_[config], 4);
	    else {
		v12 = add_c_vertex();
		add_triangle3(tiling12_2[config], 8, v12);
	    }
	} else {
	    if (test_face(test12[config][1])) {
		v12 = add_c_vertex();
		add_triangle3(tiling12_2_[config], 8, v12);
	    } else {
		if (test_interior(test12[config][2]))
		    add_triangle2(tiling12_1_1[config], 4);
		else
		    add_triangle2(tiling12_1_2[config], 8);
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
	    add_triangle2(tiling13_1[config], 4);
	    break;
	case 1:
	    add_triangle2(tiling13_2[config][0], 6);
	    break;
	case 2:
	    add_triangle2(tiling13_2[config][1], 6);
	    break;
	case 3:
	    add_triangle2(tiling13_2[config][2], 6);
	    break;
	case 4:
	    add_triangle2(tiling13_2[config][3], 6);
	    break;
	case 5:
	    add_triangle2(tiling13_2[config][4], 6);
	    break;
	case 6:
	    add_triangle2(tiling13_2[config][5], 6);
	    break;
	case 7:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][0], 10, v12);
	    break;
	case 8:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][1], 10, v12);
	    break;
	case 9:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][2], 10, v12);
	    break;
	case 10:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][3], 10, v12);
	    break;
	case 11:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][4], 10, v12);
	    break;
	case 12:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][5], 10, v12);
	    break;
	case 13:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][6], 10, v12);
	    break;
	case 14:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][7], 10, v12);
	    break;
	case 15:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][8], 10, v12);
	    break;
	case 16:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][9], 10, v12);
	    break;
	case 17:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][10], 10, v12);
	    break;
	case 18:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3[config][11], 10, v12);
	    break;
	case 19:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_4[config][0], 12, v12);
	    break;
	case 20:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_4[config][1], 12, v12);
	    break;
	case 21:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_4[config][2], 12, v12);
	    break;
	case 22:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_4[config][3], 12, v12);
	    break;
	case 23:
	    subconfig = 0;
	    if (test_interior(test13[config][6]))
		add_triangle2(tiling13_5_1[config][0], 6);
	    else
		add_triangle2(tiling13_5_2[config][0], 10);
	    break;
	case 24:
	    subconfig = 1;
	    if (test_interior(test13[config][6]))
		add_triangle2(tiling13_5_1[config][1], 6);
	    else
		add_triangle2(tiling13_5_2[config][1], 10);
	    break;
	case 25:
	    subconfig = 2;
	    if (test_interior(test13[config][6]))
		add_triangle2(tiling13_5_1[config][2], 6);
	    else
		add_triangle2(tiling13_5_2[config][2], 10);
	    break;
	case 26:
	    subconfig = 3;
	    if (test_interior(test13[config][6]))
		add_triangle2(tiling13_5_1[config][3], 6);
	    else
		add_triangle2(tiling13_5_2[config][3], 10);
	    break;
	case 27:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][0], 10, v12);
	    break;
	case 28:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][1], 10, v12);
	    break;
	case 29:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][2], 10, v12);
	    break;
	case 30:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][3], 10, v12);
	    break;
	case 31:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][4], 10, v12);
	    break;
	case 32:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][5], 10, v12);
	    break;
	case 33:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][6], 10, v12);
	    break;
	case 34:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][7], 10, v12);
	    break;
	case 35:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][8], 10, v12);
	    break;
	case 36:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][9], 10, v12);
	    break;
	case 37:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][10], 10, v12);
	    break;
	case 38:
	    v12 = add_c_vertex();
	    add_triangle3(tiling13_3_[config][11], 10, v12);
	    break;
	case 39:
	    add_triangle2(tiling13_2_[config][0], 6);
	    break;
	case 40:
	    add_triangle2(tiling13_2_[config][1], 6);
	    break;
	case 41:
	    add_triangle2(tiling13_2_[config][2], 6);
	    break;
	case 42:
	    add_triangle2(tiling13_2_[config][3], 6);
	    break;
	case 43:
	    add_triangle2(tiling13_2_[config][4], 6);
	    break;
	case 44:
	    add_triangle2(tiling13_2_[config][5], 6);
	    break;
	case 45:
	    add_triangle2(tiling13_1_[config], 4);
	    break;
	default:
	    printf("Marching Cubes: Impossible case 13?\n");
	}
	break;
    case 14:
	add_triangle2(tiling14[config], 4);
	break;
    };
}

static void
add_triangle2(const int *trig, int n)
{
    add_triangle3(trig, n, -1);
}

static void
add_triangle3(const int *trig, int n, int v12)
{
    int tv[3];
    int t;

    for (t = 0; t < 3 * n; t++) {
	switch (trig[t]) {
	case 0:
	    tv[t % 3] = get_vert(X, I, J, K);
	    break;
	case 1:
	    tv[t % 3] = get_vert(Y, I + 1, J, K);
	    break;
	case 2:
	    tv[t % 3] = get_vert(X, I, J + 1, K);
	    break;
	case 3:
	    tv[t % 3] = get_vert(Y, I, J, K);
	    break;
	case 4:
	    tv[t % 3] = get_vert(X, I, J, K + 1);
	    break;
	case 5:
	    tv[t % 3] = get_vert(Y, I + 1, J, K + 1);
	    break;
	case 6:
	    tv[t % 3] = get_vert(X, I, J + 1, K + 1);
	    break;
	case 7:
	    tv[t % 3] = get_vert(Y, I, J, K + 1);
	    break;
	case 8:
	    tv[t % 3] = get_vert(Z, I, J, K);
	    break;
	case 9:
	    tv[t % 3] = get_vert(Z, I + 1, J, K);
	    break;
	case 10:
	    tv[t % 3] = get_vert(Z, I + 1, J + 1, K);
	    break;
	case 11:
	    tv[t % 3] = get_vert(Z, I, J + 1, K);
	    break;
	case 12:
	    tv[t % 3] = v12;
	    break;
	default:
	    break;
	}
	if (tv[t % 3] == -1) {
	    printf("Marching Cubes: invalid triangle %d\n", nt + 1);
	}
	if (t % 3 == 2) {
	    if (nt >= Nt) {
		Nt *= 2;
		triangles = realloc(triangles, Nt * sizeof(*triangles));
		printf("%d allocated triangles\n", Nt);
	    }
	    triangles[3 * nt + X] = tv[0];
	    triangles[3 * nt + Y] = tv[1];
	    triangles[3 * nt + Z] = tv[2];
	    nt++;
	}
    }
}

static void
test_vertex_addition()
{
    if (nv >= Nv) {
	Nv *= 2;
	vertices = realloc(vertices, Nv * sizeof(*vertices));
	printf("%d allocated vertices\n", Nv);
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
	    rx += vertices[3 * m + X];
	    ry += vertices[3 * m + Y];
	    rz += vertices[3 * m + Z];
	}
    }
    vertices[3 * nv + X] = rx / u;
    vertices[3 * nv + Y] = ry / u;
    vertices[3 * nv + Z] = rz / u;
    return nv++;
}

void
writeObj()
{
    int t;
    int v;
    int i;
    int *tri;
    double *ver;

    t = nt;
    v = nv;
    tri = triangles;
    ver = vertices;
    printf("# File type: ASCII OBJ\n");
    for (i = 0; i < v; i++)
	printf("v %.16g %.16g %.16g\n", ver[3 * i + X], ver[3 * i + Y],
	       ver[3 * i + Z]);
    for (i = 0; i < t; i++)
	printf("f %d %d %d\n", tri[3 * i + X] + 1, tri[3 * i + Y] + 1,
	       tri[3 * i + Z] + 1);
}
