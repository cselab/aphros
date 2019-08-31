#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <float.h>
#include "mc.h"
#include "table.h"
#define ALLOC_SIZE 65536
#define SIZE(a) (sizeof(a)/sizeof(*(a)))

enum { X, Y, Z };

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
int i;
int j;
int k;
int N;
int nt;
int Nt;
int nv;
int Nv;
int subconfig;
int *verts[3];
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
    return verts[D][i + j * x + k * x * y];
}

static void
set_vert(int D, int i, int j, int k)
{
    double u;
    int m[3];

    m[X] = 1;
    m[Y] = 3;
    m[Z] = 4;

    test_vertex_addition();
    u = (cube[0]) / (cube[0] - cube[m[D]]);
    vertices[3*nv + X] = i;
    vertices[3*nv + Y] = j;
    vertices[3*nv + Z] = k;
    vertices[3*nv + D] += u;
    verts[D][i + j * x + k * x * y] = nv++;
}

void
MarchingCubes(int x0, int y0, int z0, double *data0)
{
    int i;

    x = x0;
    y = y0;
    z = z0;
    data = data0;
    N = x * y * z;
    verts[X] = malloc(N * sizeof(*verts[X]));
    verts[Y] = malloc(N * sizeof(*verts[Y]));
    verts[Z] = malloc(N * sizeof(*verts[Z]));
    for (i = 0; i < N; i++) {
	verts[X][i] = -1;
	verts[Y][i] = -1;
	verts[Z][i] = -1;
    }
    nv = nt = 0;
    Nv = Nt = ALLOC_SIZE;
    vertices = malloc(N * sizeof(*vertices));
    triangles = malloc(N * sizeof(*triangles));
}

void
run()
{
    int lut_entry;

    compute_intersection_points();
    for (k = 0; k < z - 1; k++)
	for (j = 0; j < y - 1; j++)
	    for (i = 0; i < x - 1; i++) {
		lut_entry = 0;
		for (int p = 0; p < 8; ++p) {
		    cube[p] =
			get_data(i + ((p ^ (p >> 1)) & 1),
				 j + ((p >> 1) & 1), k + ((p >> 2) & 1));
		    if (cube[p] > 0)
			lut_entry += 1 << p;
		}
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
		cube[0] = get_data(i, j, k);
		if (i < x - 1)
		    cube[1] = get_data(i + 1, j, k);
		else
		    cube[1] = cube[0];
		if (j < y - 1)
		    cube[3] = get_data(i, j + 1, k);
		else
		    cube[3] = cube[0];
		if (k < z - 1)
		    cube[4] = get_data(i, j, k + 1);
		else
		    cube[4] = cube[0];
		if (cube[0] < 0) {
		    if (cube[1] > 0)
			set_vert(X, i, j, k);
		    if (cube[3] > 0)
			set_vert(Y, i, j, k);
		    if (cube[4] > 0)
			set_vert(Z, i, j, k);
		} else {
		    if (cube[1] < 0)
			set_vert(X, i, j, k);
		    if (cube[3] < 0)
			set_vert(Y, i, j, k);
		    if (cube[4] < 0)
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
	A = cube[0];
	B = cube[4];
	C = cube[5];
	D = cube[1];
	break;
    case -2:
    case 2:
	A = cube[1];
	B = cube[5];
	C = cube[6];
	D = cube[2];
	break;
    case -3:
    case 3:
	A = cube[2];
	B = cube[6];
	C = cube[7];
	D = cube[3];
	break;
    case -4:
    case 4:
	A = cube[3];
	B = cube[7];
	C = cube[4];
	D = cube[0];
	break;
    case -5:
    case 5:
	A = cube[0];
	B = cube[3];
	C = cube[2];
	D = cube[1];
	break;
    case -6:
    case 6:
	A = cube[4];
	B = cube[7];
	C = cube[6];
	D = cube[5];
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
	a = (cube[4] - cube[0]) * (cube[6] - cube[2]) - (cube[7] -
							 cube[3]) *
	    (cube[5] - cube[1]);
	b = cube[2] * (cube[4] - cube[0]) + cube[0] * (cube[6] - cube[2])
	    - cube[1] * (cube[7] - cube[3]) - cube[3] * (cube[5] -
							 cube[1]);
	t = -b / (2 * a);
	if (t < 0 || t > 1)
	    return s > 0;
	At = cube[0] + (cube[4] - cube[0]) * t;
	Bt = cube[3] + (cube[7] - cube[3]) * t;
	Ct = cube[2] + (cube[6] - cube[2]) * t;
	Dt = cube[1] + (cube[5] - cube[1]) * t;
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
	    t = cube[0] / (cube[0] - cube[1]);
	    At = 0;
	    Bt = cube[3] + (cube[2] - cube[3]) * t;
	    Ct = cube[7] + (cube[6] - cube[7]) * t;
	    Dt = cube[4] + (cube[5] - cube[4]) * t;
	    break;
	case 1:
	    t = cube[1] / (cube[1] - cube[2]);
	    At = 0;
	    Bt = cube[0] + (cube[3] - cube[0]) * t;
	    Ct = cube[4] + (cube[7] - cube[4]) * t;
	    Dt = cube[5] + (cube[6] - cube[5]) * t;
	    break;
	case 2:
	    t = cube[2] / (cube[2] - cube[3]);
	    At = 0;
	    Bt = cube[1] + (cube[0] - cube[1]) * t;
	    Ct = cube[5] + (cube[4] - cube[5]) * t;
	    Dt = cube[6] + (cube[7] - cube[6]) * t;
	    break;
	case 3:
	    t = cube[3] / (cube[3] - cube[0]);
	    At = 0;
	    Bt = cube[2] + (cube[1] - cube[2]) * t;
	    Ct = cube[6] + (cube[5] - cube[6]) * t;
	    Dt = cube[7] + (cube[4] - cube[7]) * t;
	    break;
	case 4:
	    t = cube[4] / (cube[4] - cube[5]);
	    At = 0;
	    Bt = cube[7] + (cube[6] - cube[7]) * t;
	    Ct = cube[3] + (cube[2] - cube[3]) * t;
	    Dt = cube[0] + (cube[1] - cube[0]) * t;
	    break;
	case 5:
	    t = cube[5] / (cube[5] - cube[6]);
	    At = 0;
	    Bt = cube[4] + (cube[7] - cube[4]) * t;
	    Ct = cube[0] + (cube[3] - cube[0]) * t;
	    Dt = cube[1] + (cube[2] - cube[1]) * t;
	    break;
	case 6:
	    t = cube[6] / (cube[6] - cube[7]);
	    At = 0;
	    Bt = cube[5] + (cube[4] - cube[5]) * t;
	    Ct = cube[1] + (cube[0] - cube[1]) * t;
	    Dt = cube[2] + (cube[3] - cube[2]) * t;
	    break;
	case 7:
	    t = cube[7] / (cube[7] - cube[4]);
	    At = 0;
	    Bt = cube[6] + (cube[5] - cube[6]) * t;
	    Ct = cube[2] + (cube[1] - cube[2]) * t;
	    Dt = cube[3] + (cube[0] - cube[3]) * t;
	    break;
	case 8:
	    t = cube[0] / (cube[0] - cube[4]);
	    At = 0;
	    Bt = cube[3] + (cube[7] - cube[3]) * t;
	    Ct = cube[2] + (cube[6] - cube[2]) * t;
	    Dt = cube[1] + (cube[5] - cube[1]) * t;
	    break;
	case 9:
	    t = cube[1] / (cube[1] - cube[5]);
	    At = 0;
	    Bt = cube[0] + (cube[4] - cube[0]) * t;
	    Ct = cube[3] + (cube[7] - cube[3]) * t;
	    Dt = cube[2] + (cube[6] - cube[2]) * t;
	    break;
	case 10:
	    t = cube[2] / (cube[2] - cube[6]);
	    At = 0;
	    Bt = cube[1] + (cube[5] - cube[1]) * t;
	    Ct = cube[0] + (cube[4] - cube[0]) * t;
	    Dt = cube[3] + (cube[7] - cube[3]) * t;
	    break;
	case 11:
	    t = cube[3] / (cube[3] - cube[7]);
	    At = 0;
	    Bt = cube[2] + (cube[6] - cube[2]) * t;
	    Ct = cube[1] + (cube[5] - cube[1]) * t;
	    Dt = cube[0] + (cube[4] - cube[0]) * t;
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
	    tv[t % 3] = get_vert(X, i, j, k);
	    break;
	case 1:
	    tv[t % 3] = get_vert(Y, i + 1, j, k);
	    break;
	case 2:
	    tv[t % 3] = get_vert(X, i, j + 1, k);
	    break;
	case 3:
	    tv[t % 3] = get_vert(Y, i, j, k);
	    break;
	case 4:
	    tv[t % 3] = get_vert(X, i, j, k + 1);
	    break;
	case 5:
	    tv[t % 3] = get_vert(Y, i + 1, j, k + 1);
	    break;
	case 6:
	    tv[t % 3] = get_vert(X, i, j + 1, k + 1);
	    break;
	case 7:
	    tv[t % 3] = get_vert(Y, i, j, k + 1);
	    break;
	case 8:
	    tv[t % 3] = get_vert(Z, i, j, k);
	    break;
	case 9:
	    tv[t % 3] = get_vert(Z, i + 1, j, k);
	    break;
	case 10:
	    tv[t % 3] = get_vert(Z, i + 1, j + 1, k);
	    break;
	case 11:
	    tv[t % 3] = get_vert(Z, i, j + 1, k);
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
    int Off[][4] = {
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
    double u, rx, ry, rz;
    int g, m, *o;
    u = 0;
    test_vertex_addition();
    rx = ry = rz = 0;
    for (g = 0; g < SIZE(Off); g++) {
	o = Off[g];
	m = get_vert(o[3], i + o[X], j + o[Y], k + o[Z]);
	if (m != -1) {
	    u++;
	    rx += vertices[3*m + X];
	    ry += vertices[3*m + Y];
	    rz += vertices[3*m + Z];
	}
    }
    vertices[3*nv + X] = rx / u;
    vertices[3*nv + Y] = ry / u;
    vertices[3*nv + Z] = rz / u;
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
	printf("v %.16g %.16g %.16g\n", ver[3*i + X], ver[3*i + Y],
	       ver[3*i + Z]);
    for (i = 0; i < t; i++)
	printf("f %d %d %d\n", tri[3 * i + X] + 1, tri[3 * i + Y] + 1,
	       tri[3 * i + Z] + 1);
}
