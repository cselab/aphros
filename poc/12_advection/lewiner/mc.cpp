#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <float.h>
#include "mc.h"
#include "table.h"
#define ALLOC_SIZE 65536

double
 MarchingCubes::get_data(int i, int j, int k)
{
    return data[i + j * x + k * x * y];
}

int
 MarchingCubes::get_x_vert(int i, int j, int k)
{
    return x_verts[i + j * x + k * x * y];
}

int
 MarchingCubes::get_y_vert(int i, int j, int k)
{
    return y_verts[i + j * x + k * x * y];
}

int
 MarchingCubes::get_z_vert(int i, int j, int k)
{
    return z_verts[i + j * x + k * x * y];
}

void
 MarchingCubes::set_x_vert(int val, int i, int j, int k)
{
    x_verts[i + j * x + k * x * y] = val;
}

void
 MarchingCubes::set_y_vert(int val, int i, int j, int k)
{
    y_verts[i + j * x + k * x * y] = val;
}

void
 MarchingCubes::set_z_vert(int val, int i, int j, int k)
{
    z_verts[i + j * x + k * x * y] = val;
}


MarchingCubes::MarchingCubes(int x0, int y0, int z0, double *data0)
{
    int N, i;

    x = x0;
    y = y0;
    z = z0;
    data = data0;
    
    N = x * y * z;

    x_verts = new int[N];
    y_verts = new int[N];
    z_verts = new int[N];

    for (i = 0; i < N; i++) {
	x_verts[i] = -1;
	y_verts[i] = -1;
	z_verts[i] = -1;
    }
    nverts = ntrigs = 0;
    Nverts = Ntrigs = ALLOC_SIZE;
    vertices = new Vertex[Nverts];
    triangles = new Triangle[Ntrigs];
}

void
 MarchingCubes::run()
{
    compute_intersection_points();
    for (k = 0; k < z - 1; k++)
	for (j = 0; j < y - 1; j++)
	    for (i = 0; i < x - 1; i++) {
		lut_entry = 0;
		for (int p = 0; p < 8; ++p) {
		    cube[p] =
			get_data(i + ((p ^ (p >> 1)) & 1),
				 j + ((p >> 1) & 1), k + ((p >> 2) & 1));
		    if (fabs(cube[p]) < FLT_EPSILON)
			cube[p] = FLT_EPSILON;
		    if (cube[p] > 0)
			lut_entry += 1 << p;
		}
		process_cube();
	    }
}

void
 MarchingCubes::compute_intersection_points()
{
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
		if (fabs(cube[0]) < FLT_EPSILON)
		    cube[0] = FLT_EPSILON;
		if (fabs(cube[1]) < FLT_EPSILON)
		    cube[1] = FLT_EPSILON;
		if (fabs(cube[3]) < FLT_EPSILON)
		    cube[3] = FLT_EPSILON;
		if (fabs(cube[4]) < FLT_EPSILON)
		    cube[4] = FLT_EPSILON;
		if (cube[0] < 0) {
		    if (cube[1] > 0)
			set_x_vert(add_x_vertex(), i, j, k);
		    if (cube[3] > 0)
			set_y_vert(add_y_vertex(), i, j, k);
		    if (cube[4] > 0)
			set_z_vert(add_z_vertex(), i, j, k);
		} else {
		    if (cube[1] < 0)
			set_x_vert(add_x_vertex(), i, j, k);
		    if (cube[3] < 0)
			set_y_vert(add_y_vertex(), i, j, k);
		    if (cube[4] < 0)
			set_z_vert(add_z_vertex(), i, j, k);
		}
	    }
}

bool
 MarchingCubes::test_face(schar face)
{
    double
     A,
	B,
	C,
	D;

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

bool
 MarchingCubes::test_interior(schar s)
{
    double
     t,
	At = 0, Bt = 0, Ct = 0, Dt = 0, a, b;
    char
     test = 0;
    char
     edge = -1;

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
	return s > 0;
    case 1:
	return s > 0;
    case 2:
	return s > 0;
    case 3:
	return s > 0;
    case 4:
	return s > 0;
    case 5:
	if (At * Ct - Bt * Dt < FLT_EPSILON)
	    return s > 0;
	break;
    case 6:
	return s > 0;
    case 7:
	return s < 0;
    case 8:
	return s > 0;
    case 9:
	return s > 0;
    case 10:
	if (At * Ct - Bt * Dt >= FLT_EPSILON)
	    return s > 0;
	break;
    case 11:
	return s < 0;
    case 12:
	return s > 0;
    case 13:
	return s < 0;
    case 14:
	return s < 0;
    case 15:
	return s < 0;
    }
    return s < 0;
}

void
 MarchingCubes::process_cube()
{
    int
     v12 = -1;

    Case = cases[lut_entry][0];
    config = cases[lut_entry][1];
    subconfig = 0;
    switch (Case) {
    case 0:
	break;
    case 1:
	add_triangle(tiling1[config], 1);
	break;
    case 2:
	add_triangle(tiling2[config], 2);
	break;
    case 3:
	if (test_face(test3[config]))
	    add_triangle(tiling3_2[config], 4);
	else
	    add_triangle(tiling3_1[config], 2);
	break;
    case 4:
	if (test_interior(test4[config]))
	    add_triangle(tiling4_1[config], 2);
	else
	    add_triangle(tiling4_2[config], 6);
	break;
    case 5:
	add_triangle(tiling5[config], 3);
	break;
    case 6:
	if (test_face(test6[config][0]))
	    add_triangle(tiling6_2[config], 5);
	else {
	    if (test_interior(test6[config][1]))
		add_triangle(tiling6_1_1[config], 3);
	    else {
		v12 = add_c_vertex();
		add_triangle(tiling6_1_2[config], 9, v12);
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
	    add_triangle(tiling7_1[config], 3);
	    break;
	case 1:
	    add_triangle(tiling7_2[config][0], 5);
	    break;
	case 2:
	    add_triangle(tiling7_2[config][1], 5);
	    break;
	case 3:
	    v12 = add_c_vertex();
	    add_triangle(tiling7_3[config][0], 9, v12);
	    break;
	case 4:
	    add_triangle(tiling7_2[config][2], 5);
	    break;
	case 5:
	    v12 = add_c_vertex();
	    add_triangle(tiling7_3[config][1], 9, v12);
	    break;
	case 6:
	    v12 = add_c_vertex();
	    add_triangle(tiling7_3[config][2], 9, v12);
	    break;
	case 7:
	    if (test_interior(test7[config][3]))
		add_triangle(tiling7_4_2[config], 9);
	    else
		add_triangle(tiling7_4_1[config], 5);
	    break;
	};
	break;
    case 8:
	add_triangle(tiling8[config], 2);
	break;
    case 9:
	add_triangle(tiling9[config], 4);
	break;
    case 10:
	if (test_face(test10[config][0])) {
	    if (test_face(test10[config][1]))
		add_triangle(tiling10_1_1_[config], 4);
	    else {
		v12 = add_c_vertex();
		add_triangle(tiling10_2[config], 8, v12);
	    }
	} else {
	    if (test_face(test10[config][1])) {
		v12 = add_c_vertex();
		add_triangle(tiling10_2_[config], 8, v12);
	    } else {
		if (test_interior(test10[config][2]))
		    add_triangle(tiling10_1_1[config], 4);
		else
		    add_triangle(tiling10_1_2[config], 8);
	    }
	}
	break;
    case 11:
	add_triangle(tiling11[config], 4);
	break;
    case 12:
	if (test_face(test12[config][0])) {
	    if (test_face(test12[config][1]))
		add_triangle(tiling12_1_1_[config], 4);
	    else {
		v12 = add_c_vertex();
		add_triangle(tiling12_2[config], 8, v12);
	    }
	} else {
	    if (test_face(test12[config][1])) {
		v12 = add_c_vertex();
		add_triangle(tiling12_2_[config], 8, v12);
	    } else {
		if (test_interior(test12[config][2]))
		    add_triangle(tiling12_1_1[config], 4);
		else
		    add_triangle(tiling12_1_2[config], 8);
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
	    add_triangle(tiling13_1[config], 4);
	    break;
	case 1:
	    add_triangle(tiling13_2[config][0], 6);
	    break;
	case 2:
	    add_triangle(tiling13_2[config][1], 6);
	    break;
	case 3:
	    add_triangle(tiling13_2[config][2], 6);
	    break;
	case 4:
	    add_triangle(tiling13_2[config][3], 6);
	    break;
	case 5:
	    add_triangle(tiling13_2[config][4], 6);
	    break;
	case 6:
	    add_triangle(tiling13_2[config][5], 6);
	    break;
	case 7:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][0], 10, v12);
	    break;
	case 8:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][1], 10, v12);
	    break;
	case 9:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][2], 10, v12);
	    break;
	case 10:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][3], 10, v12);
	    break;
	case 11:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][4], 10, v12);
	    break;
	case 12:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][5], 10, v12);
	    break;
	case 13:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][6], 10, v12);
	    break;
	case 14:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][7], 10, v12);
	    break;
	case 15:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][8], 10, v12);
	    break;
	case 16:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][9], 10, v12);
	    break;
	case 17:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][10], 10, v12);
	    break;
	case 18:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3[config][11], 10, v12);
	    break;
	case 19:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_4[config][0], 12, v12);
	    break;
	case 20:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_4[config][1], 12, v12);
	    break;
	case 21:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_4[config][2], 12, v12);
	    break;
	case 22:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_4[config][3], 12, v12);
	    break;
	case 23:
	    subconfig = 0;
	    if (test_interior(test13[config][6]))
		add_triangle(tiling13_5_1[config][0], 6);
	    else
		add_triangle(tiling13_5_2[config][0], 10);
	    break;
	case 24:
	    subconfig = 1;
	    if (test_interior(test13[config][6]))
		add_triangle(tiling13_5_1[config][1], 6);
	    else
		add_triangle(tiling13_5_2[config][1], 10);
	    break;
	case 25:
	    subconfig = 2;
	    if (test_interior(test13[config][6]))
		add_triangle(tiling13_5_1[config][2], 6);
	    else
		add_triangle(tiling13_5_2[config][2], 10);
	    break;
	case 26:
	    subconfig = 3;
	    if (test_interior(test13[config][6]))
		add_triangle(tiling13_5_1[config][3], 6);
	    else
		add_triangle(tiling13_5_2[config][3], 10);
	    break;
	case 27:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][0], 10, v12);
	    break;
	case 28:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][1], 10, v12);
	    break;
	case 29:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][2], 10, v12);
	    break;
	case 30:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][3], 10, v12);
	    break;
	case 31:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][4], 10, v12);
	    break;
	case 32:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][5], 10, v12);
	    break;
	case 33:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][6], 10, v12);
	    break;
	case 34:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][7], 10, v12);
	    break;
	case 35:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][8], 10, v12);
	    break;
	case 36:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][9], 10, v12);
	    break;
	case 37:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][10], 10, v12);
	    break;
	case 38:
	    v12 = add_c_vertex();
	    add_triangle(tiling13_3_[config][11], 10, v12);
	    break;
	case 39:
	    add_triangle(tiling13_2_[config][0], 6);
	    break;
	case 40:
	    add_triangle(tiling13_2_[config][1], 6);
	    break;
	case 41:
	    add_triangle(tiling13_2_[config][2], 6);
	    break;
	case 42:
	    add_triangle(tiling13_2_[config][3], 6);
	    break;
	case 43:
	    add_triangle(tiling13_2_[config][4], 6);
	    break;
	case 44:
	    add_triangle(tiling13_2_[config][5], 6);
	    break;
	case 45:
	    add_triangle(tiling13_1_[config], 4);
	    break;
	default:
	    printf("Marching Cubes: Impossible case 13?\n");
	}
	break;
    case 14:
	add_triangle(tiling14[config], 4);
	break;
    };
}

void
 MarchingCubes::add_triangle(const char *trig, char n, int v12)
{
    int
     tv[3];

    for (int t = 0; t < 3 * n; t++) {
	switch (trig[t]) {
	case 0:
	    tv[t % 3] = get_x_vert(i, j, k);
	    break;
	case 1:
	    tv[t % 3] = get_y_vert(i + 1, j, k);
	    break;
	case 2:
	    tv[t % 3] = get_x_vert(i, j + 1, k);
	    break;
	case 3:
	    tv[t % 3] = get_y_vert(i, j, k);
	    break;
	case 4:
	    tv[t % 3] = get_x_vert(i, j, k + 1);
	    break;
	case 5:
	    tv[t % 3] = get_y_vert(i + 1, j, k + 1);
	    break;
	case 6:
	    tv[t % 3] = get_x_vert(i, j + 1, k + 1);
	    break;
	case 7:
	    tv[t % 3] = get_y_vert(i, j, k + 1);
	    break;
	case 8:
	    tv[t % 3] = get_z_vert(i, j, k);
	    break;
	case 9:
	    tv[t % 3] = get_z_vert(i + 1, j, k);
	    break;
	case 10:
	    tv[t % 3] = get_z_vert(i + 1, j + 1, k);
	    break;
	case 11:
	    tv[t % 3] = get_z_vert(i, j + 1, k);
	    break;
	case 12:
	    tv[t % 3] = v12;
	    break;
	default:
	    break;
	}
	if (tv[t % 3] == -1) {
	    printf("Marching Cubes: invalid triangle %d\n", ntrigs + 1);
	}
	if (t % 3 == 2) {
	    if (ntrigs >= Ntrigs) {
		Triangle *
		    temp = triangles;

		triangles = new Triangle[2 * Ntrigs];
		memcpy(triangles, temp, Ntrigs * sizeof(Triangle));
		delete[]temp;
		printf("%d allocated triangles\n", Ntrigs);
		Ntrigs *= 2;
	    }
	    Triangle *
		T = triangles + ntrigs++;

	    T->v1 = tv[0];
	    T->v2 = tv[1];
	    T->v3 = tv[2];
	}
    }
}

void
 MarchingCubes::test_vertex_addition()
{
    if (nverts >= Nverts) {
	Vertex *
	    temp = vertices;

	vertices = new Vertex[Nverts * 2];
	memcpy(vertices, temp, Nverts * sizeof(Vertex));
	delete[]temp;
	printf("%d allocated vertices\n", Nverts);
	Nverts *= 2;
    }
}

int
 MarchingCubes::add_x_vertex()
{
    test_vertex_addition();
    Vertex *
	vert = vertices + nverts++;

    double
     u = (cube[0]) / (cube[0] - cube[1]);

    vert->x = (double) i + u;
    vert->y = (double) j;
    vert->z = (double) k;
    return nverts - 1;
}

int
 MarchingCubes::add_y_vertex()
{
    test_vertex_addition();
    Vertex *
	vert = vertices + nverts++;

    double
     u = (cube[0]) / (cube[0] - cube[3]);

    vert->x = (double) i;
    vert->y = (double) j + u;
    vert->z = (double) k;
    return nverts - 1;
}

int
 MarchingCubes::add_z_vertex()
{
    test_vertex_addition();
    Vertex *
	vert = vertices + nverts++;

    double
     u = (cube[0]) / (cube[0] - cube[4]);

    vert->x = (double) i;
    vert->y = (double) j;
    vert->z = (double) k + u;
    return nverts - 1;
}

int
 MarchingCubes::add_c_vertex()
{
    test_vertex_addition();
    Vertex *
	vert = vertices + nverts++;

    double
     u = 0;
    int
     vid;

    vert->x = vert->y = vert->z = 0;
    vid = get_x_vert(i, j, k);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_y_vert(i + 1, j, k);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_x_vert(i, j + 1, k);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_y_vert(i, j, k);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_x_vert(i, j, k + 1);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_y_vert(i + 1, j, k + 1);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_x_vert(i, j + 1, k + 1);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_y_vert(i, j, k + 1);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_z_vert(i, j, k);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_z_vert(i + 1, j, k);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_z_vert(i + 1, j + 1, k);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vid = get_z_vert(i, j + 1, k);
    if (vid != -1) {
	++u;
	const
	 Vertex &
	    v = vertices[vid];

	vert->x += v.x;
	vert->y += v.y;
	vert->z += v.z;
    }
    vert->x /= u;
    vert->y /= u;
    vert->z /= u;
    return nverts - 1;
}

void
writeObj(MarchingCubes * q)
{
    int
     t;
    int
     v;
    int
     i;
    Triangle *
	tri;
    Vertex *
	ver;

    t = q->ntrigs;
    v = q->nverts;
    tri = q->triangles;
    ver = q->vertices;
    printf("# File type: ASCII OBJ\n");
    for (i = 0; i < v; i++)
	printf("v %.16g %.16g %.16g\n", ver[i].x, ver[i].y, ver[i].z);
    for (i = 0; i < t; i++)
	printf("f %d %d %d\n", tri[i].v1 + 1, tri[i].v2 + 1,
	       tri[i].v3 + 1);
}
