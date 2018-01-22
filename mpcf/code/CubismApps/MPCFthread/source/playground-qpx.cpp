
/*
 *		
 *		 *  TryThis
 *		  *
 *		   *  Created by Diego Rossinelli on 1/24/13.
 *		    *  Copyright 2013 ETH Zurich. All rights reserved.
 *		     *
 *		      */
#include <cstdio>

int main()
{
	printf("asd!\n");
	
	vector4double a= (vector4double)(0.,1.,2.,3.);
	vector4double b= (vector4double)(4.,5.,6.,7.);
	vector4double c= (vector4double)(8.,9.,10.,11.);
	vector4double d= (vector4double)(12.,13.,14.,15.);

	float data[4][4];
	
	vec_st(a, 0, &data[0][0]);
	vec_st(b, 0, &data[1][0]);
	vec_st(c, 0, &data[2][0]);
	vec_st(d, 0, &data[3][0]);
	
	printf("a is %f %f %f %f\n", data[0][0], data[0][1], data[0][2], data[0][3]);
	printf("b is %f %f %f %f\n", data[1][0], data[1][1], data[1][2], data[1][3]);
	printf("c is %f %f %f %f\n", data[2][0], data[2][1], data[2][2], data[2][3]);
	printf("d is %f %f %f %f\n", data[3][0], data[3][1], data[3][2], data[3][3]);
	
	const vector4double v01L = vec_perm(a, b, vec_gpci(00415));
	const vector4double v01H = vec_perm(a, b, vec_gpci(02637));
	const vector4double v23L = vec_perm(c, d, vec_gpci(00415));
	const vector4double v23H = vec_perm(c, d, vec_gpci(02637));
	
	a = vec_perm(v01L, v23L, vec_gpci(00145));
	b = vec_perm(v01L, v23L, vec_gpci(02367));
	c = vec_perm(v01H, v23H, vec_gpci(00145));
	d = vec_perm(v01H, v23H, vec_gpci(02367));
	
	vec_st(a, 0, &data[0][0]);
	vec_st(b, 0, &data[1][0]);
	vec_st(c, 0, &data[2][0]);
	vec_st(d, 0, &data[3][0]);

	printf("\n\ntransposition now!!\n\n");
	
	printf("a is %f %f %f %f\n", data[0][0], data[0][1], data[0][2], data[0][3]);
	printf("b is %f %f %f %f\n", data[1][0], data[1][1], data[1][2], data[1][3]);
	printf("c is %f %f %f %f\n", data[2][0], data[2][1], data[2][2], data[2][3]);
	printf("d is %f %f %f %f\n", data[3][0], data[3][1], data[3][2], data[3][3]);
	
	return 0;
}
