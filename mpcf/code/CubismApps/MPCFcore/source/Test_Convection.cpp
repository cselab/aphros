/*
 *  Test_Convection.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cstring>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>

using namespace std;

#include "TestTypes.h"
#include "Test_Convection.h"

void Test_Convection::_initialize_lab(TestLab& lab)
{
#if 0
	const double Simulation_Environment__GAMMA1 = 1.4, Simulation_Environment__GAMMA2 = 1.4;
	const double Simulation_Environment__PC1 = 0, Simulation_Environment__PC2 = 0;

	const double G1 = Simulation_Environment__GAMMA1-1;
	const double G2 = Simulation_Environment__GAMMA2-1;
	const double F1 = Simulation_Environment__GAMMA1*Simulation_Environment__PC1;

 	const double F2 = Simulation_Environment__GAMMA2*Simulation_Environment__PC2;

	const double bubble = 0;
	const double shock = 0;

	for(int iz = -3; iz<_BLOCKSIZE_+3; iz++)
		for(int iy = -3; iy<_BLOCKSIZE_+3; iy++)
		{
			for(int ix = -3; ix<_BLOCKSIZE_+3; ix++)
			{
				lab(ix, iy. iz).clear();

				lab(ix, iy, iz).s.r  = 1;//shock*post_shock[0] + (1-shock)*(0.01*bubble+pre_shock[0]*(1-bubble));
				lab(ix, iy, iz).s.u  = 0;//post_shock[0]*post_shock[1]*shock;
				lab(ix, iy, iz).s.v  = 0;
				lab(ix, iy, iz).s.w  = 0;

				const double pressure = 1;//post_shock[2]*shock + pre_shock[2]*(1-shock);

				const double mix_gamma = 1 + (G2*G1)/(G1*bubble+G2*(1-bubble));
				const double mix_pinf  = (mix_gamma-1)/mix_gamma * (F1/G1*(1-bubble) + F2/G2*bubble);
				lab(ix, iy, iz).s.G = 1./(mix_gamma-1);
				lab(ix, iy, iz).s.P = mix_gamma*mix_pinf/(mix_gamma-1);

				const double ke = 0.5*(pow(lab(ix, iy, iz).s.u,2)+pow(lab(ix, iy, iz).s.v,2)+pow(lab(ix, iy, iz).s.w,2))/lab(ix, iy, iz).s.r;
				lab(ix, iy, iz).s.s = pressure*lab(ix, iy, iz).s.G + lab(ix, iy, iz).s.P + ke;
			}
		}
#else
        for(int iz = -3; iz<_BLOCKSIZE_+3; iz++)
                for(int iy = -3; iy<_BLOCKSIZE_+3; iy++)
                {
                        for(int ix = -3; ix<_BLOCKSIZE_+3; ix++)
                        {
				lab(ix, iy, iz).clear();

                                const int a = iy + 3;
                                const int b = iz + 3;
                                const int c = ix + 3;
                                const double L = _BLOCKSIZE_;

                                /* lab(ix, iy, iz).s.r = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L; */
                                /* lab(ix, iy, iz).s.u = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L/(double)L; */
                                /* lab(ix, iy, iz).s.v = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L/(double)L; */
                                /* lab(ix, iy, iz).s.w = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L/(double)L; */
                                /* lab(ix, iy, iz).s.s = 100+b/(double)L; */
                                /* lab(ix, iy, iz).s.G = (1+c+b*L + a*L*L)/(double)L; */
                                /* lab(ix, iy, iz).s.P = (1+c+b*L + a*L*L)/(double)L; */

                                lab(ix, iy, iz).s.r = 1.0;
                                lab(ix, iy, iz).s.u = 1.0;
                                lab(ix, iy, iz).s.v = 0.1;
                                lab(ix, iy, iz).s.w = -0.4;
                                lab(ix, iy, iz).s.s = 0.0;
                                lab(ix, iy, iz).s.G = 100.0;
                                lab(ix, iy, iz).s.P = 0.0;

                        }
                }

#endif
}

void Test_Convection::_initialize_block(Block& block)
{
	srand48(61651);

	for(int iz = 0; iz<_BLOCKSIZE_; iz++)
		for(int iy = 0; iy<_BLOCKSIZE_; iy++)
			for(int ix = 0; ix<_BLOCKSIZE_; ix++)
			{
                block(ix, iy, iz).clear();

                /* block(ix, iy, iz).dsdt.r = 0.1+iz; */
                /* block(ix, iy, iz).dsdt.u = ix; */
                /* block(ix, iy, iz).dsdt.v = iy; */
                /* block(ix, iy, iz).dsdt.w = (ix*iy/(double)_BLOCKSIZE_); */
                /* block(ix, iy, iz).dsdt.s = 1+(iy*iz)/(double)_BLOCKSIZE_; */
                /* block(ix, iy, iz).dsdt.G = -1 +  (iz+ix+iy)/(double)_BLOCKSIZE_; */
                /* block(ix, iy, iz).dsdt.P = -1 +  (iz+ix+iy)/(double)_BLOCKSIZE_; */

                block(ix, iy, iz).dsdt.r = 1.0;
                block(ix, iy, iz).dsdt.u = 1.0;
                block(ix, iy, iz).dsdt.v = 0.1;
                block(ix, iy, iz).dsdt.w = -0.4;
                block(ix, iy, iz).dsdt.s = 0.0;
                block(ix, iy, iz).dsdt.G = 100.0;
                block(ix, iy, iz).dsdt.P = 0.0;

            }
}

void Test_Convection::_print(Block& block)
{
	for(int iz = 0; iz<_BLOCKSIZE_; iz++)
		for(int iy = 0; iy<_BLOCKSIZE_; iy++)
			for(int ix = 0; ix<_BLOCKSIZE_; ix++)
			{
				if (((iz == ix) && (iz == iy))) {
					StateVector v = block(ix, iy, iz).dsdt;
					Real s = v.r + v.u + v.v + v.w + v.s + v.G + v.P;
					//if (isnan(s)) {
						//printf("(%d, %d, %d):", iz, iy, ix);
                		printf("%d %15.15e %15.15e %15.15e %15.15e %15.15e %15.15e %15.15e\n", ix+ iy*_BLOCKSIZE_ + iz*_BLOCKSIZE_*_BLOCKSIZE_, v.r, v.u, v.v, v.w, v.s, v.G, v.P);
					//}
				}
			}
}

/*
void Test_Convection::_print_lab(TestLab& lab)
{
	for(int iz = -3; iz<_BLOCKSIZE_+3; iz++)
		for(int iy = -3; iy<_BLOCKSIZE_+3; iy++)
			for(int ix = -3; ix<_BLOCKSIZE_+3; ix++)
			{
				if (((iz == ix) && (iz == iy))) {
					StateVector v = lab(ix, iy, iz).s;
					Real s = v.r + v.u + v.v + v.w + v.s + v.G + v.P;
					//if (isnan(s)) {
						//printf("(%d, %d, %d):", iz, iy, ix);
                		//printf("%d %15.15e %15.15e %15.15e %15.15e %15.15e %15.15e %15.15e\n", ix+ iy*_BLOCKSIZE_ + iz*_BLOCKSIZE_*_BLOCKSIZE_, v.r, v.u, v.v, v.w, v.s, v.G, v.P);
						printf("(%d,%d,%d) %15.15e %15.15e %15.15e %15.15e %15.15e %15.15e %15.15e\n", iz, iy, ix, v.r, v.u, v.v, v.w, v.s, v.G, v.P);
					//}
				}
			}
}*/
