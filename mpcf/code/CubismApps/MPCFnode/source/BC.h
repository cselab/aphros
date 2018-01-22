/*
 *  BoundaryConditions.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/16/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "BlockLab.h"
#include <Matrix3D.h>

template<typename TBlock, typename TElement, template<typename X> class allocator=std::allocator>
class BoundaryCondition
{
protected:
	
	int s[3], e[3];
	int stencilStart[3], stencilEnd[3];	
	Matrix3D<TElement, true, allocator> * cacheBlock;	
    
	template<int dir, int side>
	void _setup()
	{
		s[0] =	dir==0? (side==0? stencilStart[0]: TBlock::sizeX) : 0;
		s[1] =	dir==1? (side==0? stencilStart[1]: TBlock::sizeY) : 0;
		s[2] =	dir==2? (side==0? stencilStart[2]: TBlock::sizeZ) : 0;
		
		e[0] =	dir==0? (side==0? 0: TBlock::sizeX + stencilEnd[0]-1) : TBlock::sizeX;
		e[1] =	dir==1? (side==0? 0: TBlock::sizeY + stencilEnd[1]-1) : TBlock::sizeY;
		e[2] =	dir==2? (side==0? 0: TBlock::sizeZ + stencilEnd[2]-1) : TBlock::sizeZ;
	}
    
    Real _pulse(const Real t_star, const Real p_ratio)
    {
        const Real Pa = p_ratio*2.38*10;//target peak pressure (factor 2.38 if no velocity)
        const Real omega = 2*M_PI*0.5/6e-6;//tensile part set to be 6microseconds (6e-6)
        const Real rise  = 1.03*(1-exp(-9.21e7*t_star));//50ns rise time
        const Real alpha = 9.1e5;
        
        Real p = rise*2*Pa*exp(-alpha*t_star)*cos(omega*t_star + M_PI/3.);

        return p;
    }
    
public:
	
	BoundaryCondition(const int ss[3], const int se[3], Matrix3D<TElement, true, allocator> * cacheBlock): 
    cacheBlock(cacheBlock)
	{
		s[0]=s[1]=s[2]=0;
		e[0]=e[1]=e[2]=0;
		
		stencilStart[0] = ss[0];
		stencilStart[1] = ss[1];
		stencilStart[2] = ss[2];
		
		stencilEnd[0] = se[0];
		stencilEnd[1] = se[1];
		stencilEnd[2] = se[2];
	}
	
	TElement& operator()(int ix, int iy, int iz)
	{
		return cacheBlock->Access(ix-stencilStart[0],iy-stencilStart[1],iz-stencilStart[2]);
	}
	
	template<int dir, int side>
	void applyBC_absorbing()
	{
		_setup<dir,side>();
        
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
					(*this)(ix,iy,iz) = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1):ix,
												dir==1? (side==0? 0:TBlock::sizeY-1):iy,
												dir==2? (side==0? 0:TBlock::sizeZ-1):iz);
				}
	}
    
    template<int dir, int side>
	void applyBC_absorbing_better_faces()
	{
		_setup<dir,side>();

        //Zeroth order extrapolation for faces. 
        //This fills up ALL the ghost values although corners and edges will not be right. 
        //Corners and edges will be overwritten on request by void applyBC_absorbing_better_tensorials<int,int>().
        
#if 0
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
					(*this)(ix,iy,iz) = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1):ix,
												dir==1? (side==0? 0:TBlock::sizeY-1):iy,
												dir==2? (side==0? 0:TBlock::sizeZ-1):iz);
				}
#else
//		printf("s=[%d,%d,%d] e=[%d,%d,%d]\n", s[0], s[1], s[2], e[0], e[1], e[2]);
		if (dir == 0)
		{
			int px;
			if (side == 0)	px = 0;
			else			px = TBlock::sizeX-1;
			
			for(int iz=s[2]; iz<e[2]; iz++)
				for(int iy=s[1]; iy<e[1]; iy++)
				{
					#if 1
					//TElement tmp =  (*this)(px, iy, iz);
					TElement * ptrSrc = &cacheBlock->Access(px-stencilStart[0],iy-stencilStart[1],iz-stencilStart[2]);
					TElement * ptrDest = &cacheBlock->Access(s[0]-stencilStart[0],iy-stencilStart[1],iz-stencilStart[2]);
					for (int i = 0; i < e[0]-s[0]; i++)
						//ptrDest[i] = tmp;
						memcpy(ptrDest + i, ptrSrc, sizeof(TElement));
					#else
					TElement tmp =  (*this)(px, iy, iz);
					for(int ix=s[0]; ix<e[0]; ix++)
					{
						(*this)(ix,iy,iz) = tmp;	//(*this)(px, iy, iz);
					}
					#endif
				}
		}
		else if (dir == 1)
		{
			int py;
			if (side == 0)	py = 0;
			else			py = TBlock::sizeY-1;

			int nbytes = (e[0]-s[0])*sizeof(TElement);
			const int ix = s[0];

			for(int iz=s[2]; iz<e[2]; iz++)
			{
				TElement * ptrSrc = &cacheBlock->Access(ix-stencilStart[0],py-stencilStart[1],iz-stencilStart[2]);
				for(int iy=s[1]; iy<e[1]; iy++)
				{
					TElement * ptrDest = &cacheBlock->Access(ix-stencilStart[0],iy-stencilStart[1],iz-stencilStart[2]);
					memcpy((char *)ptrDest, (char *)ptrSrc, nbytes);

					//for(int i=0; i<e[0]-s[0]; i++)
					//{
					//	ptrDest[i] = ptrSrc[i];
					//}


					//for(int ix=s[0]; ix<e[0]; ix++)
					//{
					//	(*this)(ix,iy,iz) = (*this)(ix, py, iz);
					//}
				}
			}
		}
		else /* (dir == 2) */
		{
			int pz;
			if (side == 0)	pz = 0;
			else			pz = TBlock::sizeZ-1;

			int nbytes = (e[0]-s[0])*sizeof(TElement);
			const int ix = s[0];
			for(int iz=s[2]; iz<e[2]; iz++)
				for(int iy=s[1]; iy<e[1]; iy++)
				{
					TElement * ptrSrc = &cacheBlock->Access(ix-stencilStart[0],iy-stencilStart[1],pz-stencilStart[2]);
					TElement * ptrDest = &cacheBlock->Access(ix-stencilStart[0],iy-stencilStart[1],iz-stencilStart[2]);
					memcpy((char *)ptrDest, (char *)ptrSrc, nbytes);
					
					//for(int i=0; i<e[0]-s[0]; i++)
					//{
					//	ptrDest[i] = ptrSrc[i];
					//}

					
					//for(int ix=s[0]; ix<e[0]; ix++)
					//{
					//	(*this)(ix,iy,iz) = (*this)(ix, iy, pz);
					//}
				}
		}

#endif
	}
    
	void applyBC_absorbing_better_tensorials_edges()
    {
        const int bsize[3] = {TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ};
        int s[3], e[3];
        
        //Edges
        {           
            for(int d=0; d<3; ++d)
                for(int b=0; b<2; ++b)
                    for(int a=0; a<2; ++a)
                    {
                        const int d1 = (d + 1) % 3;
                        const int d2 = (d + 2) % 3;
                        
                        s[d]  =	stencilStart[d];
                        s[d1] =	a*(bsize[d1]-stencilStart[d1])+stencilStart[d1];
                        s[d2] =	b*(bsize[d2]-stencilStart[d2])+stencilStart[d2];
                        
                        e[d]  =	bsize[d]-1+stencilEnd[d];
                        e[d1] =	a*(bsize[d1]-1+stencilEnd[d1]);
                        e[d2] =	b*(bsize[d2]-1+stencilEnd[d2]);
                        
                        for(int iz=s[2]; iz<e[2]; iz++)
                            for(int iy=s[1]; iy<e[1]; iy++)
                                for(int ix=s[0]; ix<e[0]; ix++)
                                {
                                    (*this)(ix,iy,iz) = d==0? (*this)(ix,a*(bsize[1]-1),b*(bsize[2]-1)) : (d==1? (*this)(a*(bsize[0]-1),iy,b*(bsize[2]-1)) : (*this)(a*(bsize[0]-1),b*(bsize[1]-1),iz));
                                }
                    }         
        }
    }
    
    void applyBC_absorbing_better_tensorials_corners()
    {
        const int bsize[3] = {TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ};
        int s[3], e[3];
        
        //Corners
        {
            for(int c=0; c<2; ++c)
                for(int b=0; b<2; ++b)
                    for(int a=0; a<2; ++a)
                    {                        
                        s[0]  =	a*(bsize[0]-stencilStart[0])+stencilStart[0];
                        s[1] =	b*(bsize[1]-stencilStart[1])+stencilStart[1];
                        s[2] =	c*(bsize[2]-stencilStart[2])+stencilStart[2];
                        
                        e[0]  =	a*(bsize[0]-1+stencilEnd[0]);
                        e[1] =	b*(bsize[1]-1+stencilEnd[1]);
                        e[2] =	c*(bsize[2]-1+stencilEnd[2]);
                        
                        for(int iz=s[2]; iz<e[2]; iz++)
                            for(int iy=s[1]; iy<e[1]; iy++)
                                for(int ix=s[0]; ix<e[0]; ix++)
                                {
                                    (*this)(ix,iy,iz) = (*this)(a*(bsize[0]-1),b*(bsize[1]-1),c*(bsize[2]-1));
                                }
                    }
        }
    }
    
	template<int dir, int side>
	void applyBC_reflecting()
	{
		_setup<dir,side>();
		
#if 0	
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
					TElement source = (*this)(dir==0? (side==0? -ix-1:2*TBlock::sizeX-ix-1):ix,
											  dir==1? (side==0? -iy-1:2*TBlock::sizeY-iy-1):iy,
											  dir==2? (side==0? -iz-1:2*TBlock::sizeZ-iz-1):iz);
					
					(*this)(ix,iy,iz).rho      = source.rho;
					(*this)(ix,iy,iz).u        = ((dir==0)? -1:1)*source.u;
					(*this)(ix,iy,iz).v        = ((dir==1)? -1:1)*source.v;
					(*this)(ix,iy,iz).w        = ((dir==2)? -1:1)*source.w;
					(*this)(ix,iy,iz).energy   = source.energy;
					(*this)(ix,iy,iz).G = source.G;
                    (*this)(ix,iy,iz).P = source.P;
				}
#else
/*
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
					TElement source = (*this)(dir==0? (side==0? -ix-1:2*TBlock::sizeX-ix-1):ix,
											  dir==1? (side==0? -iy-1:2*TBlock::sizeY-iy-1):iy,
											  dir==2? (side==0? -iz-1:2*TBlock::sizeZ-iz-1):iz);
					
                    (*this)(ix,iy,iz) = source;
                    if  (dir==0) (*this)(ix,iy,iz).u = -(*this)(ix,iy,iz).u;
                    if  (dir==1) (*this)(ix,iy,iz).v = -(*this)(ix,iy,iz).v;
                    if  (dir==2) (*this)(ix,iy,iz).w = -(*this)(ix,iy,iz).w;
				}
*/

		if (dir == 0)
		{
			int px;
			if (side == 0)	px = -1; //-ix;
			else			px = 2*TBlock::sizeX-1; //-ix;
			
			for(int iz=s[2]; iz<e[2]; iz++)
				for(int iy=s[1]; iy<e[1]; iy++)
					for(int ix=s[0]; ix<e[0]; ix++)
					{
						TElement source = (*this)(px-ix, iy, iz);
						source.u = -source.u;
						(*this)(ix,iy,iz) = source;
					}
		}
		else if (dir == 1)
		{
			int py;
			if (side == 0)	py = -1; //-iy;
			else			py = 2*TBlock::sizeY-1; //-iy;

			for(int iz=s[2]; iz<e[2]; iz++)
				for(int iy=s[1]; iy<e[1]; iy++)
				{
					for(int ix=s[0]; ix<e[0]; ix++)
					{
						TElement source = (*this)(ix, py-iy, iz);
						source.v = -source.v;
						(*this)(ix,iy,iz) = source;
					}
				}
		}
		else /* (dir == 2) */
		{
			int pz;
			if (side == 0)	pz = -1; //-iz;
			else			pz = 2*TBlock::sizeZ-1; //-iz;

			for(int iz=s[2]; iz<e[2]; iz++)
				for(int iy=s[1]; iy<e[1]; iy++)
					for(int ix=s[0]; ix<e[0]; ix++)
					{
						TElement source = (*this)(ix, iy, pz-iz);
						source.w = -source.w;
						(*this)(ix,iy,iz) = source;
					}

		}

#endif
	}
	
	template<int dir, int side>
	void applyBC_dirichlet(const TElement& p)
	{
		_setup<dir,side>();
		
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
					(*this)(ix,iy,iz).rho      = p.rho;
					(*this)(ix,iy,iz).u        = p.u;
					(*this)(ix,iy,iz).v        = p.v;
					(*this)(ix,iy,iz).w        = p.w;
					(*this)(ix,iy,iz).energy   = p.energy;
					(*this)(ix,iy,iz).G        = p.G;
                    (*this)(ix,iy,iz).P        = p.P;
				}
	}
    
    template<int dir, int side>
	void applyBC_spaceDirichlet(const TElement& p, const Real t, const Real h)
	{
		_setup<dir,side>();
		
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
                    Real pc = 0;
					(*this)(ix,iy,iz).rho      = p.rho;
					(*this)(ix,iy,iz).u        = p.u;
					(*this)(ix,iy,iz).v        = p.v;
					(*this)(ix,iy,iz).w        = p.w;
					(*this)(ix,iy,iz).G        = p.G;
                    pc = p.P;
                    (*this)(ix,iy,iz).P        = p.P;
                    /// time is scaled by tc = R0/c_L, R0=0.1 equivalent to 50microns, domain size: 500microns
                    Real c_L                 = sqrt((1/p.G+1)*(pc+10)/p.rho);//speed of sound
                    Real p_ratio             = 1000;
                    Real radius              = 50e-6;
                    Real time_scale          = 0.1/c_L*1500/radius; //time scale between us and real dimension
                    Real t_off               = (ix<0? abs(ix): ix-TBlock::sizeX+1)*h/c_L;//space offset
                    //printf("t: %e, t_off: %e, c_L: %f\n",t,t_off,c_L*10);
                    //(*this)(ix,iy,iz).energy   = _pulse((t+t_off)/time_scale,p_ratio)*p.G+pc;
                    (*this)(ix,iy,iz).energy = p.energy;
				}
	}
};
