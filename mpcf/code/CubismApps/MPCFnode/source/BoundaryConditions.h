/*
 *  BoundaryConditions.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/16/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <cassert>
#include <vector>
#include "BlockInfo.h"
#include "BlockLab.h"
#include "Types.h"
#include "BoundaryGhostCells.h"
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
		const int m_vSize0 = cacheBlock->getSize(0); //m_vSize[0];
		const int m_nElemsPerSlice = cacheBlock->getNumberOfElementsPerSlice(); //m_nElementsPerSlice;

		if (dir == 0)
		{
//			printf("dir0: from %d to %d\n", s[0], e[0]);

			int px;
			if (side == 0)	px = 0;
			else			px = TBlock::sizeX-1;

			for(int iz=s[2]; iz<e[2]; iz++)
			{
				const int my_zx_src = (iz-stencilStart[2])*m_nElemsPerSlice + px-stencilStart[0];
				const int my_zx_dst = (iz-stencilStart[2])*m_nElemsPerSlice + s[0]-stencilStart[0];

				for(int iy=s[1]; iy<e[1]; iy++)
				{
					const int my_zyx_src = (iy-stencilStart[1])*m_vSize0 + my_zx_src;
					const int my_zyx_dst = (iy-stencilStart[1])*m_vSize0 + my_zx_dst;

					TElement * ptrSrc = &cacheBlock->LinAccess(my_zyx_src);
					TElement * ptrDest = &cacheBlock->LinAccess(my_zyx_dst);

					#if 0
					for (int i = 0; i < e[0]-s[0]; i++)
						memcpy2((char *)(ptrDest + i), (char *)ptrSrc, sizeof(TElement));
					#else
					if ((e[0]-s[0]) == 3)
					{
						memcpy2((char *)(ptrDest + 0), (char *)ptrSrc, sizeof(TElement));
						memcpy2((char *)(ptrDest + 1), (char *)ptrSrc, sizeof(TElement));
						memcpy2((char *)(ptrDest + 2), (char *)ptrSrc, sizeof(TElement));
					}
					else
					{
						printf("does this happen? if yes, just remove this line\n"); abort();
						if ((e[0]-s[0])%4 != 0)
						{
							for (int i = 0; i < e[0]-s[0]; i++)
								memcpy2((char *)(ptrDest + i), (char *)ptrSrc, sizeof(TElement));
						}
						else
						{
							for (int i = 0; i < e[0]-s[0]; i+=4)
							{
								memcpy2((char *)(ptrDest + i + 0), (char *)ptrSrc, sizeof(TElement));
								memcpy2((char *)(ptrDest + i + 1), (char *)ptrSrc, sizeof(TElement));
								memcpy2((char *)(ptrDest + i + 2), (char *)ptrSrc, sizeof(TElement));
								memcpy2((char *)(ptrDest + i + 3), (char *)ptrSrc, sizeof(TElement));
							}
						}
					}
					#endif
				}
			}
		}
		else if (dir == 1)
		{
			int py;
			if (side == 0)	py = 0;
			else			py = TBlock::sizeY-1;

			int nbytes = (e[0]-s[0])*sizeof(TElement);
			const int ix = s[0];

//			printf("dir1: from %d to %d\n", s[1], e[1]);

			const int my_y_src = (py-stencilStart[1])*m_vSize0;
			for(int iz=s[2]; iz<e[2]; iz++)
			{
				const int my_xz = (iz-stencilStart[2])*m_nElemsPerSlice + ix-stencilStart[0];

				//TElement * ptrSrc = &cacheBlock->Access(ix-stencilStart[0],py-stencilStart[1],iz-stencilStart[2]);
				TElement * ptrSrc = &cacheBlock->LinAccess(my_xz + my_y_src);


				#if 0
				for(int iy=s[1]; iy<e[1]; iy++)
				{
					const int my_y_dst = (iy-stencilStart[1])*m_vSize0;

					TElement * ptrDest = &cacheBlock->LinAccess(my_xz + my_y_dst);
					memcpy2((char *)ptrDest, (char *)ptrSrc, nbytes);
				}
				#else
				if ((e[1]-s[1]) == 3)
				{
					const int iy = s[1];
					const int my_y_dst0 = (iy-stencilStart[1])*m_vSize0;
					//const int my_y_dst1 = my_y_dst0 + m_vSize0;
					//const int my_y_dst2 = my_y_dst1 + m_vSize0;

					TElement * ptrDest0 = &cacheBlock->LinAccess(my_xz + my_y_dst0);
//					TElement * ptrDest1 = &cacheBlock->LinAccess(my_xz + my_y_dst1);
//					TElement * ptrDest2 = &cacheBlock->LinAccess(my_xz + my_y_dst2);
					TElement * ptrDest1 = ptrDest0 + m_vSize0;
					TElement * ptrDest2 = ptrDest1 + m_vSize0;

//					printf("D1: %p %p\n", ptrDest1, ptrD1);
//					printf("D2: %p %p\n", ptrDest2, ptrD2);

					memcpy2((char *)ptrDest0, (char *)ptrSrc, nbytes);
					memcpy2((char *)ptrDest1, (char *)ptrSrc, nbytes);
					memcpy2((char *)ptrDest2, (char *)ptrSrc, nbytes);
				}
				else
				{
					printf("does this happen? if yes, just remove this line\n"); abort();
					if ((e[1]-s[1])%4 != 0)
					{
						for(int iy=s[1]; iy<e[1]; iy++)
						{
							const int my_y_dst = (iy-stencilStart[1])*m_vSize0;

							TElement * ptrDest = &cacheBlock->LinAccess(my_xz + my_y_dst);
							memcpy2((char *)ptrDest, (char *)ptrSrc, nbytes);
						}
					}
					else
					{
						for(int iy=s[1]; iy<e[1]; iy+=4)
						{
							const int my_y_dst0 = (iy-stencilStart[1])*m_vSize0;
							const int my_y_dst1 = my_y_dst0 + m_vSize0;
							const int my_y_dst2 = my_y_dst1 + m_vSize0;
							const int my_y_dst3 = my_y_dst2 + m_vSize0;

							TElement * ptrDest0 = &cacheBlock->LinAccess(my_xz + my_y_dst0);
							TElement * ptrDest1 = &cacheBlock->LinAccess(my_xz + my_y_dst1);
							TElement * ptrDest2 = &cacheBlock->LinAccess(my_xz + my_y_dst2);
							TElement * ptrDest3 = &cacheBlock->LinAccess(my_xz + my_y_dst3);

							memcpy2((char *)ptrDest0, (char *)ptrSrc, nbytes);
							memcpy2((char *)ptrDest1, (char *)ptrSrc, nbytes);
							memcpy2((char *)ptrDest2, (char *)ptrSrc, nbytes);
							memcpy2((char *)ptrDest3, (char *)ptrSrc, nbytes);
						}
					}
				}
				#endif
			}
		}
		else /* (dir == 2) */
		{
			int pz;
			if (side == 0)	pz = 0;
			else			pz = TBlock::sizeZ-1;

			int nbytes = (e[0]-s[0])*sizeof(TElement);
			const int ix = s[0];

			const int my_x = ix-stencilStart[0];
			const int my_z_src = (pz-stencilStart[2])*m_nElemsPerSlice;

//			printf("dir2: from %d to %d\n", s[1], e[1]);
			for(int iz=s[2]; iz<e[2]; iz++)
			{
				const int my_z_dst = (iz-stencilStart[2])*m_nElemsPerSlice;

				if ((e[1]-s[1]) % 4 == 0)
				{
					for(int iy=s[1]; iy<e[1]; iy+=4)
					{
						const int my_yx0 = (iy-stencilStart[1])*m_vSize0 + my_x;
						//const int my_yx1 = my_yx0 + m_vSize0;
						//const int my_yx2 = my_yx1 + m_vSize0;
						//const int my_yx3 = my_yx2 + m_vSize0;

						TElement * ptrSrc0 = &cacheBlock->LinAccess(my_yx0 + my_z_src);
						//TElement * ptrSrc1 = &cacheBlock->LinAccess(my_yx1 + my_z_src);
						//TElement * ptrSrc2 = &cacheBlock->LinAccess(my_yx2 + my_z_src);
						//TElement * ptrSrc3 = &cacheBlock->LinAccess(my_yx3 + my_z_src);
						TElement * ptrSrc1 = ptrSrc0 + m_vSize0;
						TElement * ptrSrc2 = ptrSrc1 + m_vSize0;
						TElement * ptrSrc3 = ptrSrc2 + m_vSize0;

						TElement * ptrDest0 = &cacheBlock->LinAccess(my_yx0 + my_z_dst);
						//TElement * ptrDest1 = &cacheBlock->LinAccess(my_yx1 + my_z_dst);
						//TElement * ptrDest2 = &cacheBlock->LinAccess(my_yx2 + my_z_dst);
						//TElement * ptrDest3 = &cacheBlock->LinAccess(my_yx3 + my_z_dst);
						TElement * ptrDest1 = ptrDest0 + m_vSize0;
						TElement * ptrDest2 = ptrDest1 + m_vSize0;
						TElement * ptrDest3 = ptrDest2 + m_vSize0;

						memcpy2((char *)ptrDest0, (char *)ptrSrc0, nbytes);
						memcpy2((char *)ptrDest1, (char *)ptrSrc1, nbytes);
						memcpy2((char *)ptrDest2, (char *)ptrSrc2, nbytes);
						memcpy2((char *)ptrDest3, (char *)ptrSrc3, nbytes);
					}
				}
				else
				{
					for(int iy=s[1]; iy<e[1]; iy++)
					{
						const int my_yx = (iy-stencilStart[1])*m_vSize0 + my_x;

						TElement * ptrSrc = &cacheBlock->LinAccess(my_yx + my_z_src);
						TElement * ptrDest = &cacheBlock->LinAccess(my_yx + my_z_dst);
						memcpy2((char *)ptrDest, (char *)ptrSrc, nbytes);
					}
				}
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
#if 1
        for(int iz=s[2]; iz<e[2]; iz++)
            for(int iy=s[1]; iy<e[1]; iy++)
                for(int ix=s[0]; ix<e[0]; ix++)
                {
                    TElement source = (*this)(dir==0? (side==0? -ix-1:2*TBlock::sizeX-ix-1):ix,
                            dir==1? (side==0? -iy-1:2*TBlock::sizeY-iy-1):iy,
                            dir==2? (side==0? -iz-1:2*TBlock::sizeZ-iz-1):iz);
                    (*this)(ix,iy,iz).alpha1rho1      = source.alpha1rho1;
                    (*this)(ix,iy,iz).alpha2rho2      = source.alpha2rho2;
                    (*this)(ix,iy,iz).ru              = ((dir==0)? (Real)(-1.0):(Real)(1.0))*source.ru;
                    (*this)(ix,iy,iz).rv              = ((dir==1)? (Real)(-1.0):(Real)(1.0))*source.rv;
                    (*this)(ix,iy,iz).rw              = ((dir==2)? (Real)(-1.0):(Real)(1.0))*source.rw;
                    (*this)(ix,iy,iz).energy          = source.energy;
                    (*this)(ix,iy,iz).alpha2          = source.alpha2;
                }
#else
        // TODO: [fabianw@mavt.ethz.ch; Sun Apr 30 2017 06:53:39 PM (-0700)]
        // this needs adaption to relaxed 5equation system.
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

#if 0
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
#else
		const int m_vSize0 = cacheBlock->getSize(0); //m_vSize[0];
		const int m_nElemsPerSlice = cacheBlock->getNumberOfElementsPerSlice(); //m_nElementsPerSlice;

		if (dir == 0)
		{
			int px;
			if (side == 0)	px = -1; //-ix;
			else			px = 2*TBlock::sizeX-1; //-ix;

//			printf("dir%d: y(%d,%d), x(%d,%d)\n", dir, s[1], e[1], s[0], e[0]);
			for(int iz=s[2]; iz<e[2]; iz++)
			{
				const int my_z = (iz-stencilStart[2])*m_nElemsPerSlice;
				for(int iy=s[1]; iy<e[1]; iy++)
				{
					const int my_y = (iy-stencilStart[1])*m_vSize0;
					#if 0
					for(int ix=s[0]; ix<e[0]; ix++)
					{
						TElement *ptrsrc = &cacheBlock->LinAccess((px-ix)-stencilStart[0] + my_y + my_z);
						TElement *ptrdst = &cacheBlock->LinAccess(ix-stencilStart[0] + my_y + my_z);
						memcpy2((char *)ptrdst, (char *)ptrsrc, sizeof(TElement));
						ptrdst->u = -ptrsrc->u;
					}
					#else
					TElement *ptrsrc0 = &cacheBlock->LinAccess((px-s[0])-stencilStart[0] + my_y + my_z);
					TElement *ptrdst0 = &cacheBlock->LinAccess(s[0]-stencilStart[0] + my_y + my_z);

					for(int ix=0; ix<e[0]-s[0]; ix++)
					{
						TElement *ptrsrc = ptrsrc0 + ix;
						TElement *ptrdst = ptrdst0 + ix;
						memcpy2((char *)ptrdst, (char *)ptrsrc, sizeof(TElement));
						ptrdst->ru = -ptrsrc->ru;
					}
					#endif
				}
			}
		}
		else if (dir == 1)
		{
			int py;
			if (side == 0)	py = -1; //-iy;
			else			py = 2*TBlock::sizeY-1; //-iy;

//			printf("dir%d: y(%d,%d), x(%d,%d)\n", dir, s[1], e[1], s[0], e[0]);

			for(int iz=s[2]; iz<e[2]; iz++)
			{
				const int my_z = (iz-stencilStart[2])*m_nElemsPerSlice;
				for(int iy=s[1]; iy<e[1]; iy++)
				{
					const int my_y_src = (py-iy-stencilStart[1])*m_vSize0;
					const int my_y_dst = (iy-stencilStart[1])*m_vSize0;
					#if 0
					for(int ix=s[0]; ix<e[0]; ix++)
					{
						TElement *ptrsrc = &cacheBlock->LinAccess(ix-stencilStart[0] + my_y_src + my_z);
						TElement *ptrdst = &cacheBlock->LinAccess(ix-stencilStart[0] + my_y_dst + my_z);
						memcpy2((char *)ptrdst, (char *)ptrsrc, sizeof(TElement));
						ptrdst->v = -ptrsrc->v;
					}
					#else
					TElement *ptrsrc0 = &cacheBlock->LinAccess(s[0]-stencilStart[0] + my_y_src + my_z);
					TElement *ptrdst0 = &cacheBlock->LinAccess(s[0]-stencilStart[0] + my_y_dst + my_z);
					for(int ix=0; ix<e[0]-s[0]; ix++)
					{
						TElement *ptrsrc = ptrsrc0 + ix;
						TElement *ptrdst = ptrdst0 + ix;
						memcpy2((char *)ptrdst, (char *)ptrsrc, sizeof(TElement));
						ptrdst->rv = -ptrsrc->rv;
					}
					#endif
				}
			}
		}
		else /* (dir == 2) */
		{
			int pz;
			if (side == 0)	pz = -1; //-iz;
			else			pz = 2*TBlock::sizeZ-1; //-iz;

			//printf("dir%d: y(%d,%d), x(%d,%d)\n", dir, s[1], e[1], s[0], e[0]);

			for(int iz=s[2]; iz<e[2]; iz++)
			{
				const int my_z_src = (pz-iz-stencilStart[2])*m_nElemsPerSlice;
				const int my_z_dst = (iz-stencilStart[2])*m_nElemsPerSlice;
				for(int iy=s[1]; iy<e[1]; iy++)
				{
					const int my_y = (iy-stencilStart[1])*m_vSize0;
					#if 0
					for(int ix=s[0]; ix<e[0]; ix++)
					{
						TElement *ptrsrc = &cacheBlock->LinAccess(ix-stencilStart[0] + my_y + my_z_src);
						TElement *ptrdst = &cacheBlock->LinAccess(ix-stencilStart[0] + my_y + my_z_dst);
						memcpy2((char *)ptrdst, (char *)ptrsrc, sizeof(TElement));
						ptrdst->w = -ptrsrc->w;
					}
					#else
					TElement *ptrsrc0 = &cacheBlock->LinAccess(s[0]-stencilStart[0] + my_y + my_z_src);
					TElement *ptrdst0 = &cacheBlock->LinAccess(s[0]-stencilStart[0] + my_y + my_z_dst);
					#if 0
					for(int ix=0; ix<e[0]-s[0]; ix++)
					{
						TElement *ptrsrc = ptrsrc0 + ix;
						TElement *ptrdst = ptrdst0 + ix;
						memcpy2((char *)ptrdst, (char *)ptrsrc, sizeof(TElement));
						//ptrdst->w = -ptrsrc->w;
						ptrdst->w = -ptrdst->w;
						if (ptrdst->w != -ptrsrc->w) abort();
					}
					#else
					memcpy2((char *)ptrdst0, (char *)ptrsrc0, (e[0]-s[0])*sizeof(TElement));
					for(int ix=0; ix<e[0]-s[0]; ix++)
					{
						TElement *ptrdst = ptrdst0 + ix;
						ptrdst->rw = -ptrdst->rw;
					}
					#endif
					#endif
				}
			}

		}
#endif

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
					(*this)(ix,iy,iz).alpha1rho1 = p.alpha1rho1;
					(*this)(ix,iy,iz).alpha2rho2 = p.alpha2rho2;
					(*this)(ix,iy,iz).ru         = p.ru;
					(*this)(ix,iy,iz).rv         = p.rv;
					(*this)(ix,iy,iz).rw         = p.rw;
					(*this)(ix,iy,iz).energy     = p.energy;
					(*this)(ix,iy,iz).alpha2     = p.alpha2;
				}
	}

    template<int dir, int side>
    void applyBC_dirichletNonUniform(const BlockInfo& info, const std::vector<TElement>& p)
    {
        assert(dir==0); // we only support dir = 0 for this boundary
        _setup<dir,side>();

        int dummy_ix;
        if (side==0) dummy_ix = 0;
        else         dummy_ix = TBlock::sizeX - 1;

        for(int iz=s[2]; iz<e[2]; iz++)
            for(int iy=s[1]; iy<e[1]; iy++)
                for(int ix=s[0]; ix<e[0]; ix++)
                {
                    double pos[3];
                    info.pos(pos, dummy_ix, iy, iz);

                    if (pos[1] < 1./3.*Simulation_Environment::extents[1])
                        (*this)(ix,iy,iz) = p[0];
                    else if (pos[1] < 2./3.*Simulation_Environment::extents[1])
                        (*this)(ix,iy,iz) = p[1];
                    else
                        (*this)(ix,iy,iz) = p[0];
                }
    }

  // symmetrie boundary condition (equivalent to reflecting above)
  // rasthofer January 2016
  template<int dir, int side>
  void applyBC_symmetry()
  {
    _setup<dir,side>();

    for (int iz=s[2]; iz<e[2]; iz++)
      for (int iy=s[1]; iy<e[1]; iy++)
        for (int ix=s[0]; ix<e[0]; ix++)
        {
          TElement source = (*this)(dir==0? (side==0? -ix-1:2*TBlock::sizeX-ix-1):ix,
                                    dir==1? (side==0? -iy-1:2*TBlock::sizeY-iy-1):iy,
				    dir==2? (side==0? -iz-1:2*TBlock::sizeZ-iz-1):iz);

          (*this)(ix,iy,iz).alpha1rho1      = source.alpha1rho1;
          (*this)(ix,iy,iz).alpha2rho2      = source.alpha2rho2;
          (*this)(ix,iy,iz).ru              = ((dir==0)? (Real)(-1.0):(Real)(1.0))*source.ru;
          (*this)(ix,iy,iz).rv              = ((dir==1)? (Real)(-1.0):(Real)(1.0))*source.rv;
          (*this)(ix,iy,iz).rw              = ((dir==2)? (Real)(-1.0):(Real)(1.0))*source.rw;
          (*this)(ix,iy,iz).energy          = source.energy;
          (*this)(ix,iy,iz).alpha2          = source.alpha2;
        }
  }

  // note: the boundary values for the faces have to be set first
  // rasthofer January 2016
  template<int dir, int side>
  void applyBC_symmetry_tensorials_edges_corner()
  {
    const int bsize[3] = {TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ};
    int s[3], e[3];
    bool db[3];

    // edges, corners
    {
      s[dir] = stencilStart[dir]*(1-side) + bsize[dir]*side;
      e[dir] = (bsize[dir]-1+stencilEnd[dir])*side;

      const int d1 = (dir + 1) % 3;
      const int d2 = (dir + 2) % 3;

      for (int b=0; b<2; ++b)
        for (int a=0; a<2; ++a)
        {
          s[d1] = b*stencilStart[d1] + a*b*(bsize[d1] - stencilStart[d1]);
          s[d2] = stencilStart[d2] + (a-a*b)*(bsize[d2] - stencilStart[d2]);

          e[d1] = (1-b+a*b)*(bsize[d1] - 1 + stencilEnd[d1]) + (1-b) * stencilStart[d1];
          e[d2] = (a+b-a*b)*(bsize[d2] - 1 + stencilEnd[d2]);

          for (int dim=0; dim<3; dim++)
            db[dim] = d1==dim ? (b==0 ? true : false) : (b==1 ? true : false);
          db[dir] = true;

          for (int iz=s[2]; iz<e[2]; iz++)
           for (int iy=s[1]; iy<e[1]; iy++)
             for (int ix=s[0]; ix<e[0]; ix++)
             {
               TElement source = (*this)(dir==0? ix : (db[0]? ix : (a==0? -ix-1:2*TBlock::sizeX-ix-1)),
                                         dir==1? iy : (db[1]? iy : (a==0? -iy-1:2*TBlock::sizeY-iy-1)),
                                         dir==2? iz : (db[2]? iz : (a==0? -iz-1:2*TBlock::sizeZ-iz-1)));

               (*this)(ix,iy,iz).alpha1rho1      = source.alpha1rho1;
               (*this)(ix,iy,iz).alpha2rho2      = source.alpha2rho2;
               (*this)(ix,iy,iz).ru              = (!db[0]? (Real)(-1.0):(Real)(1.0))*source.ru;
               (*this)(ix,iy,iz).rv              = (!db[1]? (Real)(-1.0):(Real)(1.0))*source.rv;
               (*this)(ix,iy,iz).rw              = (!db[2]? (Real)(-1.0):(Real)(1.0))*source.rw;
               (*this)(ix,iy,iz).energy          = source.energy;
               (*this)(ix,iy,iz).alpha2          = source.alpha2;
             }
        }
    }
  }


  // noslip boundary condition with symmetry assumption for pressure and density
  // rasthofer June 2016
  template<int dir, int side>
  void applyBC_noslip_simple()
  {
    _setup<dir,side>();

    for (int iz=s[2]; iz<e[2]; iz++)
      for (int iy=s[1]; iy<e[1]; iy++)
        for (int ix=s[0]; ix<e[0]; ix++)
        {
          TElement source = (*this)(dir==0? (side==0? -ix-1:2*TBlock::sizeX-ix-1):ix,
                                    dir==1? (side==0? -iy-1:2*TBlock::sizeY-iy-1):iy,
                                    dir==2? (side==0? -iz-1:2*TBlock::sizeZ-iz-1):iz);

          (*this)(ix,iy,iz).alpha1rho1      = source.alpha1rho1;
          (*this)(ix,iy,iz).alpha2rho2      = source.alpha2rho2;
          (*this)(ix,iy,iz).ru              = -source.ru;
          (*this)(ix,iy,iz).rv              = -source.rv;
          (*this)(ix,iy,iz).rw              = -source.rw;
          (*this)(ix,iy,iz).energy          = source.energy;
          (*this)(ix,iy,iz).alpha2          = source.alpha2;
        }
  }


  // note: the boundary values for the faces have to be set first
  // rasthofer June 2016
  template<int dir, int side>
  void applyBC_noslip_simple_tensorial_pbc()
  {
    const int bsize[3] = {TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ};
    int s[3], e[3];

    // edges, corners
    {
      s[dir] = stencilStart[dir]*(1-side) + bsize[dir]*side;
      e[dir] = (bsize[dir]-1+stencilEnd[dir])*side;

      const int d1 = (dir + 1) % 3;
      const int d2 = (dir + 2) % 3;

      for (int b=0; b<2; ++b)
        for (int a=0; a<2; ++a)
        {
          s[d1] = b*stencilStart[d1] + a*b*(bsize[d1] - stencilStart[d1]);
          s[d2] = stencilStart[d2] + (a-a*b)*(bsize[d2] - stencilStart[d2]);

          e[d1] = (1-b+a*b)*(bsize[d1] - 1 + stencilEnd[d1]) + (1-b) * stencilStart[d1];
          e[d2] = (a+b-a*b)*(bsize[d2] - 1 + stencilEnd[d2]);

          for (int iz=s[2]; iz<e[2]; iz++)
           for (int iy=s[1]; iy<e[1]; iy++)
             for (int ix=s[0]; ix<e[0]; ix++)
             {
               TElement source = (*this)(dir==0? (side==0? -ix-1:2*TBlock::sizeX-ix-1):ix,
                                         dir==1? (side==0? -iy-1:2*TBlock::sizeY-iy-1):iy,
                                         dir==2? (side==0? -iz-1:2*TBlock::sizeZ-iz-1):iz);

               (*this)(ix,iy,iz).alpha1rho1      = source.alpha1rho1;
               (*this)(ix,iy,iz).alpha2rho2      = source.alpha2rho2;
               (*this)(ix,iy,iz).ru              = -source.ru;
               (*this)(ix,iy,iz).rv              = -source.rv;
               (*this)(ix,iy,iz).rw              = -source.rw;
               (*this)(ix,iy,iz).energy          = source.energy;
               (*this)(ix,iy,iz).alpha2          = source.alpha2;
             }
        }
    }
  }


  // algebraic characteristic-based far-field boundary condition
  // see textbook "Computational Fluid Dynamics: Principles and Applications" by J. Blazek
  // note: currently only subsonic
  // rasthofer April 2016
  template<int dir, int side>
  void applyBC_farfield(const TElement& p)
  {
    _setup<dir,side>();

    // define normal on boundary
    const Real nx = (dir==0? 1 : 0) * (2*side-1);
    const Real ny = (dir==1? 1 : 0) * (2*side-1);
    const Real nz = (dir==2? 1 : 0) * (2*side-1);

    // we assume primitive variables for the reference element
    const Real pinf = p.energy;
    const Real rinf = p.alpha1rho1;
    const Real uinf = p.ru;
    const Real vinf = p.rv;
    const Real winf = p.rw;

    for (int iz=s[2]; iz<e[2]; iz++)
      for (int iy=s[1]; iy<e[1]; iy++)
        for (int ix=s[0]; ix<e[0]; ix++)
        {
          TElement source = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1):ix,
                                    dir==1? (side==0? 0:TBlock::sizeY-1):iy,
                                    dir==2? (side==0? 0:TBlock::sizeZ-1):iz);

          // convert to primitive vales
          const Real rho = source.alpha1rho1 + source.alpha2rho2;
          const Real rInv = (Real) 1.0/rho;
          const Real u = rInv * source.ru;
          const Real v = rInv * source.rv;
          const Real w = rInv * source.rw;
          // assume pure fluid 1 for pressure computation (i.e., alpha1 = 1)
          const Real p = (Simulation_Environment::GAMMA1 - (Real)1.0) * (source.energy - (Real)0.5 * rInv * (source.ru*source.ru + source.rv * source.rv + source.rw * source.rw))
                       - Simulation_Environment::GAMMA1 * Simulation_Environment::PC1;
          // as usual, the reference state is set equalt to the state at the interior point
          const Real csq = Simulation_Environment::GAMMA1 * (p + Simulation_Environment::PC1)/rho;

          // this is only valid for subsonic flow
          const bool inflow = (u*nx + v*ny + w*nz >= (Real)0.0) ? false : true;

          if (inflow)
          {
            const Real pb = (Real)0.5 * (pinf+p + rho*sqrt(csq) * (nx*(uinf - u) + ny*(vinf - v) + nz*(winf - w)));
            const Real ub = uinf - nx*(pinf-pb) / (rho*sqrt(csq));
            const Real vb = vinf - ny*(pinf-pb) / (rho*sqrt(csq));
            const Real wb = winf - nz*(pinf-pb) / (rho*sqrt(csq));
            const Real rb = rinf + (pb-pinf)/csq;

            (*this)(ix,iy,iz).alpha1rho1      = rb;
            (*this)(ix,iy,iz).alpha2rho2      = 0.0;
            (*this)(ix,iy,iz).ru              = ub * rb;
            (*this)(ix,iy,iz).rv              = vb * rb;
            (*this)(ix,iy,iz).rw              = wb * rb;
            (*this)(ix,iy,iz).energy          = (Real)0.5 * (ub*ub + vb*vb + wb*wb) * rb
                                              + (pb + Simulation_Environment::GAMMA1*Simulation_Environment::PC1) / (Simulation_Environment::GAMMA1 - (Real)1.0);
            (*this)(ix,iy,iz).alpha2          = 0.0;
          }
          else //outflow
          {
            const Real pb = pinf;
            const Real ub = u - nx*(p-pb) / (rho*sqrt(csq));
            const Real vb = v - ny*(p-pb) / (rho*sqrt(csq));
            const Real wb = w - nz*(p-pb) / (rho*sqrt(csq));
            const Real rb = rho + (pb-p)/csq;

            (*this)(ix,iy,iz).alpha1rho1      = rb;
            (*this)(ix,iy,iz).alpha2rho2      = 0.0;
            (*this)(ix,iy,iz).ru              = ((dir==0)? ub : uinf) * rb;
            (*this)(ix,iy,iz).rv              = ((dir==1)? ub : vinf) * rb;
            (*this)(ix,iy,iz).rw              = ((dir==2)? ub : winf) * rb;
            (*this)(ix,iy,iz).energy          = (Real)0.5 * (ub*ub + vb*vb + wb*wb) * rb
                                              + (pb + Simulation_Environment::GAMMA1*Simulation_Environment::PC1) / (Simulation_Environment::GAMMA1 - (Real)1.0);
            (*this)(ix,iy,iz).alpha2          = 0.0;
          }
        }
  }

  // algebraic characteristic-based far-field boundary condition with extrapolation
  // see textbook "Computational Fluid Dynamics: Principles and Applications" by J. Blazek
  // note: currently only subsonic
  // rasthofer April 2016
  template<int dir, int side>
  void applyBC_farfield_extra(const TElement& p)
  {
    _setup<dir,side>();

    // define normal on boundary
    const Real nx = (dir==0? 1 : 0) * (2*side-1);
    const Real ny = (dir==1? 1 : 0) * (2*side-1);
    const Real nz = (dir==2? 1 : 0) * (2*side-1);

    // get loop boundaries
    const int sa = (side==1)? s[dir] : (e[dir]-1);
    const int ea = (side==0)? s[dir] : (e[dir]-1);
    const int sb = (dir==0) ? s[1]   : s[0];
    const int eb = (dir==0) ? e[1]   : e[0];
    const int sc = (dir==2) ? s[1]   : s[2];
    const int ec = (dir==2) ? e[1]   : e[2];

    // we assume primitive variables for the reference element
    const Real pinf = p.energy;
    const Real rinf = p.alpha1rho1;
    const Real uinf = p.ru;
    const Real vinf = p.rv;
    const Real winf = p.rw;

    // far-field boundary loop
    int ia = sa;
    for (int ib=sb; ib<eb; ib++)
      for (int ic=sc; ic<ec; ic++)
        {
          // convert ia, ib and ic to ix, iy and iz
          const int ix = (dir==0)? ia : ib;
          const int iy = (dir==1)? ia : ((dir==0)? ib : ic);
          const int iz = (dir==2)? ia : ic;

          TElement source = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1):ix,
                                    dir==1? (side==0? 0:TBlock::sizeY-1):iy,
                                    dir==2? (side==0? 0:TBlock::sizeZ-1):iz);

          // convert to primitive vales
          const Real rho = source.alpha1rho1 + source.alpha2rho2;
          const Real rInv = (Real) 1.0/rho;
          const Real u = rInv * source.ru;
          const Real v = rInv * source.rv;
          const Real w = rInv * source.rw;
          // assume pure fluid 1 for pressure computation (i.e., alpha1 = 1)
          const Real p = (Simulation_Environment::GAMMA1 - (Real)1.0) * (source.energy - (Real)0.5 * rInv * (source.ru*source.ru + source.rv * source.rv + source.rw * source.rw))
                       - Simulation_Environment::GAMMA1 * Simulation_Environment::PC1;
          // as usual, the reference state is set equalt to the state at the interior point
          const Real csq = Simulation_Environment::GAMMA1 * (p + Simulation_Environment::PC1)/rho;

          // this is only valid for subsonic flow
          const bool inflow = (u*nx + v*ny + w*nz >= (Real)0.0) ? false : true;

          if (inflow)
          {
            const Real pb = (Real)0.5 * (pinf+p + rho*sqrt(csq) * (nx*(uinf - u) + ny*(vinf - v) + nz*(winf - w)));
            const Real ub = uinf - nx*(pinf-pb) / (rho*sqrt(csq));
            const Real vb = vinf - ny*(pinf-pb) / (rho*sqrt(csq));
            const Real wb = winf - nz*(pinf-pb) / (rho*sqrt(csq));
            const Real rb = rinf + (pb-pinf)/csq;

            (*this)(ix,iy,iz).alpha1rho1      = rb;
            (*this)(ix,iy,iz).alpha2rho2      = 0.0;
            (*this)(ix,iy,iz).ru              = ub * rb;
            (*this)(ix,iy,iz).rv              = vb * rb;
            (*this)(ix,iy,iz).rw              = wb * rb;
            (*this)(ix,iy,iz).energy          = (Real)0.5 * (ub*ub + vb*vb + wb*wb) * rb
                                              + (pb + Simulation_Environment::GAMMA1*Simulation_Environment::PC1) / (Simulation_Environment::GAMMA1 - (Real)1.0);
            (*this)(ix,iy,iz).alpha2          = 0.0;
          }
          else //outflow
          {
            const Real pb = pinf;
            const Real ub = u - nx*(p-pb) / (rho*sqrt(csq));
            const Real vb = v - ny*(p-pb) / (rho*sqrt(csq));
            const Real wb = w - nz*(p-pb) / (rho*sqrt(csq));
            const Real rb = rho + (pb-p)/csq;

            (*this)(ix,iy,iz).alpha1rho1      = rb;
            (*this)(ix,iy,iz).alpha2rho2      = 0.0;
            (*this)(ix,iy,iz).ru              = ((dir==0)? ub : uinf) * rb;
            (*this)(ix,iy,iz).rv              = ((dir==1)? ub : vinf) * rb;
            (*this)(ix,iy,iz).rw              = ((dir==2)? ub : winf) * rb;
            (*this)(ix,iy,iz).energy          = (Real)0.5 * (ub*ub + vb*vb + wb*wb) * rb
                                              + (pb + Simulation_Environment::GAMMA1*Simulation_Environment::PC1) / (Simulation_Environment::GAMMA1 - (Real)1.0);
            (*this)(ix,iy,iz).alpha2          = 0.0;
          }
        }

    const int incr = 2*side-1;
    const int slice = 1;
    for (ia=(sa+incr); (incr*ia)<=(incr*ea); ia+=incr)
      for (int ib=sb; ib<eb; ib++)
        for (int ic=sc; ic<ec; ic++)
        {
          // convert ia, ib and ic to ix, iy and iz
          const int ix = (dir==0)? ia : ib;
          const int iy = (dir==1)? ia : ((dir==0)? ib : ic);
          const int iz = (dir==2)? ia : ic;

          // last element inside
          TElement source1 = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1):ix,
                                     dir==1? (side==0? 0:TBlock::sizeY-1):iy,
                                     dir==2? (side==0? 0:TBlock::sizeZ-1):iz);
          // first element outside
          TElement source2 = (*this)(dir==0? (side==0? -1:TBlock::sizeX):ix,
                                     dir==1? (side==0? -1:TBlock::sizeY):iy,
                                     dir==2? (side==0? -1:TBlock::sizeZ):iz);

          // do extrapolation
          (*this)(ix,iy,iz).alpha1rho1      = source1.alpha1rho1 + (Real)slice * (source2.alpha1rho1 - source1.alpha1rho1);
          (*this)(ix,iy,iz).alpha2rho2      = source1.alpha2rho2 + (Real)slice * (source2.alpha2rho2 - source1.alpha2rho2);
          (*this)(ix,iy,iz).ru              = source1.ru + (Real)slice * (source2.ru - source1.ru);
          (*this)(ix,iy,iz).rv              = source1.rv + (Real)slice * (source2.rv - source1.rv);
          (*this)(ix,iy,iz).rw              = source1.rw + (Real)slice * (source2.rw - source1.rw);
          (*this)(ix,iy,iz).energy          = source1.energy + (Real)slice * (source2.energy - source1.energy);
          (*this)(ix,iy,iz).alpha2          = source1.alpha2 + (Real)slice * (source2.alpha2 - source1.alpha2);
        }
  }

  template<int dir, int side>
  void applyBC_1dchar_ode(const BGC::BoundaryBlock& p)
  {
    _setup<dir,side>();

    // get increments for access to boundary cells
    const int delta_ix = dir!=0 ? 0 : (side==0 ? 3 : -TBlock::sizeX);
    const int delta_iy = dir!=1 ? 0 : (side==0 ? 3 : -TBlock::sizeY);
    const int delta_iz = dir!=2 ? 0 : (side==0 ? 3 : -TBlock::sizeZ);
    const int bnx = dir==0 ? 3 : TBlock::sizeX;
    const int bny = dir==1 ? 3 : TBlock::sizeY;

    for(int iz=s[2]; iz<e[2]; iz++)
      for(int iy=s[1]; iy<e[1]; iy++)
        for(int ix=s[0]; ix<e[0]; ix++)
        {
          const int index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

          const Real rho = p.fe[index].rho;
          const Real u = p.fe[index].u;
          const Real e = (Real)0.5 * rho * u*u + (p.fe[index].pressure + Simulation_Environment::GAMMA1*Simulation_Environment::PC1)/(Simulation_Environment::GAMMA1 - (Real)1.0);

          assert(!isnan(rho));
          assert(!isnan(u));
          assert(!isnan(e));

          (*this)(ix,iy,iz).alpha1rho1 = rho;
          (*this)(ix,iy,iz).alpha2rho2 = 0.0;
          (*this)(ix,iy,iz).ru         = dir==0 ? u*rho : 0.0;
          (*this)(ix,iy,iz).rv         = dir==1 ? u*rho : 0.0;
          (*this)(ix,iy,iz).rw         = dir==2 ? u*rho : 0.0;
          (*this)(ix,iy,iz).energy     = e;
          (*this)(ix,iy,iz).alpha2     = 0.0;
       }
  }

  // time update for ode using Runge-Kutta
  template<int dir, int side>
  void updateBC_1dchar_ode(BGC::BoundaryBlock& vec, const Real h, const Real dt, const Real a, const Real b, const Real pamb, const Real L, const Real lambda)
  {
    // safety check
    if (ALPHAEPS != 0.0)
     cout << "### WARNING: ode bc assumes pure fluid, no mixtures!!! ###" << endl;


    _setup<dir,side>();

    // define normal on boundary (w/o sign)
    const Real nx = (dir==0? 1 : 0);
    const Real ny = (dir==1? 1 : 0);
    const Real nz = (dir==2? 1 : 0);

    const int sa = side==0 ? -3 : (dir==0? TBlock::sizeX : (dir==1 ? TBlock::sizeY : TBlock::sizeZ));
    const int ea = side==0 ? 0 : (dir==0? TBlock::sizeX+3 : (dir==1 ? TBlock::sizeY+3 : TBlock::sizeZ+3));
    const int eb = dir==0 ? TBlock::sizeY : TBlock::sizeX;
    const int ec = dir==2 ? TBlock::sizeY : TBlock::sizeZ;

    // get increments for access to boundary cells
    const int delta_ix = dir!=0 ? 0 : (side==0 ? 3 : -TBlock::sizeX);
    const int delta_iy = dir!=1 ? 0 : (side==0 ? 3 : -TBlock::sizeY);
    const int delta_iz = dir!=2 ? 0 : (side==0 ? 3 : -TBlock::sizeZ);
    const int bnx = dir==0 ? 3 : TBlock::sizeX;
    const int bny = dir==1 ? 3 : TBlock::sizeY;

    const Real fac = (side==0? 1.0: -1.0) /((Real)12.0*h);

    for(int ia=sa; ia<ea; ia++)
      for(int ib=0; ib<eb; ib++)
        for(int ic=0; ic<ec; ic++)
        {
          // convert ia, ib and ic to ix, iy and iz
          const int ix = (dir==0)? ia : ib;
          const int iy = (dir==1)? ia : ((dir==0)? ib : ic);
          const int iz = (dir==2)? ia : ic;

          const int index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

          vector<Real> r(5,0.0);
          vector<Real> u(5,0.0);
          vector<Real> p(5,0.0);
          for (size_t ii=0; ii<5; ii++)
          {
            const int bix = dir == 0? (side == 0? ix+ii : ix-ii) : ix;
            const int biy = dir == 1? (side == 0? iy+ii : iy-ii) : iy;
            const int biz = dir == 2? (side == 0? iz+ii : iz-ii) : iz;
            r[ii] = (*this)(bix,biy,biz).alpha1rho1 + (*this)(bix,biy,biz).alpha2rho2;
            u[ii] = (*this)(bix,biy,biz).ru/r[ii] *nx + (*this)(bix,biy,biz).rv/r[ii] *ny + (*this)(bix,biy,biz).rw/r[ii] *nz;
            p[ii] = ((*this)(bix,biy,biz).energy - (Real)0.5 * ((*this)(bix,biy,biz).ru * (*this)(bix,biy,biz).ru
                                                              + (*this)(bix,biy,biz).rv * (*this)(bix,biy,biz).rv
                                                              + (*this)(bix,biy,biz).rw * (*this)(bix,biy,biz).rw) / r[ii]) * (Simulation_Environment::GAMMA1 - (Real)1.0)
                    - Simulation_Environment::GAMMA1*Simulation_Environment::PC1;

            assert(!isnan(r[ii]));
            assert(!isnan(u[ii]));
            assert(!isnan(p[ii]));
          }

          const Real drdx = (-(Real)25.0*r[0] + (Real)48.0*r[1] -  (Real)36.0*r[2] + (Real)16.0*r[3] - (Real)3.0*r[4])*fac;
          const Real dudx = (-(Real)25.0*u[0] + (Real)48.0*u[1] -  (Real)36.0*u[2] + (Real)16.0*u[3] - (Real)3.0*u[4])*fac;
          const Real dpdx = (-(Real)25.0*p[0] + (Real)48.0*p[1] -  (Real)36.0*p[2] + (Real)16.0*p[3] - (Real)3.0*p[4])*fac;

          assert(!isnan(drdx));
          assert(!isnan(dudx));
          assert(!isnan(dpdx));

          const Real rho = vec.fe[index].rho;
          const Real vel = vec.fe[index].u;
          const Real pres = vec.fe[index].pressure;
          const Real c = sqrt(Simulation_Environment::GAMMA1*(pres+Simulation_Environment::PC1)/rho);

          assert(!isnan(rho));
          assert(!isnan(vel));
          assert(!isnan(pres));
          assert(!isnan(c));

          Real L1 = side==0 ? ((vel-c)>=(Real)0.0 ? 0.0 : ((vel-c)*(dpdx-rho*c*dudx))) : ((vel-c)<=(Real)0.0 ? 0.0 : ((vel-c)*(dpdx-rho*c*dudx)));
          Real L2 = side==0 ? (vel>=(Real)0.0     ? 0.0 : (vel*(dpdx-c*c*drdx)))       : (vel<=(Real)0.0     ? 0.0 : (vel*(dpdx-c*c*drdx)));
          Real L3 = side==0 ? ((vel+c)>=(Real)0.0 ? 0.0 : ((vel+c)*(dpdx+rho*c*dudx))) : ((vel+c)<=(Real)0.0 ? 0.0 : ((vel+c)*(dpdx+rho*c*dudx)));

          // K is defined as given in Poinsot & Lele JCP 1992
          // K = lambda * (1-M^2) * c / L
          // where sigma is a user-defined parameter usually between 0.25 and 0.6
          // M the maximal Mach number
          // c the speed of sound
          // and L a charcteristic length of the domain
          // recommondation for free-field collapses: lambda=0.75, M=0 and L equal to
          // the distance from the center of the bubble/cloud to the boundary
          const Real K = lambda * c / L;
          L3 += side==0 ? K*(pres-pamb) : 0.0;
          L1 += side==1 ? K*(pres-pamb) : 0.0;

          assert(!isnan(L1));
          assert(!isnan(L2));
          assert(!isnan(L3));

          vec.rhs[index].pressure = a * vec.rhs[index].pressure - (Real)0.5 * dt * (L1+L3);
          vec.rhs[index].u = a * vec.rhs[index].u - dt * (L3-L1) / ((Real)2.0*rho*c);
          vec.rhs[index].rho = a * vec.rhs[index].rho + dt * (L2 - (Real)0.5 * (L1+L3))/(c*c);

          vec.fe[index].pressure += b*vec.rhs[index].pressure;
          vec.fe[index].u += b*vec.rhs[index].u;
          vec.fe[index].rho += b*vec.rhs[index].rho;

          assert(!isnan(vec.fe[index].pressure));
          assert(!isnan(vec.fe[index].u));
          assert(!isnan(vec.fe[index].rho));
          assert(!isnan(vec.rhs[index].pressure));
          assert(!isnan(vec.rhs[index].u));
          assert(!isnan(vec.rhs[index].rho));
        }
  }


#if 0

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

    // TESTING
    template<int dir, int side>
    void applyBC_dirichlet_inflow(const TElement& p)
    {
        _setup<dir,side>();

        for(int iz=s[2]; iz<e[2]; iz++)
            for(int iy=s[1]; iy<e[1]; iy++)
                for(int ix=s[0]; ix<e[0]; ix++)
                {
                    (*this)(ix,iy,iz).rho = p.rho;
                    (*this)(ix,iy,iz).u   = p.u;
                    (*this)(ix,iy,iz).v   = p.v;
                    (*this)(ix,iy,iz).w   = p.w;
                    (*this)(ix,iy,iz).G   = p.G;
                    (*this)(ix,iy,iz).P   = p.P;

                    // get interior cell
                    const TElement b = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1):ix,
                            dir==1? (side==0? 0:TBlock::sizeY-1):iy,
                            dir==2? (side==0? 0:TBlock::sizeZ-1):iz);
                    const Real pressure = (b.energy - b.P - 0.5*(b.u*b.u + b.v*b.v + b.w*b.w)/b.rho) / b.G;
                    (*this)(ix,iy,iz).energy = p.G*pressure + p.P + 0.5*(p.u*p.u + p.v*p.v + p.w*p.w)/p.rho;
                }
    }

    template<int dir, int side>
    void applyBC_absorbing_constPressure(const Real pinf)
    {
        _setup<dir,side>();

        for(int iz=s[2]; iz<e[2]; iz++)
            for(int iy=s[1]; iy<e[1]; iy++)
                for(int ix=s[0]; ix<e[0]; ix++)
                {
                    TElement b = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1):ix,
                            dir==1? (side==0? 0:TBlock::sizeY-1):iy,
                            dir==2? (side==0? 0:TBlock::sizeZ-1):iz);
                    b.energy = b.G*pinf + b.P + 0.5*(b.u*b.u + b.v*b.v + b.w*b.w)/b.rho;
                    (*this)(ix,iy,iz) = b;
                }
    }

#endif /* 0 */

    ///////////////////////////////////////////////////////////////////////////
    /* TODO: (fabianw; Mon 12 Oct 2015 05:22:25 PM CEST) The following
     * boundaries are tensorial, where corners and edges are just put as it is.
     * Apply these boundaries only where you have a 1D mean flow along one of
     * the principal coordinates, e.g. inflow/outflow in z, preiodic x,y */

    /* TODO: (fabianw; Mon 12 Oct 2015 04:59:49 PM CEST) The adjacent edges and
     * faces must be previously computed, like it is done for periodic BC */
    template <int dir, int side>
    void applyBC_absorbing_tensorials_edges_corner()
    {
        const int bsize[3] = {TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ};
        int s[3], e[3];

        // edges, corners
        {
            s[dir] = stencilStart[dir]*(1-side) + bsize[dir]*side;
            e[dir] = (bsize[dir]-1+stencilEnd[dir])*side;

            const int d1 = (dir + 1) % 3;
            const int d2 = (dir + 2) % 3;

            for(int b=0; b<2; ++b)
                for(int a=0; a<2; ++a)
                {
                    s[d1] = stencilStart[d1] + a*b*(bsize[d1] - stencilStart[d1]);
                    s[d2] = stencilStart[d2] + (a-a*b)*(bsize[d2] - stencilStart[d2]);

                    e[d1] = (1-b+a*b)*(bsize[d1] - 1 + stencilEnd[d1]);
                    e[d2] = (a+b-a*b)*(bsize[d2] - 1 + stencilEnd[d2]);

                    for(int iz=s[2]; iz<e[2]; iz++)
                        for(int iy=s[1]; iy<e[1]; iy++)
                            for(int ix=s[0]; ix<e[0]; ix++)
                            {
                                (*this)(ix,iy,iz) = dir==0? (*this)(side*(bsize[0]-1),iy,iz) : (dir==1? (*this)(ix,side*(bsize[1]-1),iz) : (*this)(ix,iy,side*(bsize[2]-1)));
                            }
                }
        }
    }

    template<int dir, int side>
    void applyBC_dirichlet_tensorial_all(const TElement& p)
    {
        _setup<dir,side>();

        // faces
        for(int iz=s[2]; iz<e[2]; iz++)
            for(int iy=s[1]; iy<e[1]; iy++)
                for(int ix=s[0]; ix<e[0]; ix++)
                    (*this)(ix,iy,iz) = p;

        // edges, corners
        {
            const int bsize[3] = {TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ};
            int s[3], e[3];

            s[dir] = stencilStart[dir]*(1-side) + bsize[dir]*side;
            e[dir] = (bsize[dir]-1+stencilEnd[dir])*side;

            const int d1 = (dir + 1) % 3;
            const int d2 = (dir + 2) % 3;

            for(int b=0; b<2; ++b)
                for(int a=0; a<2; ++a)
                {
                    s[d1] = stencilStart[d1] + a*b*(bsize[d1] - stencilStart[d1]);
                    s[d2] = stencilStart[d2] + (a-a*b)*(bsize[d2] - stencilStart[d2]);

                    e[d1] = (1-b+a*b)*(bsize[d1] - 1 + stencilEnd[d1]);
                    e[d2] = (a+b-a*b)*(bsize[d2] - 1 + stencilEnd[d2]);

                    for(int iz=s[2]; iz<e[2]; iz++)
                        for(int iy=s[1]; iy<e[1]; iy++)
                            for(int ix=s[0]; ix<e[0]; ix++)
                                (*this)(ix,iy,iz) = p;
                }
        }
    }

    template<int dir, int side>
    void applyBC_movingWall(const double wallMom_U, const double wallMom_V=0, const double wallMom_W=0)
    {
        _setup<dir,side>();

        // faces
        for (int iz=s[2]; iz<e[2]; iz++)
            for (int iy=s[1]; iy<e[1]; iy++)
                for (int ix=s[0]; ix<e[0]; ix++)
                {
                    TElement source_reflect = (*this)(dir==0? (side==0? -ix-1:2*TBlock::sizeX-ix-1):ix,
                            dir==1? (side==0? -iy-1:2*TBlock::sizeY-iy-1):iy,
                            dir==2? (side==0? -iz-1:2*TBlock::sizeZ-iz-1):iz);

                    (*this)(ix,iy,iz).alpha1rho1 = source_reflect.alpha1rho1;
                    (*this)(ix,iy,iz).alpha2rho2 = source_reflect.alpha2rho2;
                    (*this)(ix,iy,iz).ru         = -source_reflect.ru + 2.0*wallMom_U;
                    (*this)(ix,iy,iz).rv         = -source_reflect.rv + 2.0*wallMom_V;
                    (*this)(ix,iy,iz).rw         = -source_reflect.rw + 2.0*wallMom_W;
                    (*this)(ix,iy,iz).energy     = source_reflect.energy;
                    (*this)(ix,iy,iz).alpha2     = source_reflect.alpha2;
                }
    }

    template<int dir, int side>
    void applyBC_movingWall_tensorial_all(const double wallMom_U, const double wallMom_V=0, const double wallMom_W=0)
    {
        // faces
        this->template applyBC_movingWall<dir,side>(wallMom_U, wallMom_V, wallMom_W);

        // edges, corners
        {
            const int bsize[3] = {TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ};
            int s[3], e[3];

            s[dir] = stencilStart[dir]*(1-side) + bsize[dir]*side;
            e[dir] = (bsize[dir]-1+stencilEnd[dir])*side;

            const int d1 = (dir + 1) % 3;
            const int d2 = (dir + 2) % 3;

            for(int b=0; b<2; ++b)
                for(int a=0; a<2; ++a)
                {
                    s[d1] = stencilStart[d1] + a*b*(bsize[d1] - stencilStart[d1]);
                    s[d2] = stencilStart[d2] + (a-a*b)*(bsize[d2] - stencilStart[d2]);

                    e[d1] = (1-b+a*b)*(bsize[d1] - 1 + stencilEnd[d1]);
                    e[d2] = (a+b-a*b)*(bsize[d2] - 1 + stencilEnd[d2]);

                    for (int iz=s[2]; iz<e[2]; iz++)
                        for (int iy=s[1]; iy<e[1]; iy++)
                            for (int ix=s[0]; ix<e[0]; ix++)
                            {
                                TElement source_reflect = (*this)(dir==0? (side==0? -ix-1:2*TBlock::sizeX-ix-1):ix,
                                        dir==1? (side==0? -iy-1:2*TBlock::sizeY-iy-1):iy,
                                        dir==2? (side==0? -iz-1:2*TBlock::sizeZ-iz-1):iz);

                                (*this)(ix,iy,iz).alpha1rho1 = source_reflect.alpha1rho1;
                                (*this)(ix,iy,iz).alpha2rho2 = source_reflect.alpha2rho2;
                                (*this)(ix,iy,iz).ru         = -source_reflect.ru + 2.0*wallMom_U;
                                (*this)(ix,iy,iz).rv         = -source_reflect.rv + 2.0*wallMom_V;
                                (*this)(ix,iy,iz).rw         = -source_reflect.rw + 2.0*wallMom_W;
                                (*this)(ix,iy,iz).energy     = source_reflect.energy;
                                (*this)(ix,iy,iz).alpha2     = source_reflect.alpha2;
                            }
                }
        }
    }
};

#endif // BOUNDARYCONDITIONS_H
