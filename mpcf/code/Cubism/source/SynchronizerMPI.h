/*
 *  SynchronizerMPI.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 10/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include "mpi.h"

#include "PUPkernelsMPI.h"
#include "DependencyCubeMPI.h"

using namespace std;

class SynchronizerMPI
{
	struct I3
	{
		int ix, iy, iz;

		I3(int ix, int iy, int iz):ix(ix), iy(iy), iz(iz){}
		I3(const I3& c): ix(c.ix), iy(c.iy), iz(c.iz){}

		bool operator<(const I3& a) const
		{
			return ix<a.ix || ix==a.ix && iy<a.iy || ix==a.ix && iy==a.iy && iz<a.iz;
		}
	};

	struct PackInfo { Real * block, * pack; int sx, sy, sz, ex, ey, ez; };
	struct SubpackInfo { Real * block, * pack; int sx, sy, sz, ex, ey, ez; int x0, y0, z0, xpacklenght, ypacklenght; };

	DependencyCubeMPI<MPI_Request> cube;
	bool isroot;
	const int synchID;
	int send_thickness[3][2], recv_thickness[3][2];
	int blockinfo_counter;
	StencilInfo stencil;
	vector<PackInfo> send_packinfos;
	map<Real *, vector<PackInfo> > recv_packinfos;
	map<Real *, vector<SubpackInfo> > recv_subpackinfos;
	vector<Real *> all_mallocs;

	vector<BlockInfo> globalinfos;

    map<Region, vector<BlockInfo> > region2infos;

	//?static?
	MPI_Comm cartcomm;
	int blocksize[3];
	int mypeindex[3], pesize[3], mybpd[3];
	int periodic[3];
	int neighborsrank[3][3][3];

	map<I3,int> c2i;

	struct CommData {
		Real * faces[3][2], * edges[3][2][2], * corners[2][2][2];
		set<MPI_Request> pending;
	} send, recv;

	bool _face_needed(const int d) const
	{
		return periodic[d] || mypeindex[d] > 0 && mypeindex[d] < pesize[d]-1;
	}

	bool _myself(const int indx[3])
	{
		return (indx[0]+pesize[0]) % pesize[0] == mypeindex[0] &&
		(indx[1]+pesize[1]) % pesize[1] == mypeindex[1] &&
		(indx[2]+pesize[2]) % pesize[2] == mypeindex[2];
	}

	int _rank(const int indx[3])
	{
		int indx_final[3]={indx[0],indx[1],indx[2]};

		for(int i=0; i<3; ++i)
		{
			if (pesize[i]==1) continue;
			const int d=indx[i]- mypeindex[i];
			indx_final[i]=d-pesize[i]*(int)((double)d/(pesize[i]-1))+mypeindex[i];
		}

#if !defined(NDEBUG)
		for(int i=0;i<3;++i)
			assert(indx_final[i]>=-1+mypeindex[i] && indx_final[i]<2+mypeindex[i]);
#endif
		return neighborsrank[indx_final[2]+1-mypeindex[2]][indx_final[1]+1-mypeindex[1]][indx_final[0]+1-mypeindex[0]];
	}

	template <bool computesubregions>
	map<Real *, vector<SubpackInfo> > _setup(CommData& data, const int thickness[3][2], const int blockstart[3], const int blockend[3], const int origin[3], vector<PackInfo>& packinfos)
	{
		map<Real *, vector<SubpackInfo> > retval;

		const int NC = stencil.selcomponents.size();
		const int bpd[3] = {
			mybpd[0],
			mybpd[1],
			mybpd[2]
		};

		//faces
		for(int d=0; d<3; ++d)
		{
			const int dim_other1 = (d+1)%3;
			const int dim_other2 = (d+2)%3;

			for(int s=0; s<2; ++s)
			{
				const int NFACEBLOCK = NC * thickness[d][s] * blocksize[dim_other1] * blocksize[dim_other2];
				const int NFACE = NFACEBLOCK * mybpd[dim_other1] * mybpd[dim_other2];

				const bool needed = _face_needed(d) || NFACE == 0;
				data.faces[d][s] = needed ? _myalloc(sizeof(Real)*NFACE, 16) : NULL;

				if (!needed) continue;

				int neighbor_index[3];
				neighbor_index[d] = (mypeindex[d] + 2*s-1 + pesize[d])%pesize[d];
				neighbor_index[dim_other1] = mypeindex[dim_other1];
				neighbor_index[dim_other2] = mypeindex[dim_other2];

				if (_myself(neighbor_index)) continue;

				int start[3];
				start[d] = (1-s)*blockstart[d] + s*(blockend[d]-thickness[d][s]);
				start[dim_other1] = 0;
				start[dim_other2] = 0;

				int end[3];
				end[d] = (1-s)*(blockstart[d] + thickness[d][s]) + s*blockend[d];
				end[dim_other1] = blocksize[dim_other1];
				end[dim_other2] = blocksize[dim_other2];

				const int n1 = bpd[dim_other1];
				const int n2 = bpd[dim_other2];

				for(int b=0; b<n2; ++b)
					for(int a=0; a<n1; ++a)
					{
						int index[3];
						index[d] = s*(bpd[d]-1);
						index[dim_other1] = a;
						index[dim_other2] = b;

						assert(c2i.find(I3(origin[0] + index[0], origin[1] + index[1], origin[2] + index[2]))!=c2i.end());
						const int blockid = c2i[I3(origin[0] + index[0], origin[1] + index[1], origin[2] + index[2])];

						PackInfo info = {(Real *)globalinfos[blockid].ptrBlock, data.faces[d][s] + NFACEBLOCK*(a + n1*b), start[0], start[1], start[2], end[0], end[1], end[2]};

						const bool nonempty = end[0]>start[0] && end[1]>start[1] && end[2]>start[2];
						if (nonempty) packinfos.push_back(info);
					}
			}
		}

		if (!stencil.tensorial) return retval;

		//edges
		for(int d=0; d<3; ++d)
		{
			const int dim_other1 = (d+1)%3;
			const int dim_other2 = (d+2)%3;

			for(int b=0; b<2; ++b)
				for(int a=0; a<2; ++a)
				{
					const int NEDGEBLOCK = NC * blocksize[d] * thickness[dim_other2][b] * thickness[dim_other1][a];
					const int NEDGE = NEDGEBLOCK * mybpd[d];

					const bool needed = NEDGE > 0;
					data.edges[d][b][a] = needed ? _myalloc(sizeof(Real)*NEDGE, 16) : NULL;

					if (!needed) continue;

					int neighbor_index[3];
					neighbor_index[d] = mypeindex[d];
					neighbor_index[dim_other1] = (mypeindex[dim_other1] + 2*a-1 + pesize[dim_other1])%pesize[dim_other1];
					neighbor_index[dim_other2] = (mypeindex[dim_other2] + 2*b-1 + pesize[dim_other2])%pesize[dim_other2];

					if (_myself(neighbor_index)) continue;

					int start[3];
					start[d] = 0;
					start[dim_other1] = blockstart[dim_other1]*(1-a) + a*(blockend[dim_other1]-thickness[dim_other1][1]);
					start[dim_other2] = blockstart[dim_other2]*(1-b) + b*(blockend[dim_other2]-thickness[dim_other2][1]);

					int end[3];
					end[d] = blocksize[d];
					end[dim_other1] = a*blockend[dim_other1] + (1-a)*(blockstart[dim_other1] + thickness[dim_other1][0]);
					end[dim_other2] = b*blockend[dim_other2] + (1-b)*(blockstart[dim_other2] + thickness[dim_other2][0]);

					const int n = bpd[d];
					for(int c=0; c<n; ++c)
					{
						int index[3];
						index[d] = c;
						index[dim_other1] = a*(bpd[dim_other1]-1);
						index[dim_other2] = b*(bpd[dim_other2]-1);

						assert(c2i.find(I3(origin[0] + index[0], origin[1] + index[1], origin[2] + index[2]))!=c2i.end());
						const int blockid = c2i[I3(origin[0] + index[0], origin[1] + index[1], origin[2] + index[2])];

						PackInfo info = {(Real *)globalinfos[blockid].ptrBlock, data.edges[d][b][a] + NEDGEBLOCK*c, start[0], start[1], start[2],  end[0], end[1], end[2]};

						const bool nonempty = end[0]>start[0] && end[1]>start[1] && end[2]>start[2];
						if (nonempty) packinfos.push_back(info);
					}
				}
		}

		//new part
		if (computesubregions)
		{
			for(int dface=0; dface<3; ++dface)
			{
				const int dim_other1face = (dface+1)%3;
				const int dim_other2face = (dface+2)%3;

				for(int s=0; s<2; ++s)
				{
					{
						int neighbor_pe[3];
						neighbor_pe[dface] = (mypeindex[dface] + 2*s-1 + pesize[dface])%pesize[dface];
						neighbor_pe[dim_other1face] = mypeindex[dim_other1face];
						neighbor_pe[dim_other2face] = mypeindex[dim_other2face];

						if (_myself(neighbor_pe)) continue;
					}

					const int n1 = mybpd[dim_other1face];
					const int n2 = mybpd[dim_other2face];

					const int NFACEBLOCK = NC * thickness[dface][s] * blocksize[dim_other1face] * blocksize[dim_other2face];

					int face_start[3];
					face_start[dface] = (1-s)*blockstart[dface] + s*(blockend[dface]-thickness[dface][s]);
					face_start[dim_other1face] = 0;
					face_start[dim_other2face] = 0;

					int face_end[3];
					face_end[dface] = (1-s)*(blockstart[dface] + thickness[dface][s]) + s*blockend[dface];
					face_end[dim_other1face] = blocksize[dim_other1face];
					face_end[dim_other2face] = blocksize[dim_other2face];

					assert(NFACEBLOCK == NC*(face_end[0]-face_start[0])*(face_end[1]-face_start[1])*(face_end[2]-face_start[2]));

					for(int p2=0; p2<n2; ++p2)
						for(int p1=0; p1<n1; ++p1) //iterate over inner face blocks
						{
							int index[3];
							index[dface] = s*(mybpd[dface]-1);
							index[dim_other1face] = p1 ;
							index[dim_other2face] = p2;

							assert(c2i.find(I3(origin[0] + index[0], origin[1] + index[1], origin[2] + index[2]))!=c2i.end());
							const int blockID = c2i[I3(origin[0] + index[0], origin[1] + index[1], origin[2] + index[2])];
							Real * const ptrBlock = (Real*)globalinfos[blockID].ptrBlock;

							for(int dedge=0; dedge<3; ++dedge) //iterate over edges
							{
								const int dim_other1edge = (dedge+1)%3;
								const int dim_other2edge = (dedge+2)%3;

								for(int b=0; b<2; ++b)
									for(int a=0; a<2; ++a)
								{
									{
										int asd[3];
										asd[dedge] = 0;
										asd[dim_other1edge] = a;
										asd[dim_other2edge] = b;

										if (dedge==dface || asd[dface] != s) continue;
									}

									int start[3];
									start[dedge] = 0;
									start[dim_other1edge] = blockstart[dim_other1edge]*(1-a) + a*(blockend[dim_other1edge]-thickness[dim_other1edge][1]);
									start[dim_other2edge] = blockstart[dim_other2edge]*(1-b) + b*(blockend[dim_other2edge]-thickness[dim_other2edge][1]);

									int end[3];
									end[dedge] = blocksize[dedge];
									end[dim_other1edge] = a*blockend[dim_other1edge] + (1-a)*(blockstart[dim_other1edge] + thickness[dim_other1edge][0]);
									end[dim_other2edge] = b*blockend[dim_other2edge] + (1-b)*(blockstart[dim_other2edge] + thickness[dim_other2edge][0]);

									const int vol = max(0, end[2]-start[2])*max(0, end[1]-start[1])*max(0, end[0]-start[0]);
									if (vol == 0) continue;

									int xxx[3];
									xxx[dedge] = 0;
									xxx[dim_other1edge] = 2*a-1;
									xxx[dim_other2edge] = 2*b-1;

									int neighbor[3];
									neighbor[dface] = index[dface];
									neighbor[dedge] = index[dedge];
									neighbor[3-dface-dedge] = index[3-dface-dedge] +  xxx[3-dface-dedge];

									if(c2i.find(I3(origin[0] + neighbor[0], origin[1] + neighbor[1], origin[2] + neighbor[2]))==c2i.end()) continue;

									assert(n1 > neighbor[dim_other1face]);
									assert(n2 > neighbor[dim_other2face]);
									assert(0 <= neighbor[dim_other1face]);
									assert(0 <= neighbor[dim_other2face]);

									{
										const int sregion[3] = {
											start[0]  + (index[0] - neighbor[0])*blocksize[0] - face_start[0],
											start[1]  + (index[1] - neighbor[1])*blocksize[1] - face_start[1],
											start[2]  + (index[2] - neighbor[2])*blocksize[2] - face_start[2]
										};

										const int L[3] = {
											face_end[0] - face_start[0],
											face_end[1] - face_start[1],
											face_end[2] - face_start[2]
										};

										//if (isroot)
//										{
//											printf("-----EDGE --------------->  index: %d %d %d\n", index[0], index[1], index[2]);
//											printf("neighbor: %d %d %d\n", neighbor[0], neighbor[1], neighbor[2]);
//											printf("face: %d %d\n", dface, s);
//											printf("edge: %d %d %d\n", dedge, a, b);
//											printf("facestart: %d %d %d\n", face_start[0], face_start[1], face_start[2]);
//											printf("mystart-end: %d %d %d , %d %d %d\n", start[0], start[1], start[2],  end[0], end[1], end[2]);
//											printf("s: %d %d %d\n",sregion[0], sregion[1], sregion[2]);
//											printf("L: %d %d %d\n",L[0], L[1], L[2]);
//											printf("neighbor p1, p2: %d %d\n", neighbor[dim_other1face], neighbor[dim_other2face]);
//										}

										assert(sregion[0]>= 0);
										assert(sregion[1]>= 0);
										assert(sregion[2]>= 0);

										assert(sregion[0]< L[0]);
										assert(sregion[1]< L[1]);
										assert(sregion[2]< L[2]);

										Real * src_base = data.faces[dface][s] + NFACEBLOCK*(neighbor[dim_other1face] + n1*neighbor[dim_other2face]);

										SubpackInfo subinfo = { ptrBlock, src_base,
												start[0], start[1], start[2], end[0], end[1], end[2],
												sregion[0], sregion[1], sregion[2], L[0], L[1]};

										retval[ptrBlock].push_back(subinfo);
									}
								}
							}

							//iterate over corners
							for(int z=0; z<2; ++z)
							for(int y=0; y<2; ++y)
							for(int x=0; x<2; ++x)
							{
								int xxx[3] = {x,y,z};
								if (xxx[dface] != s) continue;

								const int start[3] = {
									x*(blockend[0] - thickness[0][1]) + (1-x)*blockstart[0],
									y*(blockend[1] - thickness[1][1]) + (1-y)*blockstart[1],
									z*(blockend[2] - thickness[2][1]) + (1-z)*blockstart[2]
								};

								const int end[3] = {
									x*blockend[0] + (1-x)*(thickness[0][0] + blockstart[0]),
									y*blockend[1] + (1-y)*(thickness[1][0] + blockstart[1]),
									z*blockend[2] + (1-z)*(thickness[2][0] + blockstart[2])
								};

								const int vol = max(0, end[2]-start[2])*max(0, end[1]-start[1])*max(0, end[0]-start[0]);
								if (vol == 0) continue;

								int neighbor[3];
								neighbor[0] = index[0] + 2*x-1;
								neighbor[1] = index[1] + 2*y-1;
								neighbor[2] = index[2] + 2*z-1;
								neighbor[dface] = index[dface];

								if(c2i.find(I3(origin[0] + neighbor[0], origin[1] + neighbor[1], origin[2] + neighbor[2]))==c2i.end()) continue;

								assert(n1 > neighbor[dim_other1face]);
								assert(n2 > neighbor[dim_other2face]);
								assert(0 <= neighbor[dim_other1face]);
								assert(0 <= neighbor[dim_other2face]);

								{
									const int sregion[3] = {
										start[0]  + (index[0] - neighbor[0])*blocksize[0] - face_start[0],
										start[1]  + (index[1] - neighbor[1])*blocksize[1] - face_start[1],
										start[2]  + (index[2] - neighbor[2])*blocksize[2] - face_start[2]
									};

									const int L[3] = {
										face_end[0] - face_start[0],
										face_end[1] - face_start[1],
										face_end[2] - face_start[2]
									};

								//	if (isroot)
//									{
//										printf("---CORNER ----------------->  index: %d %d %d\n", index[0], index[1], index[2]);
//										printf("neighbor: %d %d %d\n", neighbor[0], neighbor[1], neighbor[2]);
//										printf("face: %d %d\n", dface, s);
//										printf("corner: %d %d %d\n", x, y, z);
//										printf("facestart: %d %d %d\n", face_start[0], face_start[1], face_start[2]);
//										printf("mystart: %d %d %d\n", start[0], start[1], start[2]);
//										printf("s: %d %d %d\n",sregion[0], sregion[1], sregion[2]);
//										printf("L: %d %d %d\n",L[0], L[1], L[2]);
//										printf("neighbor p1, p2: %d %d\n", neighbor[dim_other1face], neighbor[dim_other2face]);
//									}
									assert(c2i.find(I3(origin[0] + neighbor[0], origin[1] + neighbor[1], origin[2] + neighbor[2]))!=c2i.end());
									assert(sregion[0]>= 0);
									assert(sregion[1]>= 0);
									assert(sregion[2]>= 0);

									assert(sregion[0]< L[0]);
									assert(sregion[1]< L[1]);
									assert(sregion[2]< L[2]);

									Real * src_base = data.faces[dface][s] + NFACEBLOCK*(neighbor[dim_other1face] + n1*neighbor[dim_other2face]);

									SubpackInfo subinfo = { ptrBlock, src_base,
										start[0], start[1], start[2], end[0], end[1], end[2],
										sregion[0], sregion[1], sregion[2], L[0], L[1]};

									retval[ptrBlock].push_back(subinfo);
								}
							}
						}
				}
			}

			for(int d=0; d<3; ++d)
			{
				const int dim_other1 = (d+1)%3;
				const int dim_other2 = (d+2)%3;

				for(int b=0; b<2; ++b)
					for(int a=0; a<2; ++a)
					{
						{
							int neighbor_pe[3];
							neighbor_pe[d] = mypeindex[d];
							neighbor_pe[dim_other1] = (mypeindex[dim_other1] + 2*a-1 + pesize[dim_other1])%pesize[dim_other1];
							neighbor_pe[dim_other2] = (mypeindex[dim_other2] + 2*b-1 + pesize[dim_other2])%pesize[dim_other2];

							if (_myself(neighbor_pe)) continue;
						}

						const int n = bpd[d];
						const int NEDGEBLOCK = NC * blocksize[d] * thickness[dim_other2][b] * thickness[dim_other1][a];

						int edge_start[3];
						edge_start[d] = 0;
						edge_start[dim_other1] = blockstart[dim_other1]*(1-a) + a*(blockend[dim_other1]-thickness[dim_other1][1]);
						edge_start[dim_other2] = blockstart[dim_other2]*(1-b) + b*(blockend[dim_other2]-thickness[dim_other2][1]);

						int edge_end[3];
						edge_end[d] = blocksize[d];
						edge_end[dim_other1] = a*blockend[dim_other1] + (1-a)*(blockstart[dim_other1] + thickness[dim_other1][0]);
						edge_end[dim_other2] = b*blockend[dim_other2] + (1-b)*(blockstart[dim_other2] + thickness[dim_other2][0]);

						assert(NEDGEBLOCK == NC*(edge_end[0]-edge_start[0])*(edge_end[1]-edge_start[1])*(edge_end[2]-edge_start[2]));

						for(int p1=0; p1<n; ++p1) //iterate over inner edge blocks
						{
							int index[3];
							index[d] = p1;
							index[dim_other1] = a*(bpd[dim_other1]-1);
							index[dim_other2] = b*(bpd[dim_other2]-1);

							assert(c2i.find(I3(origin[0] + index[0], origin[1] + index[1], origin[2] + index[2]))!=c2i.end());
							const int blockID = c2i[I3(origin[0] + index[0], origin[1] + index[1], origin[2] + index[2])];
							Real * const ptrBlock = (Real*)globalinfos[blockID].ptrBlock;

							for(int z=0; z<2; ++z) //iterate over corners
								for(int y=0; y<2; ++y)
									for(int x=0; x<2; ++x)
									{
										int xxx[3] = {x,y,z};
										if (xxx[dim_other1] != a || xxx[dim_other2] != b) continue;

										const int start[3] = {
											x*(blockend[0] - thickness[0][1]) + (1-x)*blockstart[0],
											y*(blockend[1] - thickness[1][1]) + (1-y)*blockstart[1],
											z*(blockend[2] - thickness[2][1]) + (1-z)*blockstart[2]
										};

										const int end[3] = {
											x*blockend[0] + (1-x)*(thickness[0][0] + blockstart[0]),
											y*blockend[1] + (1-y)*(thickness[1][0] + blockstart[1]),
											z*blockend[2] + (1-z)*(thickness[2][0] + blockstart[2])
										};

										const int vol = max(0, end[2]-start[2])*max(0, end[1]-start[1])*max(0, end[0]-start[0]);
										if (vol == 0) continue;

										int neighbor[3];
										neighbor[0] = index[0];
										neighbor[1] = index[1];
										neighbor[2] = index[2];
										neighbor[d] = index[d] + xxx[d]*2-1;
										if(c2i.find(I3(origin[0] + neighbor[0], origin[1] + neighbor[1], origin[2] + neighbor[2]))==c2i.end()) continue;

										assert(n > neighbor[d]);
										assert(0 <= neighbor[d]);

										{
											const int sregion[3] = {
												start[0]  + (index[0] - neighbor[0])*blocksize[0] - edge_start[0],
												start[1]  + (index[1] - neighbor[1])*blocksize[1] - edge_start[1],
												start[2]  + (index[2] - neighbor[2])*blocksize[2] - edge_start[2]
											};

											const int L[3] = {
												edge_end[0] - edge_start[0],
												edge_end[1] - edge_start[1],
												edge_end[2] - edge_start[2]
											};

										//	if (isroot)
//											{
//												printf("---CORNER (from edge) ----------------->  index: %d %d %d\n", index[0], index[1], index[2]);
//												printf("neighbor: %d %d %d\n", neighbor[0], neighbor[1], neighbor[2]);
//												printf("edge: %d %d %d\n", d, a, b);
//												printf("corner: %d %d %d\n", x, y, z);
//												printf("edgestart: %d %d %d\n", edge_start[0], edge_start[1], edge_start[2]);
//												printf("mystart: %d %d %d\n", start[0], start[1], start[2]);
//												printf("s: %d %d %d\n",sregion[0], sregion[1], sregion[2]);
//												printf("L: %d %d %d\n",L[0], L[1], L[2]);
//												printf("neighbor p1: %d\n", neighbor[d]);
//											}
											assert(c2i.find(I3(origin[0] + neighbor[0], origin[1] + neighbor[1], origin[2] + neighbor[2]))!=c2i.end());
											assert(sregion[0]>= 0);
											assert(sregion[1]>= 0);
											assert(sregion[2]>= 0);

											assert(sregion[0]< L[0]);
											assert(sregion[1]< L[1]);
											assert(sregion[2]< L[2]);
											assert(vol <NEDGEBLOCK);

											//Real * src_base = data.faces[dface][s] + NFACEBLOCK*(neighbor[dim_other1face] + n1*neighbor[dim_other2face]);
											Real * src_base = data.edges[d][b][a] + NEDGEBLOCK*neighbor[d];

											SubpackInfo subinfo = { ptrBlock, src_base,
												start[0], start[1], start[2], end[0], end[1], end[2],
												sregion[0], sregion[1], sregion[2], L[0], L[1]};

											retval[ptrBlock].push_back(subinfo);
										}
									}
						}
					}
			}
		}

		//corners
		for(int z=0; z<2; ++z)
			for(int y=0; y<2; ++y)
				for(int x=0; x<2; ++x)
				{
					const int NCORNERBLOCK = NC * thickness[0][x]*thickness[1][y]*thickness[2][z];

					const bool needed = NCORNERBLOCK > 0;
					data.corners[z][y][x] = needed ? _myalloc(sizeof(Real)*NCORNERBLOCK, 16) : NULL;

					if (!needed) continue;

					int neighbor_index[3];
					neighbor_index[0] = (mypeindex[0] + 2*x-1 + pesize[0])%pesize[0];
					neighbor_index[1] = (mypeindex[1] + 2*y-1 + pesize[1])%pesize[1];
					neighbor_index[2] = (mypeindex[2] + 2*z-1 + pesize[2])%pesize[2];

					if (_myself(neighbor_index)) continue;

					const int start[3] = {
						x*(blockend[0] - thickness[0][1]) + (1-x)*blockstart[0],
						y*(blockend[1] - thickness[1][1]) + (1-y)*blockstart[1],
						z*(blockend[2] - thickness[2][1]) + (1-z)*blockstart[2]
					};

					const int end[3] = {
						x*blockend[0] + (1-x)*(thickness[0][0] + blockstart[0]),
						y*blockend[1] + (1-y)*(thickness[1][0] + blockstart[1]),
						z*blockend[2] + (1-z)*(thickness[2][0] + blockstart[2])
					};

					const int index[3] = {
						x*(bpd[0]-1),
						y*(bpd[1]-1),
						z*(bpd[2]-1),
					};

					assert(c2i.find(I3(origin[0] + index[0], origin[1] + index[1], origin[2] + index[2]))!=c2i.end());
					const int blockid = c2i[I3(origin[0] + index[0], origin[1] + index[1], origin[2] + index[2])];

					PackInfo info = {(Real *)globalinfos[blockid].ptrBlock, data.corners[z][y][x], start[0], start[1], start[2],  end[0], end[1], end[2]};

					const bool nonempty = end[0]>start[0] && end[1]>start[1] && end[2]>start[2];
					if (nonempty) packinfos.push_back(info);
				}

		return retval;
	}

	Real * _myalloc(const int NBYTES, const int ALIGN)
	{
		if (NBYTES>0)
		{
			//Real * ret_val = (Real *)_mm_malloc(NBYTES, ALIGN);
			Real * ret_val = NULL;

			int error = posix_memalign((void**)&ret_val, std::max(8, ALIGN), NBYTES);

			assert(error == 0);

			all_mallocs.push_back(ret_val);

			return ret_val;
		}

		return NULL;
	}

	void _myfree(Real *& ptr) {if (ptr!=NULL) { free(ptr); ptr=NULL;} }

	//forbidden methods
	SynchronizerMPI(const SynchronizerMPI& c):cube(-1,-1,-1), synchID(-1), isroot(true){ abort(); }

	void operator=(const SynchronizerMPI& c){ abort(); }

public:

	SynchronizerMPI(const int synchID, StencilInfo stencil, vector<BlockInfo> globalinfos, MPI_Comm cartcomm, const int mybpd[3], const int blocksize[3]):
	synchID(synchID), stencil(stencil), globalinfos(globalinfos), cube(mybpd[0], mybpd[1], mybpd[2]), cartcomm(cartcomm)
	{
		int myrank;
        MPI_Comm_rank(cartcomm, &myrank);
        isroot = (myrank == 0);

        MPI_Cart_get(cartcomm, 3, pesize, periodic, mypeindex);
        MPI_Cart_coords(cartcomm, myrank, 3, mypeindex);

		for(int iz=0; iz<3; iz++)
			for(int iy=0; iy<3; iy++)
				for(int ix=0; ix<3; ix++)
				{
					int s[3] = { ix-1+mypeindex[0], iy-1+mypeindex[1], iz-1+mypeindex[2]};
                    int nbrRank;
                    MPI_Cart_rank(cartcomm, s, &nbrRank);
					neighborsrank[iz][iy][ix] = nbrRank;
				}

		for(int i=0; i<3; ++i) this->mybpd[i]=mybpd[i];
		for(int i=0; i<3; ++i) this->blocksize[i]=blocksize[i];

		for(int i=0; i< globalinfos.size(); ++i)
		{
			I3 coord(globalinfos[i].index[0], globalinfos[i].index[1], globalinfos[i].index[2]);
			c2i[coord] = i;
		}

		const int origin[3] = {
			mypeindex[0]*mybpd[0],
			mypeindex[1]*mybpd[1],
			mypeindex[2]*mybpd[2]
		};

		const int s[3] = {stencil.sx, stencil.sy, stencil.sz};
		const int e[3] = {stencil.ex, stencil.ey, stencil.ez};
		const int z[3] = {0, 0, 0};

		send_thickness[0][0] = e[0] - 1;  send_thickness[0][1] = -s[0];
		send_thickness[1][0] = e[1] - 1;  send_thickness[1][1] = -s[1];
		send_thickness[2][0] = e[2] - 1;  send_thickness[2][1] = -s[2];

		_setup<false>(send, send_thickness, z, blocksize, origin, send_packinfos);

		recv_thickness[0][0] = -s[0]; recv_thickness[0][1] = e[0] - 1;
		recv_thickness[1][0] = -s[1]; recv_thickness[1][1] = e[1] - 1;
		recv_thickness[2][0] = -s[2]; recv_thickness[2][1] = e[2] - 1;

		{
			const int blockstart[3] = {
				stencil.sx ,
				stencil.sy ,
				stencil.sz
			};

			const int blockend[3] = {
				stencil.ex + blocksize[0]-1,
				stencil.ey + blocksize[1]-1,
				stencil.ez + blocksize[2]-1
			};

			vector<PackInfo> packinfos;
			recv_subpackinfos = _setup<true>(recv, recv_thickness, blockstart, blockend, origin, packinfos);

			for(vector<PackInfo>::const_iterator it = packinfos.begin(); it<packinfos.end(); ++it)
				recv_packinfos[it->block].push_back(*it);
		}

		assert(recv.pending.size() == 0);
		assert(send.pending.size() == 0);
	}

	~SynchronizerMPI()
	{
		for(int i=0;i<all_mallocs.size();++i)
			_myfree(all_mallocs[i]);
	}

	virtual void sync(unsigned int gptfloats, MPI_Datatype MPIREAL, const int timestamp)
	{
		//0. wait for pending sends, couple of checks
		//1. pack all stuff
		//2. perform send/receive requests
		//3. setup the dependency

		//0.
		{
			const int NPENDINGSENDS = send.pending.size();
			if (NPENDINGSENDS > 0)
			{
				vector<MPI_Request> pending(NPENDINGSENDS);
				copy(send.pending.begin(), send.pending.end(), pending.begin());
#if 1
                MPI_Waitall(NPENDINGSENDS, &pending.front(), MPI_STATUSES_IGNORE);
#else
				int done = false;
				while (1)
				{
					MPI_Testall(NPENDINGSENDS, &pending.front(), &done, MPI_STATUSES_IGNORE);
					if (done) break;
					sched_yield();
				};
#endif

				send.pending.clear();
			}
		}

		assert(recv.pending.size() == 0);
		assert(send.pending.size() == 0);

		cube.prepare();
		blockinfo_counter = globalinfos.size();
		const int NC = stencil.selcomponents.size();

		//1. pack
		{
			const int N = send_packinfos.size();

			vector<int> selcomponents = stencil.selcomponents;
			sort(selcomponents.begin(), selcomponents.end());

			const bool contiguous = false;//selcomponents.back()+1-selcomponents.front() == selcomponents.size();

			if (!contiguous)
			{
#pragma omp parallel for schedule(runtime)
				for(int i=0; i<N; ++i)
				{
					PackInfo info = send_packinfos[i];
					pack(info.block, info.pack, gptfloats, &selcomponents.front(), NC, info.sx, info.sy, info.sz, info.ex, info.ey, info.ez);
				}
			}
			else
			{
				const int selstart = selcomponents.front();
				const int selend = selcomponents.back()+1;

#pragma omp parallel for schedule(runtime)
				for(int i=0; i<N; ++i)
				{
					PackInfo info = send_packinfos[i];
					pack_stripes(info.block, info.pack, gptfloats, selstart, selend, info.sx, info.sy, info.sz, info.ex, info.ey, info.ez);
				}
			}
		}

		//2. send requests
		{
			//faces
			for(int d=0; d<3; ++d)
			{
				if (!_face_needed(d)) continue;

				const int dim_other1 = (d+1)%3;
				const int dim_other2 = (d+2)%3;

				for(int s=0; s<2; ++s)
				{
					const int NFACEBLOCK_SEND = NC * send_thickness[d][s] * blocksize[dim_other1] * blocksize[dim_other2];
					const int NFACEBLOCK_RECV = NC * recv_thickness[d][s] * blocksize[dim_other1] * blocksize[dim_other2];
					const int NFACE_SEND = NFACEBLOCK_SEND * mybpd[dim_other1] * mybpd[dim_other2];
					const int NFACE_RECV = NFACEBLOCK_RECV * mybpd[dim_other1] * mybpd[dim_other2];

					int neighbor_index[3];
					neighbor_index[d] = (mypeindex[d] + 2*s-1 + pesize[d])%pesize[d];
					neighbor_index[dim_other1] = mypeindex[dim_other1];
					neighbor_index[dim_other2] = mypeindex[dim_other2];

					if (_myself(neighbor_index)) continue;

                    if (NFACE_SEND > 0)
                    {
                        MPI_Request req;
                        MPI_Isend(send.faces[d][s], NFACE_SEND, MPIREAL, _rank(neighbor_index), 6*timestamp + 2*d + 1-s, cartcomm, &req);
                        send.pending.insert( req );
                    }

					if (NFACE_RECV > 0)
					{
                        MPI_Request rc;
						MPI_Irecv(recv.faces[d][s], NFACE_RECV, MPIREAL, _rank(neighbor_index), 6*timestamp + 2*d + s, cartcomm, &rc);
						recv.pending.insert(rc);
						cube.face(rc, d, s);
					}
				}
			}

			if (stencil.tensorial)
			{
				//edges
				for(int d=0; d<3; ++d)
				{
					const int dim_other1 = (d+1)%3;
					const int dim_other2 = (d+2)%3;

					for(int b=0; b<2; ++b)
						for(int a=0; a<2; ++a)
						{
							const int NEDGEBLOCK_SEND = NC * blocksize[d] * send_thickness[dim_other2][b] * send_thickness[dim_other1][a];
							const int NEDGEBLOCK_RECV = NC * blocksize[d] * recv_thickness[dim_other2][b] * recv_thickness[dim_other1][a];
							const int NEDGE_SEND = NEDGEBLOCK_SEND * mybpd[d];
							const int NEDGE_RECV = NEDGEBLOCK_RECV * mybpd[d];

							int neighbor_index[3];
							neighbor_index[d] = mypeindex[d];
							neighbor_index[dim_other1] = (mypeindex[dim_other1] + 2*a-1 + pesize[dim_other1])%pesize[dim_other1];
							neighbor_index[dim_other2] = (mypeindex[dim_other2] + 2*b-1 + pesize[dim_other2])%pesize[dim_other2];

							if (_myself(neighbor_index)) continue;

							if (NEDGE_RECV > 0)
							{
                                MPI_Request rc;
								MPI_Irecv(recv.edges[d][b][a], NEDGE_RECV, MPIREAL, _rank(neighbor_index), 12*timestamp + 4*d + 2*b + a, cartcomm, &rc);

								recv.pending.insert(rc);

								cube.edge(rc, d, a, b);
							}

                            if (NEDGE_SEND > 0)
                            {
                                MPI_Request req;
                                MPI_Isend(send.edges[d][b][a], NEDGE_SEND, MPIREAL, _rank(neighbor_index), 12*timestamp + 4*d + 2*(1-b) + (1-a), cartcomm, &req);
                                send.pending.insert(req);
                            }
						}
				}

				//corners
				{
					for(int z=0; z<2; ++z)
						for(int y=0; y<2; ++y)
							for(int x=0; x<2; ++x)
								{
									const int NCORNERBLOCK_SEND = NC * send_thickness[0][x]*send_thickness[1][y]*send_thickness[2][z];
									const int NCORNERBLOCK_RECV = NC * recv_thickness[0][x]*recv_thickness[1][y]*recv_thickness[2][z];

									int neighbor_index[3];
									neighbor_index[0] = (mypeindex[0] + 2*x-1 + pesize[0])%pesize[0];
									neighbor_index[1] = (mypeindex[1] + 2*y-1 + pesize[1])%pesize[1];
									neighbor_index[2] = (mypeindex[2] + 2*z-1 + pesize[2])%pesize[2];

									if (_myself(neighbor_index)) continue;

									if (NCORNERBLOCK_RECV)
									{
                                        MPI_Request rc;
										MPI_Irecv(recv.corners[z][y][x], NCORNERBLOCK_RECV, MPIREAL, _rank(neighbor_index), 8*timestamp + 4*z + 2*y + x, cartcomm, &rc);

										recv.pending.insert(rc);

										cube.corner(rc, x, y, z);
									}

									if (NCORNERBLOCK_SEND)
                                    {
                                        MPI_Request req;
                                        MPI_Isend(send.corners[z][y][x], NCORNERBLOCK_SEND, MPIREAL, _rank(neighbor_index), 8*timestamp + 4*(1-z) + 2*(1-y) + (1-x), cartcomm, &req);
                                        send.pending.insert(req);
                                    }
								}
				}
            }
		}

		//3.
		cube.make_dependencies(isroot);
	}

    //peh
	virtual void sync0(unsigned int gptfloats, MPI_Datatype MPIREAL, const int timestamp)
	{
		double t0, t1;

		//0. wait for pending sends, couple of checks
		//1. pack all stuff
		//2. perform send/receive requests
		//3. setup the dependency

		//0.
		{
			const int NPENDINGSENDS = send.pending.size();
			if (NPENDINGSENDS > 0)
			{
				vector<MPI_Request> pending(NPENDINGSENDS);
				copy(send.pending.begin(), send.pending.end(), pending.begin());
#if 1
                MPI_Waitall(NPENDINGSENDS, &pending.front(), MPI_STATUSES_IGNORE);
#else
				int done = false;
				while (1)
				{
					MPI_Testall(NPENDINGSENDS, &pending.front(), &done, MPI_STATUSES_IGNORE);
					if (done) break;
					pthread_yield();
				};
#endif

				send.pending.clear();
			}
		}

		assert(recv.pending.size() == 0);
		assert(send.pending.size() == 0);

		cube.prepare();

		blockinfo_counter = globalinfos.size();
		const int NC = stencil.selcomponents.size();

		//1. pack
		{
			double t0, t1;


			const int N = send_packinfos.size();

			vector<int> selcomponents = stencil.selcomponents;
			sort(selcomponents.begin(), selcomponents.end());

            const bool contiguous = false;//true; //false; //true;//selcomponents.back()+1-selcomponents.front() == selcomponents.size();

			if (!contiguous)
			{
#pragma omp parallel for schedule(runtime)
				for(int i=0; i<N; ++i)
				{
					PackInfo info = send_packinfos[i];
					pack(info.block, info.pack, gptfloats, &selcomponents.front(), NC, info.sx, info.sy, info.sz, info.ex, info.ey, info.ez);
				}
			}
			else
			{
				const int selstart = selcomponents.front();
				const int selend = selcomponents.back()+1;

#pragma omp parallel for //schedule(runtime)
				for(int i=0; i<N; ++i)
				{
					PackInfo info = send_packinfos[i];
					pack_stripes(info.block, info.pack, gptfloats, selstart, selend, info.sx, info.sy, info.sz, info.ex, info.ey, info.ez);
				}
			}

		}

		// recvs
		for (int pass = 0; pass <= 1; pass++)
		//2. send requests
		{
			//faces
			for(int d=0; d<3; ++d)
			{
				if (!_face_needed(d)) continue;

				const int dim_other1 = (d+1)%3;
				const int dim_other2 = (d+2)%3;

				for(int s=0; s<2; ++s)
				{
					int neighbor_index[3];
					neighbor_index[d] = (mypeindex[d] + 2*s-1 + pesize[d])%pesize[d];
					neighbor_index[dim_other1] = mypeindex[dim_other1];
					neighbor_index[dim_other2] = mypeindex[dim_other2];

					if (_myself(neighbor_index)) continue;

					if (pass == 0)
					{
						const int NFACEBLOCK_RECV = NC * recv_thickness[d][s] * blocksize[dim_other1] * blocksize[dim_other2];
						const int NFACE_RECV = NFACEBLOCK_RECV * mybpd[dim_other1] * mybpd[dim_other2];

						if (NFACE_RECV > 0)
						{
                            MPI_Request rc;
							MPI_Irecv(recv.faces[d][s], NFACE_RECV, MPIREAL, _rank(neighbor_index), 6*timestamp + 2*d + s, cartcomm, &rc);
							recv.pending.insert(rc);
							cube.face(rc, d, s);
						}
					}
					else
					{
						const int NFACEBLOCK_SEND = NC * send_thickness[d][s] * blocksize[dim_other1] * blocksize[dim_other2];
						const int NFACE_SEND = NFACEBLOCK_SEND * mybpd[dim_other1] * mybpd[dim_other2];

						if (NFACE_SEND > 0)
                        {
                            MPI_Request req;
                            MPI_Isend(send.faces[d][s], NFACE_SEND, MPIREAL, _rank(neighbor_index), 6*timestamp + 2*d + 1-s, cartcomm, &req);
                            send.pending.insert(req);
                        }
					}
				}
			}

			if (stencil.tensorial)
			{
				//edges
				for(int d=0; d<3; ++d)
				{
					const int dim_other1 = (d+1)%3;
					const int dim_other2 = (d+2)%3;

					for(int b=0; b<2; ++b)
						for(int a=0; a<2; ++a)
						{
							const int NEDGEBLOCK_SEND = NC * blocksize[d] * send_thickness[dim_other2][b] * send_thickness[dim_other1][a];
							const int NEDGEBLOCK_RECV = NC * blocksize[d] * recv_thickness[dim_other2][b] * recv_thickness[dim_other1][a];
							const int NEDGE_SEND = NEDGEBLOCK_SEND * mybpd[d];
							const int NEDGE_RECV = NEDGEBLOCK_RECV * mybpd[d];

							int neighbor_index[3];
							neighbor_index[d] = mypeindex[d];
							neighbor_index[dim_other1] = (mypeindex[dim_other1] + 2*a-1 + pesize[dim_other1])%pesize[dim_other1];
							neighbor_index[dim_other2] = (mypeindex[dim_other2] + 2*b-1 + pesize[dim_other2])%pesize[dim_other2];

							if (_myself(neighbor_index)) continue;

							if (pass == 0)
							{
								if (NEDGE_RECV > 0)
								{
                                    MPI_Request rc;
									MPI_Irecv(recv.edges[d][b][a], NEDGE_RECV, MPIREAL, _rank(neighbor_index), 12*timestamp + 4*d + 2*b + a, cartcomm, &rc);

									recv.pending.insert(rc);

									cube.edge(rc, d, a, b);
								}
							}
							else
							{
                                if (NEDGE_SEND > 0)
                                {
                                    MPI_Request req;
                                    MPI_Isend(send.edges[d][b][a], NEDGE_SEND, MPIREAL, _rank(neighbor_index), 12*timestamp + 4*d + 2*(1-b) + (1-a), cartcomm, &req);
                                    send.pending.insert(req);
                                }
							}
						}
				}

				//corners
				{
					for(int z=0; z<2; ++z)
						for(int y=0; y<2; ++y)
							for(int x=0; x<2; ++x)
								{
									const int NCORNERBLOCK_SEND = NC * send_thickness[0][x]*send_thickness[1][y]*send_thickness[2][z];
									const int NCORNERBLOCK_RECV = NC * recv_thickness[0][x]*recv_thickness[1][y]*recv_thickness[2][z];

									int neighbor_index[3];
									neighbor_index[0] = (mypeindex[0] + 2*x-1 + pesize[0])%pesize[0];
									neighbor_index[1] = (mypeindex[1] + 2*y-1 + pesize[1])%pesize[1];
									neighbor_index[2] = (mypeindex[2] + 2*z-1 + pesize[2])%pesize[2];

									if (_myself(neighbor_index)) continue;

									if (pass == 0)
									{
										if (NCORNERBLOCK_RECV)
										{
                                            MPI_Request rc;
											MPI_Irecv(recv.corners[z][y][x], NCORNERBLOCK_RECV, MPIREAL, _rank(neighbor_index), 8*timestamp + 4*z + 2*y + x, cartcomm, &rc);

											recv.pending.insert(rc);

											cube.corner(rc, x, y, z);
										}
									}
									else
									{
										if (NCORNERBLOCK_SEND)
                                        {
                                            MPI_Request req;
											MPI_Isend(send.corners[z][y][x], NCORNERBLOCK_SEND, MPIREAL, _rank(neighbor_index), 8*timestamp + 4*(1-z) + 2*(1-y) + (1-x), cartcomm, &req);
											send.pending.insert(req);
                                        }
									}
								}
				}
			}
		}

		//3.
		cube.make_dependencies(isroot);
	}

	vector<BlockInfo> avail_inner()
	{
		vector<BlockInfo> retval;

		const int xorigin = mypeindex[0]*mybpd[0];
		const int yorigin =	mypeindex[1]*mybpd[1];
		const int zorigin =	mypeindex[2]*mybpd[2];

		vector<Region> regions = cube.avail();

		for(vector<Region>::const_iterator it=regions.begin(); it!=regions.end(); ++it)
		{
            map<Region, vector<BlockInfo> >::const_iterator r2v = region2infos.find(*it);

            if(r2v!=region2infos.end())
            {
                retval.insert(retval.end(), r2v->second.begin(), r2v->second.end());
                blockinfo_counter -=  r2v->second.size();
            }
            else
            {
                vector<BlockInfo> entry;

                const int sx = it->s[0];
                const int sy = it->s[1];
                const int sz = it->s[2];
                const int ex = it->e[0];
                const int ey = it->e[1];
                const int ez = it->e[2];

                for(int iz=sz; iz<ez; ++iz)
                    for(int iy=sy; iy<ey; ++iy)
                        for(int ix=sx; ix<ex; ++ix, blockinfo_counter--)
                        {
                            assert(c2i.find(I3(ix + xorigin, iy + yorigin, iz + zorigin)) != c2i.end());
                            entry.push_back(globalinfos[ c2i[I3(ix + xorigin, iy + yorigin, iz + zorigin)] ]);
                        }

                retval.insert(retval.end(), entry.begin(), entry.end());

                region2infos[*it] = entry;
            }
		}

		assert(cube.pendingcount() != 0 || blockinfo_counter == cube.pendingcount());
		assert(blockinfo_counter != 0 || blockinfo_counter == cube.pendingcount());
		assert(blockinfo_counter != 0 || recv.pending.size() == 0);

		return retval;
	}

	vector<BlockInfo> avail_halo()
	{
		vector<BlockInfo> retval;

		const int NPENDING = recv.pending.size();

		vector<MPI_Request> pending(NPENDING);

		copy(recv.pending.begin(), recv.pending.end(), pending.begin());

		vector<MPI_Request> old = pending;

#if 1
        MPI_Waitall(NPENDING, &pending.front(), MPI_STATUSES_IGNORE);
#else
		int done = false;
		while (1)
		{
            MPI_Testall(NPENDING, &pending.front(), &done, MPI_STATUSES_IGNORE);
			if (done) break;
			pthread_yield();
		};
#endif
		for(int i=0; i<NPENDING; ++i)
		{
			cube.received(old[i]);
			recv.pending.erase(old[i]);
		}

		const int xorigin = mypeindex[0]*mybpd[0];
		const int yorigin =	mypeindex[1]*mybpd[1];
		const int zorigin =	mypeindex[2]*mybpd[2];

		vector<Region> regions = cube.avail();

		for(vector<Region>::const_iterator it=regions.begin(); it!=regions.end(); ++it)
		{
            map<Region, vector<BlockInfo> >::const_iterator r2v = region2infos.find(*it);

            if(r2v!=region2infos.end())
            {
                retval.insert(retval.end(), r2v->second.begin(), r2v->second.end());
                blockinfo_counter -=  r2v->second.size();
            }
            else
            {
                vector<BlockInfo> entry;

                const int sx = it->s[0];
                const int sy = it->s[1];
                const int sz = it->s[2];
                const int ex = it->e[0];
                const int ey = it->e[1];
                const int ez = it->e[2];

                for(int iz=sz; iz<ez; ++iz)
                    for(int iy=sy; iy<ey; ++iy)
                        for(int ix=sx; ix<ex; ++ix, blockinfo_counter--)
                        {
                            assert(c2i.find(I3(ix + xorigin, iy + yorigin, iz + zorigin)) != c2i.end());
                            entry.push_back(globalinfos[ c2i[I3(ix + xorigin, iy + yorigin, iz + zorigin)] ]);
                        }

                retval.insert(retval.end(), entry.begin(), entry.end());

                region2infos[*it] = entry;
            }
		}

		assert(cube.pendingcount() != 0 || blockinfo_counter == cube.pendingcount());
		assert(blockinfo_counter != 0 || blockinfo_counter == cube.pendingcount());
		assert(blockinfo_counter != 0 || recv.pending.size() == 0);

		return retval;
	}


	bool test_halo()
	{
		vector<BlockInfo> retval;

		const int NPENDING = recv.pending.size();

		if (NPENDING == 0) return true;

		vector<MPI_Request> pending(NPENDING);

		copy(recv.pending.begin(), recv.pending.end(), pending.begin());

		int done = false;
        MPI_Testall(NPENDING, &pending.front(), &done, MPI_STATUSES_IGNORE);

		return done;
	}

	vector<BlockInfo> avail()
	{
		vector<BlockInfo> retval;

		const int NPENDING = recv.pending.size();

		vector<MPI_Request> pending(NPENDING);

		copy(recv.pending.begin(), recv.pending.end(), pending.begin());

		vector<MPI_Request> old = pending;

		if(NPENDING > 0)
		{
			if(mybpd[0]==1 || mybpd[1]==1 || mybpd[2] == 1) //IS THERE SOMETHING MORE INTELLIGENT?!
			{
                MPI_Waitall(NPENDING, &pending.front(), MPI_STATUSES_IGNORE);
				for(int i=0; i<NPENDING; ++i)
				{
					cube.received(old[i]);
					recv.pending.erase(old[i]);
				}
			}
			else
			{
				vector<int> indices(NPENDING);
				int NSOLVED = 0;
				if (blockinfo_counter == globalinfos.size())
                    MPI_Testsome(NPENDING, &pending.front(), &NSOLVED, &indices.front(), MPI_STATUSES_IGNORE);
				else
				{
					MPI_Waitsome(NPENDING, &pending.front(), &NSOLVED, &indices.front(), MPI_STATUSES_IGNORE);
					assert(NSOLVED > 0);
				}

				for(int i=0; i<NSOLVED; ++i)
				{
					cube.received(old[indices[i]]);
					recv.pending.erase(old[indices[i]]);
				}
			}
		}

		const int xorigin = mypeindex[0]*mybpd[0];
		const int yorigin =	mypeindex[1]*mybpd[1];
		const int zorigin =	mypeindex[2]*mybpd[2];

		vector<Region> regions = cube.avail();

		for(vector<Region>::const_iterator it=regions.begin(); it!=regions.end(); ++it)
		{
            map<Region, vector<BlockInfo> >::const_iterator r2v = region2infos.find(*it);

            if(r2v!=region2infos.end())
            {
                retval.insert(retval.end(), r2v->second.begin(), r2v->second.end());
                blockinfo_counter -=  r2v->second.size();
            }
            else
            {
                vector<BlockInfo> entry;

                const int sx = it->s[0];
                const int sy = it->s[1];
                const int sz = it->s[2];
                const int ex = it->e[0];
                const int ey = it->e[1];
                const int ez = it->e[2];

                for(int iz=sz; iz<ez; ++iz)
                    for(int iy=sy; iy<ey; ++iy)
                        for(int ix=sx; ix<ex; ++ix, blockinfo_counter--)
                        {
                            assert(c2i.find(I3(ix + xorigin, iy + yorigin, iz + zorigin)) != c2i.end());
                            entry.push_back(globalinfos[ c2i[I3(ix + xorigin, iy + yorigin, iz + zorigin)] ]);
                        }

                retval.insert(retval.end(), entry.begin(), entry.end());

                region2infos[*it] = entry;
            }
		}

		assert(cube.pendingcount() != 0 || blockinfo_counter == cube.pendingcount());
		assert(blockinfo_counter != 0 || blockinfo_counter == cube.pendingcount());
		assert(blockinfo_counter != 0 || recv.pending.size() == 0);

		return retval;
	}

	vector<BlockInfo> avail(const int smallest)
	{
		vector<BlockInfo> accumulator;

		while(accumulator.size()<smallest && !done())
		{
			const vector<BlockInfo> r = avail();

			accumulator.insert(accumulator.end(), r.begin(), r.end());
		}

		return accumulator;
	}

	bool done() const
	{
		assert(!(blockinfo_counter == 0) || recv.pending.size() == 0);

		return blockinfo_counter == 0;
	}

	StencilInfo getstencil() const
	{
		return stencil;
	}

	void getpedata(int mypeindex[3], int pesize[3], int mybpd[3]) const
	{
		for(int i=0; i<3; ++i) mypeindex[i] = this->mypeindex[i];
		for(int i=0; i<3; ++i) pesize[i] = this->pesize[i];
		for(int i=0; i<3; ++i) mybpd[i] = this->mybpd[i];
	}

class MyRange
{
  const int sx, sy, sz, ex, ey, ez;

 public:

 MyRange(const int sx, const int ex, const int sy, const int ey, const int sz, const int ez):
  sx(sx), sy(sy), sz(sz), ex(ex), ey(ey), ez(ez) { }


  bool outside(MyRange range) const
  {
    const int x0 = max(sx, range.sx);
    const int y0 = max(sy, range.sy);
    const int z0 = max(sz, range.sz);
    const int x1 = min(ex, range.ex);
    const int y1 = min(ey, range.ey);
    const int z1 = min(ez, range.ez);

    return (x0 >= x1) || (y0 >= y1) || (z0 >= z1);
  }
};

	void fetch(const Real * const ptrBlock, Real * const ptrLab, const int x0, const int y0, const int z0,
		   const int xsize, const int ysize, const int zsize, const int gptfloats, const int rsx, const int rex, const int rsy, const int rey, const int rsz, const int rez) const
	{
	  //build range
	  MyRange myrange(rsx, rex, rsy, rey, rsz, rez);

		//packs
		{
			map<Real *, vector<PackInfo> >::const_iterator it = recv_packinfos.find(const_cast<Real *>(ptrBlock));

			if( it!=recv_packinfos.end() )
			{
				vector<PackInfo> packs = it->second;

				//assert(!stencil.tensorial || packs.size() <= 7 || mybpd[0]*mybpd[1]*mybpd[2] == 1);
				//assert(stencil.tensorial || packs.size()<=3 || mybpd[0]*mybpd[1]*mybpd[2] == 1);

				for(vector<PackInfo>::const_iterator itpack=packs.begin(); itpack!=packs.end(); ++itpack)
				{
				  MyRange packrange(itpack->sx, itpack->ex, itpack->sy, itpack->ey, itpack->sz, itpack->ez);

				  if (myrange.outside(packrange)) continue;

					const int nsrc = (itpack->ex-itpack->sx)*(itpack->ey-itpack->sy)*(itpack->ez-itpack->sz);

					unpack(itpack->pack, ptrLab, gptfloats, &stencil.selcomponents.front(), stencil.selcomponents.size(), nsrc,
						   itpack->sx-x0, itpack->sy-y0, itpack->sz-z0,
						   itpack->ex-x0, itpack->ey-y0, itpack->ez-z0,
						   xsize, ysize, zsize);
				}
			}
		}

		//subregions inside packs
		if (stencil.tensorial)
		{
			map<Real *, vector<SubpackInfo> >::const_iterator it = recv_subpackinfos.find(const_cast<Real *>(ptrBlock));

			assert(stencil.tensorial || it==recv_subpackinfos.end());

			if( it!=recv_subpackinfos.end() )
			{
				vector<SubpackInfo> subpacks = it->second;

			//	assert(subpacks.size()<=12+8);

				for(vector<SubpackInfo>::const_iterator itsubpack=subpacks.begin(); itsubpack!=subpacks.end(); ++itsubpack)
				  {
				    MyRange packrange(itsubpack->sx, itsubpack->ex, itsubpack->sy, itsubpack->ey, itsubpack->sz, itsubpack->ez);

				    if (myrange.outside(packrange)) continue;

					unpack_subregion(itsubpack->pack, ptrLab, gptfloats, &stencil.selcomponents.front(), stencil.selcomponents.size(),
									 itsubpack->x0, itsubpack->y0, itsubpack->z0,
									 itsubpack->xpacklenght, itsubpack->ypacklenght,
									 itsubpack->sx-x0, itsubpack->sy-y0, itsubpack->sz-z0,
									 itsubpack->ex-x0, itsubpack->ey-y0, itsubpack->ez-z0,
									 xsize, ysize, zsize);
				  }
			}
		}
	}
};
