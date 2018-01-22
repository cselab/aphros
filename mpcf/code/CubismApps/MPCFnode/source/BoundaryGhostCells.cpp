/*
 *  BoundaryGhostCells.cpp
 *  MPCFnode
 *
 *  Ursula Rasthofer, May 2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include "BoundaryGhostCells.h"
#include "PlainBinDumper_MPI.h"

namespace BGC
{
  int BPDX, BPDY, BPDZ;
  BGBlocks bgblock_list_dir0_side0;
  BGBlocks bgblock_list_dir0_side1;
  BGBlocks bgblock_list_dir1_side0;
  BGBlocks bgblock_list_dir1_side1;
  BGBlocks bgblock_list_dir2_side0;
  BGBlocks bgblock_list_dir2_side1;
  Real pamb, L, lambda;


  void bgc_save(int restart_id, MPI_Comm comm)
  {
#if DBG
    cout << "INSIDE bgc_save!!!!!\n";
#endif

//  set-up boundary ghost cells
    int bgc_nblocks =
	BGC::bgblock_list_dir0_side0.size() +
	BGC::bgblock_list_dir0_side1.size() +
	BGC::bgblock_list_dir1_side0.size() +
	BGC::bgblock_list_dir1_side1.size() +
	BGC::bgblock_list_dir2_side0.size() +
	BGC::bgblock_list_dir2_side1.size();

    long data_size = bgc_nblocks * sizeof(BGC::BoundaryFluidElement)* (3*_BLOCKSIZE_*_BLOCKSIZE_);  // bgc_nblocks * sizeof(BGC::BoundaryBlock) / 2;

    long header_size =
	3 * sizeof(Real) + 	// BPDX, BPDY, BPDZ (saved as Real)
	6 * sizeof(Real) +	// # elements in each of the six lists (saved as Real)
	3 * sizeof(Real);	// pamb, L, lambda

    long bgc_size = header_size + data_size;
    Real *bgc_buffer = (Real *)malloc(bgc_size);

#define doserialize(a) \
{ \
  Real v; \
  v = (Real)a; bgc_buffer[bufidx] = v; bufidx++; \
}

#define doserializev() \
{ \
  Real v; \
  for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++) { v = it->fe[i].rho; bgc_buffer[bufidx] = v; bufidx++; } \
  for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++) { v = it->fe[i].u; bgc_buffer[bufidx] = v; bufidx++; } \
  for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++) { v = it->fe[i].pressure; bgc_buffer[bufidx] = v; bufidx++; } \
}

#if DBG
    printf("BGC::BPDX = %d\n", BGC::BPDX);
    printf("BGC::BPDY = %d\n", BGC::BPDY);
    printf("BGC::BPDZ = %d\n", BGC::BPDZ);

    printf("BGC::s00 = %d\n", BGC::bgblock_list_dir0_side0.size());
    printf("BGC::s01 = %d\n", BGC::bgblock_list_dir0_side1.size());
    printf("BGC::s10 = %d\n", BGC::bgblock_list_dir1_side0.size());
    printf("BGC::s11 = %d\n", BGC::bgblock_list_dir1_side1.size());
    printf("BGC::s20 = %d\n", BGC::bgblock_list_dir2_side0.size());
    printf("BGC::s21 = %d\n", BGC::bgblock_list_dir2_side1.size());

    printf("BGC::pamb   = %f\n", BGC::pamb);
    printf("BGC::L      = %f\n", BGC::L);
    printf("BGC::lambda = %f\n", BGC::lambda);
#endif

    int bufidx = 0;
    doserialize(BGC::BPDX);
    doserialize(BGC::BPDY);
    doserialize(BGC::BPDZ);

    doserialize(BGC::bgblock_list_dir0_side0.size());
    doserialize(BGC::bgblock_list_dir0_side1.size());
    doserialize(BGC::bgblock_list_dir1_side0.size());
    doserialize(BGC::bgblock_list_dir1_side1.size());
    doserialize(BGC::bgblock_list_dir2_side0.size());
    doserialize(BGC::bgblock_list_dir2_side1.size());

    doserialize(BGC::pamb);
    doserialize(BGC::L);
    doserialize(BGC::lambda);

    vector<BGC::BoundaryBlock>::iterator it;
    for(it=BGC::bgblock_list_dir0_side0.begin(); it != BGC::bgblock_list_dir0_side0.end(); it++) doserializev();
    for(it=BGC::bgblock_list_dir0_side1.begin(); it != BGC::bgblock_list_dir0_side1.end(); it++) doserializev();
    for(it=BGC::bgblock_list_dir1_side0.begin(); it != BGC::bgblock_list_dir1_side0.end(); it++) doserializev();
    for(it=BGC::bgblock_list_dir1_side1.begin(); it != BGC::bgblock_list_dir1_side1.end(); it++) doserializev();
    for(it=BGC::bgblock_list_dir2_side0.begin(); it != BGC::bgblock_list_dir2_side0.end(); it++) doserializev();
    for(it=BGC::bgblock_list_dir2_side1.begin(); it != BGC::bgblock_list_dir2_side1.end(); it++) doserializev();

#if DBG
    cout << "bufidx = " << bufidx << "\n";
#endif

    //store in the bgc_restart_file the buffer size and the buffer
#if 0
    char datafile[256];
    sprintf(datafile,"restart_ic_%d.bin", restart_id);
    FILE *fp = fopen(datafile, "wb");
    fwrite(&bgc_size, 1, sizeof(bgc_size), fp);
    fwrite(bgc_buffer, bgc_size, sizeof(char), fp);
    fclose(fp);
#else
    const string path = "."; //parser("-fpath").asString(".");
    stringstream datafile;
    datafile << "restart_ic_" << restart_id;
    PlainDumpBin_MPI(comm, bgc_buffer, bgc_size, datafile.str().c_str(), path.c_str());
    if (restart_id>0)
    {
      datafile << "restart_ic_" << restart_id-1;
      remove(datafile.str().c_str());
    }
#endif

    free(bgc_buffer);
// #endif
 }

  void bgc_restart(int restart_id, MPI_Comm comm)
  {
#if DBG
    cout << "INSIDE bgc_restart!!!!!\n";
#endif

//  set-up boundary ghost cells
    long bgc_size;
    Real *bgc_buffer = NULL;

//  open bgc_restart_file, read the buffer size, allocate the buffer, fill the buffer, copy contents to the right place

#if 0
    char datafile[256];
    sprintf(datafile,"restart_ic_%d.bin", restart_id);
    printf("reading restart file %s\n", datafile);
    FILE *fp = fopen(datafile, "rb");
    fread(&bgc_size, 1, sizeof(bgc_size), fp);
    bgc_buffer = (Real *)malloc(bgc_size);
    fread(bgc_buffer, bgc_size, sizeof(char), fp);
    fclose(fp);
#else
    const string path = "."; //parser("-fpath").asString(".");
    stringstream datafile;
    datafile << "restart_ic_" << restart_id;
    PlainReadBin_MPI(comm, &bgc_buffer, &bgc_size, datafile.str().c_str(), path.c_str());
#endif

#define dodeserialize(a) \
{ \
  Real v; \
  v = bgc_buffer[bufidx]; a = v; bufidx++; \
}

#define dodeserialize2int(a) \
{ \
  Real v; \
  v = bgc_buffer[bufidx]; a = (int)v; bufidx++; \
}

#define dodeserializev() \
{ \
  Real v; \
  for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++) { v = bgc_buffer[bufidx]; it->fe[i].rho = v; bufidx++; } \
  for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++) { v = bgc_buffer[bufidx]; it->fe[i].u = v; bufidx++; } \
  for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++) { v = bgc_buffer[bufidx]; it->fe[i].pressure = v; bufidx++; } \
}

    int bufidx = 0;
    dodeserialize2int(BGC::BPDX);
    dodeserialize2int(BGC::BPDY);
    dodeserialize2int(BGC::BPDZ);

    int bgblock_list_dir0_side0_size;
    int bgblock_list_dir0_side1_size;
    int bgblock_list_dir1_side0_size;
    int bgblock_list_dir1_side1_size;
    int bgblock_list_dir2_side0_size;
    int bgblock_list_dir2_side1_size;

    dodeserialize2int(bgblock_list_dir0_side0_size);
    dodeserialize2int(bgblock_list_dir0_side1_size);
    dodeserialize2int(bgblock_list_dir1_side0_size);
    dodeserialize2int(bgblock_list_dir1_side1_size);
    dodeserialize2int(bgblock_list_dir2_side0_size);
    dodeserialize2int(bgblock_list_dir2_side1_size);

    dodeserialize(BGC::pamb);
    dodeserialize(BGC::L);
    dodeserialize(BGC::lambda);

    BGC::BoundaryBlock dummy;
    for (int i=0; i < bgblock_list_dir0_side0_size; i++) BGC::bgblock_list_dir0_side0.push_back(dummy);
    for (int i=0; i < bgblock_list_dir0_side1_size; i++) BGC::bgblock_list_dir0_side1.push_back(dummy);
    for (int i=0; i < bgblock_list_dir1_side0_size; i++) BGC::bgblock_list_dir1_side0.push_back(dummy);
    for (int i=0; i < bgblock_list_dir1_side1_size; i++) BGC::bgblock_list_dir1_side1.push_back(dummy);
    for (int i=0; i < bgblock_list_dir2_side0_size; i++) BGC::bgblock_list_dir2_side0.push_back(dummy);
    for (int i=0; i < bgblock_list_dir2_side1_size; i++) BGC::bgblock_list_dir2_side1.push_back(dummy);

    vector<BGC::BoundaryBlock>::iterator it;
    for(it=BGC::bgblock_list_dir0_side0.begin(); it != BGC::bgblock_list_dir0_side0.end(); it++) dodeserializev();
    for(it=BGC::bgblock_list_dir0_side1.begin(); it != BGC::bgblock_list_dir0_side1.end(); it++) dodeserializev();
    for(it=BGC::bgblock_list_dir1_side0.begin(); it != BGC::bgblock_list_dir1_side0.end(); it++) dodeserializev();
    for(it=BGC::bgblock_list_dir1_side1.begin(); it != BGC::bgblock_list_dir1_side1.end(); it++) dodeserializev();
    for(it=BGC::bgblock_list_dir2_side0.begin(); it != BGC::bgblock_list_dir2_side0.end(); it++) dodeserializev();
    for(it=BGC::bgblock_list_dir2_side1.begin(); it != BGC::bgblock_list_dir2_side1.end(); it++) dodeserializev();

    free(bgc_buffer);

#if DBG
    cout << "bufidx = " << bufidx << "\n";

    printf("BGC::BPDX = %d\n", BGC::BPDX);
    printf("BGC::BPDY = %d\n", BGC::BPDY);
    printf("BGC::BPDZ = %d\n", BGC::BPDZ);

    printf("BGC::s00 = %d\n", BGC::bgblock_list_dir0_side0.size());
    printf("BGC::s01 = %d\n", BGC::bgblock_list_dir0_side1.size());
    printf("BGC::s10 = %d\n", BGC::bgblock_list_dir1_side0.size());
    printf("BGC::s11 = %d\n", BGC::bgblock_list_dir1_side1.size());
    printf("BGC::s20 = %d\n", BGC::bgblock_list_dir2_side0.size());
    printf("BGC::s21 = %d\n", BGC::bgblock_list_dir2_side1.size());

    printf("BGC::pamb   = %f\n", BGC::pamb);
    printf("BGC::L      = %f\n", BGC::L);
    printf("BGC::lambda = %f\n", BGC::lambda);
#endif
 }

}

