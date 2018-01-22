/*
 *  BoundaryGhostCells.h
 *  MPCFnode
 *
 *  Ursula Rasthofer, May 2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */

#ifndef BOUNDARYGHOSTCELLS_H
#define BOUNDARYGHOSTCELLS_H

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <mpi.h>
using namespace std;

namespace BGC
{
  extern int BPDX, BPDY, BPDZ;

  struct BoundaryFluidElement
  {
    Real rho, u, pressure;

    BoundaryFluidElement() { rho = u = pressure = 0.0; }

    BoundaryFluidElement(const BoundaryFluidElement& a) {rho=a.rho; u=a.u; pressure=a.pressure;}

    void clear() { rho = u = pressure = 0.0; }

    void init(const Real val) { rho = u = pressure = val; }

    BoundaryFluidElement& operator = (const BoundaryFluidElement & gp)
    {
        this->rho = gp.rho;
        this->u = gp.u;
        this->pressure=gp.pressure;

        return *this;
    }
  };

  struct BoundaryBlock
  {
    BoundaryFluidElement fe[3*_BLOCKSIZE_*_BLOCKSIZE_];
    BoundaryFluidElement rhs[3*_BLOCKSIZE_*_BLOCKSIZE_];

    void print()
    {
      cout << " density: " << endl;
      for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++)
       cout << i << " " << fe[i].rho << endl;
      cout << " velocity: " << endl;
      for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++)
       cout << i << "  " << fe[i].u << endl;
      cout << " pressure: " << endl;
      for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++)
       cout << i << " " << fe[i].pressure << endl;
    }

    void printrhs()
    {
      cout << " density: " << endl;
      for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++)
       cout << i << " " << rhs[i].rho << endl;
      cout << " velocity: " << endl;
      for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++)
       cout << i << "  " << rhs[i].u << endl;
      cout << " pressure: " << endl;
      for (int i=0; i < (3*_BLOCKSIZE_*_BLOCKSIZE_); i++)
       cout << i << " " << rhs[i].pressure << endl;
    }
  };

  typedef std::vector<BoundaryBlock> BGBlocks;

  extern BGBlocks bgblock_list_dir0_side0;
  extern BGBlocks bgblock_list_dir0_side1;
  extern BGBlocks bgblock_list_dir1_side0;
  extern BGBlocks bgblock_list_dir1_side1;
  extern BGBlocks bgblock_list_dir2_side0;
  extern BGBlocks bgblock_list_dir2_side1;

  extern Real pamb; // ambient pressure
  extern Real L; // characteristic domain size
  extern Real lambda; // penalty parameter


  extern void bgc_save(int restart_id, MPI_Comm comm);
  extern void bgc_restart(int restart_id, MPI_Comm comm);

  // template <typename TBlockList>
  // void bgc_init_blocklist()
  // {
  //       BGC::BPDX = this->BPDX;
  //       BGC::BPDY = this->BPDY;
  //       BGC::BPDZ = this->BPDZ;

  //       const int xpesize = this->parser("-xpesize").asInt(1);
  //       const int ypesize = this->parser("-ypesize").asInt(1);
  //       const int zpesize = this->parser("-zpesize").asInt(1);

  //       const int NBX = this->BPDX * xpesize;
  //       const int NBY = this->BPDY * ypesize;
  //       const int NBZ = this->BPDZ * zpesize;

  //       // loop all block and check whether they are close to the boundary
  //       for(int i=0; i<(int)vInfo.size(); i++)
  //       {
  //           BlockInfo info = vInfo[i];

  //           BGC::BoundaryBlock dummy;
  //           if (info.index[0] == 0)
  //               BGC::bgblock_list_dir0_side0.push_back(dummy);
  //           // if (info.index[0] == (NBX-1))
  //           //     BGC::bgblock_list_dir0_side1.push_back(dummy);
  //           // if (info.index[1] == 0)
  //           //     BGC::bgblock_list_dir1_side0.push_back(dummy);
  //           // if (info.index[1] == (NBY-1))
  //           //     BGC::bgblock_list_dir1_side1.push_back(dummy);
  //           // if (info.index[2] == 0)
  //           //     BGC::bgblock_list_dir2_side0.push_back(dummy);
  //           // if (info.index[2] == (NBZ-1))
  //           //     BGC::bgblock_list_dir2_side1.push_back(dummy);
  //       }
  // }

  template <typename TBlock, typename TInfo, typename TElement>
  void bgc_set_block(const int ix, const int iy, const int iz, const TInfo& info, const TElement& element)
  {
      BGC::BoundaryFluidElement mybgc;
      mybgc.rho = element.alpha1rho1+element.alpha2rho2;
      // IMPORTANT REMARK: if we need the dummy for something else
      // we can recompute the pressure from energy
      mybgc.pressure = element.dummy;
      const double myrhoinv = 1.0/mybgc.rho;
      assert(!std::isnan(mybgc.rho));
      assert(!std::isnan(mybgc.pressure));
      assert(!std::isnan(myrhoinv));

      if (iz<0 && iy>=0 && iy<TBlock::sizeY && ix>=0 && ix<TBlock::sizeX)
      {
          const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY);
          const int delta_ix = 0;
          const int delta_iy = 0;
          const int delta_iz = 3;
          const int bnx = TBlock::sizeX;
          const int bny = TBlock::sizeY;
          const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

          BGC::bgblock_list_dir2_side0[bgb_index].fe[bgc_index]=mybgc;
          BGC::bgblock_list_dir2_side0[bgb_index].fe[bgc_index].u=element.rw*myrhoinv;
      }
      else if  (iz>=TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix>=0 && ix<TBlock::sizeX)
      {
          const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY);
          const int delta_ix = 0;
          const int delta_iy = 0;
          const int delta_iz = -TBlock::sizeZ;
          const int bnx = TBlock::sizeX;
          const int bny = TBlock::sizeY;
          const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

          BGC::bgblock_list_dir2_side1[bgb_index].fe[bgc_index]=mybgc;
          BGC::bgblock_list_dir2_side1[bgb_index].fe[bgc_index].u=element.rw*myrhoinv;
      }
      else if (iz>=0 && iz<TBlock::sizeZ && iy<0 && ix>=0 && ix<TBlock::sizeX)
      {
          const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ);
          const int delta_iy = 3;
          const int delta_ix = 0;
          const int delta_iz = 0;
          const int bnx = TBlock::sizeX;
          const int bny = 3;
          const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

          BGC::bgblock_list_dir1_side0[bgb_index].fe[bgc_index]=mybgc;
          BGC::bgblock_list_dir1_side0[bgb_index].fe[bgc_index].u=element.rv*myrhoinv;
      }
      else if (iz>=0 && iz<TBlock::sizeZ && iy>=TBlock::sizeY && ix>=0 && ix<TBlock::sizeX)
      {
          const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ);
          const int delta_iy = -TBlock::sizeY;
          const int delta_ix = 0;
          const int delta_iz = 0;
          const int bnx = TBlock::sizeX;
          const int bny = 3;
          const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

          BGC::bgblock_list_dir1_side1[bgb_index].fe[bgc_index]=mybgc;
          BGC::bgblock_list_dir1_side1[bgb_index].fe[bgc_index].u=element.rv*myrhoinv;
      }
      else if (iz>=0 && iz<TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix<0)
      {
          const int bgb_index = info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ);
          const int delta_ix = 3;
          const int delta_iy = 0;
          const int delta_iz = 0;
          const int bnx = 3;
          const int bny = TBlock::sizeY;
          const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

          BGC::bgblock_list_dir0_side0[bgb_index].fe[bgc_index]=mybgc;
          BGC::bgblock_list_dir0_side0[bgb_index].fe[bgc_index].u=element.ru*myrhoinv;
      }
      else if (iz>=0 && iz<TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix>=TBlock::sizeX)
      {
          const int bgb_index = info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ);
          const int delta_ix = -TBlock::sizeX;
          const int delta_iy = 0;
          const int delta_iz = 0;
          const int bny = TBlock::sizeY;
          const int bnx = 3;
          const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

          BGC::bgblock_list_dir0_side1[bgb_index].fe[bgc_index]=mybgc;
          BGC::bgblock_list_dir0_side1[bgb_index].fe[bgc_index].u=element.ru*myrhoinv;
      }
  }
}

#endif /*BOUNDARYGHOSTCELLS_H*/
