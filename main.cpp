// File       : test-cubismo.cpp
// Date       : Fri 01 Apr 2016 05:50:11 PM CEST
// Author     : Fabian Wermelinger
// Description: Short Cubism MPI test
// Copyright 2016 ETH Zurich. All Rights Reserved.

#include <mpi.h>
#include <cmath>

#include "Types.h"
#include "BlockProcessorMPI.h"
#include "HDF5Dumper_MPI.h"

using namespace std;

const double twopi = atan(1.0)*8.0;

void initGrid(MPIGrid& grid, const double Lx, const double Ly, const double Lz)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();

#pragma omp parallel for
    for(int i=0; i<(int)vInfo.size(); i++)
    {
        BlockInfo info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;

        for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
                for(int ix=0; ix<FluidBlock::sizeX; ix++)
                {
                    double pos[3];
                    info.pos(pos, ix, iy, iz);
                    const double x = pos[0]/Lx;
                    const double y = pos[1]/Ly;
                    const double z = pos[2]/Lz;
                    b(ix, iy, iz).u = x;
                    b(ix, iy, iz).volume = 1.;
                }
    }
}

struct Evaluate123Divergence_CPP
{
    // second order FD divergence of data[0], data[1], data[2] (of
    // FluidElement)

    StencilInfo stencil;

    // stencil: first 3 are stencil to the left from zero (inclusive)
    //          next 3 are stencil to the right from zero (exclusive)
    //          next 1 says tensorial true/false
    //          next 1 specifies the n fields you need to evaluate the operator
    //          next n are indices of the fields to be communicated by MPI
    Evaluate123Divergence_CPP(): stencil(-1,-1,-1,2,2,2, false, 3, 0,1,2) {}

    // operator interface
    template<typename LabType, typename BlockType>
    void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        typedef BlockType B;
        const double h = info.h_gridpoint;

        for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const Real dd0_dx = lab(ix+1,iy,iz).data[0] - lab(ix-1,iy,iz).data[0];
                    const Real dd1_dy = lab(ix,iy+1,iz).data[1] - lab(ix,iy-1,iz).data[1];
                    const Real dd2_dz = lab(ix,iy,iz+1).data[2] - lab(ix,iy,iz-1).data[2];

                    // we stuff it into data[3]
                    o(ix, iy, iz).data[3] = 0.5*(dd0_dx + dd1_dy + dd2_dz)/h;
                }
    }
};

struct OpAdd
{
    StencilInfo stencil;
    OpAdd(): stencil(-1,-1,-1,2,2,2, false, 2, 0,1) {}

    template<typename LabType, typename BlockType>
    void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        typedef BlockType B;
        const double h = info.h_gridpoint;

        for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    o(ix, iy, iz).u += 1.;
                }
    }
};

template <int _ID>
struct StreamerPickOne_HDF5
{
    static const int NCHANNELS = 1; // we dump scalar fields with this streamer
    struct AssumedType { Real team[8]; };

    FluidBlock& ref;

    StreamerPickOne_HDF5(FluidBlock& b): ref(b) {}

    // write
    void operate(const int ix, const int iy, const int iz, Real output[0]) const
    {
        const AssumedType& input = *((AssumedType*)&ref.data[iz][iy][ix].u);
        output[0] = input.team[_ID];
    }

    // read
    void operate(const Real output[0], const int ix, const int iy, const int iz) const
    {
        AssumedType& input = *((AssumedType*)&ref.data[iz][iy][ix].u);
        input.team[_ID] = output[0];
    }

    static const char * getAttributeName() { return "Scalar"; }
};


int main(int argc, char* argv[])
{
    int provided;
    MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &provided);

    ArgumentParser parser(argc, (const char**)argv);

    const int px = parser("-xpesize").asInt(1);
    const int py = parser("-ypesize").asInt(1);
    const int pz = parser("-zpesize").asInt(2);
    const int bx    = parser("-bpdx").asInt(2);
    const int by    = parser("-bpdy").asInt(2);
    const int bz    = parser("-bpdz").asInt(1);

    const double e = parser("-maxextent").asDouble(1.0);
    const int bmax = max(max(bx, by), bz);
    const double lx = e * bx / (double)bmax;
    const double ly = e * by / (double)bmax;
    const double lz = e * bz / (double)bmax;

    MPIGrid * const g = new MPIGrid(px, py, pz, bx, by, bz, e);

    initGrid(*g, lx, ly, lz);

    Real dt = 0.1;
    OpAdd a;
    for (int i = 0; i < 10; ++i) {
        auto suff = "_" + std::to_string(i);
        DumpHDF5_MPI<MPIGrid, StreamerPickOne_HDF5<0>>(*g, i, i*dt, "p" + suff);
        DumpHDF5_MPI<MPIGrid, StreamerPickOne_HDF5<1>>(*g, i, i*dt, "volume" + suff);

        BlockProcessorMPI<AMPILab>(a, *g);  
    }


    MPI_Finalize();
    return 0;
}
