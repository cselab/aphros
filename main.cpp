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

#ifndef _2DSTENCIL_
        for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
                for(int ix=0; ix<FluidBlock::sizeX; ix++)
                {
                    double pos[3];
                    info.pos(pos, ix, iy, iz);
                    const double x = pos[0]/Lx;
                    const double y = pos[1]/Ly;
                    const double z = pos[2]/Lz;
                    const double q1 = x;
                    const double q2 = y;
                    const double q3 = z;
                    b(ix, iy, iz).data[0] = q1;
                    b(ix, iy, iz).data[1] = q2;
                    b(ix, iy, iz).data[2] = q3;
                    b(ix, iy, iz).data[3] = 5.0*q1;
                    b(ix, iy, iz).data[4] = 5.0*q2;
                    b(ix, iy, iz).data[5] = 5.0*q3;
                    b(ix, iy, iz).data[6] = 1.0;
                    b(ix, iy, iz).data[7] = 0.0;
                }
#else
        for(int iy=0; iy<FluidBlock::sizeY; iy++)
            for(int ix=0; ix<FluidBlock::sizeX; ix++)
            {
                double pos[2];
                info.pos(pos, ix, iy);
                const double x = pos[0]/Lx;
                const double y = pos[1]/Ly;
                const double q1 = sin(twopi*x);
                const double q2 = cos(twopi*y);
                b(ix, iy).data[0] = q1;
                b(ix, iy).data[1] = q2;
                b(ix, iy).data[2] = q1;
                b(ix, iy).data[3] = 0.5*q1;
                b(ix, iy).data[4] = 0.5*q2;
                b(ix, iy).data[5] = 0.5*q1;
                b(ix, iy).data[6] = 1.0;
                b(ix, iy).data[7] = 0.0;
            }
#endif /* _2DSTENCIL_ */
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
#ifndef _2DSTENCIL_
    Evaluate123Divergence_CPP(): stencil(-1,-1,-1,2,2,2, false, 3, 0,1,2) {}
#else
    Evaluate123Divergence_CPP(): stencil(-1,-1,0,2,2,1, false, 2, 0,1) {}
#endif /* _2DSTENCIL_ */

    // operator interface
    template<typename LabType, typename BlockType>
    void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        typedef BlockType B;
        const double h = info.h_gridpoint;

#ifndef _2DSTENCIL_
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
#else
        for(int iy=0; iy<BlockType::sizeY; iy++)
            for(int ix=0; ix<BlockType::sizeX; ix++)
            {
                const Real dd0_dx = lab(ix+1,iy).data[0] - lab(ix-1,iy).data[0];
                const Real dd1_dy = lab(ix,iy+1).data[1] - lab(ix,iy-1).data[1];

                // we stuff it into data[3]
                o(ix, iy).data[3] = 0.5*(dd0_dx + dd1_dy)/h;
            }
#endif /* _2DSTENCIL_ */
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
        const AssumedType& input = *((AssumedType*)&ref.data[iz][iy][ix].data[0]);
        output[0] = input.team[_ID];
    }

    // read
    void operate(const Real output[0], const int ix, const int iy, const int iz) const
    {
        AssumedType& input = *((AssumedType*)&ref.data[iz][iy][ix].data[0]);
        input.team[_ID] = output[0];
    }

    static const char * getAttributeName() { return "Scalar"; }
};


int main(int argc, char* argv[])
{
    int provided;
    MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &provided);

    ArgumentParser parser(argc, (const char**)argv);

    const int xpesize = parser("-xpesize").asInt(1);
    const int ypesize = parser("-ypesize").asInt(xpesize);
    const int zpesize = parser("-zpesize").asInt(xpesize);
    const int bpdx    = parser("-bpdx").asInt(1);
    const int bpdy    = parser("-bpdy").asInt(bpdx);
    const int bpdz    = parser("-bpdz").asInt(bpdx);

    const double maxextent = parser("-maxextent").asDouble(1.0);
    const int bpd_max = max(max(bpdx, bpdy), bpdz);
    const double Lx = maxextent * bpdx / (double)bpd_max;
    const double Ly = maxextent * bpdy / (double)bpd_max;
    const double Lz = maxextent * bpdz / (double)bpd_max;

    // allocate the grid
    MPIGrid * const mygrid = new MPIGrid(xpesize, ypesize, zpesize, bpdx, bpdy, bpdz, maxextent);

    // initialize values
    initGrid(*mygrid, Lx, Ly, Lz);
    DumpHDF5_MPI<MPIGrid, StreamerPickOne_HDF5<0> >(*mygrid, 0, 0, "data0");
    DumpHDF5_MPI<MPIGrid, StreamerPickOne_HDF5<1> >(*mygrid, 0, 0, "data1");
#ifndef _2DSTENCIL_
    DumpHDF5_MPI<MPIGrid, StreamerPickOne_HDF5<2> >(*mygrid, 0, 0, "data2");
#endif

    // evaluate div([data[0], data[1], data[2]]')
    Evaluate123Divergence_CPP myOperator;
#ifndef _2DSTENCIL_
    BlockProcessorMPI<AMPILab>(myOperator, *mygrid);  // absorbing BC
#else
    BlockProcessorMPI<PMPILab>(myOperator, *mygrid);  // periodic BC
#endif /* _2DSTENCIL_ */

    // dump result
    DumpHDF5_MPI<MPIGrid, StreamerPickOne_HDF5<3> >(*mygrid, 0, 0, "div123");

    MPI_Finalize();
    return 0;
}
