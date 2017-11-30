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

#include "hydro/vect.hpp"


using namespace std;

void initGrid(MPIGrid& grid, const double lx, const double ly, const double lz)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();

#pragma omp parallel for
    for(int i=0; i<(int)vInfo.size(); i++)
    {
        using B = FluidBlock;
        BlockInfo info = vInfo[i];
        B& b = *(B*)info.ptrBlock;
        MIdx size(B::sizeX, B::sizeY, B::sizeZ);
        double pos[3];
        info.pos(pos, 0, 0, 0);
        Vect d0(pos[0], pos[1], pos[2]);
        info.pos(pos, size[0] - 1, size[1] - 1, size[2] - 1);
        Vect d1(pos[0], pos[1], pos[2]);
        Rect domain(d0, d1);
        
        b.mesh = std::unique_ptr<Mesh>(new Mesh());
        Mesh& mesh = (*b.mesh);
        geom::InitUniformMesh(mesh, domain, size);

        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double pos[3];
                    MIdx midx(ix, iy, iz);
                    IdxCell idxcell(mesh.GetBlockCells().GetIdx(midx));
                    info.pos(pos, ix, iy, iz);
                    auto& e = b(ix, iy, iz);
                    const double x = pos[0]/lx;
                    const double y = pos[1]/ly;
                    const double z = pos[2]/lz;
                    const double k = 10.;
                    e.p = std::sin(k*x) * std::sin(k*y) * std::sin(k*z);
                    e.v = mesh.GetCenter(idxcell);
                    e.volume = mesh.GetCenter(idxcell)[0];
                }
    }
}

struct OpAdd
{
    StencilInfo stencil;
    OpAdd(): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}

    template<typename LabType, typename BlockType>
    void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        using B = BlockType;
        const double h = info.h_gridpoint;

        auto& mesh = (*o.mesh);
        FieldCell<Scal> fc(mesh);
        

        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++) {
                    o(ix, iy, iz).p += 0.1;
                }

        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++) {
                    o(ix, iy, iz).p += 0.1;
                }
    }
};

template <int _ID>
struct StreamerPickOne_HDF5
{
    static const int NCHANNELS = 1; // we dump scalar fields with this streamer
    struct AssumedType { Real team[5]; };

    FluidBlock& ref;

    StreamerPickOne_HDF5(FluidBlock& b): ref(b) {}

    // write
    void operate(const int ix, const int iy, const int iz, Real output[0]) const
    {
        const AssumedType& input = *((AssumedType*)&ref.data[iz][iy][ix].p);
        output[0] = input.team[_ID];
    }

    // read
    void operate(const Real output[0], const int ix, const int iy, const int iz) const
    {
        AssumedType& input = *((AssumedType*)&ref.data[iz][iy][ix].p);
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
    const double lx = e * px * bx / (double)bmax;
    const double ly = e * py * by / (double)bmax;
    const double lz = e * pz * bz / (double)bmax;

    MPIGrid * const g = new MPIGrid(px, py, pz, bx, by, bz, e);

    initGrid(*g, lx, ly, lz);

    Real dt = 0.1;
    OpAdd a;
    for (int i = 0; i < 10; ++i) {
        auto suff = "_" + std::to_string(i);
        DumpHDF5_MPI<MPIGrid, StreamerPickOne_HDF5<0>>(*g, i, i*dt, "p" + suff);
        DumpHDF5_MPI<MPIGrid, StreamerPickOne_HDF5<4>>(*g, i, i*dt, "volume" + suff);

        BlockProcessorMPI<AMPILab>(a, *g);  
    }


    MPI_Finalize();
    return 0;
}
