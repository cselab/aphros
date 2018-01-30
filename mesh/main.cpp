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
#include "hydro/solver.hpp"
#include "hydro/advection.hpp"


using namespace std;

const int w = 1;

void initGrid(MPIGrid& grid, const double lx, const double ly, const double lz)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();

#pragma omp parallel for
    for(int i=0; i<(int)vInfo.size(); i++)
    {
        using B = FluidBlock;
        BlockInfo info = vInfo[i];
        B& b = *(B*)info.ptrBlock;
        MIdx size(B::sizeX+2, B::sizeY+2, B::sizeZ+2);
        double pos[3];
        info.pos(pos, -1, -1, -1);
        Vect d0(pos[0]/lx, pos[1]/ly, pos[2]/lz);
        info.pos(pos, B::sizeX, B::sizeY, B::sizeZ);
        Vect d1(pos[0]/lx, pos[1]/ly, pos[2]/lz);
        Rect domain(d0, d1);
        std::cout << d0 << " " << d1 << std::endl;
        
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
    Scal dt;
    Vect vel;
    OpAdd(Scal dt, Vect vel): stencil(-w,-w,-w,w+1,w+1,w+1, true, 5, 0,1,2,3,4), dt(dt), vel(vel) {}

    template<typename LabType, typename BlockType>
    void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        using B = BlockType;
        const double h = info.h_gridpoint;

        auto& mesh = (*o.mesh);
        FieldCell<Scal> fc(mesh);

        // copy from block
        //for(int iz=0; iz<B::sizeZ; iz++)
        //    for(int iy=0; iy<B::sizeY; iy++)
        //        for(int ix=0; ix<B::sizeX; ix++) {
        for(int iz=-w; iz<B::sizeZ+w; iz++)
            for(int iy=-w; iy<B::sizeY+w; iy++)
                for(int ix=-w; ix<B::sizeX+w; ix++) {
                    MIdx midx(ix+w, iy+w, iz+w);
                    IdxCell idxcell(mesh.GetBlockCells().GetIdx(midx));
                    //fc[idxcell] = o(ix, iy, iz).p;
                    fc[idxcell] = lab(ix, iy, iz).p;
                }

        // velocity and flux
        FieldFace<Scal> ff_flux(mesh);
        for (auto idxface : mesh.Faces()) {
          ff_flux[idxface] = vel.dot(mesh.GetSurface(idxface));
        }

        // zero-derivative boundary conditions
        geom::MapFace<std::shared_ptr<solver::ConditionFace>> mf_cond;
        for (auto idxface : mesh.Faces()) {
          if (!mesh.IsExcluded(idxface) && !mesh.IsInner(idxface)) {
            mf_cond[idxface] =
                std::make_shared<solver::ConditionFaceDerivativeFixed<Scal>>(Scal(0));
          }
        }

        FieldCell<Scal> fc_src(mesh, 0.);

        //Scal dt = 0.5 * h / vel.norm();

        solver::AdvectionSolverExplicit<Mesh, FieldFace<Scal>> 
        a(mesh, fc, mf_cond, &ff_flux, &fc_src, 0., dt);
        a.StartStep();
        a.MakeIteration();
        a.FinishStep();

        fc = a.GetField();

        // copy to block
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++) {
                    MIdx midx(ix+w, iy+w, iz+w);
                    IdxCell idxcell(mesh.GetBlockCells().GetIdx(midx));
                    o(ix, iy, iz).p = fc[idxcell];
                }
    }
};

/*
struct StreamerAlpha2
{

    const Block_t& ref;

    StreamerAlpha2(const Block_t& b): ref(b) {}

    inline void operate(const int ix, const int iy, const int iz, Real output[NCHANNELS]) const
    {
        const FluidElement& el = ref.data[iz][iy][ix];
        output[0] = el.alpha2;
    }

    inline Real operate(const int ix, const int iy, const int iz) const
    {
        const FluidElement& el = ref.data[iz][iy][ix];
        return el.alpha2;
    }

    static const char * getAttributeName() { return "Scalar"; }
};
*/


template <int _ID>
struct StreamerPickOne_HDF5
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;
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

template <int i>
const std::string StreamerPickOne_HDF5<i>::NAME = "alpha";
template <int i>
const std::string StreamerPickOne_HDF5<i>::EXT = "0";


int main(int argc, char* argv[])
{
    int provided;
    MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &provided);

    ArgumentParser parser(argc, (const char**)argv);

    const int px = parser("-xpesize").asInt(2);
    const int py = parser("-ypesize").asInt(2);
    const int pz = parser("-zpesize").asInt(2);
    const int bx    = parser("-bpdx").asInt(1);
    const int by    = parser("-bpdy").asInt(1);
    const int bz    = parser("-bpdz").asInt(1);

    const double e = parser("-maxextent").asDouble(1.0);
    const int bmax = max(max(px * bx, py * by), pz * bz);
    const double lx = e * px * bx / bmax;
    const double ly = e * py * by / bmax;
    const double lz = e * pz * bz / bmax;

    MPIGrid * const g = new MPIGrid(px, py, pz, bx, by, bz, e);

    initGrid(*g, lx, ly, lz);

    Scal h = g->getBlocksInfo()[0].h_gridpoint;
    Vect vel(0.2, 0.2, 0.2);
    Real dt = 0.1;
    OpAdd a(0.25 * h / vel.norm(), vel);
    Scal t = 0.;
    for (int i = 0; i < 10; ++i) {
        Scal t0 = t;
        auto suff = "_" + std::to_string(i);
        DumpHDF5_MPI<MPIGrid, StreamerPickOne_HDF5<0>>(*g, i, i*dt, "p" + suff);
        DumpHDF5_MPI<MPIGrid, StreamerPickOne_HDF5<4>>(*g, i, i*dt, "volume" + suff);

        while (t < t0 + dt) {
          BlockProcessorMPI<AMPILab>(a, *g);  
          printf("t=%f, a.dt=%f\n", t, a.dt);
          t += a.dt;
        }
    }


    MPI_Finalize();
    return 0;
}
