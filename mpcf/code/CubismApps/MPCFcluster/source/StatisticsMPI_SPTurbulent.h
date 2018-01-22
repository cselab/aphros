/*
 *  StatisticsMPI_SPTurbulent.h
 *  MPCFcluster
 *  Single Phase Turbulent Statistics
 *
 *  Created by Fabian Wermelinegr 04/14/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef STATISTICSMPI_SPTURBULENT_H_VAIDNWTG
#define STATISTICSMPI_SPTURBULENT_H_VAIDNWTG

#include <cassert>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <mpi.h>

#include "Types.h"
#include "BlockInfo.h"
#include "ArgumentParser.h"

using namespace std;

#define USE_OPENMP 1

#ifdef _FLOAT_PRECISION_
#define POS_MPI_REAL MPI_FLOAT
#else
#define POS_MPI_REAL MPI_DOUBLE
#endif /* _FLOAT_PRECISION_ */

namespace SinglePhaseTurbulenceStatistics
{
    struct AverageQoI
    {
        size_t updateCount;
        std::vector<double> quantity;
        AverageQoI() : quantity(24) { clear(); }

        void clear()
        {
            updateCount = 0;
            for (size_t i = 0; i < quantity.size(); ++i)
                quantity[i] = 0.0;
        }
        void update(const double rho,
                const double rhoU, const double rhoV, const double rhoW,
                const double U, const double V, const double W,
                const double p)
        {
            // mean velocity products
            quantity[0] += U;   quantity[1] += V;   quantity[2] += W;
            quantity[3] += U*U; quantity[4] += U*V; quantity[5] += U*W;
            quantity[6] += V*V; quantity[7] += V*W; quantity[8] += W*W;

            // mean pressure products
            quantity[9] += p;
            quantity[10] += p*U; quantity[11] += p*V; quantity[12] += p*W;
            quantity[13] += rho*p;

            // mean density products
            quantity[14] += rho;
            quantity[15] += rhoU;   quantity[16] += rhoV;   quantity[17] += rhoW;
            quantity[18] += rhoU*U; quantity[19] += rhoU*V; quantity[20] += rhoU*W;
            quantity[21] += rhoV*V; quantity[22] += rhoV*W; quantity[23] += rhoW*W;
        }
        void update(const AverageQoI& rhs)
        {
            for (size_t i = 0; i < quantity.size(); ++i)
                quantity[i] += rhs.quantity[i];
        }
    };


    template <typename TGrid>
    void updateStatistics(TGrid& grid, std::vector<AverageQoI>& hist)
    {
        typedef typename TGrid::BlockType TBlock;

        MPI_Comm world_comm = grid.getCartComm();
        int myrank;
        MPI_Comm_rank(world_comm, &myrank);
        const bool isroot = (0 == myrank);

        // Node wise update, no MPI communication here
        const int BPDY = grid.getResidentBlocksPerDimension(1);
        assert(hist.size() == BPDY*TBlock::sizeY);
        double t0, t1;
        std::vector<BlockInfo> vInfo = grid.getBlocksInfo();

        // single phase, material properties are assumed in GAMMA1 and PC1
        const double G = 1.0/(Simulation_Environment::GAMMA1-1.0);
        const double F = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;

        t0 = omp_get_wtime();

#ifdef USE_OPENMP
#pragma omp parallel
#endif
        {
            std::vector<AverageQoI> myshare(hist.size());

#pragma omp for schedule(dynamic,1)
            for(int i=0; i<(int)vInfo.size(); i++)
            {
                BlockInfo info = vInfo[i];
                TBlock& b = *(TBlock*)info.ptrBlock;

                const int by = (info.index[1] % BPDY)*TBlock::sizeY;

                // process xz-planes
                for(int iy=0; iy<TBlock::sizeY; iy++)
                {
                    AverageQoI& layer = myshare[by+iy];
                    for(int iz=0; iz<TBlock::sizeZ; iz++)
                        for(int ix=0; ix<TBlock::sizeX; ix++)
                        {
                            const double rho  = b(ix,iy,iz).alpha1rho1 + b(ix,iy,iz).alpha2rho2;
                            const double rhoInv = 1.0/rho;
                            const double rhoU = b(ix,iy,iz).ru;
                            const double rhoV = b(ix,iy,iz).rv;
                            const double rhoW = b(ix,iy,iz).rw;
                            const double U    = rhoU*rhoInv;
                            const double V    = rhoV*rhoInv;
                            const double W    = rhoW*rhoInv;

                            const double ke = 0.5*(rhoU*rhoU + rhoV*rhoV + rhoW*rhoW)*rhoInv;
                            const double p  = (b(ix,iy,iz).energy - F*G - ke)/G;

                            layer.update(rho, rhoU, rhoV, rhoW, U, V, W, p);
                        }
                }
            }

            // update history (OMP reduction)
#pragma omp critical
            {
                for (size_t i = 0; i < hist.size(); ++i)
                    hist[i].update(myshare[i]);
            }
        }

        ++(hist.front().updateCount); // we only update the count of the first element (reduce redundancy)

        t1 = omp_get_wtime();
        if (isroot) printf("statistics turbulent update: Took %lf seconds\n", t1-t0);
    }


    template <typename TGrid>
    void dumpStatistics(TGrid& grid, std::vector<AverageQoI>& hist, const int step_id, const Real t, const Real dt, ArgumentParser& parser)
    {
        typedef typename TGrid::BlockType TBlock;

        const std::string path = parser("-fpath").asString(".");
        std::vector<BlockInfo> vInfo = grid.getBlocksInfo();

        const double h  = vInfo[0].h_gridpoint;

        int peidx[3];
        grid.peindex(peidx);

        double t0, t1;

        t0 = omp_get_wtime();

        MPI_Comm world_comm = grid.getCartComm();
        int worldRank, worldSize;
        MPI_Comm_size(world_comm, &worldSize);
        MPI_Comm_rank(world_comm, &worldRank);
        const bool isroot = (worldRank == 0);

        MPI_Comm comm_subgroup;
        MPI_Comm_split(world_comm, peidx[1], worldRank, &comm_subgroup); // color is y-index

        int subrank;
        MPI_Comm_rank(comm_subgroup, &subrank);

        // reduction on comm_subgroups
        for (size_t i = 0; i < hist.size(); ++i)
        {
            if (0==subrank)
                MPI_Reduce(MPI_IN_PLACE, hist[i].quantity.data(), hist[i].quantity.size(), MPI_DOUBLE, MPI_SUM, 0, comm_subgroup);
            else
                MPI_Reduce(hist[i].quantity.data(), hist[i].quantity.data(), hist[i].quantity.size(), MPI_DOUBLE, MPI_SUM, 0, comm_subgroup);
        }

        // gather to world root
        const int packSize = hist[0].quantity.size();
        int localCount = packSize*hist.size();
        std::vector<double> nAverages(localCount);
        int nSerial = 0;
        if (0==subrank) nSerial = nAverages.size();

        // copy data into serial array (for gatherv)
        for (size_t i = 0; i < hist.size(); ++i)
            memcpy(nAverages.data() + i*packSize, hist[i].quantity.data(), packSize*sizeof(double));

        // compute offsets and counts (for gatherv)
        std::vector<int> howMany(worldSize, 0);
        std::vector<int> offsets(worldSize, 0);
        MPI_Gather(&nSerial, 1, MPI_INT, howMany.data(), 1, MPI_INT, 0, world_comm);
        for (size_t i = 1; i < offsets.size(); ++i)
            offsets[i] = offsets[i-1] + howMany[i-1];

        // gather all data to root
        const size_t yCells = grid.getBlocksPerDimension(1)*TBlock::sizeY;
        std::vector<double> allAverages(packSize*yCells);
        MPI_Gatherv(nAverages.data(), nSerial, MPI_DOUBLE, allAverages.data(), howMany.data(), offsets.data(), MPI_DOUBLE, 0, world_comm);

        MPI_Comm_free(&comm_subgroup);

        t1 = omp_get_wtime();
        if (isroot) printf("statistics turbulent: A took %lf seconds\n", t1-t0);

        t0 = omp_get_wtime();

        if (isroot)
        {
            const double xCells = grid.getBlocksPerDimension(0)*TBlock::sizeX;
            const double zCells = grid.getBlocksPerDimension(2)*TBlock::sizeZ;
            const size_t updateCount = hist.front().updateCount; // must use hist.front(), see updateStatistics function
            const double factor = 1.0/(updateCount*xCells*zCells);
            std::ostringstream outname;
            outname << "statistics_channel_y_" << setw(6) << setfill('0') << step_id << ".dat";
            std::ofstream out(outname.str().c_str(), std::ios::out);
            out << "# Turbulent Channel Flow Y-statistics" << std::endl;
            out << "# Step ID      = " << step_id << std::endl;
            out.setf(std::ios::scientific, std::ios::floatfield);
            out.precision(12);
            out << "# Time         = " << t << std::endl;
            out << "# Iter Samples = " << updateCount << std::endl;
            out << "# Grid Spacing = " << h << std::endl;
            out << "#" << std::endl;
            out << "# Mean quantities:" << std::endl;
            out << "# y-pos\tU\tV\tW\tU*U\tU*V\tU*W\tV*V\tV*W\tW*W\tp\tp*U\tp*V\tp*W\trho*p\trho\trho*U\trho*V\trho*W\trho*U*U\trho*U*V\trho*U*W\trho*V*V\trho*V*W\trho*W*W" << std::endl;
            for (size_t cell = 0; cell < yCells; ++cell)
            {
                out << h*(cell+0.5);
                for (size_t item = 0; item < packSize; ++item)
                    out << '\t' << factor*allAverages[item + cell*packSize];
                out << std::endl;
            }
            out.close();
        }

        // clear
        for (size_t i = 0; i < hist.size(); ++i)
            hist[i].clear();

        t1 = omp_get_wtime();
        if (isroot) printf("statistics turbulent: B took %lf seconds\n", t1-t0);
    }
}

#undef USE_OPENMP

#endif /* STATISTICSMPI_SPTURBULENT_H_VAIDNWTG */
