/*
 *  StatisticsMPI.h
 *  MPCFcluster
 *
 *  Created by Fabian Wermelinegr 04/14/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef STATISTICSMPI_H_MXLD8JPJ
#define STATISTICSMPI_H_MXLD8JPJ

#include <mpi.h>
#include <omp.h>
#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "ArgumentParser.h"
#include "Types.h"
#include "FlowStep_LSRK3MPI.h"
#include "BlockInfo.h"
#include "VectorCalculusMPI.h"

namespace SinglePhaseStatistics
{
    template <typename TGrid>
    void dumpStatistics(TGrid& grid, const int step_id, const Real t, const Real dt, ArgumentParser& parser);
}

namespace MultiPhaseStatistics
{
    template <typename TGrid>
    void dumpStatistics(TGrid& grid, const int step_id, const Real t, const Real dt, ArgumentParser& parser);
}

// namespace SinglePhaseTurbulenceStatistics
// {
//     struct AverageQoI
//     {
//         size_t updateCount;
//         std::vector<double> quantity;
//         AverageQoI() : quantity(24) { clear(); }

//         void clear()
//         {
//             updateCount = 0;
//             for (size_t i = 0; i < quantity.size(); ++i)
//                 quantity[i] = 0.0;
//         }
//         void update(const double rho,
//                 const double rhoU, const double rhoV, const double rhoW,
//                 const double U, const double V, const double W,
//                 const double p)
//         {
//             // mean velocity products
//             quantity[0] += U;   quantity[1] += V;   quantity[2] += W;
//             quantity[3] += U*U; quantity[4] += U*V; quantity[5] += U*W;
//             quantity[6] += V*V; quantity[7] += V*W; quantity[8] += W*W;

//             // mean pressure products
//             quantity[9] += p;
//             quantity[10] += p*U; quantity[11] += p*V; quantity[12] += p*W;
//             quantity[13] += rho*p;

//             // mean density products
//             quantity[14] += rho;
//             quantity[15] += rhoU;   quantity[16] += rhoV;   quantity[17] += rhoW;
//             quantity[18] += rhoU*U; quantity[19] += rhoU*V; quantity[20] += rhoU*W;
//             quantity[21] += rhoV*V; quantity[22] += rhoV*W; quantity[23] += rhoW*W;
//         }
//         void update(const AverageQoI& rhs)
//         {
//             for (size_t i = 0; i < quantity.size(); ++i)
//                 quantity[i] += rhs.quantity[i];
//         }
//     };

//     template <typename TGrid>
//     void updateStatistics(TGrid& grid, std::vector<AverageQoI>& hist);

//     template <typename TGrid>
//     void dumpStatistics(TGrid& grid, std::vector<AverageQoI>& hist, const int step_id, const Real t, const Real dt, ArgumentParser& parser);
// }

// Basic types
struct double_int {
    double value;
    int    rank;
};

struct GlobalData
{
    double h;
    double V_block;
    double center[3];
    int nBlocks;
    unsigned long long int wallCells;
    ArgumentParser& parser;
    GlobalData(ArgumentParser& p) :
        h(0.), V_block(0.), nBlocks(0), wallCells(0), parser(p)
    { center[0] = center[1] = center[2] = 0.; }
};

struct BlockData
{
    double r, r1, r2, u, v, w, e, ke, p, a1, a2, W, W0, W1, W2, IuI, c, M, K, divu, gradp;
    double energyDissipationRate, kenstrophy, dummy;
    double V;
    double V2;
    double V2_01;
    double V2_05;
    double V2_10;
    double V2_50;
    double V2_90;
    double V2_95;
    double V2_99;

    unsigned long long int wallCells;

    BlockData() :
        r(0.), r1(0.), r2(0.), u(0.), v(0.), w(0.), e(0.), ke(0.), p(0.),
        a1(0.), a2(0.), W(0.), W0(0.), W1(0.), W2(0.), IuI(0.), c(0.), M(0.),
        K(0.), divu(0.), gradp(0.), energyDissipationRate(0.), kenstrophy(0.),
        dummy(0.), V(0.), V2(0.), V2_01(0.), V2_05(0.), V2_10(0.), V2_50(0.),
        V2_90(0.), V2_95(0.), V2_99(0.),
        wallCells(0)
    {}
};


struct LightQoI
{
    std::string name;
    double value;

    LightQoI(const std::string name="default", const double value=0) : name(name), value(value) {}

    inline void update(const double up) { value += up; }
};

struct FatQoI : public LightQoI
{
    double alpha2; // associated volume fraction
    double lower_a2thresh;
    double upper_a2thresh;
    Real pos[3];
#ifdef _JONAS_STATS_
    double distance;
    double distance_x, distance_y, distance_z;
#endif /* _JONAS_STATS_ */

#ifdef _JONAS_STATS_
    FatQoI(const std::string name="default", const double value=0, const double lthresh=-HUGE_VAL, const double uthresh=HUGE_VAL) :
        LightQoI(name,value), alpha2(-1), lower_a2thresh(lthresh), upper_a2thresh(uthresh), distance(0), distance_x(0), distance_y(0), distance_z(0)
#else
    FatQoI(const std::string name="default", const double value=0, const double lthresh=-HUGE_VAL, const double uthresh=HUGE_VAL) :
        LightQoI(name,value), alpha2(-1), lower_a2thresh(lthresh), upper_a2thresh(uthresh)
#endif /* _JONAS_STATS_ */
    {
        for (int i = 0; i < 3; ++i)
            pos[i] = 0.0;
    }

    inline void update(const double new_value, const double new_alpha2, const Real new_pos[3])
    {
        value = new_value;
        alpha2 = new_alpha2;
        pos[0] = new_pos[0];
        pos[1] = new_pos[1];
        pos[2] = new_pos[2];
    }
};

template <bool _MIN>
struct ConditionalMinMaxQoI : public FatQoI
{
    ConditionalMinMaxQoI<_MIN>(const std::string name="default", const double value=0, const double lthresh=-HUGE_VAL, const double uthresh=HUGE_VAL) :
        FatQoI(name,value,lthresh,uthresh) {}

    inline void update(const double new_value, const double new_alpha2, const Real new_pos[3])
    {
        if (_MIN)
        {
            // conditional min
            if (new_value < value && (lower_a2thresh <= new_alpha2 && new_alpha2 <= upper_a2thresh))
            {
                value = new_value;
                alpha2 = new_alpha2;
                pos[0] = new_pos[0];
                pos[1] = new_pos[1];
                pos[2] = new_pos[2];
            }
        }
        else
        {
            // conditional max
            if (new_value > value && (lower_a2thresh <= new_alpha2 && new_alpha2 <= upper_a2thresh))
            {
                value = new_value;
                alpha2 = new_alpha2;
                pos[0] = new_pos[0];
                pos[1] = new_pos[1];
                pos[2] = new_pos[2];
            }
        }
    }
};

struct SphericalObserverQoI : public ConditionalMinMaxQoI<true> // we don't care whether true/false here, as update member is overloaded below
{
    double radius;
    double value_min;
    double value_max;
    unsigned long long int count;

    SphericalObserverQoI(const std::string name, const double value, const Real pos_[3], const double radius) :
        ConditionalMinMaxQoI<true>(name, value), radius(radius), count(0), value_min(HUGE_VAL), value_max(-HUGE_VAL)
    {
        pos[0] = pos_[0];
        pos[1] = pos_[1];
        pos[2] = pos_[2];
    }

    inline void update(const double new_value, const Real new_pos[3])
    {
        double new_radius2 = 0;
        for (int d = 0; d < 3; d++)
            new_radius2 += (pos[d] - new_pos[d]) * (pos[d] - new_pos[d]);

        if (new_radius2 <= radius * radius)
        {
            value += new_value;
            value_min = (new_value < value_min ? new_value : value_min);
            value_max = (new_value > value_max ? new_value : value_max);
            ++count;
        }
    }

    inline void update(const double new_value, const double new_min, const double new_max, const unsigned long long int lcount)
    {
        value += new_value;
        value_min = (new_min < value_min ? new_min : value_min);
        value_max = (new_max > value_max ? new_max : value_max);
        count += lcount;
    }
};

#endif /* STATISTICSMPI_H_MXLD8JPJ */
