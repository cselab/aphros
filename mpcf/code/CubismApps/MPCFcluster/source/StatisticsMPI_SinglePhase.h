/*
 *  StatisticsMPI_SinglePhase.h
 *  MPCFcluster (Single-Phase Statistics)
 *
 *  Created by Fabian Wermelinegr 04/14/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef STATISTICSMPI_SINGLEPHASE_H_UMSPG42P
#define STATISTICSMPI_SINGLEPHASE_H_UMSPG42P

#include "FlowStep_LSRK3.h"
#include "StatisticsMPI.h"
#include "VectorOperator.h"

using namespace std;

#define USE_OPENMP 1

#ifdef _FLOAT_PRECISION_
#define POS_MPI_REAL MPI_FLOAT
#else
#define POS_MPI_REAL MPI_DOUBLE
#endif /* _FLOAT_PRECISION_ */

namespace SinglePhaseStatistics
{
    // ============================================================================
    // definition of sets of QoI's
    class IntegralSet
    {
        std::vector<LightQoI*> m_pelements;
        void _push_elements()
        {
            m_pelements.push_back (&r_avg);
            m_pelements.push_back (&u_avg);
            m_pelements.push_back (&v_avg);
            m_pelements.push_back (&w_avg);
            m_pelements.push_back (&IuI_avg);
            m_pelements.push_back (&W_avg);
            m_pelements.push_back (&ke_avg);
            m_pelements.push_back (&e_avg);
            m_pelements.push_back (&p_avg);
            m_pelements.push_back (&pw_avg);
            m_pelements.push_back (&c_avg);
            m_pelements.push_back (&M_avg);
            m_pelements.push_back (&rInt);
            m_pelements.push_back (&uInt);
            m_pelements.push_back (&vInt);
            m_pelements.push_back (&wInt);
            m_pelements.push_back (&EInt);
            m_pelements.push_back (&EkInt);
            m_pelements.push_back (&enstropy);
            m_pelements.push_back (&epsilon1);
            m_pelements.push_back (&epsilon2);
            m_pelements.push_back (&epsilon3);
        }

    public:
        LightQoI r_avg, u_avg, v_avg, w_avg, IuI_avg, W_avg, ke_avg, e_avg,
                 p_avg, pw_avg, c_avg, M_avg, rInt, uInt, vInt, wInt, EInt,
                 EkInt, enstropy, epsilon1, epsilon2, epsilon3;

        IntegralSet(const double val=0.) :
            r_avg  ("r_avg",  val), // density
            u_avg  ("u_avg",  val),
            v_avg  ("v_avg",  val),
            w_avg  ("w_avg",  val),
            IuI_avg("m_avg",  val), // magnitude of velocity vector
            W_avg  ("W_avg",  val), // magnitude of vorticity vector
            ke_avg ("ke_avg", val),
            e_avg  ("e_avg",  val),
            p_avg  ("p_avg",  val),
            pw_avg ("pw_avg", val),
            c_avg  ("c_avg",  val),
            M_avg  ("M_avg",  val),
            rInt   ("rInt",   val), // mass integral
            uInt   ("uInt",   val), // u-mom integral
            vInt   ("vInt",   val), // v-mom integral
            wInt   ("wInt",   val), // w-mom integral
            EInt   ("EInt",   val), // energy integral
            EkInt  ("EkInt",  val), // kinetic energy integral
            enstropy("enstrophy", val),
            epsilon1("epsilon1", val),
            epsilon2("epsilon2", val),
            epsilon3("epsilon3", val)
        {
            _push_elements();
        }

        IntegralSet(const IntegralSet& rhs) :
            r_avg(rhs.r_avg), u_avg(rhs.u_avg), v_avg(rhs.v_avg),
            w_avg(rhs.w_avg), IuI_avg(rhs.IuI_avg), W_avg(rhs.W_avg),
            ke_avg(rhs.ke_avg), e_avg(rhs.e_avg), p_avg(rhs.p_avg),
            pw_avg(rhs.pw_avg), c_avg(rhs.c_avg), M_avg(rhs.M_avg),
            rInt(rhs.rInt), uInt(rhs.uInt), vInt(rhs.vInt), wInt(rhs.wInt),
            EInt(rhs.EInt), EkInt(rhs.EkInt), enstropy(rhs.enstropy),
            epsilon1(rhs.epsilon1), epsilon2(rhs.epsilon2),
            epsilon3(rhs.epsilon3)
        {
            _push_elements();
        }

        inline void update(const IntegralSet& rhs)
        {
            for (int i = 0; i < m_pelements.size(); ++i)
                m_pelements[i]->update(rhs.m_pelements[i]->value);
        }

        inline void updateCell(const BlockData& data, const int ix, const int iy, const int iz, const BlockInfo& info)
        {
            r_avg.update(data.r);
            u_avg.update(data.u);
            v_avg.update(data.v);
            w_avg.update(data.w);
            IuI_avg.update(data.IuI);
            W_avg.update(data.W);
            ke_avg.update(data.ke);
            e_avg.update(data.e);
            p_avg.update(data.p);
            c_avg.update(data.c);
            M_avg.update(data.M);
            rInt.update(data.r);
            uInt.update(data.u*data.r);
            vInt.update(data.v*data.r);
            wInt.update(data.w*data.r);
            EInt.update(data.e);
            EkInt.update(data.ke);
            enstropy.update(data.kenstrophy);
            epsilon1.update(data.energyDissipationRate);
            epsilon2.update(data.divu*data.divu);
            epsilon3.update(data.p*data.divu);
            if (info.index[2]==0 && iz==0)
                pw_avg.update(data.p);
        }

        inline void updateBlock(const BlockData& data, const BlockInfo& info)
        {
        }

        inline void reduce(MPI_Comm comm, MPI_Op operation)
        {
            int myrank;
            MPI_Comm_rank(comm, &myrank);
            const bool isroot = (0 == myrank);

            std::vector<double> values(m_pelements.size());
            for (int i = 0; i < m_pelements.size(); i++)
                values[i] = m_pelements[i]->value;

            MPI_Reduce(isroot ? MPI_IN_PLACE : &values[0], &values[0], m_pelements.size(), MPI_DOUBLE, operation, 0, comm);

            for (int i = 0; i < m_pelements.size(); i++)
                m_pelements[i]->value = values[i];
        }

        inline void postMPI(const GlobalData& data)
        {
            const double factor = 1.0 / (data.nBlocks*data.V_block);
            const double h3 = data.h * data.h * data.h;
            // used for non-dimensionalization (this certainly depends on how
            // your case is set up and might result in wrong statistics results
            // if you use the same identifiers!
            const double rho0 = data.parser("rho0").asDouble(1.0);
            const double u0 = data.parser("u0").asDouble(1.0);
            const double L  = data.parser("L").asDouble(1.0);
            const double mu = Simulation_Environment::MU_MAX; // single phase
            const double tc_inv = u0/L;
            const double tc_inv2 = tc_inv*tc_inv;
            const double u0_inv2 = 1.0/(u0*u0);
            const double eps_normfac = (u0*u0*u0)/L;

            r_avg.value   *= factor;
            u_avg.value   *= factor;
            v_avg.value   *= factor;
            w_avg.value   *= factor;
            IuI_avg.value *= factor;
            W_avg.value   *= factor;
            ke_avg.value  *= factor;
            e_avg.value   *= factor;
            p_avg.value   *= factor;
            c_avg.value   *= factor;
            M_avg.value   *= factor;
            enstropy.value *= factor/rho0*tc_inv2;
            epsilon1.value *= factor*2.0*mu/rho0*eps_normfac;
            epsilon2.value *= factor/rho0*eps_normfac; // misses the factor mu_v (bulk viscosity)
            epsilon3.value *= -factor/rho0*eps_normfac;
            rInt.value    *= h3;
            uInt.value    *= h3;
            vInt.value    *= h3;
            wInt.value    *= h3;
            EInt.value    *= h3;
            EkInt.value   *= h3;
            if (data.wallCells > 0)
                pw_avg.value /= data.wallCells;
        }

        inline void writeHeader(ofstream& stream) const
        {
            for (int i = 0; i < m_pelements.size(); ++i)
                stream << " " << m_pelements[i]->name;
        }

        inline void writeData(ofstream& stream) const
        {
            for (int i = 0; i < m_pelements.size(); ++i)
                stream << " " << m_pelements[i]->value;
        }
    };

    template <typename T>
    class GenericFatSet
    {
        std::vector<T*> m_pelements;
        void _push_elements()
        {
            m_pelements.push_back (&r);
            m_pelements.push_back (&u);
            m_pelements.push_back (&v);
            m_pelements.push_back (&w);
            m_pelements.push_back (&IuI);
            m_pelements.push_back (&W);
            m_pelements.push_back (&ke);
            m_pelements.push_back (&e);
            m_pelements.push_back (&p);
            m_pelements.push_back (&c);
            m_pelements.push_back (&M);
            m_pelements.push_back (&pw);
            m_pelements.push_back (&divu);     // minimum value of divergence of velocity
            m_pelements.push_back (&re); // Reynolds number
        }

    public:
        T r, u, v, w, IuI, W, ke, e, p, c, M, pw, divu, re;

        GenericFatSet(const std::string suffix, const double val, const double lower=-HUGE_VAL, const double upper=HUGE_VAL) :
            r    ("r_" + suffix,    val, lower, upper),
            u    ("u_" + suffix,    val, lower, upper),
            v    ("v_" + suffix,    val, lower, upper),
            w    ("w_" + suffix,    val, lower, upper),
            IuI  ("m_" + suffix,    val, lower, upper),
            W    ("W_" + suffix,    val, lower, upper),
            ke   ("ke_" + suffix,   val, lower, upper),
            e    ("e_" + suffix,    val, lower, upper),
            p    ("p_" + suffix,    val, lower, upper),
            c    ("c_" + suffix,    val, lower, upper),
            M    ("M_" + suffix,    val, lower, upper),
            pw   ("pw_" + suffix,   val, lower, upper),
            divu ("divu_" + suffix, val, lower, upper),
            re   ("re_" + suffix, val, lower, upper)
        {
            _push_elements();
        }

        GenericFatSet(const GenericFatSet& rhs) :
            r(rhs.r), u(rhs.u), v(rhs.v), w(rhs.w), IuI(rhs.IuI), W(rhs.W),
            ke(rhs.ke), e(rhs.e), p(rhs.p), c(rhs.c), M(rhs.M),
            pw(rhs.pw), divu (rhs.divu), re(rhs.re)
        {
            _push_elements();
        }

        inline void update(const GenericFatSet<T>& rhs)
        {
            for (int i = 0; i < m_pelements.size(); ++i)
                m_pelements[i]->update(rhs.m_pelements[i]->value, rhs.m_pelements[i]->alpha2, rhs.m_pelements[i]->pos);
        }

        inline void updateCell(const BlockData& data, const int ix, const int iy, const int iz, const BlockInfo& info)
        {
            Real pos[3];
            info.pos(pos, ix, iy, iz);
            r.update(data.r, data.a2, pos);
            u.update(data.u, data.a2, pos);
            v.update(data.v, data.a2, pos);
            w.update(data.w, data.a2, pos);
            ke.update(data.ke, data.a2, pos);
            IuI.update(data.IuI, data.a2, pos);
            W.update(data.W, data.a2, pos);
            e.update(data.e, data.a2, pos);
            p.update(data.p, data.a2, pos);
            c.update(data.c, data.a2, pos);
            M.update(data.M, data.a2, pos);
            divu.update(data.divu, data.a2, pos);
            re.update(data.IuI, data.a2, pos);
            if (info.index[2]==0 && iz==0)
                pw.update(data.p, data.a2, pos);
        }

        inline void updateBlock(const BlockData& data, const BlockInfo& info) {}

        inline void reduce(MPI_Comm comm, MPI_Op operation)
        {
            int myrank;
            MPI_Comm_rank(comm, &myrank);
            const bool isroot = (0 == myrank);

            std::vector<double_int> values(m_pelements.size());
            for (int i = 0; i < m_pelements.size(); i++)
            {
                values[i].value = m_pelements[i]->value;
                values[i].rank  = myrank;
            }

            MPI_Allreduce(MPI_IN_PLACE, &values[0], values.size(), MPI_DOUBLE_INT, operation, comm);

            std::vector<MPI_Request> request_pos(m_pelements.size(), MPI_REQUEST_NULL);
            std::vector<MPI_Request> request_alpha2(m_pelements.size(), MPI_REQUEST_NULL);
            for (int i = 0; i < m_pelements.size(); i++)
            {
                int targetRank = values[i].rank;
                m_pelements[i]->value = values[i].value;
                if (targetRank != 0)
                {
                    if (myrank == targetRank)
                    {
                        MPI_Isend(&(m_pelements[i]->pos[0]), 3, POS_MPI_REAL, 0, 101, comm, &request_pos[i]);
                        MPI_Isend(&(m_pelements[i]->alpha2), 1, MPI_DOUBLE, 0, 102, comm, &request_alpha2[i]);
                    }
                    else if (isroot)
                    {
                        MPI_Irecv(&(m_pelements[i]->pos[0]), 3, POS_MPI_REAL, targetRank, 101, comm, &request_pos[i]);
                        MPI_Irecv(&(m_pelements[i]->alpha2), 1, MPI_DOUBLE, targetRank, 102, comm, &request_alpha2[i]);
                    }
                }
            }
            MPI_Waitall(m_pelements.size(), request_pos.data(), MPI_STATUSES_IGNORE);
            MPI_Waitall(m_pelements.size(), request_alpha2.data(), MPI_STATUSES_IGNORE);
        }

        inline void postMPI(const GlobalData& data)
        {
            // used for non-dimensionalization (this certainly depends on how
            // your case is set up and might result in wrong statistics results
            // if you use the same identifiers!
            const double rho0 = data.parser("rho0").asDouble(1.0);
            const double L = data.parser("L").asDouble(1.0);
            const double reFac = L*rho0/Simulation_Environment::MU_MAX;
            re.value *= reFac;

#ifdef _JONAS_STATS_
            for (int i = 0; i < m_pelements.size(); i++)
            {
                m_pelements[i]->distance =   sqrt((data.center[0] - m_pelements[i]->pos[0]) * (data.center[0] - m_pelements[i]->pos[0]) + (data.center[1] - m_pelements[i]->pos[1]) * (data.center[1] - m_pelements[i]->pos[1]) + (data.center[2] - m_pelements[i]->pos[2]) * (data.center[2] - m_pelements[i]->pos[2]) );
                m_pelements[i]->distance_x = sqrt((data.center[1] - m_pelements[i]->pos[1]) * (data.center[1] - m_pelements[i]->pos[1]) + (data.center[2] - m_pelements[i]->pos[2]) * (data.center[2] - m_pelements[i]->pos[2]) );
                m_pelements[i]->distance_y = sqrt((data.center[0] - m_pelements[i]->pos[0]) * (data.center[0] - m_pelements[i]->pos[0]) + (data.center[2] - m_pelements[i]->pos[2]) * (data.center[2] - m_pelements[i]->pos[2]) );
                m_pelements[i]->distance_z = sqrt((data.center[0] - m_pelements[i]->pos[0]) * (data.center[0] - m_pelements[i]->pos[0]) + (data.center[1] - m_pelements[i]->pos[1]) * (data.center[1] - m_pelements[i]->pos[1]) );
            }
#endif /* _JONAS_STATS_ */
        }

        inline void writeHeader(ofstream& stream) const
        {
            for (int i = 0; i < m_pelements.size(); ++i)
            {
                stream << " " << m_pelements[i]->name;
                stream << " " << m_pelements[i]->name << "_pos_x";
                stream << " " << m_pelements[i]->name << "_pos_y";
                stream << " " << m_pelements[i]->name << "_pos_z";
#ifdef _JONAS_STATS_
                stream << " " << m_pelements[i]->name << "_pos_d";
                stream << " " << m_pelements[i]->name << "_pos_d_x";
                stream << " " << m_pelements[i]->name << "_pos_d_y";
                stream << " " << m_pelements[i]->name << "_pos_d_z";
#endif /* _JONAS_STATS_ */
            }
        }

        inline void writeData(ofstream& stream) const
        {
            for (int i = 0; i < m_pelements.size(); ++i)
            {
                stream << " " << m_pelements[i]->value;
                stream << " " << m_pelements[i]->pos[0];
                stream << " " << m_pelements[i]->pos[1];
                stream << " " << m_pelements[i]->pos[2];
#ifdef _JONAS_STATS_
                stream << " " << m_pelements[i]->distance;
                stream << " " << m_pelements[i]->distance_x;
                stream << " " << m_pelements[i]->distance_y;
                stream << " " << m_pelements[i]->distance_z;
#endif /* _JONAS_STATS_ */
            }
        }
    };

    class SphericalObserverSet
    {
        std::vector<SphericalObserverQoI*> m_pelements;
        void _push_elements()
        {
            m_pelements.push_back(&r);
            m_pelements.push_back(&u);
            m_pelements.push_back(&v);
            m_pelements.push_back(&w);
            m_pelements.push_back(&IuI);
            m_pelements.push_back(&W);
            m_pelements.push_back(&ke);
            m_pelements.push_back(&e);
            m_pelements.push_back(&p);
            m_pelements.push_back(&c);
            m_pelements.push_back(&M);
            m_pelements.push_back(&divu);
            m_pelements.push_back(&gradp);
        }

    public:
        SphericalObserverQoI r, u, v, w, IuI, W, ke, e, p, c, M, divu, gradp;

        SphericalObserverSet(const std::string suffix, const double val, const Real pos[3], const Real rad) :
            r    ("r_"     + suffix, val, pos, rad),
            u    ("u_"     + suffix, val, pos, rad),
            v    ("v_"     + suffix, val, pos, rad),
            w    ("w_"     + suffix, val, pos, rad),
            IuI  ("IuI_"   + suffix, val, pos, rad),
            W    ("W_"     + suffix, val, pos, rad),
            ke   ("ke_"    + suffix, val, pos, rad),
            e    ("e_"     + suffix, val, pos, rad),
            p    ("p_"     + suffix, val, pos, rad),
            c    ("c_"     + suffix, val, pos, rad),
            M    ("M_"     + suffix, val, pos, rad),
            divu ("divu_"  + suffix, val, pos, rad),
            gradp("gradp_" + suffix, val, pos, rad)
        {
            _push_elements();
        }

        SphericalObserverSet(const SphericalObserverSet& rhs) :
            r(rhs.r), u(rhs.u), v(rhs.v), w(rhs.w), IuI(rhs.IuI),
            W(rhs.W), ke(rhs.ke), e(rhs.e), p(rhs.p), c(rhs.c), M(rhs.M),
            divu(rhs.divu), gradp(rhs.gradp)
        {
            _push_elements();
        }

        inline void update(const SphericalObserverSet& rhs)
        {
            for (int i = 0; i < m_pelements.size(); ++i)
                m_pelements[i]->update(rhs.m_pelements[i]->value, rhs.m_pelements[i]->value_min, rhs.m_pelements[i]->value_max, rhs.m_pelements[i]->count);
        }

        inline void updateCell(const BlockData& data, const int ix, const int iy, const int iz, const BlockInfo& info)
        {
            Real pos[3];
            info.pos(pos, ix, iy, iz);
            r.update(data.r, pos);
            u.update(data.u, pos);
            v.update(data.v, pos);
            w.update(data.w, pos);
            IuI.update(data.IuI, pos);
            W.update(data.W, pos);
            ke.update(data.ke, pos);
            e.update(data.e, pos);
            p.update(data.p, pos);
            c.update(data.c, pos);
            M.update(data.M, pos);
            divu.update(data.divu, pos);
            gradp.update(data.gradp, pos);
        }

        inline void updateBlock(const BlockData& data, const BlockInfo& info) {}

        inline void reduce(MPI_Comm comm, MPI_Op operation)
        {
            int myrank;
            MPI_Comm_rank(comm, &myrank);
            const bool isroot = (0 == myrank);

            std::vector<double> values(m_pelements.size());
            std::vector<double> values_min(m_pelements.size());
            std::vector<double> values_max(m_pelements.size());
            std::vector<unsigned long long int> counts(m_pelements.size());
            for (int j = 0; j < m_pelements.size(); j++)
            {
                values[j]     = m_pelements[j]->value;
                values_min[j] = m_pelements[j]->value_min;
                values_max[j] = m_pelements[j]->value_max;
                counts[j]     = m_pelements[j]->count;
            }

            MPI_Reduce(isroot ? MPI_IN_PLACE : &values[0], &values[0], values.size(), MPI_DOUBLE, operation, 0, comm);
            MPI_Reduce(isroot ? MPI_IN_PLACE : &values_min[0], &values_min[0], values_min.size(), MPI_DOUBLE, MPI_MIN, 0, comm);
            MPI_Reduce(isroot ? MPI_IN_PLACE : &values_max[0], &values_max[0], values_max.size(), MPI_DOUBLE, MPI_MAX, 0, comm);
            MPI_Reduce(isroot ? MPI_IN_PLACE : &counts[0], &counts[0], counts.size(), MPI_UNSIGNED_LONG_LONG, operation, 0, comm);

            for (int j = 0; j < m_pelements.size(); j++)
            {
                m_pelements[j]->value     = values[j];
                m_pelements[j]->value_min = values_min[j];
                m_pelements[j]->value_max = values_max[j];
                m_pelements[j]->count     = counts[j];
            }
        }

        inline void postMPI(const GlobalData& data)
        {
            for (int j = 0; j < m_pelements.size(); j++)
                if (m_pelements[j]->count != 0)
                    m_pelements[j]->value /= m_pelements[j]->count;
        }

        inline void postMPI(const GlobalData& data, const SphericalObserverSet& other)
        {
            for (int j = 0; j < m_pelements.size(); j++)
            {
                const double sum1 = m_pelements[j]->value;
                const double c1   = static_cast<double>(m_pelements[j]->count);
                const double sum0 = other.m_pelements[j]->value;
                const double c0   = static_cast<double>(other.m_pelements[j]->count);
                m_pelements[j]->value = ((c1-c0) > 0. ? (sum1-sum0)/(c1-c0) : 0.0);
            }
        }

        inline void writeHeader(ofstream& stream) const
        {
            for (int i = 0; i < m_pelements.size(); ++i)
            {
                stream << " " << m_pelements[i]->name;
                stream << " " << m_pelements[i]->name << "_min";
                stream << " " << m_pelements[i]->name << "_max";
            }
        }

        inline void writeData(ofstream& stream) const
        {
            for (int i = 0; i < m_pelements.size(); ++i)
            {
                stream << " " << m_pelements[i]->value;
                stream << " " << m_pelements[i]->value_min;
                stream << " " << m_pelements[i]->value_max;
            }
        }
    };

    typedef ConditionalMinMaxQoI<true>  ConditionalMinQoI;
    typedef ConditionalMinMaxQoI<false> ConditionalMaxQoI;
    typedef GenericFatSet<ConditionalMinQoI> ConditionalMinimumSet;
    typedef GenericFatSet<ConditionalMaxQoI> ConditionalMaximumSet;
    // ============================================================================
    // helper
    vector<SphericalObserverSet> getSensors(ArgumentParser& parser, const GlobalData& data)
    {
        const int sensors = parser("-sensors").asInt(0);
        vector<SphericalObserverSet> ret;
        for (int i = 0; i < sensors; i++)
        {
            std::stringstream id;
            id << i + 1;

            // position of the sensor
            Real sen_pos [3];
            sen_pos [0] = parser ( "sensor" + id.str() + "_pos_x") .asDouble (data.center[0]);
            sen_pos [1] = parser ( "sensor" + id.str() + "_pos_y") .asDouble (data.center[1]);
            sen_pos [2] = parser ( "sensor" + id.str() + "_pos_z") .asDouble (data.center[2]);

            // size of the sensor
            double fallback = 0.01 * parser("-extent").asDouble(1.0);
            double radius   = parser ( "sensor" + id.str() + "_rad").asDouble(fallback);

            ret.push_back( SphericalObserverSet("sensor" + id.str(), 0.0, sen_pos, radius) );
        }
        return ret;
    }

    vector<SphericalObserverSet> getShells(ArgumentParser& parser, const GlobalData& data)
    {
        int n_shells = parser("-shells").asInt(0);
        vector<SphericalObserverSet> ret;
        if (n_shells > 0)
        {
            // get radius of cloud
            const double fallback  = 0.1 * parser("-extent").asDouble(1.0);
            const double cloud_rad = parser("Rc").asDouble(fallback);
            const double delta = cloud_rad/n_shells; // increment

            Real center[3];
            for (int i = 0; i < 3; ++i)
                center[i] = data.center[i];

            for (int i = 0; i < n_shells; i++)
            {
                std::stringstream id;
                id << i + 1;

                const double radius = (i+1)*delta;
                ret.push_back( SphericalObserverSet("shell_avg" + id.str(), 0.0, center, radius) );
            }
        }
        return ret;
    }

    // ============================================================================
    template <typename TGrid>
    void dumpStatistics(TGrid& grid, const int step_id, const Real t, const Real dt, ArgumentParser& parser)
    {
        typedef typename TGrid::BlockType TBlock;

        double t0, t1;

        MPI_Comm world_comm = grid.getCartComm();
        int myrank;
        MPI_Comm_rank(world_comm, &myrank);
        const bool isroot = (0 == myrank);

        // global stuff
        GlobalData gdata(parser);
        std::vector<BlockInfo> g_vInfo = grid.getBlocksInfo();
        gdata.h = g_vInfo.front().h_gridpoint;
        gdata.V_block = TBlock::sizeX * TBlock::sizeY * TBlock::sizeZ;
        gdata.nBlocks = 1;
        for (int i = 0; i < 3; ++i)
            gdata.nBlocks *= grid.getBlocksPerDimension(i);
        gdata.center[0] = 0.5 * Simulation_Environment::extents[0];
        gdata.center[1] = 0.5 * Simulation_Environment::extents[1];
        gdata.center[2] = 0.5 * Simulation_Environment::extents[2];

        t0 = omp_get_wtime();

        // integral quantities
        IntegralSet integralSet;

        // global minimum quantities
        ConditionalMinimumSet globalMinSet("global_min", HUGE_VAL);

        // global maximum quantities
        ConditionalMaximumSet globalMaxSet("global_max", -HUGE_VAL);

        // sensors
        std::vector<SphericalObserverSet> sensors = getSensors(parser, gdata);

        // shells
        std::vector<SphericalObserverSet> shells = getShells(parser, gdata);


        const double inv_gam1m1 =  1.0 /(Simulation_Environment::GAMMA1-1.0);
        const double inv_gam2m1 =  1.0 /(Simulation_Environment::GAMMA2-1.0);

        const int sponge = parser("-sponge").asInt(0);

        const bool bAll = parser("-statistics_computeAll").asBool(false);

        // compute vorticity
        const bool bVort = parser("-statistics_computeVorticity").asBool(false);
        if (bVort || bAll) evaluateOperator_MPI<LabMPI, OVort_4>(grid); // this writes into tmp.[0], tmp.[1] and tmp.[2]
        // compute div(u)
        const bool bDivu = parser("-statistics_computeDivu").asBool(false);
        if (bDivu || bAll) evaluateOperator_MPI<LabMPI, OdivU_4>(grid); // this writes into tmp.[3]
        // S^d : S^d
        const bool bSS = parser("-statistics_computeS:S").asBool(false);
        if (bSS || bAll) evaluateOperator_MPI<LabMPI, Operator_SSDeviatoric_4th>(grid); // @dummy
        // |grap(p)|
        const bool bGradp = false;
        // const bool bGradp = parser("-statistics_ComputeGradp").asBool(false);
        // if (bGradp || bAll) evaluateOperator_MPI<LabMPI, EvaluateMagGradientPressureFourthOrder_CPP>(grid); // @tmp[4]

        t1 = omp_get_wtime();
        if (isroot) printf("statistics: A took %lf seconds\n", t1-t0);
        t0 = omp_get_wtime();

#ifdef USE_OPENMP
#pragma omp parallel
#endif
        {

#if 1
            // thread local integral quantities
            IntegralSet myintegralSet;

            // global minimum quantities
            ConditionalMinimumSet myglobalMinSet("global_min", HUGE_VAL);

            // global maximum quantities
            ConditionalMaximumSet myglobalMaxSet("global_max", -HUGE_VAL);

            // sensors
            std::vector<SphericalObserverSet> mysensors = getSensors(parser, gdata);

            // shells
            std::vector<SphericalObserverSet> myshells = getShells(parser, gdata);

#endif
            unsigned long long int myWallCells = 0;


#pragma omp for schedule(dynamic,1)
            for(int i=0; i<(int)g_vInfo.size(); i++)
            {
                BlockInfo info = g_vInfo[i];

                bool bBoundary = 0;
                for (int d = 0; d < 3; d++)
                    if ((info.index[d]==0 && ! Simulation_Environment::BC_PERIODIC[d]) || (info.index[d]==grid.getBlocksPerDimension(d)-1 && ! Simulation_Environment::BC_PERIODIC[d]))
                        bBoundary = 1;

                // to handle symmetric domains with sponge correctly
#ifdef _BCLABCLOUDSYMABSORB_
                for (int d = 0; d < 3; d++)
                    if (info.index[d]==0)
                        bBoundary = 0;
#endif

                // ignore sponge blocks
                if (sponge && bBoundary)
                    continue;

                TBlock& b = *(TBlock*)info.ptrBlock;

                BlockData bdata;
                bdata.V = gdata.V_block;

                for(int iz=0; iz<TBlock::sizeZ; iz++)
                    for(int iy=0; iy<TBlock::sizeY; iy++)
                        for(int ix=0; ix<TBlock::sizeX; ix++)
                        {
#ifndef _CONVERTCLIP_
                            bdata.r1 = (double)b(ix, iy, iz).alpha1rho1;
                            bdata.r2 = (double)b(ix, iy, iz).alpha2rho2;
#else
                            bdata.r1 = max(0.0, (double)b(ix, iy, iz).alpha1rho1);
                            bdata.r2 = max(0.0, (double)b(ix, iy, iz).alpha2rho2);
#endif
                            bdata.r = bdata.r1 + bdata.r2;

                            bdata.u = (double)b(ix, iy, iz).ru / bdata.r;
                            bdata.v = (double)b(ix, iy, iz).rv / bdata.r;
                            bdata.w = (double)b(ix, iy, iz).rw / bdata.r;
                            bdata.e = (double)b(ix, iy, iz).energy;
#ifndef _CONVERTCLIP_
#ifdef _ALPHACLIP_
                            bdata.a2 = max(ALPHAEPS, min(1.0-ALPHAEPS,(double)b(ix, iy, iz).alpha2));
                            bdata.a1 = 1.0-bdata.a2;
#else
                            bdata.a2 = (double)b(ix, iy, iz).alpha2;
                            bdata.a1 = 1.0-bdata.a2;
#endif /* _ALPHACLIP_ */
#else
                            bdata.a2 = max(ALPHAEPS, min(1.0-ALPHAEPS,(double)b(ix, iy, iz).alpha2));
                            bdata.a1 = 1.0-bdata.a2;
#endif
                            // vorticity components
                            bdata.W0 = ((bVort || bAll) ? b.tmp[ix][iy][iz][0] : 0.0);
                            bdata.W1 = ((bVort || bAll) ? b.tmp[ix][iy][iz][1] : 0.0);
                            bdata.W2 = ((bVort || bAll) ? b.tmp[ix][iy][iz][2] : 0.0);

                            const double uvw2 = bdata.u*bdata.u + bdata.v*bdata.v + bdata.w*bdata.w;
                            bdata.IuI = sqrt (uvw2);

                            const double IWI2 = bdata.W0*bdata.W0 + bdata.W1*bdata.W1 + bdata.W2*bdata.W2;
                            bdata.W  = sqrt(IWI2);
                            bdata.ke = 0.5 * bdata.r * (uvw2);
                            bdata.kenstrophy = 0.5 * bdata.r * IWI2;

                            // pressure
                            bdata.p = Simulation_Environment::pressure_5eq(bdata.r, bdata.e, b(ix,iy,iz).ru, b(ix,iy,iz).rv, b(ix,iy,iz).rw, bdata.a2);

                            // speed of sound
#ifdef _NOK_
                            const double gamma_mix_m1 = 1.0 / (bdata.a1 * inv_gam1m1 + bdata.a2 * inv_gam2m1);
                            const double gamPc_mix = bdata.a1 * Simulation_Environment::GAMMA1 * Simulation_Environment::PC1 * inv_gam1m1 + bdata.a2 * Simulation_Environment::GAMMA2 * Simulation_Environment::PC2 * inv_gam2m1;
                            bdata.c = sqrt(((gamma_mix_m1 + 1.0) * bdata.p + gamPc_mix * gamma_mix_m1)/bdata.r);
#else
                            const double _g1 = Simulation_Environment::GAMMA1;
                            const double _g2 = Simulation_Environment::GAMMA2;
                            const double _pc1 = Simulation_Environment::PC1;
                            const double _pc2 = Simulation_Environment::PC2;
                            bdata.c = sqrt(_g1*_g2*(bdata.p+_pc1)*(bdata.p+_pc2) / ((bdata.a1*_g2*(bdata.p+_pc2) + bdata.a2*_g1*(bdata.p+_pc1))*bdata.r));
#endif

                            bdata.dummy = b(ix,iy,iz).dummy;
                            bdata.M = bdata.IuI / bdata.c;
                            bdata.divu = ((bDivu || bAll) ? b.tmp[iz][iy][ix][3] : 0.0);
                            bdata.energyDissipationRate = ((bSS || bAll) ? bdata.dummy : 0.0);
                            bdata.gradp = ((bGradp || bAll) ? b.tmp[iz][iy][ix][4] : 0.0);

                            if (info.index[2]==0 && iz==0)
                                ++bdata.wallCells;

                            // === update integral quantities
                            myintegralSet.updateCell(bdata, ix, iy, iz, info);

                            // === update extremal quantities
                            myglobalMinSet.updateCell(bdata, ix, iy, iz, info);
                            myglobalMaxSet.updateCell(bdata, ix, iy, iz, info);

                            // === update sensors
                            for (int i = 0; i < mysensors.size(); i++)
                            {
                                SphericalObserverSet& thisSensor = mysensors[i];
                                thisSensor.updateCell(bdata, ix, iy, iz, info);
                            }

                            // === update shells
                            for (int i = 0; i < myshells.size(); i++)
                            {
                                SphericalObserverSet& thisShell = myshells[i];
                                thisShell.updateCell(bdata, ix, iy, iz, info);
                            }
                        }

                myWallCells += bdata.wallCells;
            }  // omp for

#if 1

            // update global statistics
#pragma omp critical
            {
                gdata.wallCells += myWallCells;

                integralSet.update(myintegralSet);

                globalMinSet.update(myglobalMinSet);
                globalMaxSet.update(myglobalMaxSet);

                // === update sensors
                for (int i = 0; i < sensors.size(); i++)
                    sensors[i].update(mysensors[i]);

                // === update shells
                for (int i = 0; i < shells.size(); i++)
                    shells[i].update(myshells[i]);
            }  // critical

#endif

        }  // parallel

        t1 = omp_get_wtime();
        if (isroot) printf("statistics: B took %lf seconds\n", t1-t0);
        t0 = omp_get_wtime();

        integralSet.reduce(world_comm, MPI_SUM);

        globalMinSet.reduce(world_comm, MPI_MINLOC);
        globalMaxSet.reduce(world_comm, MPI_MAXLOC);

        for (int i = 0; i < sensors.size(); i++)
        {
            SphericalObserverSet& thisSensor = sensors[i];
            thisSensor.reduce(world_comm, MPI_SUM);
        }

        for (int i = 0; i < shells.size(); i++)
        {
            SphericalObserverSet& thisShell = shells[i];
            thisShell.reduce(world_comm, MPI_SUM);
        }

        // Non-Set related reductions
        // reduce the total number of wall cells
        MPI_Reduce(isroot ? MPI_IN_PLACE : &gdata.wallCells, &gdata.wallCells, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, world_comm);

        MPI_Barrier(world_comm);

        // Post MPI computations
        integralSet.postMPI(gdata); // computes averages

        globalMinSet.postMPI(gdata);
        globalMaxSet.postMPI(gdata);

        for (int i = 0; i < sensors.size(); i++)
        {
            SphericalObserverSet& thisSensor = sensors[i];
            thisSensor.postMPI(gdata); // computes averages
        }

        for (int i = (int)shells.size()-1; i > 0; --i)
        {
            SphericalObserverSet& outerObserver = shells[i];
            const SphericalObserverSet& innerObserver = shells[i-1];
            // compute shell averages by taking difference of two spherical
            // spherical observer from outside inwards (skips innermost
            // spherical observer which is treated as a sphere instead)
            outerObserver.postMPI(gdata, innerObserver);
        }
        if (shells.size() > 0)
            shells[0].postMPI(gdata); // innermost observer

        t1 = omp_get_wtime();
        if (isroot) printf("statistics: C took %lf seconds\n", t1-t0);
        t0 = omp_get_wtime();

        // write output file
        if (isroot)
        {
            std::string filename = "statistics.dat";

            // write header
            std::ofstream f;
            f.open(filename.c_str(), std::ofstream::in);
            if ( f.good() )
                f.close();
            else {
                f.close();
                f.open(filename.c_str(), std::ofstream::out);
                f << "step t dt";

                integralSet.writeHeader(f);

                globalMinSet.writeHeader(f);
                globalMaxSet.writeHeader(f);

                for (int i = 0; i < sensors.size(); i++)
                {
                    const SphericalObserverSet& thisSensor = sensors[i];
                    thisSensor.writeHeader(f);
                }

                for (int i = 0; i < shells.size(); i++)
                {
                    const SphericalObserverSet& thisShell = shells[i];
                    thisShell.writeHeader(f);
                }

                f << std::endl;
                f.close();
            }

            f.open (filename.c_str(), std::ofstream::app);
            f << step_id << " " << t << " " << dt;
            f << std::scientific;

            integralSet.writeData(f);
            globalMinSet.writeData(f);
            globalMaxSet.writeData(f);

            for (int i = 0; i < sensors.size(); i++)
            {
                const SphericalObserverSet& thisSensor = sensors[i];
                thisSensor.writeData(f);
            }

            for (int i = 0; i < shells.size(); i++)
            {
                const SphericalObserverSet& thisShell = shells[i];
                thisShell.writeData(f);
            }

            f << std::endl;
            f.close();
        }

        t1 = omp_get_wtime();
        if (isroot) printf("statistics: D took %lf seconds\n", t1-t0);
    }
}

#endif /* STATISTICSMPI_SINGLEPHASE_H_UMSPG42P */
