/*
 *  OutputProcessingMPI.h
 *  MPCFcluster
 *
 *  Created by Fabian Wermelinger 09/18/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef OUTPUTPROCESSINGMPI_H_TQSCO9A3
#define OUTPUTPROCESSINGMPI_H_TQSCO9A3

#include <mpi.h>

#include "OutputProcessing.h"
#include "GridMPI.h"
#include "BlockLabMPI.h"
#include "BlockProcessor_MPI.h"

// MPI based dumper
#include "ZBinDumper_MPI.h"
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#include "HDF5Dumper_MPI.h"
#include "HDF5SliceDumperMPI.h"


#define __REGISTER_VP_ENTITY__(KEY, INFO, TYPE, HEAVY, TENTITY, TGRID, GRID, PARS, PROC, LAB) \
    this->_register(Item(#KEY"_rho",     new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Density", new SerializerWavelet<TGRID, StreamerDensity>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_Ux",      new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Velocity Ux", new SerializerWavelet<TGRID, StreamerVelocity<0> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_Uy",      new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Velocity Uy", new SerializerWavelet<TGRID, StreamerVelocity<1> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_Uz",      new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Velocity Uz", new SerializerWavelet<TGRID, StreamerVelocity<2> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_IUI",     new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Velocity magnitude", new SerializerWavelet<TGRID, StreamerVelocityMagnitude>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_divU",    new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Velocity divergence", new SerializerWavelet<TGRID, StreamerDivU>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY, StreamerDivU::CLASS, new OperatorType<TGRID,OdivU_4>(&PROC<LAB,OdivU_4,TGRID>)))); \
    this->_register(Item(#KEY"_KdivU",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Kdiv(U)", new SerializerWavelet<TGRID, StreamerKDivU>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY, StreamerKDivU::CLASS, new OperatorType<TGRID,OdivU_4>(&PROC<LAB,OdivU_4,TGRID>)))); \
    this->_register(Item(#KEY"_E",       new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Total energy", new SerializerWavelet<TGRID, StreamerEnergy>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_a2",      new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Alpha 2", new SerializerWavelet<TGRID, StreamerAlpha2>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_K",       new ProcessingElement<TENTITY,void*,TGRID>(#INFO": K", new SerializerWavelet<TGRID, StreamerK>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_p",       new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Pressure", new SerializerWavelet<TGRID, StreamerPressure>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_M",       new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Mach", new SerializerWavelet<TGRID, StreamerMach>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_c",       new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Speed of sound", new SerializerWavelet<TGRID, StreamerSpeedOfSound>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_Omegax",  new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Vorticity Omegax", new SerializerWavelet<TGRID, StreamerVorticity<0> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY, StreamerVorticity<0>::CLASS, new OperatorType<TGRID,OVort_4>(&PROC<LAB,OVort_4,TGRID>)))); \
    this->_register(Item(#KEY"_Omegay",  new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Vorticity Omegay", new SerializerWavelet<TGRID, StreamerVorticity<1> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY, StreamerVorticity<1>::CLASS, new OperatorType<TGRID,OVort_4>(&PROC<LAB,OVort_4,TGRID>)))); \
    this->_register(Item(#KEY"_Omegaz",  new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Vorticity Omegaz", new SerializerWavelet<TGRID, StreamerVorticity<2> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY, StreamerVorticity<2>::CLASS, new OperatorType<TGRID,OVort_4>(&PROC<LAB,OVort_4,TGRID>)))); \
    this->_register(Item(#KEY"_IOomegaI",new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Vorticity magnitude", new SerializerWavelet<TGRID, StreamerVorticityMagnitude>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY, StreamerVorticityMagnitude::CLASS, new OperatorType<TGRID,OVort_4>(&PROC<LAB,OVort_4,TGRID>)))); \
    this->_register(Item(#KEY"_Qcrit",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Q-criterion", new SerializerWavelet<TGRID, StreamerQcriterion>(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY, StreamerQcriterion::CLASS, new OperatorType<TGRID,OQcrit_4>(&PROC<LAB,OQcrit_4,TGRID>)))); \
    this->_register(Item(#KEY"_comp0",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 0", new SerializerWavelet<TGRID, StreamerAoScomponent<0> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp1",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 1", new SerializerWavelet<TGRID, StreamerAoScomponent<1> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp2",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 2", new SerializerWavelet<TGRID, StreamerAoScomponent<2> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp3",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 3", new SerializerWavelet<TGRID, StreamerAoScomponent<3> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp4",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 4", new SerializerWavelet<TGRID, StreamerAoScomponent<4> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp5",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 5", new SerializerWavelet<TGRID, StreamerAoScomponent<5> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp6",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 6", new SerializerWavelet<TGRID, StreamerAoScomponent<6> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp7",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 7", new SerializerWavelet<TGRID, StreamerAoScomponent<7> >(&PARS), &GRID, &PARS, NULL, TYPE, HEAVY))); \


template <typename TGrid>
class SerializerWaveletBase
{
public:
    virtual void write(const TGrid& grid, const int step_id, const Real t, const std::string basename, const std::string path) = 0;
};

template <typename TGrid, typename TStreamer, template <typename,typename> class TCompressor=SerializerIO_WaveletCompression_MPI_SimpleBlocking>
class SerializerWavelet : public SerializerWaveletBase<TGrid>
{
public:
    SerializerWavelet(ArgumentParser* const p)
    {
        m_serializer.verbose();
        // TODO: (fabianw@mavt.ethz.ch; Wed 19 Oct 2016 12:09:47 PM CEST) We
        // need to transition to wtype=3 as default.  Need to test difference
        // between QPX/SSE of type=1 and non-vectorized type=3.
        m_serializer.set_wtype_write((*p)("wtype").asInt(1));

        // streamer specific settings
        const Real vpeps = (*p)("vpeps").asDouble(1.0e-3);
        if (TStreamer::NAME == "Alpha2")
            m_serializer.set_threshold((*p)("vpeps_a2").asDouble(vpeps));
        else if (TStreamer::NAME == "Pressure")
            m_serializer.set_threshold((*p)("vpeps_p").asDouble(1000.0*vpeps));
        else if (TStreamer::NAME == "Density")
            m_serializer.set_threshold((*p)("vpeps_rho").asDouble(vpeps));
        else if (TStreamer::NAME == "Energy")
            m_serializer.set_threshold((*p)("vpeps_E").asDouble(vpeps));
        else if (TStreamer::NAME == "Qcriterion")
            m_serializer.set_threshold((*p)("vpeps_qcrit").asDouble(vpeps));
        else if (TStreamer::NAME == "Velocity Ux")
            m_serializer.set_threshold((*p)("vpeps_Ux").asDouble(vpeps));
        else if (TStreamer::NAME == "Velocity Uy")
            m_serializer.set_threshold((*p)("vpeps_Uy").asDouble(vpeps));
        else if (TStreamer::NAME == "Velocity Uz")
            m_serializer.set_threshold((*p)("vpeps_Uz").asDouble(vpeps));
        else if (TStreamer::NAME == "Velocity magnitude")
            m_serializer.set_threshold((*p)("vpeps_IUI").asDouble(vpeps));
        else if (TStreamer::NAME == "Vorticity Ox")
            m_serializer.set_threshold((*p)("vpeps_Ox").asDouble(vpeps));
        else if (TStreamer::NAME == "Vorticity Oy")
            m_serializer.set_threshold((*p)("vpeps_Oy").asDouble(vpeps));
        else if (TStreamer::NAME == "Vorticity Oz")
            m_serializer.set_threshold((*p)("vpeps_Oz").asDouble(vpeps));
        else if (TStreamer::NAME == "Vorticity magnitude")
            m_serializer.set_threshold((*p)("vpeps_IOI").asDouble(vpeps));
        else if (TStreamer::NAME == "K")
            m_serializer.set_threshold((*p)("vpeps_K").asDouble(vpeps));
        else if (TStreamer::NAME == "divU")
            m_serializer.set_threshold((*p)("vpeps_divU").asDouble(vpeps));
        else if (TStreamer::NAME == "Mach")
            m_serializer.set_threshold((*p)("vpeps_M").asDouble(vpeps));
        else if (TStreamer::NAME == "Speed of sound")
            m_serializer.set_threshold((*p)("vpeps_c").asDouble(vpeps));
        else if (TStreamer::NAME == "KdivU")
            m_serializer.set_threshold((*p)("vpeps_KdivU").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 0")
            m_serializer.set_threshold((*p)("vpeps_comp0").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 1")
            m_serializer.set_threshold((*p)("vpeps_comp1").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 2")
            m_serializer.set_threshold((*p)("vpeps_comp2").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 3")
            m_serializer.set_threshold((*p)("vpeps_comp3").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 4")
            m_serializer.set_threshold((*p)("vpeps_comp4").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 5")
            m_serializer.set_threshold((*p)("vpeps_comp5").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 6")
            m_serializer.set_threshold((*p)("vpeps_comp6").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 7")
            m_serializer.set_threshold((*p)("vpeps_comp7").asDouble(vpeps));
        else
            m_serializer.set_threshold(vpeps);
    }

    virtual void write(const TGrid& grid, const int step_id, const Real t, const std::string basename, const std::string path)
    {
        m_serializer.Write(grid, basename, path);
    }

private:
    TCompressor<TGrid,TStreamer> m_serializer;
};


template <typename TGrid>
struct SliceMPI : public Slice<TGrid>
{
    typedef TGrid GridType;

    MPI_Comm sliceComm;
    int localWidth, localHeight;
    int offsetWidth, offsetHeight;
    SliceMPI() : sliceComm(MPI_COMM_NULL), localWidth(-1), localHeight(-1), offsetWidth(-1), offsetHeight(-1) {}

    template <typename TSlice>
    static std::vector<TSlice> getSlices(ArgumentParser& parser, TGrid& grid)
    {
        std::vector<TSlice> slices = Slice<TGrid>::template getSlices<TSlice>(parser, grid);

        typedef typename TGrid::BlockType B;
        int Dim[3];
        Dim[0] = grid.getResidentBlocksPerDimension(0)*B::sizeX;
        Dim[1] = grid.getResidentBlocksPerDimension(1)*B::sizeY;
        Dim[2] = grid.getResidentBlocksPerDimension(2)*B::sizeZ;

        // get slice communicators
        int myRank;
        MPI_Comm_rank(grid.getCartComm(), &myRank);
        int peIdx[3];
        grid.peindex(peIdx);
        int myStart[3], myEnd[3];
        for (int i = 0; i < 3; ++i)
        {
            myStart[i] = Dim[i]*peIdx[i];
            myEnd[i]   = myStart[i] + Dim[i];
        }
        for (size_t i = 0; i < slices.size(); ++i)
        {
            TSlice& s = slices[i];
            const int sIdx = s.idx;
            const int dir  = s.dir;
            int color = 0;
            if (myStart[dir] <= sIdx && sIdx < myEnd[dir])
                color = 1; // gotcha!
            else
                s.valid = false;

            MPI_Comm_split(grid.getCartComm(), color, myRank, &s.sliceComm);

            // scale index to process local index
            s.idx = s.idx % Dim[s.dir];

            if (s.dir == 0)
            {
                s.localWidth  = Dim[2];
                s.localHeight = Dim[1];
                s.offsetWidth = peIdx[2]*Dim[2];
                s.offsetHeight= peIdx[1]*Dim[1];
            }
            else if (s.dir == 1)
            {
                s.localWidth  = Dim[2];
                s.localHeight = Dim[0];
                s.offsetWidth = peIdx[2]*Dim[2];
                s.offsetHeight= peIdx[0]*Dim[0];
            }
            else if (s.dir == 2)
            {
                s.localWidth  = Dim[0];
                s.localHeight = Dim[1];
                s.offsetWidth = peIdx[0]*Dim[0];
                s.offsetHeight= peIdx[1]*Dim[1];
            }
        }
        return slices;
    }
};


typedef BlockLabMPI<Lab> LabMPI;
typedef GridMPI<Grid_t> GridMPI_t;
typedef SliceMPI<GridMPI_t> SliceGridMPI_t;

typedef void (*Dumper_h5_full_grid_mpi)(const GridMPI_t&, const int, const Real, const std::string, const std::string, const bool);
typedef void (*Dumper_h5_slice_grid_mpi)(const SliceGridMPI_t&, const int, const Real, const std::string, const std::string, const bool);

// specialization
template <>
void ProcessingElement<GridMPI_t, Dumper_h5_full_grid_mpi, GridMPI_t>::operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool verbose)
{
    assert(m_outputFunctor != NULL);
    m_outputFunctor(*m_grid, step_id, t, basename, path, true);
}

template <>
void ProcessingElement<SliceProcessor<SliceGridMPI_t>, Dumper_h5_slice_grid_mpi, GridMPI_t>::operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool verbose)
{
    assert(m_outputFunctor != NULL);
    m_entity->setFunctor(m_outputFunctor);
    m_entity->process_all(step_id, t, basename, path);
}

template <>
void ProcessingElement<SerializerWaveletBase<GridMPI_t>, void*, GridMPI_t>::operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool verbose)
{
    assert(m_entity != NULL);
    m_entity->write(*m_grid, step_id, t, basename, path);
}

template <>
void ProcessingElement<SerializerWaveletBase<GridMPI_t>, void*, GridMPI_t>::dispose()
{
    if (m_entity)
        delete m_entity;
}


template <typename TGrid, template <typename> class TSlice=SliceMPI>
class OutputProcessingMPI : public OutputProcessing<TGrid,TSlice>
{
public:
    OutputProcessingMPI(ArgumentParser& p, TGrid& grid, const bool verbose=true) :
        OutputProcessing<TGrid,TSlice>(p,grid,verbose)
    {}


protected:
    typedef typename OutputProcessing<TGrid,TSlice>::TFunc_h5_grid  TFunc_h5_gridMPI;
    typedef typename OutputProcessing<TGrid,TSlice>::TFunc_h5_slice TFunc_h5_sliceMPI;
    typedef typename OutputProcessing<TGrid,TSlice>::TSP TSP;
    typedef typename OutputProcessing<TGrid,TSlice>::USlice USlice;
    typedef SerializerWaveletBase<TGrid> TVP;

    virtual void _register_h5_grid(TGrid& grid)
    {
        __REGISTER_H5_ENTITY__(h5, HDF5_FullMPI, PElement::H5, true, TGrid, TFunc_h5_gridMPI, TGrid, grid, grid, (this->m_parser), DumpHDF5_MPI, TGrid, process, LabMPI)
    }

    virtual void _register_h5_slice(TGrid& grid)
    {
        __REGISTER_H5_ENTITY__(h5s, HDF5_SliceMPI, PElement::H5SLICE, false, TSP, TFunc_h5_sliceMPI, TGrid, (this->m_sliceProcessor), grid, (this->m_parser), DumpSliceHDF5MPI, USlice, process, LabMPI)
    }

    void _register_vp(TGrid& grid)
    {
        __REGISTER_VP_ENTITY__(vp, VP_Wavelet, PElement::VP, true, TVP, TGrid, grid, (this->m_parser), process, LabMPI)
    }

    virtual void _register_all(TGrid& grid)
    {
        _register_h5_grid(grid);
        _register_h5_slice(grid);
        _register_vp(grid);
    }
};

#endif /* OUTPUTPROCESSINGMPI_H_TQSCO9A3 */
