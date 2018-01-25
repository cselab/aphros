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
#include "HDF5Dumper_MPI.h"
#include "HDF5SliceDumperMPI.h"


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

    virtual void _register_h5_grid(TGrid& grid)
    {
        __REGISTER_H5_ENTITY__(h5, HDF5_FullMPI, PElement::H5, true, TGrid, TFunc_h5_gridMPI, TGrid, grid, grid, (this->m_parser), DumpHDF5_MPI, TGrid, process, LabMPI)
    }

    virtual void _register_h5_slice(TGrid& grid)
    {
        __REGISTER_H5_ENTITY__(h5s, HDF5_SliceMPI, PElement::H5SLICE, false, TSP, TFunc_h5_sliceMPI, TGrid, (this->m_sliceProcessor), grid, (this->m_parser), DumpSliceHDF5MPI, USlice, process, LabMPI)
    }

    virtual void _register_all(TGrid& grid)
    {
        _register_h5_grid(grid);
        _register_h5_slice(grid);
    }
};

#endif /* OUTPUTPROCESSINGMPI_H_TQSCO9A3 */
