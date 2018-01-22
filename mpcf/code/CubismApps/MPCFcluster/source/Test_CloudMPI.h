/*
 *  Test_CloudMPI.h
 *  MPCFcluster
 *
 *  Created by Babak Hejazialhosseini on 2/25/13.
 *  Completly revised by Ursula Rasthofer in 2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_CLOUDMPI_H_IBF8M2BV
#define TEST_CLOUDMPI_H_IBF8M2BV

#include "Test_Cloud.h"
#include "Test_SteadyStateMPI.h"
// #include "StatisticsMPI_MultiPhase.h"

template <typename TGrid, typename TStepper, template <typename> class TSlice=SliceMPI>
class Test_CloudMPI :
    public Test_Cloud<TGrid,TStepper,TSlice>,
    public Test_SteadyStateMPI<TGrid,TStepper,TSlice>
{
public:
    Test_CloudMPI(MPI_Comm comm, ArgumentParser& _parser) :
        Test_SteadyState<TGrid,TStepper,TSlice>(_parser),
        Test_Cloud<TGrid,TStepper,TSlice>(_parser),
        Test_SteadyStateMPI<TGrid,TStepper,TSlice>(comm, _parser)
    { }
    virtual ~Test_CloudMPI() {}


protected:
    typedef typename TGrid::BlockType TBlock;
    typedef typename TBlock::ElementType TElement;

    // TODO: (fabianw@mavt.ethz.ch; Tue 25 Oct 2016 11:27:03 AM CEST) The BGC
    // needs more refactoring.  For the time being, it is the same as it was
    // before (only in Test_CloudMPI)

    // case specific
    virtual void _post_save()
    {
#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
        BGC::bgc_save(this->restart_id, this->m_comm_world);
#endif
    }

    virtual void _post_restart()
    {
#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
        BGC::bgc_restart(this->restart_id, this->m_comm_world);
#endif
    }

    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////         TEST CLOUD MPI          ///////////////\n");
        printf("////////////////////////////////////////////////////////////\n");
        typedef typename TGrid::BlockType B;
        std::cout << "Domain size:   [" << this->BPDX*B::sizeX*this->XPESIZE;
        std::cout << " x " << this->BPDY*B::sizeY*this->YPESIZE;
        std::cout << " x " <<  this->BPDZ*B::sizeZ*this->ZPESIZE << "]" << std::endl;

        std::cout << "Domain extent: [" << Simulation_Environment::extents[0];
        std::cout << " x " << Simulation_Environment::extents[1];
        std::cout << " x " <<  Simulation_Environment::extents[2] << "]" << std::endl;
    }
    virtual void _setup_parameter()
    {
        Test_Cloud<TGrid,TStepper,TSlice>::_setup_parameter();
        Test_SteadyStateMPI<TGrid,TStepper,TSlice>::_setup_parameter();
    }

    virtual void _analysis()
    {
        if (this->ANALYSISPERIOD != 0 && this->step_id%this->ANALYSISPERIOD == 0 && !this->bRESTART)
        {
            this->profiler.push_start("ANALYSIS");
            // MultiPhaseStatistics::dumpStatistics(*this->grid, this->step_id, this->t, this->dt, this->parser);
            this->profiler.pop_stop();
        }
    }

    virtual Seed<shape> _make_shapes()
    {
        Seed<shape> base_seed = Test_Cloud<TGrid,TStepper,TSlice>::_make_shapes();

        // setting maxlayer-rank=0 is fine for all ic, except for tiwari-cloud
        // in this case, each rank as to see its closest bubbles
        // recommended value for this case is 3, but you should carefully check
        // your ic before subitting a production run
        Seed<shape> new_seed;
        const int maxlayer = this->parser("maxlayer-rank").asInt(0);

        int bbox = 1;
        int layerplus = 0;
        while (true)
        {
            double mystart[3], myextent[3];
            _get_rank_coords(mystart, myextent, bbox);

            new_seed = base_seed.retain_shapes(mystart,myextent);

            if (new_seed.get_shapes().size() == 0)
                bbox++;
            else
            {
                if (layerplus<maxlayer)
                {
                    bbox++;
                    layerplus++;
                }
                else
                    break;
            }
        }

        MPI_Barrier(this->m_comm_world);
        return new_seed;
    }

    virtual void _relax()
    {
        if(this->parser("-relax").asInt(0)>0)
        {
            if (this->isroot)
                std::cout << "relaxing pressure a little bit..." << std::endl;

            std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();
            const int Niter = this->parser("-relax").asInt(0);

            for(int i = 0; i <Niter; ++i)
            {
                if (this->isroot && i%100==0)
                    printf("gonna do iteration %i out of %i\n", i, Niter);

                // TODO: disabled pressure relaxation for the time being
                // evaluateOperator_MPI<LabMPI, GaussSeidel>(grid);

#pragma omp parallel for schedule(dynamic,1)
                for(int i=0; i<(int)vInfo.size(); i++)
                {
                    BlockInfo info = vInfo[i];
                    TBlock& b = *(TBlock*)info.ptrBlock;

                    for(int iz=0; iz<TBlock::sizeZ; iz++)
                        for(int iy=0; iy<TBlock::sizeY; iy++)
                            for(int ix=0; ix<TBlock::sizeX; ix++)
                                b(ix,iy,iz).energy = b(ix,iy,iz).dummy;
                }
            }

            if (this->isroot)
                std::cout << "done!"<< std::endl;
        }
    }


private:
    void _get_rank_coords(double mystart[3], double myextent[3], int bbox=1)
    {
        int peidx[3];
        this->grid->peindex(peidx);
        const double spacing = this->grid->getH()*_BLOCKSIZE_;

        const int BPD[3] = {this->BPDX, this->BPDY, this->BPDZ};

        for(int i=0; i<3; i++)
        {
            mystart[i] = (peidx[i]-bbox)*BPD[i]*spacing;
            myextent[i] = (1+2*bbox)*BPD[i]*spacing;
        }
    }
};

#endif /* TEST_CLOUDMPI_H_IBF8M2BV */
