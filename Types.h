// File       : Types.h
// Date       : Fri 01 Apr 2016 05:54:03 PM CEST
// Author     : Fabian Wermelinger
// Description: Some required types
// Copyright 2016 ETH Zurich. All Rights Reserved.
#ifndef TYPES_H_NTYWLJPY
#define TYPES_H_NTYWLJPY

#include "ArgumentParser.h"
#include "GridMPI.h"
#include "BlockLab.h"
#include "BlockLabMPI.h"
#include "BlockInfo.h"

#include "BoundaryConditions.h"

#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"

#include <memory>

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

using Scal = Real;

using Mesh = geom::geom3d::MeshStructured<Scal>;
using Vect = typename Mesh::Vect;
using MIdx = typename Mesh::MIdx;
using Rect = geom::Rect<Vect>;
using IdxCell = geom::IdxCell;
using IdxFace = geom::IdxFace;
using IdxNode = geom::IdxNode;

template <class T>
using FieldCell = geom::FieldCell<T>;
template <class T>
using FieldFace = geom::FieldFace<T>;
template <class T>
using FieldNode = geom::FieldNode<T>;


struct FluidElement
{
    Scal p;
    Vect v;
    Scal volume;
    FluidElement() = default;
    void clear() { p = 0.; v = Vect(0); volume = 0.; }
};

struct FluidBlock
{
    static const int sizeX = _BLOCKSIZEX_;
    static const int sizeY = _BLOCKSIZEY_;
    static const int sizeZ = _BLOCKSIZEZ_;

    static const int gptfloats = sizeof(FluidElement)/sizeof(Real);

    typedef FluidElement ElementType;
    typedef FluidElement element_type;

    //Mesh mesh;
    std::unique_ptr<Mesh> mesh;


    FluidElement __attribute__((__aligned__(_ALIGNBYTES_))) data[_BLOCKSIZEZ_][_BLOCKSIZEY_][_BLOCKSIZEX_];

    Real __attribute__((__aligned__(_ALIGNBYTES_))) tmp[_BLOCKSIZEZ_][_BLOCKSIZEY_][_BLOCKSIZEX_][gptfloats];

    void clear_data()
    {
        const int N = sizeX*sizeY*sizeZ;
        FluidElement * const e = &data[0][0][0];
        for(int i=0; i<N; ++i) e[i].clear();
    }

    void clear_tmp()
    {
        const int N = sizeX * sizeY * sizeZ * gptfloats;

        Real * const e = &tmp[0][0][0][0];
        for(int i=0; i<N; ++i) e[i] = 0;
    }

    void clear()
    {
        clear_data();
        clear_tmp();
    }

    inline FluidElement& operator()(int ix, int iy=0, int iz=0)
    {
        assert(ix>=0 && ix<sizeX);
        assert(iy>=0 && iy<sizeY);
        assert(iz>=0 && iz<sizeZ);

        return data[iz][iy][ix];
    }

    template <typename Streamer>
    inline void Write(ofstream& output, Streamer streamer) const
    {
        for(int iz=0; iz<sizeZ; iz++)
            for(int iy=0; iy<sizeY; iy++)
                for(int ix=0; ix<sizeX; ix++)
                    streamer.operate(data[iz][iy][ix], output);
    }

    template <typename Streamer>
    inline void Read(ifstream& input, Streamer streamer)
    {
        for(int iz=0; iz<sizeZ; iz++)
            for(int iy=0; iy<sizeY; iy++)
                for(int ix=0; ix<sizeX; ix++)
                    streamer.operate(input, data[iz][iy][ix]);
    }
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class MyLabAbsorbing: public BlockLab<BlockType,allocator>
{
    typedef typename BlockType::ElementType ElementTypeBlock;

public:
    MyLabAbsorbing(): BlockLab<BlockType,allocator>(){}

    virtual inline std::string name() const { return "MyTestLab"; }
    bool is_xperiodic() {return false;}
    bool is_yperiodic() {return false;}
    bool is_zperiodic() {return false;}

    void _apply_bc(const BlockInfo& info, const Real t=0)
    {
        // "this->" references attributes inherited from BlockLab
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> 
				bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);

        if (info.index[0]==0)           bc.template applyBC_absorbing<0,0>();
        if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing<0,1>();
        if (info.index[1]==0)           bc.template applyBC_absorbing<1,0>();
        if (info.index[1]==this->NY-1)  bc.template applyBC_absorbing<1,1>();
        if (info.index[2]==0)           bc.template applyBC_absorbing<2,0>();
        if (info.index[2]==this->NZ-1)  bc.template applyBC_absorbing<2,1>();
    }
};

typedef Grid<FluidBlock, std::allocator> NodeGrid;
typedef GridMPI<NodeGrid> MPIGrid;
typedef MyLabAbsorbing<FluidBlock> ALab; // absorbing BC (Lab)
typedef BlockLab<FluidBlock> PLab;       // periodic BC (Lab)
typedef BlockLabMPI<ALab> AMPILab;
typedef BlockLabMPI<PLab> PMPILab;

#endif /* TYPES_H_NTYWLJPY */
