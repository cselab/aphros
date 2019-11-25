/*
 *  MeshMap.h
 *  Cubism
 *
 *  Created by Fabian Wermelinger on 05/03/17.
 *  Copyright 2017 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#ifndef MESHMAP_H_UAYWTJDH
#define MESHMAP_H_UAYWTJDH

#include <cassert>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

class MeshDensity
{
public:
    const bool uniform;
    MeshDensity(const bool _uniform) : uniform(_uniform) {}

    virtual void compute_spacing(const double xS, const double xE, const size_t ncells, double* const ary,
            const size_t ghostS=0, const size_t ghostE=0, double* const ghost_spacing=NULL) const = 0;
};


class UniformDensity : public MeshDensity
{
public:
    UniformDensity() : MeshDensity(true) {}

    virtual void compute_spacing(const double xS, const double xE, const size_t ncells, double* const ary,
            const size_t ghostS=0, const size_t ghostE=0, double* const ghost_spacing=NULL) const
    {
        const double h = (xE - xS) / ncells;
        for (size_t i = 0; i < ncells; ++i)
            ary[i] = h;

        // ghost cells are given by ghost start (ghostS) and ghost end
        // (ghostE) and count the number of ghosts on either side (inclusive).
        // For example, for a symmetric 6-point stencil -> ghostS = 3 and
        // ghostE = 3.  ghost_spacing must provide valid memory for it.
        if (ghost_spacing)
            for (size_t i = 0; i < ghostS+ghostE; ++i)
                ghost_spacing[i] = h;
    }
};


template <typename TBlock>
class MeshMap
{
public:
    MeshMap(const double xS, const double xE, const size_t Nblocks,
            const size_t bs) :
        m_xS(xS), m_xE(xE), m_extent(xE-xS), m_Nblocks(Nblocks),
        m_Ncells(Nblocks*bs), bs(bs),
        m_uniform(true), m_initialized(false)
    {}

    ~MeshMap()
    {
        if (m_initialized)
        {
            delete[] m_grid_spacing;
            delete[] m_block_spacing;
        }
    }

    void init(const MeshDensity* const kernel, const size_t ghostS=0, const size_t ghostE=0, double* const ghost_spacing=NULL)
    {
        _alloc();

        kernel->compute_spacing(m_xS, m_xE, m_Ncells, m_grid_spacing, ghostS, ghostE, ghost_spacing);

        assert(m_Nblocks > 0);
        for (size_t i = 0; i < m_Nblocks; ++i)
        {
            double delta_block = 0.0;
            for (size_t j = 0; j < bs; ++j)
                delta_block += m_grid_spacing[i*bs + j];
            m_block_spacing[i] = delta_block;
        }

        m_uniform = kernel->uniform;
        m_initialized = true;
    }

    inline double start() const { return m_xS; }
    inline double end() const { return m_xE; }
    inline double extent() const { return m_extent; }
    inline size_t nblocks() const { return m_Nblocks; }
    inline size_t ncells() const { return m_Ncells; }
    inline bool uniform() const { return m_uniform; }

    inline double cell_width(const int ix) const
    {
        assert(m_initialized && ix >= 0 && ix < int(m_Ncells));
        return m_grid_spacing[ix];
    }

    inline double block_width(const int bix) const
    {
        assert(m_initialized && bix >= 0 && bix < int(m_Nblocks));
        return m_block_spacing[bix];
    }

    inline double block_origin(const int bix) const
    {
        assert(m_initialized && bix >= 0 && bix < int(m_Nblocks));
        double offset = m_xS;
        for (int i = 0; i < bix; ++i)
            offset += m_block_spacing[i];
        return offset;
    }

    inline double* get_grid_spacing(const int bix)
    {
        assert(m_initialized && bix >= 0 && bix < int(m_Nblocks));
        return &m_grid_spacing[bix*bs];
    }

    inline double* data_grid_spacing() { return m_grid_spacing; }

private:
    const double m_xS;
    const double m_xE;
    const double m_extent;
    const size_t m_Nblocks;
    const size_t m_Ncells;
    const size_t bs; // block size (current direction)

    bool m_uniform;
    bool m_initialized;
    double* m_grid_spacing;
    double* m_block_spacing;

    inline void _alloc()
    {
        m_grid_spacing = new double[m_Ncells];
        m_block_spacing= new double[m_Nblocks];
    }
};

#endif /* MESHMAP_H_UAYWTJDH */
