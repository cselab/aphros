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
#include "ArgumentParser.h"

class MeshDensity
{
public:
    const bool uniform;
    MeshDensity(const bool _uniform) : uniform(_uniform) {}

    virtual void compute_spacing(const double xS, const double xE, const unsigned int ncells, double* const ary,
            const unsigned int ghostS=0, const unsigned int ghostE=0, double* const ghost_spacing=NULL) const = 0;
};


class UniformDensity : public MeshDensity
{
public:
    UniformDensity() : MeshDensity(true) {}

    virtual void compute_spacing(const double xS, const double xE, const unsigned int ncells, double* const ary,
            const unsigned int ghostS=0, const unsigned int ghostE=0, double* const ghost_spacing=NULL) const
    {
        const double h = (xE - xS) / ncells;
        for (int i = 0; i < ncells; ++i)
            ary[i] = h;

        // ghost cells are given by ghost start (ghostS) and ghost end
        // (ghostE) and count the number of ghosts on either side (inclusive).
        // For example, for a symmetric 6-point stencil -> ghostS = 3 and
        // ghostE = 3.  ghost_spacing must provide valid memory for it.
        if (ghost_spacing)
            for (int i = 0; i < ghostS+ghostE; ++i)
                ghost_spacing[i] = h;
    }
};

class GaussianDensity : public MeshDensity
{
    const double A;
    const double B;

public:
    struct DefaultParameter
    {
        double A, B;
        DefaultParameter() : A(1.0), B(0.25) {}
    };

    GaussianDensity(const DefaultParameter d=DefaultParameter()) : MeshDensity(false), A(d.A), B(d.B) {}

    virtual void compute_spacing(const double xS, const double xE, const unsigned int ncells, double* const ary,
            const unsigned int ghostS=0, const unsigned int ghostE=0, double* const ghost_spacing=NULL) const
    {
        const unsigned int total_cells = ncells + ghostE + ghostS;
        double* const buf = new double[total_cells];

        const double y = 1.0/(B*(total_cells+1));
        double ducky = 0.0;
        for (int i = 0; i < total_cells; ++i)
        {
            const double x = i - (total_cells+1)*0.5;
            buf[i] = 1.0/(A*std::exp(-0.5*x*x*y*y) + 1.0);

            if (i >= ghostS && i < ncells + ghostS)
                ducky += buf[i];
        }

        const double scale = (xE-xS)/ducky;
        for (int i = 0; i < total_cells; ++i)
            buf[i] *= scale;

        for (int i = 0; i < ncells; ++i)
            ary[i] = buf[i+ghostS];

        if (ghost_spacing)
        {
            for (int i = 0; i < ghostS; ++i)
                ghost_spacing[i] = buf[i];
            for (int i = 0; i < ghostE; ++i)
                ghost_spacing[i+ghostS] = buf[i+ncells+ghostS];
        }

        // clean up
        delete[] buf;
    }
};

class SmoothHeavisideDensity : public MeshDensity
{
protected:
    const double A;
    const double B;
    const double c;
    const double eps;

    inline double _smooth_heaviside(const double r, const double A, const double B, const double c, const double eps, const bool mirror=false) const
    {
        const double theta = M_PI*std::max(0.0, std::min(1.0, (1.0-2.0*mirror)/eps*(r - (c-0.5*eps)) + 1.0*mirror));
        const double q = 0.5*(std::cos(theta) + 1.0);
        return B + (A-B)*q;
    }

public:
    struct DefaultParameter
    {
        double A, B, c1, eps1;
        DefaultParameter() : A(2.0), B(1.0), c1(0.5), eps1(0.2) {}
    };

    SmoothHeavisideDensity(const double A, const double B, const double c, const double eps) : MeshDensity(false), A(A), B(B), c(c), eps(eps) {}
    SmoothHeavisideDensity(const DefaultParameter d=DefaultParameter()) : MeshDensity(false), A(d.A), B(d.B), c(d.c1), eps(d.eps1) {}

    virtual void compute_spacing(const double xS, const double xE, const unsigned int ncells, double* const ary,
            const unsigned int ghostS=0, const unsigned int ghostE=0, double* const ghost_spacing=NULL) const
    {
        assert(std::min(A,B) > 0.0);
        assert(0.5*eps <= c && c <= 1.0-0.5*eps);

        const unsigned int total_cells = ncells + ghostE + ghostS;
        double* const buf = new double[total_cells];

        double ducky = 0.0;
        for (int i = 0; i < ncells; ++i)
        {
            const double r = (static_cast<double>(i)+0.5)/ncells;
            double rho;
            if (A >= B)
                rho = _smooth_heaviside(r, A, B, c, eps);
            else
                rho = _smooth_heaviside(r, B, A, c, eps, true);
            buf[i+ghostS] = 1.0/rho;
            ducky += buf[i+ghostS];
        }
        for (int i = 0; i < ghostS; ++i)
            buf[i] = buf[ghostS];
        for (int i = 0; i < ghostE; ++i)
            buf[i+ncells+ghostS] = buf[ncells+ghostS-1];

        const double scale = (xE-xS)/ducky;
        for (int i = 0; i < total_cells; ++i)
            buf[i] *= scale;

        for (int i = 0; i < ncells; ++i)
            ary[i] = buf[i+ghostS];

        if (ghost_spacing)
        {
            for (int i = 0; i < ghostS; ++i)
                ghost_spacing[i] = buf[i];
            for (int i = 0; i < ghostE; ++i)
                ghost_spacing[i+ghostS] = buf[i+ncells+ghostS];
        }

        // clean up
        delete[] buf;
    }
};

class SmoothHatDensity : public SmoothHeavisideDensity
{
protected:
    const double c2;
    const double eps2;

public:
    struct DefaultParameter
    {
        double A, B, c1, eps1, c2, eps2;
        DefaultParameter() : A(2.0), B(1.0), c1(0.25), eps1(0.15), c2(0.75), eps2(0.15) {}
    };

    SmoothHatDensity(const DefaultParameter d=DefaultParameter()) : SmoothHeavisideDensity(d.A,d.B,d.c1,d.eps1), c2(d.c2), eps2(d.eps2) {}

    virtual void compute_spacing(const double xS, const double xE, const unsigned int ncells, double* const ary,
            const unsigned int ghostS=0, const unsigned int ghostE=0, double* const ghost_spacing=NULL) const
    {
        assert(std::min(A,B) > 0.0);
        assert(0.5*eps <= c && c <= c2-0.5*eps2);
        assert(c2 <= 1.0-0.5*eps2);

        const unsigned int total_cells = ncells + ghostE + ghostS;
        double* const buf = new double[total_cells];

        double ducky = 0.0;
        for (int i = 0; i < ncells; ++i)
        {
            const double r = (static_cast<double>(i)+0.5)/ncells;
            double rho;
            if (r <= c+0.5*eps)
                if (A >= B)
                    rho = _smooth_heaviside(r, A, B, c, eps);
                else
                    rho = _smooth_heaviside(r, B, A, c, eps, true);
            else if (r >= c2-0.5*eps2)
                if (A >= B)
                    rho = _smooth_heaviside(r, A, B, c2, eps2, true);
                else
                    rho = _smooth_heaviside(r, B, A, c2, eps2);
            else
                rho = B;
            buf[i+ghostS] = 1.0/rho;
            ducky += buf[i+ghostS];
        }
        for (int i = 0; i < ghostS; ++i)
            buf[i] = buf[ghostS];
        for (int i = 0; i < ghostE; ++i)
            buf[i+ncells+ghostS] = buf[ncells+ghostS-1];

        const double scale = (xE-xS)/ducky;
        for (int i = 0; i < total_cells; ++i)
            buf[i] *= scale;

        for (int i = 0; i < ncells; ++i)
            ary[i] = buf[i+ghostS];

        if (ghost_spacing)
        {
            for (int i = 0; i < ghostS; ++i)
                ghost_spacing[i] = buf[i];
            for (int i = 0; i < ghostE; ++i)
                ghost_spacing[i+ghostS] = buf[i+ncells+ghostS];
        }

        // clean up
        delete[] buf;
    }
};


class BoxFadeDensity : public MeshDensity
{
protected:
    const double box_left;
    const double box_right;
    const double abs_tol;
    const int box_ncells;

public:
    struct DefaultParameter
    {
        double box_left, box_right;
        double abs_tol;
        int box_ncells;
        DefaultParameter() : box_left(0.0), box_right(1.0), abs_tol(1.0e-8),
        box_ncells(256) {}
    };

    BoxFadeDensity(const DefaultParameter d=DefaultParameter()) : MeshDensity(false),
    box_left(d.box_left), box_right(d.box_right), abs_tol(d.abs_tol),
    box_ncells(d.box_ncells) {}

    virtual void compute_spacing(const double xS, const double xE, const unsigned int ncells, double* const ary,
            const unsigned int ghostS=0, const unsigned int ghostE=0, double* const ghost_spacing=NULL) const
    {
        const double box_width = box_right - box_left;
        const double length_left = box_left - xS;
        const double length_right= xE - box_right;
        assert(box_width > 0.0);
        assert(length_left >= 0.0);
        assert(length_right>= 0.0);

        int cells_left = 0;
        int cells_right = 0;
        if (length_left > length_right)
        {
            const double distribution_ratio = length_right/length_left;
            cells_right = 0.5*(static_cast<int>(ncells) - box_ncells)*distribution_ratio;
            cells_left = static_cast<int>(ncells) - box_ncells - cells_right;
        }
        else
        {
            const double distribution_ratio = length_left/length_right;
            cells_left = 0.5*(static_cast<int>(ncells) - box_ncells)*distribution_ratio;
            cells_right = static_cast<int>(ncells) - box_ncells - cells_left;
        }
        assert(cells_left >= 0);
        assert(cells_right >= 0);

        const unsigned int total_cells = ncells + ghostE + ghostS;
        double* const buf = new double[total_cells];
        double* const buf_left  = &buf[ghostS];
        double* const buf_middle= &buf[ghostS+cells_left];
        double* const buf_right = &buf[ghostS+cells_left+box_ncells];

        const double h_window = box_width/box_ncells;

        _compute_fade(buf_left, cells_left, length_left, h_window, 6);
        for (int i = 0; i < box_ncells; ++i) buf_middle[i] = h_window;
        _compute_fade(buf_right, cells_right, length_right, h_window, 6);

        // flip left side
        for (int i = 0; i < cells_left/2; ++i)
        {
            const double tmp = buf_left[i];
            buf_left[i] = buf_left[cells_left-1-i];
            buf_left[cells_left-1-i] = tmp;
        }

        // assemble ghosts
        for (int i = 0; i < ghostS; ++i)
            buf[i] = buf[ghostS];
        for (int i = 0; i < ghostE; ++i)
            buf[i+ncells+ghostS] = buf[ghostS+ncells-1];

        // distribute
        for (int i = 0; i < ncells; ++i)
            ary[i] = buf[i+ghostS];

        if (ghost_spacing)
        {
            for (int i = 0; i < ghostS; ++i)
                ghost_spacing[i] = buf[i];
            for (int i = 0; i < ghostE; ++i)
                ghost_spacing[i+ghostS] = buf[i+ncells+ghostS];
        }

        // clean up
        delete[] buf;
    }

private:

    // helpers
    inline double _heaviside(const double x) const { return (x>0.0) ? 1.0 : 0.0; }
    inline double _theta(const double xi) const { return M_PI*(xi*_heaviside(xi) - (xi-1.0)*_heaviside(xi-1.0)); }
    inline double _f(const double xi, const double alpha) const { return 0.5*alpha*(1.0 - std::cos(_theta(xi))) + 1.0; }
    inline double _dfdalpha(const double xi) const { return 0.5*(1.0 - std::cos(_theta(xi))); }

    inline double _F(const std::vector<double>& xi, const double alpha, const double L) const
    {
        double sum = 0.0;
        for (int i = 0; i < (int)xi.size(); ++i)
            sum += _f(xi[i], alpha);
        return L - sum;
    }

    inline double _dFdalpha(const std::vector<double>& xi) const
    {
        double sum = 0.0;
        for (int i = 0; i < (int)xi.size(); ++i)
            sum += _dfdalpha(xi[i]);
        return -sum;
    }

    double _newton(const std::vector<double>& xi, const double L, const double tol=1.0e-8) const
    {
        double alpha = 0.0; // initial guess
        unsigned int steps = 0;
        while (true)
        {
            const double F = _F(xi,alpha,L);
            const double dF= _dFdalpha(xi);
            const double tmp = alpha;
            alpha = alpha - F/dF;
            ++steps;
            if (std::abs(alpha-tmp) <= tol)
                break;
        }
        // std::cout << "MeshMap: Newton solver: Found solution in " << steps << " steps" << std::endl;
        return alpha;
    }

    void _compute_fade(double* const buf, const unsigned int N, const double L, const double h0, const int nLast=6) const
    {
        if (N > 0)
        {
            assert(N > 20); // not a strict limit but should prefereably be at least this
            const double h_xi = 1.0/N;
            std::vector<double> xi(N);
            for (unsigned int i = 0; i < N; ++i)
                xi[i] = h_xi*(i+0.5);

            // find growth factor
            const double alpha = _newton(xi, L/h0, abs_tol);

            // adjust the last 6 cells to be of equal size (for FD scheme used in
            // characteristic 1D boundaries, otherwise need non-uniform scheme
            // there)
            assert(nLast < N);
            for (unsigned int i = 0; i < N; ++i)
                buf[i] = h0*_f(xi[i], alpha);

            const unsigned int k = N-nLast;
            for (unsigned int i = k; i < N; ++i)
                buf[i] = buf[k-1];

            double sum = 0.0;
            for (unsigned int i = 0; i < N; ++i)
                sum += buf[i];

            const double scale = L/sum;
            for (unsigned int i = 0; i < N; ++i)
                buf[i] *= scale;
        }
    }
};



class MeshDensityFactory
{
public:
    MeshDensityFactory(ArgumentParser& parser) : m_parser(parser) { _make_mesh_kernels(); }
    ~MeshDensityFactory() { _dealloc(); }

    inline MeshDensity* get_mesh_kernel(const int i) { return m_mesh_kernels[i]; }

private:
    ArgumentParser& m_parser;
    std::vector<MeshDensity*> m_mesh_kernels;

    void _dealloc()
    {
        for (int i = 0; i < (int)m_mesh_kernels.size(); ++i)
            delete m_mesh_kernels[i];
    }

    void _make_mesh_kernels()
    {
        std::vector<std::string> suffix;
        suffix.push_back("_x");
        suffix.push_back("_y");
        suffix.push_back("_z");
        for (int i = 0; i < (int)suffix.size(); ++i)
        {
            const std::string mesh_density("mesh_density" + suffix[i]);
            if (m_parser.exist(mesh_density))
            {
                const std::string density_function(m_parser(mesh_density).asString());
                const std::string A("A" + suffix[i]);
                const std::string B("B" + suffix[i]);
                const std::string c1("c1" + suffix[i]);
                const std::string eps1("eps1" + suffix[i]);
                const std::string c2("c2" + suffix[i]);
                const std::string eps2("eps2" + suffix[i]);
                if (density_function == "GaussianDensity")
                {
                    typename GaussianDensity::DefaultParameter p;
                    if (m_parser.exist(A)) p.A = m_parser(A).asDouble();
                    if (m_parser.exist(B)) p.B = m_parser(B).asDouble();
                    m_mesh_kernels.push_back(new GaussianDensity(p));
                }
                else if (density_function == "SmoothHeavisideDensity")
                {
                    typename SmoothHeavisideDensity::DefaultParameter p;
                    if (m_parser.exist(A)) p.A = m_parser(A).asDouble();
                    if (m_parser.exist(B)) p.B = m_parser(B).asDouble();
                    if (m_parser.exist(c1)) p.c1 = m_parser(c1).asDouble();
                    if (m_parser.exist(eps1)) p.eps1 = m_parser(eps1).asDouble();
                    m_mesh_kernels.push_back(new SmoothHeavisideDensity(p));
                }
                else if (density_function == "SmoothHatDensity")
                {
                    typename SmoothHatDensity::DefaultParameter p;
                    if (m_parser.exist(A)) p.A = m_parser(A).asDouble();
                    if (m_parser.exist(B)) p.B = m_parser(B).asDouble();
                    if (m_parser.exist(c1)) p.c1 = m_parser(c1).asDouble();
                    if (m_parser.exist(eps1)) p.eps1 = m_parser(eps1).asDouble();
                    if (m_parser.exist(c2)) p.c2 = m_parser(c2).asDouble();
                    if (m_parser.exist(eps2)) p.eps2 = m_parser(eps2).asDouble();
                    m_mesh_kernels.push_back(new SmoothHatDensity(p));
                }
                else if (density_function == "BoxFadeDensity")
                {
                    typename BoxFadeDensity::DefaultParameter p;
                    const std::string box_left("box_left" + suffix[i]);
                    const std::string box_right("box_right" + suffix[i]);
                    const std::string box_ncells("box_ncells" + suffix[i]);
                    const std::string abs_tol("abs_tol" + suffix[i]);
                    if (m_parser.exist(box_left)) p.box_left = m_parser(box_left).asDouble();
                    if (m_parser.exist(box_right)) p.box_right = m_parser(box_right).asDouble();
                    if (m_parser.exist(box_ncells)) p.box_ncells= m_parser(box_ncells).asDouble();
                    if (m_parser.exist(abs_tol)) p.abs_tol = m_parser(abs_tol).asDouble();
                    m_mesh_kernels.push_back(new BoxFadeDensity(p));
                }
                else
                {
                    std::cerr << "ERROR: MeshMap.h: Undefined mesh density function." << std::endl;
                    abort();
                }
            }
            else
                m_mesh_kernels.push_back(new UniformDensity);
        }
    }
};


template <typename TBlock>
class MeshMap
{
public:
    MeshMap(const double xS, const double xE, const unsigned int Nblocks) :
        m_xS(xS), m_xE(xE), m_extent(xE-xS), m_Nblocks(Nblocks),
        m_Ncells(Nblocks*TBlock::sizeX), // assumes uniform cells in all directions!
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

    void init(const MeshDensity* const kernel, const unsigned int ghostS=0, const unsigned int ghostE=0, double* const ghost_spacing=NULL)
    {
        _alloc();

        kernel->compute_spacing(m_xS, m_xE, m_Ncells, m_grid_spacing, ghostS, ghostE, ghost_spacing);

        assert(m_Nblocks > 0);
        for (int i = 0; i < m_Nblocks; ++i)
        {
            double delta_block = 0.0;
            for (int j = 0; j < TBlock::sizeX; ++j)
                delta_block += m_grid_spacing[i*TBlock::sizeX + j];
            m_block_spacing[i] = delta_block;
        }

        m_uniform = kernel->uniform;
        m_initialized = true;
    }

    inline double start() const { return m_xS; }
    inline double end() const { return m_xE; }
    inline double extent() const { return m_extent; }
    inline unsigned int nblocks() const { return m_Nblocks; }
    inline unsigned int ncells() const { return m_Ncells; }
    inline bool uniform() const { return m_uniform; }

    inline double cell_width(const int ix) const
    {
        assert(m_initialized && ix >= 0 && ix < m_Ncells);
        return m_grid_spacing[ix];
    }

    inline double block_width(const int bix) const
    {
        assert(m_initialized && bix >= 0 && bix < m_Nblocks);
        return m_block_spacing[bix];
    }

    inline double block_origin(const int bix) const
    {
        assert(m_initialized && bix >= 0 && bix < m_Nblocks);
        double offset = m_xS;
        for (int i = 0; i < bix; ++i)
            offset += m_block_spacing[i];
        return offset;
    }

    inline double* get_grid_spacing(const int bix)
    {
        assert(m_initialized && bix >= 0 && bix < m_Nblocks);
        return &m_grid_spacing[bix*TBlock::sizeX];
    }

    inline double* data_grid_spacing() { return m_grid_spacing; }

private:
    const double m_xS;
    const double m_xE;
    const double m_extent;
    const unsigned int m_Nblocks;
    const unsigned int m_Ncells;

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
