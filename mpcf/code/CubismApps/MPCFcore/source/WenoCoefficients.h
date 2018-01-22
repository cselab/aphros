/*
 *  WenoCoefficients.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 05/07/2017
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef WENOCOEFFICIENTS_H_QX69ZBST
#define WENOCOEFFICIENTS_H_QX69ZBST

#include <iostream>
#include <cassert>
#include <cstdlib>
#include "common.h"
#include "MeshMap.h"

// per block
struct Coefficients_t
{
    // polynomials
    Real c[9][_BLOCKSIZE_+1];

    // ideal weights
    Real d[3][_BLOCKSIZE_+1];

    // smoothness indicators
    Real b[9][_BLOCKSIZE_+1];
};


// low level
struct Polynomial_coeffs
{
    Real* c00;
    Real* c01;
    Real* c02;

    Real* c10;
    Real* c11;
    Real* c12;

    Real* c20;
    Real* c21;
    Real* c22;
};

struct Ideal_weights
{
    Real* d0;
    Real* d1;
    Real* d2;
};

struct Smoothness_coeffs
{
    Real* b00;
    Real* b01;
    Real* b02;

    Real* b10;
    Real* b11;
    Real* b12;

    Real* b20;
    Real* b21;
    Real* b22;
};

struct WenoCoefficients_t
{
    Smoothness_coeffs b;
    Polynomial_coeffs c;
    Ideal_weights d;
};

class WenoCoefficients
{
public:
    WenoCoefficients(const unsigned int N) : m_initialized(false), m_N(N) {}
    virtual ~WenoCoefficients() {}

    virtual void setup(const double* const grid_spacing, const unsigned int ncells,
            const unsigned int ghostS=0, const unsigned int ghostE=0) = 0;

    inline WenoCoefficients_t coefficients_minus() const { return m_minus; }
    inline WenoCoefficients_t coefficients_plus() const { return m_plus; }

protected:
    bool m_initialized;
    const unsigned int m_N;
    WenoCoefficients_t m_minus;
    WenoCoefficients_t m_plus;
    Smoothness_coeffs m_b;

    void _alloc()
    {
        // MINUS
        posix_memalign((void**)&m_minus.c.c00, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_minus.c.c01, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_minus.c.c02, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_minus.c.c10, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_minus.c.c11, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_minus.c.c12, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_minus.c.c20, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_minus.c.c21, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_minus.c.c22, _ALIGNBYTES_, m_N*sizeof(Real));

        posix_memalign((void**)&m_minus.d.d0, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_minus.d.d1, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_minus.d.d2, _ALIGNBYTES_, m_N*sizeof(Real));

        // PLUS
        posix_memalign((void**)&m_plus.c.c00, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_plus.c.c01, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_plus.c.c02, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_plus.c.c10, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_plus.c.c11, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_plus.c.c12, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_plus.c.c20, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_plus.c.c21, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_plus.c.c22, _ALIGNBYTES_, m_N*sizeof(Real));

        posix_memalign((void**)&m_plus.d.d0, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_plus.d.d1, _ALIGNBYTES_, m_N*sizeof(Real));
        posix_memalign((void**)&m_plus.d.d2, _ALIGNBYTES_, m_N*sizeof(Real));

        // SMOOTHNESS INDICATORS
        posix_memalign((void**)&m_b.b00, _ALIGNBYTES_, (m_N+1)*sizeof(Real));
        posix_memalign((void**)&m_b.b01, _ALIGNBYTES_, (m_N+1)*sizeof(Real));
        posix_memalign((void**)&m_b.b02, _ALIGNBYTES_, (m_N+1)*sizeof(Real));
        posix_memalign((void**)&m_b.b10, _ALIGNBYTES_, (m_N+1)*sizeof(Real));
        posix_memalign((void**)&m_b.b11, _ALIGNBYTES_, (m_N+1)*sizeof(Real));
        posix_memalign((void**)&m_b.b12, _ALIGNBYTES_, (m_N+1)*sizeof(Real));
        posix_memalign((void**)&m_b.b20, _ALIGNBYTES_, (m_N+1)*sizeof(Real));
        posix_memalign((void**)&m_b.b21, _ALIGNBYTES_, (m_N+1)*sizeof(Real));
        posix_memalign((void**)&m_b.b22, _ALIGNBYTES_, (m_N+1)*sizeof(Real));

        m_initialized = true;
    }

    void _dealloc()
    {
        // MINUS
        free(m_minus.c.c00);
        free(m_minus.c.c01);
        free(m_minus.c.c02);
        free(m_minus.c.c10);
        free(m_minus.c.c11);
        free(m_minus.c.c12);
        free(m_minus.c.c20);
        free(m_minus.c.c21);
        free(m_minus.c.c22);

        free(m_minus.d.d0);
        free(m_minus.d.d1);
        free(m_minus.d.d2);

        // PLUS
        free(m_plus.c.c00);
        free(m_plus.c.c01);
        free(m_plus.c.c02);
        free(m_plus.c.c10);
        free(m_plus.c.c11);
        free(m_plus.c.c12);
        free(m_plus.c.c20);
        free(m_plus.c.c21);
        free(m_plus.c.c22);

        free(m_plus.d.d0);
        free(m_plus.d.d1);
        free(m_plus.d.d2);

        // SMOOTHNESS INDICATORS
        free(m_b.b00);
        free(m_b.b01);
        free(m_b.b02);
        free(m_b.b10);
        free(m_b.b11);
        free(m_b.b12);
        free(m_b.b20);
        free(m_b.b21);
        free(m_b.b22);
    }
};


///////////////////////////////////////////////////////////////////////////////
// Weno5 coefficients based on Coralic & Colonius, JCP 274 (2014) 95-121
class Weno5Coefficients_Coralic : public WenoCoefficients
{
public:
    Weno5Coefficients_Coralic(const unsigned int N) : WenoCoefficients(N) {}
    virtual ~Weno5Coefficients_Coralic() { if (m_initialized) _dealloc(); }

    static const unsigned int HALO_S = 3;
    static const unsigned int HALO_E = 3;

    virtual void setup(const double* const grid_spacing, const unsigned int ncells,
            const unsigned int ghostS=0, const unsigned int ghostE=0)
    {
        _alloc();

        // compute coefficients
        assert(ncells+1 == m_N);
        assert(ghostS == HALO_S && ghostE == HALO_E);

        // minus: coefficients correspond to faces located at x_{i - 1/2}
        // plus:  coefficients correspond to faces located at x_{i + 1/2}
        //
        // NOTE: The WENO kernels in the solver are implemented using the
        // notation  -|+, where | corresponds to a face located at x_{i + 1/2}.
        // Therefore, when the reconstruction is carried out for -|, this
        // corresponds to the PLUS coefficients in cell i, whereas |+
        // corresponds to the MINUS coefficients in cell i+1.  The coefficients
        // computed below are stored for all faces n_N = ncells+1 along some
        // direction with respect to the cell i, that is, +|-.

        for (int i = 0; i < (int)m_N; ++i)
        {
            // minus
            ///////////////////////////////////////////////////////////////////
            {
                const double* const si = grid_spacing + ghostS + i;
                const double hm2 = *(si-2);
                const double hm1 = *(si-1);
                const double h00 = *si;
                const double hp1 = *(si+1);
                const double hp2 = *(si+2);

                m_minus.c.c00[i] = 1.0;
                m_minus.c.c01[i] = -(h00*(h00+hp1 + h00+hp1+hp2)) / ((h00+hp1)*(h00+hp1+hp2));
                m_minus.c.c02[i] =  (h00*(h00+hp1)) / ((h00+hp1+hp2)*(hp1+hp2));
                m_minus.c.c10[i] = 1.0;
                m_minus.c.c11[i] = -(h00*(h00+hp1)) / ((hm1+h00)*(hm1+h00+hp1));
                m_minus.c.c12[i] = -(hm1*h00) / ((hm1+h00+hp1)*(h00+hp1));
                m_minus.c.c20[i] = 1.0;
                m_minus.c.c21[i] =  (hm1*h00) / ((hm2+hm1)*(hm2+hm1+h00));
                m_minus.c.c22[i] = -(h00*(hm2+hm1 + hm1+h00)) / ((hm2+hm1+h00)*(hm1+h00));

                m_minus.d.d0[i] = ((hm2+hm1)*hm1) / ((hm2+hm1+h00+hp1+hp2)*(hm1+h00+hp1+hp2));
                m_minus.d.d2[i] = ((h00+hp1)*(h00+hp1+hp2)) / ((hm2+hm1+h00+hp1)*(hm2+hm1+h00+hp1+hp2));
                m_minus.d.d1[i] = 1.0 - m_minus.d.d0[i] - m_minus.d.d2[i];
            }

            // plus
            ///////////////////////////////////////////////////////////////////
            {
                const double* const si = grid_spacing + ghostS + (i-1);
                const double hm2 = *(si-2);
                const double hm1 = *(si-1);
                const double h00 = *si;
                const double hp1 = *(si+1);
                const double hp2 = *(si+2);

                m_plus.c.c00[i] =  1.0;
                m_plus.c.c01[i] =  (h00*(h00+hp1 + hp1+hp2)) / ((h00+hp1)*(h00+hp1+hp2));
                m_plus.c.c02[i] = -(h00*hp1) / ((h00+hp1+hp2)*(hp1+hp2));
                m_plus.c.c10[i] =  1.0;
                m_plus.c.c11[i] =  (h00*hp1) / ((hm1+h00)*(hm1+h00+hp1));
                m_plus.c.c12[i] =  ((hm1+h00)*h00) / ((hm1+h00+hp1)*(h00+hp1));
                m_plus.c.c20[i] =  1.0;
                m_plus.c.c21[i] = -((hm1+h00)*h00) / ((hm2+hm1)*(hm2+hm1+h00));
                m_plus.c.c22[i] =  (h00*(hm2+hm1+h00 + hm1+h00)) / ((hm2+hm1+h00)*(hm1+h00));

                m_plus.d.d0[i] = ((hm2+hm1+h00)*(hm1+h00)) / ((hm2+hm1+h00+hp1+hp2)*(hm1+h00+hp1+hp2));
                m_plus.d.d2[i] = (hp1*(hp1+hp2)) / ((hm2+hm1+h00+hp1)*(hm2+hm1+h00+hp1+hp2));
                m_plus.d.d1[i] = 1.0 - m_plus.d.d0[i] - m_plus.d.d2[i];
            }
        }

        // smoothness indicator coefficients
        for (int i = 0; i < (int)(m_N+1); ++i)
        {
            const double* const si = grid_spacing + ghostS + (i-1);
            const double hm2 = *(si-2);
            const double hm1 = *(si-1);
            const double h00 = *si;
            const double hp1 = *(si+1);
            const double hp2 = *(si+2);

            const double fac = 4.0*h00*h00;
            {
                const double c0 = h00+hp1 + hp1+hp2;
                const double c1 = (h00+hp1)*(h00+hp1+hp2);
                m_b.b00[i] =  fac*( (10.0*h00*h00 + h00*(h00+hp1 + hp1+hp2) + c0*c0) / (c1*c1) );
            }
            {
                const double c0 = h00+hp1+hp2;
                m_b.b01[i] = -fac*( (19.0*h00*h00 - h00*(hp1+hp2) + 2.0*(h00+hp1)*(h00+hp1 + hp1+hp2)) / ((h00+hp1)*(hp1+hp2)*c0*c0) );
            }
            {
                const double c0 = (h00+hp1+hp2)*(hp1+hp2);
                m_b.b02[i] =  fac*( (10.0*h00*h00 + h00*hp1 + hp1*hp1) / (c0*c0) );
            }

            {
                const double c0 = hm1+h00+hp1;
                m_b.b10[i] = -fac*( (h00*(hm1 + 20.0*h00) - (h00+hp1)*(2.0*hm1 + h00)) / ((hm1+h00)*(h00+hp1)*c0*c0) );
            }
            {
                const double c0 = (hm1+h00)*(hm1+h00+hp1);
                m_b.b11[i] =  fac*( (10.0*h00*h00 + h00*hp1 + hp1*hp1) / (c0*c0) );
            }
            {
                const double c0 = (hm1+h00+hp1)*(h00+hp1);
                m_b.b12[i] =  fac*( (10.0*h00*h00 + hm1*h00 + hm1*hm1) / (c0*c0) );
            }

            {
                const double c0 = hm2+hm1 + hm1;
                const double c1 = (hm2+hm1+h00)*(hm1+h00);
                m_b.b20[i] =  fac*( (12.0*h00*h00 + 3.0*h00*(hm2+hm1 + hm1) + c0*c0) / (c1*c1) );
            }
            {
                const double c0 = hm2+hm1+h00;
                m_b.b21[i] = -fac*( (19.0*h00*h00 - (hm2+hm1)*h00 + 2.0*(hm1+h00)*(hm2+hm1 + hm1+h00)) / ((hm2+hm1)*(hm1+h00)*c0*c0) );
            }
            {
                const double c0 = (hm2+hm1)*(hm2+hm1+h00);
                m_b.b22[i] =  fac*( (10.0*h00*h00 + hm1*h00 + hm1*hm1) / (c0*c0) );
            }
        }

        m_minus.b.b00 = &m_b.b00[1];
        m_minus.b.b01 = &m_b.b01[1];
        m_minus.b.b02 = &m_b.b02[1];
        m_minus.b.b10 = &m_b.b10[1];
        m_minus.b.b11 = &m_b.b11[1];
        m_minus.b.b12 = &m_b.b12[1];
        m_minus.b.b20 = &m_b.b20[1];
        m_minus.b.b21 = &m_b.b21[1];
        m_minus.b.b22 = &m_b.b22[1];

        m_plus.b.b00 = &m_b.b00[0];
        m_plus.b.b01 = &m_b.b01[0];
        m_plus.b.b02 = &m_b.b02[0];
        m_plus.b.b10 = &m_b.b10[0];
        m_plus.b.b11 = &m_b.b11[0];
        m_plus.b.b12 = &m_b.b12[0];
        m_plus.b.b20 = &m_b.b20[0];
        m_plus.b.b21 = &m_b.b21[0];
        m_plus.b.b22 = &m_b.b22[0];
    }
};
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Weno3 coefficients: Derived by hand Wed Jun 21 2017 05:16:33 PM (-0700)
//                     (Fabian notebook)
class Weno3Coefficients : public WenoCoefficients
{
public:
    Weno3Coefficients(const unsigned int N) : WenoCoefficients(N) {}
    virtual ~Weno3Coefficients() { if (m_initialized) _dealloc(); }

    static const unsigned int HALO_S = 2;
    static const unsigned int HALO_E = 2;

    virtual void setup(const double* const grid_spacing, const unsigned int ncells,
            const unsigned int ghostS=0, const unsigned int ghostE=0)
    {
        _alloc();

        // compute coefficients
        assert(ncells+1 == m_N);
        assert(ghostS == HALO_S && ghostE == HALO_E);

        // minus: coefficients correspond to faces located at x_{i - 1/2}
        // plus:  coefficients correspond to faces located at x_{i + 1/2}
        //
        // NOTE: The WENO kernels in the solver are implemented using the
        // notation  -|+, where | corresponds to a face located at x_{i + 1/2}.
        // Therefore, when the reconstruction is carried out for -|, this
        // corresponds to the PLUS coefficients in cell i, whereas |+
        // corresponds to the MINUS coefficients in cell i+1.  The coefficients
        // computed below are stored for all faces n_N = ncells+1 along some
        // direction with respect to the cell i, that is, +|-.

        for (int i = 0; i < (int)m_N; ++i)
        {
            // minus
            ///////////////////////////////////////////////////////////////////
            {
                const double* const si = grid_spacing + ghostS + i;
                const double hm1 = *(si-1);
                const double h00 = *si;
                const double hp1 = *(si+1);

                m_minus.c.c00[i] = (h00+hp1)/hp1 - h00*h00/((h00+hp1)*hp1);
                m_minus.c.c01[i] = -h00/(h00+hp1);
                m_minus.c.c10[i] = 1.0 - hm1/h00 + hm1*hm1/((hm1+h00)*h00);
                m_minus.c.c11[i] = hm1/(hm1+h00);

                m_minus.d.d0[i] = hm1/(hm1+h00+hp1);
                m_minus.d.d1[i] = 1.0 - m_minus.d.d0[i];
            }

            // plus
            ///////////////////////////////////////////////////////////////////
            {
                const double* const si = grid_spacing + ghostS + (i-1);
                const double hm1 = *(si-1);
                const double h00 = *si;
                const double hp1 = *(si+1);

                m_plus.c.c00[i] = 1.0 - h00/hp1 + h00*h00/((h00+hp1)*hp1);
                m_plus.c.c01[i] = h00/(h00+hp1);
                m_plus.c.c10[i] = hm1/(hm1+h00) - 1.0;
                m_plus.c.c11[i] = h00/(hm1+h00) + 1.0;

                m_plus.d.d0[i] = (hm1+h00)/(hm1+h00+hp1);
                m_plus.d.d1[i] = 1.0 - m_plus.d.d0[i];
            }

            m_minus.c.c02[i] = 0.0;
            m_minus.c.c12[i] = 0.0;
            m_minus.c.c20[i] = 0.0;
            m_minus.c.c21[i] = 0.0;
            m_minus.c.c22[i] = 0.0;
            m_minus.d.d2[i] = 0.0;

            m_plus.c.c02[i] = 0.0;
            m_plus.c.c12[i] = 0.0;
            m_plus.c.c20[i] = 0.0;
            m_plus.c.c21[i] = 0.0;
            m_plus.c.c22[i] = 0.0;
            m_plus.d.d2[i] = 0.0;
        }

        // smoothness indicator coefficients
        for (int i = 0; i < (int)(m_N+1); ++i)
        {
            const double* const si = grid_spacing + ghostS + (i-1);
            const double hm1 = *(si-1);
            const double h00 = *si;
            const double hp1 = *(si+1);

            const double fac = 4.0*h00*h00;
            {
                const double c0 = (h00/(h00+hp1) - 1.0)/hp1;
                const double c1 = 1.0/(h00+hp1);
                m_b.b00[i] = fac*c0*c0;
                m_b.b01[i] = fac*2.0*c0*c1;
                m_b.b02[i] = fac*c1*c1;
            }

            {
                const double c0 = (hm1/(hm1+h00) - 1.0)/h00;
                const double c1 = 1.0/(hm1+h00);
                m_b.b10[i] = fac*c0*c0;
                m_b.b11[i] = fac*2.0*c0*c1;
                m_b.b12[i] = fac*c1*c1;
            }

            m_b.b20[i] = 0.0;
            m_b.b21[i] = 0.0;
            m_b.b22[i] = 0.0;
        }

        m_minus.b.b00 = &m_b.b00[1];
        m_minus.b.b01 = &m_b.b01[1];
        m_minus.b.b02 = &m_b.b02[1];
        m_minus.b.b10 = &m_b.b10[1];
        m_minus.b.b11 = &m_b.b11[1];
        m_minus.b.b12 = &m_b.b12[1];
        m_minus.b.b20 = &m_b.b20[1];
        m_minus.b.b21 = &m_b.b21[1];
        m_minus.b.b22 = &m_b.b22[1];

        m_plus.b.b00 = &m_b.b00[0];
        m_plus.b.b01 = &m_b.b01[0];
        m_plus.b.b02 = &m_b.b02[0];
        m_plus.b.b10 = &m_b.b10[0];
        m_plus.b.b11 = &m_b.b11[0];
        m_plus.b.b12 = &m_b.b12[0];
        m_plus.b.b20 = &m_b.b20[0];
        m_plus.b.b21 = &m_b.b21[0];
        m_plus.b.b22 = &m_b.b22[0];
    }
};
///////////////////////////////////////////////////////////////////////////////

#endif /* WENOCOEFFICIENTS_H_QX69ZBST */
