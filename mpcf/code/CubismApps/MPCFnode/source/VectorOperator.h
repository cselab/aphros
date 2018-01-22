/*
 *  VectorOperator.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 03/18/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef VECTOROPERATOR_H_UBJ0NRPP
#define VECTOROPERATOR_H_UBJ0NRPP

#include <string>
#include "common.h"
#include "Types.h"
#include "BlockInfo.h"
#include "StencilInfo.h"

// Operator classes:
// Some of the streamers depend on a specific operator, which must be evaluated
// in advance and writes the result into a specific location.  When evaluating
// a set of operators, they must be evaluated in _increasing_ class order.  The
// following classes for operators with writes are defined:
// class 0:   no operator (default)
// class 1:   writes into dummy (scalar)              (MASK = 0x1;   INVALID = 0x0)
// class 2:   writes into tmp[iz][iy][ix][0] (scalar) (MASK = 0x2;   INVALID = 0x1)
// class 3:   writes into tmp[iz][iy][ix][1] (scalar) (MASK = 0x4;   INVALID = 0x3)
// class 4:   writes into tmp[iz][iy][ix][2] (scalar) (MASK = 0x8;   INVALID = 0x7)
// class 5:   writes into tmp[iz][iy][ix][3] (scalar) (MASK = 0x10;  INVALID = 0xf)
// class 6:   writes into tmp[iz][iy][ix][4] (scalar) (MASK = 0x20;  INVALID = 0x1f)
// class 7:   writes into tmp[iz][iy][ix][5] (scalar) (MASK = 0x40;  INVALID = 0x3f)
// class 8:   writes into tmp[iz][iy][ix][6] (scalar) (MASK = 0x80;  INVALID = 0x7f)
// class 9:   writes into tmp[iz][iy][ix][7] (scalar) (MASK = 0x100; INVALID = 0xff)
// class 10:  writes into tmp[iz][iy][ix][0-1] (two-component vector) (MASK = 0x200; INVALID = 0x1f9)
// class 11:  writes into tmp[iz][iy][ix][2-3] (two-component vector) (MASK = 0x400; INVALID = 0x3e7)
// class 12:  writes into tmp[iz][iy][ix][4-5] (two-component vector) (MASK = 0x800; INVALID = 0x79f)
// class 13:  writes into tmp[iz][iy][ix][6-7] (two-component vector) (MASK = 0x1000; INVALID = 0xe7f)
// class 14:  writes into tmp[iz][iy][ix][0-2] (three-component vector) (MASK = 0x2000; INVALID = 0x19f1)
// class 15:  writes into tmp[iz][iy][ix][3-5] (three-component vector) (MASK = 0x4000; INVALID = 0x338f)
// class 16:  writes into dummy + tmp[iz][iy][ix][0-7] (nine-component rank-2 tensor) (MASK = 0x8000; INVALID = 0x0)
#define OPMAXCLASS 17

// 4th-order finite differences (non-vectorized)
///////////////////////////////////////////////////////////////////////////////
struct Operator_Vort_4th
{
    static const std::string NAME;
    static const int CLASS = 14;
    static const unsigned int MASK  = 0x2000;
    static const unsigned int INVALID = 0x19f1;
    StencilInfo stencil;

    Operator_Vort_4th(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_Vort_4th(const Operator_Vort_4th& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    const Real uy = -lab(ix,iy+2,iz).ru/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).ru/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).ru/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).ru/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                    const Real vx = -lab(ix+2,iy,iz).rv/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rv/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rv/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rv/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real uz = -lab(ix,iy,iz+2).ru/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).ru/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).ru/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).ru/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                    const Real wx = -lab(ix+2,iy,iz).rw/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rw/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rw/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rw/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real vz = -lab(ix,iy,iz+2).rv/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rv/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rv/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rv/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                    const Real wy = -lab(ix,iy+2,iz).rw/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rw/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rw/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rw/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);

                    o.tmp[iz][iy][ix][0] = (wy-vz)/(12*h);
                    o.tmp[iz][iy][ix][1] = (uz-wx)/(12*h);
                    o.tmp[iz][iy][ix][2] = (vx-uy)/(12*h);
                }
    }
};


struct Operator_IgradA2I_4th
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

    Operator_IgradA2I_4th(): stencil(-2,-2,-2,3,3,3, false, 1,6) {}
    Operator_IgradA2I_4th(const Operator_IgradA2I_4th& c): stencil(-2,-2,-2,3,3,3, false, 1, 6) {}

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    const Real dadx = -lab(ix+2,iy,iz).alpha2 + 8*lab(ix+1,iy,iz).alpha2 - 8*lab(ix-1,iy,iz).alpha2 + lab(ix-2,iy,iz).alpha2;
                    const Real dady = -lab(ix,iy+2,iz).alpha2 + 8*lab(ix,iy+1,iz).alpha2 - 8*lab(ix,iy-1,iz).alpha2 + lab(ix,iy-2,iz).alpha2;
                    const Real dadz = -lab(ix,iy,iz+2).alpha2 + 8*lab(ix,iy,iz+1).alpha2 - 8*lab(ix,iy,iz-1).alpha2 + lab(ix,iy,iz-2).alpha2;

                    o(ix,iy,iz).dummy = mysqrt((dadx*dadx + dady*dady + dadz*dadz)/(144.0*h*h));
                }
    }
};

struct Operator_divU_4th
{
    static const std::string NAME;
    static const int CLASS = 5;
    static const unsigned int MASK  = 0x10;
    static const unsigned int INVALID = 0xf;
    StencilInfo stencil;

    Operator_divU_4th(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_divU_4th(const Operator_divU_4th& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    const Real dudx = -lab(ix+2,iy,iz).ru/(lab(ix+2,iy,iz).alpha1rho1+lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).ru/(lab(ix+1,iy,iz).alpha1rho1+lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).ru/(lab(ix-1,iy,iz).alpha1rho1+lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).ru/(lab(ix-2,iy,iz).alpha1rho1+lab(ix-2,iy,iz).alpha2rho2);
                    const Real dvdy = -lab(ix,iy+2,iz).rv/(lab(ix,iy+2,iz).alpha1rho1+lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rv/(lab(ix,iy+1,iz).alpha1rho1+lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rv/(lab(ix,iy-1,iz).alpha1rho1+lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rv/(lab(ix,iy-2,iz).alpha1rho1+lab(ix,iy-2,iz).alpha2rho2);
                    const Real dwdz = -lab(ix,iy,iz+2).rw/(lab(ix,iy,iz+2).alpha1rho1+lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rw/(lab(ix,iy,iz+1).alpha1rho1+lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rw/(lab(ix,iy,iz-1).alpha1rho1+lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rw/(lab(ix,iy,iz-2).alpha1rho1+lab(ix,iy,iz-2).alpha2rho2);

                    o.tmp[iz][iy][ix][3] = (dudx + dvdy + dwdz)/(12.0*h);
                }
    }
};


struct Operator_SSDeviatoric_4th
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

    Operator_SSDeviatoric_4th(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_SSDeviatoric_4th(const Operator_SSDeviatoric_4th& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    const Real dudx = -lab(ix+2,iy,iz).ru/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).ru/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).ru/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).ru/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real dudy = -lab(ix,iy+2,iz).ru/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).ru/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).ru/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).ru/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                    const Real dudz = -lab(ix,iy,iz+2).ru/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).ru/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).ru/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).ru/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                    const Real dvdx = -lab(ix+2,iy,iz).rv/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rv/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rv/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rv/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real dvdy = -lab(ix,iy+2,iz).rv/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rv/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rv/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rv/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                    const Real dvdz = -lab(ix,iy,iz+2).rv/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rv/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rv/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rv/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                    const Real dwdx = -lab(ix+2,iy,iz).rw/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rw/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rw/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rw/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real dwdy = -lab(ix,iy+2,iz).rw/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rw/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rw/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rw/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                    const Real dwdz = -lab(ix,iy,iz+2).rw/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rw/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rw/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rw/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);

                    const Real t1 = 2.0*(dudx*dudx + dvdy*dvdy + dwdz*dwdz)
                    + (dvdx + dudy)*dudy + (dwdx + dudz)*dudz
                    + (dudy + dvdx)*dvdx + (dwdy + dvdz)*dvdz
                    + (dudz + dwdx)*dwdx + (dvdz + dwdy)*dwdy;

                    const Real t2 = (dudx + dvdy + dwdz)*(dudx + dvdy + dwdz);

                    o(ix,iy,iz).dummy = (0.5*t1 - 1.0/3.0*t2)/(144*h*h);
                }
    }
};


struct Operator_Qcrit_4th
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

    Operator_Qcrit_4th(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_Qcrit_4th(const Operator_Qcrit_4th& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;
        const Real inv12h = 1.0/(12.0*h);

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    // velocity gradient
                    const Real ux = inv12h*(-lab(ix+2,iy,iz).ru/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).ru/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).ru/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).ru/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2));
                    const Real uy = inv12h*(-lab(ix,iy+2,iz).ru/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).ru/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).ru/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).ru/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2));
                    const Real uz = inv12h*(-lab(ix,iy,iz+2).ru/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).ru/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).ru/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).ru/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2));
                    const Real vx = inv12h*(-lab(ix+2,iy,iz).rv/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rv/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rv/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rv/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2));
                    const Real vy = inv12h*(-lab(ix,iy+2,iz).rv/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rv/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rv/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rv/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2));
                    const Real vz = inv12h*(-lab(ix,iy,iz+2).rv/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rv/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rv/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rv/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2));
                    const Real wx = inv12h*(-lab(ix+2,iy,iz).rw/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rw/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rw/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rw/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2));
                    const Real wy = inv12h*(-lab(ix,iy+2,iz).rw/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rw/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rw/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rw/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2));
                    const Real wz = inv12h*(-lab(ix,iy,iz+2).rw/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rw/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rw/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rw/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2));

                    // rate of rotation / rate of strain
                    const Real O[9] = {
                        0.0,         0.5*(uy-vx), 0.5*(uz-wx),
                        0.5*(vx-uy), 0.0,         0.5*(vz-wy),
                        0.5*(wx-uz), 0.5*(wy-vz), 0.0};

                    const Real S[9] = {
                        ux,          0.5*(uy+vx), 0.5*(uz+wx),
                        0.5*(vx+uy), vy,          0.5*(vz+wy),
                        0.5*(wx+uz), 0.5*(wy+vz), wz};

                    // frobenius norm
                    Real IOI2 = 0.0;
                    Real ISI2 = 0.0;
                    for (int i=0; i<9; ++i)
                    {
                        IOI2 += O[i]*O[i];
                        ISI2 += S[i]*S[i];
                    }

                    o(ix,iy,iz).dummy = 0.5*(IOI2 - ISI2);
                }
    }
};


struct Operator_gradU_4th
{
    static const std::string NAME;
    static const int CLASS = 16;
    static const unsigned int MASK  = 0x8000;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

    Operator_gradU_4th(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_gradU_4th(const Operator_gradU_4th& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double fac = 1.0/(12.0*info.h_gridpoint);

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    const Real dudx = -lab(ix+2,iy,iz).ru/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).ru/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).ru/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).ru/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real dudy = -lab(ix,iy+2,iz).ru/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).ru/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).ru/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).ru/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                    const Real dudz = -lab(ix,iy,iz+2).ru/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).ru/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).ru/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).ru/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                    const Real dvdx = -lab(ix+2,iy,iz).rv/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rv/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rv/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rv/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real dvdy = -lab(ix,iy+2,iz).rv/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rv/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rv/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rv/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                    const Real dvdz = -lab(ix,iy,iz+2).rv/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rv/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rv/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rv/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                    const Real dwdx = -lab(ix+2,iy,iz).rw/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rw/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rw/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rw/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real dwdy = -lab(ix,iy+2,iz).rw/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rw/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rw/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rw/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                    const Real dwdz = -lab(ix,iy,iz+2).rw/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rw/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rw/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rw/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);

                    o(ix,iy,iz).dummy = dudx*fac;
                    o.tmp[iz][iy][ix][0] = dudy*fac;
                    o.tmp[iz][iy][ix][1] = dudz*fac;
                    o.tmp[iz][iy][ix][2] = dvdx*fac;
                    o.tmp[iz][iy][ix][3] = dvdy*fac;
                    o.tmp[iz][iy][ix][4] = dvdz*fac;
                    o.tmp[iz][iy][ix][5] = dwdx*fac;
                    o.tmp[iz][iy][ix][6] = dwdy*fac;
                    o.tmp[iz][iy][ix][7] = dwdz*fac;
                }
    }
};


// special
///////////////////////////////////////////////////////////////////////////////
template <int comp>
struct Operator_gradUij_4th
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

    Operator_gradUij_4th(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_gradUij_4th(const Operator_gradUij_4th& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;
        const double fac = 1.0/(12.0*h);

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    if (0 == comp)
                    {
                        const Real dudx = -lab(ix+2,iy,iz).ru/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).ru/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).ru/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).ru/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dudx*fac;
                    }
                    else if (1 == comp)
                    {
                        const Real dudy = -lab(ix,iy+2,iz).ru/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).ru/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).ru/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).ru/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dudy*fac;
                    }
                    else if (2 == comp)
                    {
                        const Real dudz = -lab(ix,iy,iz+2).ru/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).ru/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).ru/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).ru/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                        o(ix,iy,iz).dummy = dudz*fac;
                    }
                    else if (3 == comp)
                    {
                        const Real dvdx = -lab(ix+2,iy,iz).rv/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rv/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rv/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rv/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dvdx*fac;
                    }
                    else if (4 == comp)
                    {
                        const Real dvdy = -lab(ix,iy+2,iz).rv/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rv/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rv/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rv/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dvdy*fac;
                    }
                    else if (5 == comp)
                    {
                        const Real dvdz = -lab(ix,iy,iz+2).rv/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rv/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rv/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rv/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                        o(ix,iy,iz).dummy = dvdz*fac;
                    }
                    else if (6 == comp)
                    {
                        const Real dwdx = -lab(ix+2,iy,iz).rw/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rw/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rw/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rw/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dwdx*fac;
                    }
                    else if (7 == comp)
                    {
                        const Real dwdy = -lab(ix,iy+2,iz).rw/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rw/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rw/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rw/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dwdy*fac;
                    }
                    else if (8 == comp)
                    {
                        const Real dwdz = -lab(ix,iy,iz+2).rw/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rw/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rw/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rw/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                        o(ix,iy,iz).dummy = dwdz*fac;
                    }
                }
    }
};


struct Operator_Ucontraction_4th
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

    Operator_Ucontraction_4th(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_Ucontraction_4th(const Operator_Ucontraction_4th& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {

                    const Real dudx = -lab(ix+2,iy,iz).ru/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).ru/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).ru/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).ru/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real dudy = -lab(ix,iy+2,iz).ru/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).ru/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).ru/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).ru/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                    const Real dudz = -lab(ix,iy,iz+2).ru/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).ru/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).ru/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).ru/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                    const Real dvdx = -lab(ix+2,iy,iz).rv/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rv/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rv/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rv/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real dvdy = -lab(ix,iy+2,iz).rv/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rv/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rv/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rv/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                    const Real dvdz = -lab(ix,iy,iz+2).rv/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rv/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rv/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rv/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);
                    const Real dwdx = -lab(ix+2,iy,iz).rw/(lab(ix+2,iy,iz).alpha1rho1 + lab(ix+2,iy,iz).alpha2rho2) + 8*lab(ix+1,iy,iz).rw/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - 8*lab(ix-1,iy,iz).rw/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2) + lab(ix-2,iy,iz).rw/(lab(ix-2,iy,iz).alpha1rho1 + lab(ix-2,iy,iz).alpha2rho2);
                    const Real dwdy = -lab(ix,iy+2,iz).rw/(lab(ix,iy+2,iz).alpha1rho1 + lab(ix,iy+2,iz).alpha2rho2) + 8*lab(ix,iy+1,iz).rw/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - 8*lab(ix,iy-1,iz).rw/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2) + lab(ix,iy-2,iz).rw/(lab(ix,iy-2,iz).alpha1rho1 + lab(ix,iy-2,iz).alpha2rho2);
                    const Real dwdz = -lab(ix,iy,iz+2).rw/(lab(ix,iy,iz+2).alpha1rho1 + lab(ix,iy,iz+2).alpha2rho2) + 8*lab(ix,iy,iz+1).rw/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - 8*lab(ix,iy,iz-1).rw/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2) + lab(ix,iy,iz-2).rw/(lab(ix,iy,iz-2).alpha1rho1 + lab(ix,iy,iz-2).alpha2rho2);

                    o(ix,iy,iz).dummy = -(dudx*dudx+dvdy*dvdy+dwdz*dwdz + 2.0*(dudy*dvdx+dudz*dwdx+dvdz*dwdy))/(144*h*h);
                }
    }
};


struct Operator_PIC_4th
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

    Operator_PIC_4th(): stencil(-2,-2,-2,3,3,3, false, 3, 2,3,4) {}
    Operator_PIC_4th(const Operator_PIC_4th& c): stencil(-2,-2,-2,3,3,3, false, 3, 2,3,4) {}

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double gamma = Simulation_Environment::GAMMA1;
        const double factor = (gamma-1.0)/gamma;

        const double h = info.h_gridpoint;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {

                    const Real dudx = -lab(ix+2,iy,iz).ru + 8*lab(ix+1,iy,iz).ru - 8*lab(ix-1,iy,iz).ru + lab(ix-2,iy,iz).ru;
                    const Real dudy = -lab(ix,iy+2,iz).ru + 8*lab(ix,iy+1,iz).ru - 8*lab(ix,iy-1,iz).ru + lab(ix,iy-2,iz).ru;
                    const Real dudz = -lab(ix,iy,iz+2).ru + 8*lab(ix,iy,iz+1).ru - 8*lab(ix,iy,iz-1).ru + lab(ix,iy,iz-2).ru;
                    const Real dvdx = -lab(ix+2,iy,iz).rv + 8*lab(ix+1,iy,iz).rv - 8*lab(ix-1,iy,iz).rv + lab(ix-2,iy,iz).rv;
                    const Real dvdy = -lab(ix,iy+2,iz).rv + 8*lab(ix,iy+1,iz).rv - 8*lab(ix,iy-1,iz).rv + lab(ix,iy-2,iz).rv;
                    const Real dvdz = -lab(ix,iy,iz+2).rv + 8*lab(ix,iy,iz+1).rv - 8*lab(ix,iy,iz-1).rv + lab(ix,iy,iz-2).rv;
                    const Real dwdx = -lab(ix+2,iy,iz).rw + 8*lab(ix+1,iy,iz).rw - 8*lab(ix-1,iy,iz).rw + lab(ix-2,iy,iz).rw;
                    const Real dwdy = -lab(ix,iy+2,iz).rw + 8*lab(ix,iy+1,iz).rw - 8*lab(ix,iy-1,iz).rw + lab(ix,iy-2,iz).rw;
                    const Real dwdz = -lab(ix,iy,iz+2).rw + 8*lab(ix,iy,iz+1).rw - 8*lab(ix,iy,iz-1).rw + lab(ix,iy,iz-2).rw;

                    o(ix,iy,iz).dummy = -factor*(dudx*dudx+dvdy*dvdy+dwdz*dwdz + 2.0*(dudy*dvdx+dudz*dwdx+dvdz*dwdy))/(144*h*h);
                }
    }
};

// abbreviations
///////////////////////////////////////////////////////////////////////////////
typedef Operator_Vort_4th         OVort_4;
typedef Operator_IgradA2I_4th     OIgradA2I_4;
typedef Operator_divU_4th         OdivU_4;
typedef Operator_SSDeviatoric_4th OSSd_4;
typedef Operator_Qcrit_4th        OQcrit_4;
typedef Operator_gradU_4th        OgradU_4;



// TODO: (fabianw@mavt.ethz.ch; Fri Oct 14 07:44:37 2016) TESTING
#if 0

struct EvaluateVorticity_CPP
{
    StencilInfo stencil;
    int stencil_start[3], stencil_end[3];

#if 0 /* TODO: (fabianw; Wed 24 Jun 2015 10:45:16 AM CEST) disabled for now */
    EvaluateVorticity_CPP(): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
    }

    EvaluateVorticity_CPP(const EvaluateVorticity_CPP& c): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
    }
#else
    EvaluateVorticity_CPP(): stencil(-3,-3,-3,4,4,4, false, 5, 0,1,2,3,4)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -3;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
    }

    EvaluateVorticity_CPP(const EvaluateVorticity_CPP& c): stencil(-3,-3,-3,4,4,4, false, 5, 0,1,2,3,4)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -3;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
    }
#endif

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    const Real uy = lab(ix,iy+1,iz).ru/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - lab(ix,iy-1,iz).ru/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2);
                    const Real vx = lab(ix+1,iy,iz).rv/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - lab(ix-1,iy,iz).rv/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2);
                    const Real uz = lab(ix,iy,iz+1).ru/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - lab(ix,iy,iz-1).ru/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2);
                    const Real wx = lab(ix+1,iy,iz).rw/(lab(ix+1,iy,iz).alpha1rho1 + lab(ix+1,iy,iz).alpha2rho2) - lab(ix-1,iy,iz).rw/(lab(ix-1,iy,iz).alpha1rho1 + lab(ix-1,iy,iz).alpha2rho2);
                    const Real vz = lab(ix,iy,iz+1).rv/(lab(ix,iy,iz+1).alpha1rho1 + lab(ix,iy,iz+1).alpha2rho2) - lab(ix,iy,iz-1).rv/(lab(ix,iy,iz-1).alpha1rho1 + lab(ix,iy,iz-1).alpha2rho2);
                    const Real wy = lab(ix,iy+1,iz).rw/(lab(ix,iy+1,iz).alpha1rho1 + lab(ix,iy+1,iz).alpha2rho2) - lab(ix,iy-1,iz).rw/(lab(ix,iy-1,iz).alpha1rho1 + lab(ix,iy-1,iz).alpha2rho2);

                    o.tmp[iz][iy][ix][0] = 0.5*(wy-vz)/h;
                    o.tmp[iz][iy][ix][1] = 0.5*(uz-wx)/h;
                    o.tmp[iz][iy][ix][2] = 0.5*(vx-uy)/h;
                }
    }
};



struct EvaluateVelocityDivergence_CPP
{
    StencilInfo stencil;
    int stencil_start[3], stencil_end[3];

#if 1
    EvaluateVelocityDivergence_CPP(): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
    }

    EvaluateVelocityDivergence_CPP(const EvaluateVelocityDivergence_CPP& c): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
    }
#else
    EvaluateVelocityDivergence_CPP(): stencil(-3,-3,-3,4,4,4, false, 5, 0,1,2,3,4)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -3;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
    }

    EvaluateVelocityDivergence_CPP(const EvaluateVelocityDivergence_CPP& c): stencil(-3,-3,-3,4,4,4, false, 5, 0,1,2,3,4)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -3;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
    }
#endif

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    const Real ux= lab(ix+1,iy,iz).ru/(lab(ix+1,iy,iz).alpha1rho1+lab(ix+1,iy,iz).alpha2rho2)-lab(ix-1,iy,iz).ru/(lab(ix-1,iy,iz).alpha1rho1+lab(ix-1,iy,iz).alpha2rho2);
                    const Real vy= lab(ix,iy+1,iz).rv/(lab(ix,iy+1,iz).alpha1rho1+lab(ix,iy+1,iz).alpha2rho2)-lab(ix,iy-1,iz).rv/(lab(ix,iy-1,iz).alpha1rho1+lab(ix,iy-1,iz).alpha2rho2);
                    const Real wz= lab(ix,iy,iz+1).rw/(lab(ix,iy,iz+1).alpha1rho1+lab(ix,iy,iz+1).alpha2rho2)-lab(ix,iy,iz-1).rw/(lab(ix,iy,iz-1).alpha1rho1+lab(ix,iy,iz-1).alpha2rho2);

                    o.tmp[iz][iy][ix][3] = 0.5*(ux+vy+wz)/h;
                }
    }
};







template <int comp>
struct GradU_CPP
{
    StencilInfo stencil;
    int stencil_start[3], stencil_end[3];

    GradU_CPP(): stencil(-3,-3,-3,4,4,4, false, 5, 0,1,2,3,4)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -3;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
    }

    GradU_CPP(const GradU_CPP& c): stencil(-3,-3,-3,4,4,4, false, 5, 0,1,2,3,4)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -3;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
    }

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;
        const double fac = 0.5/h;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    if (0 == comp)
                    {
                        const Real dudx= lab(ix+1,iy,iz).ru/(lab(ix+1,iy,iz).alpha1rho1+lab(ix+1,iy,iz).alpha2rho2)-lab(ix-1,iy,iz).ru/(lab(ix-1,iy,iz).alpha1rho1+lab(ix-1,iy,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dudx*fac;
                    }
                    else if (1 == comp)
                    {
                        const Real dudy= lab(ix,iy+1,iz).ru/(lab(ix,iy+1,iz).alpha1rho1+lab(ix,iy+1,iz).alpha2rho2)-lab(ix,iy-1,iz).ru/(lab(ix,iy-1,iz).alpha1rho1+lab(ix,iy-1,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dudy*fac;
                    }
                    else if (2 == comp)
                    {
                        const Real dudz= lab(ix,iy,iz+1).ru/(lab(ix,iy,iz+1).alpha1rho1+lab(ix,iy,iz+1).alpha2rho2)-lab(ix,iy,iz-1).ru/(lab(ix,iy,iz-1).alpha1rho1+lab(ix,iy,iz-1).alpha2rho2);
                        o(ix,iy,iz).dummy = dudz*fac;
                    }
                    else if (3 == comp)
                    {
                        const Real dvdx= lab(ix+1,iy,iz).rv/(lab(ix+1,iy,iz).alpha1rho1+lab(ix+1,iy,iz).alpha2rho2)-lab(ix-1,iy,iz).rv/(lab(ix-1,iy,iz).alpha1rho1+lab(ix-1,iy,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dvdx*fac;
                    }
                    else if (4 == comp)
                    {
                        const Real dvdy= lab(ix,iy+1,iz).rv/(lab(ix,iy+1,iz).alpha1rho1+lab(ix,iy+1,iz).alpha2rho2)-lab(ix,iy-1,iz).rv/(lab(ix,iy-1,iz).alpha1rho1+lab(ix,iy-1,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dvdy*fac;
                    }
                    else if (5 == comp)
                    {
                        const Real dvdz= lab(ix,iy,iz+1).rv/(lab(ix,iy,iz+1).alpha1rho1+lab(ix,iy,iz+1).alpha2rho2)-lab(ix,iy,iz-1).rv/(lab(ix,iy,iz-1).alpha1rho1+lab(ix,iy,iz-1).alpha2rho2);
                        o(ix,iy,iz).dummy = dvdz*fac;
                    }
                    else if (6 == comp)
                    {
                        const Real dwdx= lab(ix+1,iy,iz).rw/(lab(ix+1,iy,iz).alpha1rho1+lab(ix+1,iy,iz).alpha2rho2)-lab(ix-1,iy,iz).rw/(lab(ix-1,iy,iz).alpha1rho1+lab(ix-1,iy,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dwdx*fac;
                    }
                    else if (7 == comp)
                    {
                        const Real dwdy= lab(ix,iy+1,iz).rw/(lab(ix,iy+1,iz).alpha1rho1+lab(ix,iy+1,iz).alpha2rho2)-lab(ix,iy-1,iz).rw/(lab(ix,iy-1,iz).alpha1rho1+lab(ix,iy-1,iz).alpha2rho2);
                        o(ix,iy,iz).dummy = dwdy*fac;
                    }
                    else if (8 == comp)
                    {
                        const Real dwdz= lab(ix,iy,iz+1).rw/(lab(ix,iy,iz+1).alpha1rho1+lab(ix,iy,iz+1).alpha2rho2)-lab(ix,iy,iz-1).rw/(lab(ix,iy,iz-1).alpha1rho1+lab(ix,iy,iz-1).alpha2rho2);
                        o(ix,iy,iz).dummy = dwdz*fac;
                    }
                }
    }
};



struct EvaluateMagGradientPressureFourthOrder_CPP
{
    StencilInfo stencil;
    int stencil_start[3], stencil_end[3];

    EvaluateMagGradientPressureFourthOrder_CPP(): stencil(-2,-2,-2,3,3,3, false, 4, 0,1,2,3)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -2;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 3;
    }

    EvaluateMagGradientPressureFourthOrder_CPP(const EvaluateMagGradientPressureFourthOrder_CPP& c): stencil(-2,-2,-2,3,3,3, false, 1, 6)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -2;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 3;
    }

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double h = info.h_gridpoint;

        const Real g1m1Inv = (Real)1.0/(Simulation_Environment::GAMMA1 - (Real)1.0);
        const Real g2m1Inv = (Real)1.0/(Simulation_Environment::GAMMA2 - (Real)1.0);

        const Real pc1_g1m1Inv = g1m1Inv * Simulation_Environment::GAMMA1 * Simulation_Environment::PC1;
        const Real pc2_g2m1Inv = g2m1Inv * Simulation_Environment::GAMMA2 * Simulation_Environment::PC2;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    // convert energy to pressure
#ifndef _CONVERTCLIP_
                    const Real r1a1_ixppiyiz = lab(ix+2,iy,iz).alpha1rho1;
                    const Real r1a1_ixpiyiz = lab(ix+1,iy,iz).alpha1rho1;
                    const Real r1a1_ixmmiyiz = lab(ix-2,iy,iz).alpha1rho1;
                    const Real r1a1_ixmiyiz = lab(ix-1,iy,iz).alpha1rho1;

                    const Real r2a2_ixppiyiz = lab(ix+2,iy,iz).alpha2rho2;
                    const Real r2a2_ixpiyiz = lab(ix+1,iy,iz).alpha2rho2;
                    const Real r2a2_ixmmiyiz = lab(ix-2,iy,iz).alpha2rho2;
                    const Real r2a2_ixmiyiz = lab(ix-1,iy,iz).alpha2rho2;

                    const Real a2_ixppiyiz = lab(ix+2,iy,iz).alpha2;
                    const Real a2_ixpiyiz = lab(ix+1,iy,iz).alpha2;
                    const Real a2_ixmmiyiz = lab(ix-2,iy,iz).alpha2;
                    const Real a2_ixmiyiz = lab(ix-1,iy,iz).alpha2;

                    const Real a1_ixppiyiz = (Real)1.0 - a2_ixppiyiz;
                    const Real a1_ixpiyiz = (Real)1.0 - a2_ixpiyiz;
                    const Real a1_ixmmiyiz = (Real)1.0 - a2_ixmmiyiz;
                    const Real a1_ixmiyiz = (Real)1.0 - a2_ixmiyiz;
#else
                    const Real r1a1_ixppiyiz = max(Real(0.0),lab(ix+2,iy,iz).alpha1rho1);
                    const Real r1a1_ixpiyiz = max(Real(0.0),lab(ix+1,iy,iz).alpha1rho1);
                    const Real r1a1_ixmmiyiz = max(Real(0.0),lab(ix-2,iy,iz).alpha1rho1);
                    const Real r1a1_ixmiyiz = max(Real(0.0),lab(ix-1,iy,iz).alpha1rho1);

                    const Real r2a2_ixppiyiz = max(Real(0.0),lab(ix+2,iy,iz).alpha2rho2);
                    const Real r2a2_ixpiyiz = max(Real(0.0),lab(ix+1,iy,iz).alpha2rho2);
                    const Real r2a2_ixmmiyiz = max(Real(0.0),lab(ix-2,iy,iz).alpha2rho2);
                    const Real r2a2_ixmiyiz = max(Real(0.0),lab(ix-1,iy,iz).alpha2rho2);

                    const Real a2_ixppiyiz = max(Real(0.0),min(Real(1.0),lab(ix+2,iy,iz).alpha2));
                    const Real a2_ixpiyiz = max(Real(0.0),min(Real(1.0),lab(ix+1,iy,iz).alpha2));
                    const Real a2_ixmmiyiz = max(Real(0.0),min(Real(1.0),lab(ix-2,iy,iz).alpha2));
                    const Real a2_ixmiyiz = max(Real(0.0),min(Real(1.0),lab(ix-1,iy,iz).alpha2));

                    const Real a1_ixppiyiz = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixppiyiz));
                    const Real a1_ixpiyiz = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixpiyiz));
                    const Real a1_ixmmiyiz = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixmmiyiz));
                    const Real a1_ixmiyiz = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixmiyiz));
#endif

                    const Real rinv_ixppiyiz = (Real)1.0/(r1a1_ixppiyiz + r2a2_ixppiyiz);
                    const Real rinv_ixpiyiz = (Real)1.0/(r1a1_ixpiyiz + r2a2_ixpiyiz);
                    const Real rinv_ixmmiyiz = (Real)1.0/(r1a1_ixmmiyiz + r2a2_ixmmiyiz);
                    const Real rinv_ixmiyiz = (Real)1.0/(r1a1_ixmiyiz + r2a2_ixmiyiz);

                    const Real ke_ixppiyiz = (Real)0.5 * rinv_ixppiyiz * (lab(ix+2,iy,iz).ru*lab(ix+2,iy,iz).ru + lab(ix+2,iy,iz).rv*lab(ix+2,iy,iz).rv + lab(ix+2,iy,iz).rw*lab(ix+2,iy,iz).rw);
                    const Real ke_ixpiyiz = (Real)0.5 * rinv_ixpiyiz * (lab(ix+1,iy,iz).ru*lab(ix+1,iy,iz).ru + lab(ix+1,iy,iz).rv*lab(ix+1,iy,iz).rv + lab(ix+1,iy,iz).rw*lab(ix+1,iy,iz).rw);
                    const Real ke_ixmmiyiz = (Real)0.5 * rinv_ixmmiyiz * (lab(ix-2,iy,iz).ru*lab(ix-2,iy,iz).ru + lab(ix-2,iy,iz).rv*lab(ix-2,iy,iz).rv + lab(ix-2,iy,iz).rw*lab(ix-2,iy,iz).rw);
                    const Real ke_ixmiyiz = (Real)0.5 * rinv_ixmiyiz * (lab(ix-1,iy,iz).ru*lab(ix-1,iy,iz).ru + lab(ix-1,iy,iz).rv*lab(ix-1,iy,iz).rv + lab(ix-1,iy,iz).rw*lab(ix-1,iy,iz).rw);

                    const Real GmixInv_ixppiyiz = (Real)1.0 / (a1_ixppiyiz * g1m1Inv + a2_ixppiyiz * g2m1Inv);
                    const Real GmixInv_ixpiyiz = (Real)1.0 / (a1_ixpiyiz * g1m1Inv + a2_ixpiyiz * g2m1Inv);
                    const Real GmixInv_ixmmiyiz = (Real)1.0 / (a1_ixmmiyiz * g1m1Inv + a2_ixmmiyiz * g2m1Inv);
                    const Real GmixInv_ixmiyiz = (Real)1.0 / (a1_ixmiyiz * g1m1Inv + a2_ixmiyiz * g2m1Inv);

                    const Real Pmix_ixppiyiz = a1_ixppiyiz * pc1_g1m1Inv + a2_ixppiyiz * pc2_g2m1Inv;
                    const Real Pmix_ixpiyiz = a1_ixpiyiz * pc1_g1m1Inv + a2_ixpiyiz * pc2_g2m1Inv;
                    const Real Pmix_ixmmiyiz = a1_ixmmiyiz * pc1_g1m1Inv + a2_ixmmiyiz * pc2_g2m1Inv;
                    const Real Pmix_ixmiyiz = a1_ixmiyiz * pc1_g1m1Inv + a2_ixmiyiz * pc2_g2m1Inv;

                    Real p_ixppiyiz = GmixInv_ixppiyiz * lab(ix+2,iy,iz).energy - GmixInv_ixppiyiz * ke_ixppiyiz - GmixInv_ixppiyiz * Pmix_ixppiyiz;
                    Real p_ixpiyiz = GmixInv_ixpiyiz * lab(ix+1,iy,iz).energy - GmixInv_ixpiyiz * ke_ixpiyiz - GmixInv_ixpiyiz * Pmix_ixpiyiz;
                    Real p_ixmmiyiz = GmixInv_ixmmiyiz * lab(ix-2,iy,iz).energy - GmixInv_ixmmiyiz * ke_ixmmiyiz - GmixInv_ixmmiyiz * Pmix_ixmmiyiz;
                    Real p_ixmiyiz = GmixInv_ixmiyiz * lab(ix-1,iy,iz).energy - GmixInv_ixmiyiz * ke_ixmiyiz - GmixInv_ixmiyiz * Pmix_ixmiyiz;
#ifdef _CONVERTCLIP_
#ifdef _NOK_
                    const Real pth_ixppiyiz = Pmix_ixppiyiz*GmixInv_ixppiyiz/(GmixInv_ixppiyiz+1.0);
                    const Real pth_ixpiyiz = Pmix_ixpiyiz*GmixInv_ixpiyiz/(GmixInv_ixpiyiz+1.0);
                    const Real pth_ixmmiyiz = Pmix_ixmmiyiz*GmixInv_ixmmiyiz/(GmixInv_ixmmiyiz+1.0);
                    const Real pth_ixmiyiz = Pmix_ixmiyiz*GmixInv_ixmiyiz/(GmixInv_ixmiyiz+1.0);
#else
                    const Real pth_ixppiyiz = a2_ixppiyiz == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixppiyiz == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
                    const Real pth_ixpiyiz = a2_ixpiyiz == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixpiyiz == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
                    const Real pth_ixmmiyiz = a2_ixmmiyiz == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixmmiyiz == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
                    const Real pth_ixmiyiz = a2_ixmiyiz == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixmiyiz == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
#endif
                    //                   const Real deltap = pthresh < (-static_cast<Real>(2.0)*actpres) ? (-static_cast<Real>(4.0)*actpres - pthresh) : static_cast<Real>(0.0);
                    const Real deltap_ixppiyiz =  pth_ixppiyiz <= (-p_ixppiyiz + static_cast<Real>(PRESEPS)) ? (-p_ixppiyiz - pth_ixppiyiz + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
                    const Real deltap_ixpiyiz =  pth_ixpiyiz <= (-p_ixpiyiz + static_cast<Real>(PRESEPS)) ? (-p_ixpiyiz - pth_ixpiyiz + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
                    const Real deltap_ixmmiyiz =  pth_ixmmiyiz <= (-p_ixmmiyiz + static_cast<Real>(PRESEPS)) ? (-p_ixmmiyiz - pth_ixmmiyiz + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
                    const Real deltap_ixmiyiz =  pth_ixmiyiz <= (-p_ixmiyiz + static_cast<Real>(PRESEPS)) ? (-p_ixmiyiz - pth_ixmiyiz + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);

                    p_ixppiyiz += deltap_ixppiyiz;
                    p_ixpiyiz += deltap_ixpiyiz;
                    p_ixmmiyiz += deltap_ixmmiyiz;
                    p_ixmiyiz += deltap_ixmiyiz;
#endif

                    const Real dadx = -p_ixppiyiz + 8*p_ixpiyiz - 8*p_ixmiyiz + p_ixmmiyiz;

                    // convert energy to pressure
#ifndef _CONVERTCLIP_
                    const Real r1a1_ixiyppiz = lab(ix+2,iy,iz).alpha1rho1;
                    const Real r1a1_ixiypiz = lab(ix+1,iy,iz).alpha1rho1;
                    const Real r1a1_ixiymmiz = lab(ix-2,iy,iz).alpha1rho1;
                    const Real r1a1_ixiymiz = lab(ix-1,iy,iz).alpha1rho1;

                    const Real r2a2_ixiyppiz = lab(ix+2,iy,iz).alpha2rho2;
                    const Real r2a2_ixiypiz = lab(ix+1,iy,iz).alpha2rho2;
                    const Real r2a2_ixiymmiz = lab(ix-2,iy,iz).alpha2rho2;
                    const Real r2a2_ixiymiz = lab(ix-1,iy,iz).alpha2rho2;

                    const Real a2_ixiyppiz = lab(ix+2,iy,iz).alpha2;
                    const Real a2_ixiypiz = lab(ix+1,iy,iz).alpha2;
                    const Real a2_ixiymmiz = lab(ix-2,iy,iz).alpha2;
                    const Real a2_ixiymiz = lab(ix-1,iy,iz).alpha2;

                    const Real a1_ixiyppiz = (Real)1.0 - a2_ixiyppiz;
                    const Real a1_ixiypiz = (Real)1.0 - a2_ixiypiz;
                    const Real a1_ixiymmiz = (Real)1.0 - a2_ixiymmiz;
                    const Real a1_ixiymiz = (Real)1.0 - a2_ixiymiz;
#else
                    const Real r1a1_ixiyppiz = max(Real(0.0),lab(ix+2,iy,iz).alpha1rho1);
                    const Real r1a1_ixiypiz = max(Real(0.0),lab(ix+1,iy,iz).alpha1rho1);
                    const Real r1a1_ixiymmiz = max(Real(0.0),lab(ix-2,iy,iz).alpha1rho1);
                    const Real r1a1_ixiymiz = max(Real(0.0),lab(ix-1,iy,iz).alpha1rho1);

                    const Real r2a2_ixiyppiz = max(Real(0.0),lab(ix+2,iy,iz).alpha2rho2);
                    const Real r2a2_ixiypiz = max(Real(0.0),lab(ix+1,iy,iz).alpha2rho2);
                    const Real r2a2_ixiymmiz = max(Real(0.0),lab(ix-2,iy,iz).alpha2rho2);
                    const Real r2a2_ixiymiz = max(Real(0.0),lab(ix-1,iy,iz).alpha2rho2);

                    const Real a2_ixiyppiz = max(Real(0.0),min(Real(1.0),lab(ix+2,iy,iz).alpha2));
                    const Real a2_ixiypiz = max(Real(0.0),min(Real(1.0),lab(ix+1,iy,iz).alpha2));
                    const Real a2_ixiymmiz = max(Real(0.0),min(Real(1.0),lab(ix-2,iy,iz).alpha2));
                    const Real a2_ixiymiz = max(Real(0.0),min(Real(1.0),lab(ix-1,iy,iz).alpha2));

                    const Real a1_ixiyppiz = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixiyppiz));
                    const Real a1_ixiypiz = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixiypiz));
                    const Real a1_ixiymmiz = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixiymmiz));
                    const Real a1_ixiymiz = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixiymiz));
#endif

                    const Real rinv_ixiyppiz = (Real)1.0/(r1a1_ixiyppiz + r2a2_ixiyppiz);
                    const Real rinv_ixiypiz = (Real)1.0/(r1a1_ixiypiz + r2a2_ixiypiz);
                    const Real rinv_ixiymmiz = (Real)1.0/(r1a1_ixiymmiz + r2a2_ixiymmiz);
                    const Real rinv_ixiymiz = (Real)1.0/(r1a1_ixiymiz + r2a2_ixiymiz);

                    const Real ke_ixiyppiz = (Real)0.5 * rinv_ixiyppiz * (lab(ix+2,iy,iz).ru*lab(ix+2,iy,iz).ru + lab(ix+2,iy,iz).rv*lab(ix+2,iy,iz).rv + lab(ix+2,iy,iz).rw*lab(ix+2,iy,iz).rw);
                    const Real ke_ixiypiz = (Real)0.5 * rinv_ixiypiz * (lab(ix+1,iy,iz).ru*lab(ix+1,iy,iz).ru + lab(ix+1,iy,iz).rv*lab(ix+1,iy,iz).rv + lab(ix+1,iy,iz).rw*lab(ix+1,iy,iz).rw);
                    const Real ke_ixiymmiz = (Real)0.5 * rinv_ixiymmiz * (lab(ix-2,iy,iz).ru*lab(ix-2,iy,iz).ru + lab(ix-2,iy,iz).rv*lab(ix-2,iy,iz).rv + lab(ix-2,iy,iz).rw*lab(ix-2,iy,iz).rw);
                    const Real ke_ixiymiz = (Real)0.5 * rinv_ixiymiz * (lab(ix-1,iy,iz).ru*lab(ix-1,iy,iz).ru + lab(ix-1,iy,iz).rv*lab(ix-1,iy,iz).rv + lab(ix-1,iy,iz).rw*lab(ix-1,iy,iz).rw);

                    const Real GmixInv_ixiyppiz = (Real)1.0 / (a1_ixiyppiz * g1m1Inv + a2_ixiyppiz * g2m1Inv);
                    const Real GmixInv_ixiypiz = (Real)1.0 / (a1_ixiypiz * g1m1Inv + a2_ixiypiz * g2m1Inv);
                    const Real GmixInv_ixiymmiz = (Real)1.0 / (a1_ixiymmiz * g1m1Inv + a2_ixiymmiz * g2m1Inv);
                    const Real GmixInv_ixiymiz = (Real)1.0 / (a1_ixiymiz * g1m1Inv + a2_ixiymiz * g2m1Inv);

                    const Real Pmix_ixiyppiz = a1_ixiyppiz * pc1_g1m1Inv + a2_ixiyppiz * pc2_g2m1Inv;
                    const Real Pmix_ixiypiz = a1_ixiypiz * pc1_g1m1Inv + a2_ixiypiz * pc2_g2m1Inv;
                    const Real Pmix_ixiymmiz = a1_ixiymmiz * pc1_g1m1Inv + a2_ixiymmiz * pc2_g2m1Inv;
                    const Real Pmix_ixiymiz = a1_ixiymiz * pc1_g1m1Inv + a2_ixiymiz * pc2_g2m1Inv;

                    Real p_ixiyppiz = GmixInv_ixiyppiz * lab(ix+2,iy,iz).energy - GmixInv_ixiyppiz * ke_ixiyppiz - GmixInv_ixiyppiz * Pmix_ixiyppiz;
                    Real p_ixiypiz = GmixInv_ixiypiz * lab(ix+1,iy,iz).energy - GmixInv_ixiypiz * ke_ixiypiz - GmixInv_ixiypiz * Pmix_ixiypiz;
                    Real p_ixiymmiz = GmixInv_ixiymmiz * lab(ix-2,iy,iz).energy - GmixInv_ixiymmiz * ke_ixiymmiz - GmixInv_ixiymmiz * Pmix_ixiymmiz;
                    Real p_ixiymiz = GmixInv_ixiymiz * lab(ix-1,iy,iz).energy - GmixInv_ixiymiz * ke_ixiymiz - GmixInv_ixiymiz * Pmix_ixiymiz;
#ifdef _CONVERTCLIP_
#ifdef _NOK_
                    const Real pth_ixiyppiz = Pmix_ixiyppiz*GmixInv_ixiyppiz/(GmixInv_ixiyppiz+1.0);
                    const Real pth_ixiypiz = Pmix_ixiypiz*GmixInv_ixiypiz/(GmixInv_ixiypiz+1.0);
                    const Real pth_ixiymmiz = Pmix_ixiymmiz*GmixInv_ixiymmiz/(GmixInv_ixiymmiz+1.0);
                    const Real pth_ixiymiz = Pmix_ixiymiz*GmixInv_ixiymiz/(GmixInv_ixiymiz+1.0);
#else
                    const Real pth_ixiyppiz = a2_ixiyppiz == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixiyppiz == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
                    const Real pth_ixiypiz = a2_ixiypiz == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixiypiz == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
                    const Real pth_ixiymmiz = a2_ixiymmiz == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixiymmiz == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
                    const Real pth_ixiymiz = a2_ixiymiz == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixiymiz == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
#endif
                    //                   const Real deltap = pthresh < (-static_cast<Real>(2.0)*actpres) ? (-static_cast<Real>(4.0)*actpres - pthresh) : static_cast<Real>(0.0);
                    const Real deltap_ixiyppiz =  pth_ixiyppiz <= (-p_ixiyppiz + static_cast<Real>(PRESEPS)) ? (-p_ixiyppiz - pth_ixiyppiz + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
                    const Real deltap_ixiypiz =  pth_ixiypiz <= (-p_ixiypiz + static_cast<Real>(PRESEPS)) ? (-p_ixiypiz - pth_ixiypiz + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
                    const Real deltap_ixiymmiz =  pth_ixiymmiz <= (-p_ixiymmiz + static_cast<Real>(PRESEPS)) ? (-p_ixiymmiz - pth_ixiymmiz + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
                    const Real deltap_ixiymiz =  pth_ixiymiz <= (-p_ixiymiz + static_cast<Real>(PRESEPS)) ? (-p_ixiymiz - pth_ixiymiz + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);

                    p_ixiyppiz += deltap_ixiyppiz;
                    p_ixiypiz += deltap_ixiypiz;
                    p_ixiymmiz += deltap_ixiymmiz;
                    p_ixiymiz += deltap_ixiymiz;
#endif

                    const Real dady = -p_ixiyppiz + 8*p_ixiypiz - 8*p_ixiymiz + p_ixiymmiz;

                    // convert energy to pressure
#ifndef _CONVERTCLIP_
                    const Real r1a1_ixiyizpp = lab(ix+2,iy,iz).alpha1rho1;
                    const Real r1a1_ixiyizp = lab(ix+1,iy,iz).alpha1rho1;
                    const Real r1a1_ixiyizmm = lab(ix-2,iy,iz).alpha1rho1;
                    const Real r1a1_ixiyizm = lab(ix-1,iy,iz).alpha1rho1;

                    const Real r2a2_ixiyizpp = lab(ix+2,iy,iz).alpha2rho2;
                    const Real r2a2_ixiyizp = lab(ix+1,iy,iz).alpha2rho2;
                    const Real r2a2_ixiyizmm = lab(ix-2,iy,iz).alpha2rho2;
                    const Real r2a2_ixiyizm = lab(ix-1,iy,iz).alpha2rho2;

                    const Real a2_ixiyizpp = lab(ix+2,iy,iz).alpha2;
                    const Real a2_ixiyizp = lab(ix+1,iy,iz).alpha2;
                    const Real a2_ixiyizmm = lab(ix-2,iy,iz).alpha2;
                    const Real a2_ixiyizm = lab(ix-1,iy,iz).alpha2;

                    const Real a1_ixiyizpp = (Real)1.0 - a2_ixiyizpp;
                    const Real a1_ixiyizp = (Real)1.0 - a2_ixiyizp;
                    const Real a1_ixiyizmm = (Real)1.0 - a2_ixiyizmm;
                    const Real a1_ixiyizm = (Real)1.0 - a2_ixiyizm;
#else
                    const Real r1a1_ixiyizpp = max(Real(0.0),lab(ix+2,iy,iz).alpha1rho1);
                    const Real r1a1_ixiyizp = max(Real(0.0),lab(ix+1,iy,iz).alpha1rho1);
                    const Real r1a1_ixiyizmm = max(Real(0.0),lab(ix-2,iy,iz).alpha1rho1);
                    const Real r1a1_ixiyizm = max(Real(0.0),lab(ix-1,iy,iz).alpha1rho1);

                    const Real r2a2_ixiyizpp = max(Real(0.0),lab(ix+2,iy,iz).alpha2rho2);
                    const Real r2a2_ixiyizp = max(Real(0.0),lab(ix+1,iy,iz).alpha2rho2);
                    const Real r2a2_ixiyizmm = max(Real(0.0),lab(ix-2,iy,iz).alpha2rho2);
                    const Real r2a2_ixiyizm = max(Real(0.0),lab(ix-1,iy,iz).alpha2rho2);

                    const Real a2_ixiyizpp = max(Real(0.0),min(Real(1.0),lab(ix+2,iy,iz).alpha2));
                    const Real a2_ixiyizp = max(Real(0.0),min(Real(1.0),lab(ix+1,iy,iz).alpha2));
                    const Real a2_ixiyizmm = max(Real(0.0),min(Real(1.0),lab(ix-2,iy,iz).alpha2));
                    const Real a2_ixiyizm = max(Real(0.0),min(Real(1.0),lab(ix-1,iy,iz).alpha2));

                    const Real a1_ixiyizpp = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixiyizpp));
                    const Real a1_ixiyizp = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixiyizp));
                    const Real a1_ixiyizmm = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixiyizmm));
                    const Real a1_ixiyizm = max(Real(0.0),min(Real(1.0),(Real)1.0 - a2_ixiyizm));
#endif

                    const Real rinv_ixiyizpp = (Real)1.0/(r1a1_ixiyizpp + r2a2_ixiyizpp);
                    const Real rinv_ixiyizp = (Real)1.0/(r1a1_ixiyizp + r2a2_ixiyizp);
                    const Real rinv_ixiyizmm = (Real)1.0/(r1a1_ixiyizmm + r2a2_ixiyizmm);
                    const Real rinv_ixiyizm = (Real)1.0/(r1a1_ixiyizm + r2a2_ixiyizm);

                    const Real ke_ixiyizpp = (Real)0.5 * rinv_ixiyizpp * (lab(ix+2,iy,iz).ru*lab(ix+2,iy,iz).ru + lab(ix+2,iy,iz).rv*lab(ix+2,iy,iz).rv + lab(ix+2,iy,iz).rw*lab(ix+2,iy,iz).rw);
                    const Real ke_ixiyizp = (Real)0.5 * rinv_ixiyizp * (lab(ix+1,iy,iz).ru*lab(ix+1,iy,iz).ru + lab(ix+1,iy,iz).rv*lab(ix+1,iy,iz).rv + lab(ix+1,iy,iz).rw*lab(ix+1,iy,iz).rw);
                    const Real ke_ixiyizmm = (Real)0.5 * rinv_ixiyizmm * (lab(ix-2,iy,iz).ru*lab(ix-2,iy,iz).ru + lab(ix-2,iy,iz).rv*lab(ix-2,iy,iz).rv + lab(ix-2,iy,iz).rw*lab(ix-2,iy,iz).rw);
                    const Real ke_ixiyizm = (Real)0.5 * rinv_ixiyizm * (lab(ix-1,iy,iz).ru*lab(ix-1,iy,iz).ru + lab(ix-1,iy,iz).rv*lab(ix-1,iy,iz).rv + lab(ix-1,iy,iz).rw*lab(ix-1,iy,iz).rw);

                    const Real GmixInv_ixiyizpp = (Real)1.0 / (a1_ixiyizpp * g1m1Inv + a2_ixiyizpp * g2m1Inv);
                    const Real GmixInv_ixiyizp = (Real)1.0 / (a1_ixiyizp * g1m1Inv + a2_ixiyizp * g2m1Inv);
                    const Real GmixInv_ixiyizmm = (Real)1.0 / (a1_ixiyizmm * g1m1Inv + a2_ixiyizmm * g2m1Inv);
                    const Real GmixInv_ixiyizm = (Real)1.0 / (a1_ixiyizm * g1m1Inv + a2_ixiyizm * g2m1Inv);

                    const Real Pmix_ixiyizpp = a1_ixiyizpp * pc1_g1m1Inv + a2_ixiyizpp * pc2_g2m1Inv;
                    const Real Pmix_ixiyizp = a1_ixiyizp * pc1_g1m1Inv + a2_ixiyizp * pc2_g2m1Inv;
                    const Real Pmix_ixiyizmm = a1_ixiyizmm * pc1_g1m1Inv + a2_ixiyizmm * pc2_g2m1Inv;
                    const Real Pmix_ixiyizm = a1_ixiyizm * pc1_g1m1Inv + a2_ixiyizm * pc2_g2m1Inv;

                    Real p_ixiyizpp = GmixInv_ixiyizpp * lab(ix+2,iy,iz).energy - GmixInv_ixiyizpp * ke_ixiyizpp - GmixInv_ixiyizpp * Pmix_ixiyizpp;
                    Real p_ixiyizp = GmixInv_ixiyizp * lab(ix+1,iy,iz).energy - GmixInv_ixiyizp * ke_ixiyizp - GmixInv_ixiyizp * Pmix_ixiyizp;
                    Real p_ixiyizmm = GmixInv_ixiyizmm * lab(ix-2,iy,iz).energy - GmixInv_ixiyizmm * ke_ixiyizmm - GmixInv_ixiyizmm * Pmix_ixiyizmm;
                    Real p_ixiyizm = GmixInv_ixiyizm * lab(ix-1,iy,iz).energy - GmixInv_ixiyizm * ke_ixiyizm - GmixInv_ixiyizm * Pmix_ixiyizm;
#ifdef _CONVERTCLIP_
#ifdef _NOK_
                    const Real pth_ixiyizpp = Pmix_ixiyizpp*GmixInv_ixiyizpp/(GmixInv_ixiyizpp+1.0);
                    const Real pth_ixiyizp = Pmix_ixiyizp*GmixInv_ixiyizp/(GmixInv_ixiyizp+1.0);
                    const Real pth_ixiyizmm = Pmix_ixiyizmm*GmixInv_ixiyizmm/(GmixInv_ixiyizmm+1.0);
                    const Real pth_ixiyizm = Pmix_ixiyizm*GmixInv_ixiyizm/(GmixInv_ixiyizm+1.0);
#else
                    const Real pth_ixiyizpp = a2_ixiyizpp == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixiyizpp == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
                    const Real pth_ixiyizp = a2_ixiyizp == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixiyizp == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
                    const Real pth_ixiyizmm = a2_ixiyizmm == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixiyizmm == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
                    const Real pth_ixiyizm = a2_ixiyizm == static_cast<Real>(1.0) ? Simulation_Environment::PC2 : (a2_ixiyizm == static_cast<Real>(0.0) ? Simulation_Environment::PC1 : min(Simulation_Environment::PC1, Simulation_Environment::PC2));
#endif
                    //                   const Real deltap = pthresh < (-static_cast<Real>(2.0)*actpres) ? (-static_cast<Real>(4.0)*actpres - pthresh) : static_cast<Real>(0.0);
                    const Real deltap_ixiyizpp =  pth_ixiyizpp <= (-p_ixiyizpp + static_cast<Real>(PRESEPS)) ? (-p_ixiyizpp - pth_ixiyizpp + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
                    const Real deltap_ixiyizp =  pth_ixiyizp <= (-p_ixiyizp + static_cast<Real>(PRESEPS)) ? (-p_ixiyizp - pth_ixiyizp + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
                    const Real deltap_ixiyizmm =  pth_ixiyizmm <= (-p_ixiyizmm + static_cast<Real>(PRESEPS)) ? (-p_ixiyizmm - pth_ixiyizmm + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
                    const Real deltap_ixiyizm =  pth_ixiyizm <= (-p_ixiyizm + static_cast<Real>(PRESEPS)) ? (-p_ixiyizm - pth_ixiyizm + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);

                    p_ixiyizpp += deltap_ixiyizpp;
                    p_ixiyizp += deltap_ixiyizp;
                    p_ixiyizmm += deltap_ixiyizmm;
                    p_ixiyizm += deltap_ixiyizm;
#endif

                    const Real dadz = -p_ixiyizpp + 8*p_ixiyizp - 8*p_ixiyizm + p_ixiyizmm;

                    o.tmp[iz][iy][ix][4] = mysqrt((dadx*dadx + dady*dady + dadz*dadz)/(144.0*h*h));
                }
    }
};




struct GaussSeidel
{
    StencilInfo stencil;

    int stencil_start[3];
    int stencil_end[3];

#if 1	// peh: disabled for now
    GaussSeidel(): stencil(-2, -2, -2, +3, +3, +3, false, 1, 4)
    {
        stencil_start[0] = stencil_start[1] = stencil_start[2] = -2;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 3;
    }

    GaussSeidel(const GaussSeidel& c): stencil(-2, -2, -2, +3, +3, +3, false, 1, 4)
    {
        stencil_start[0] = stencil_start[1] = stencil_start[2] = -2;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 3;
    }
#else
    GaussSeidel(): stencil(-3, -3, -3, +4, +4, +4, false, 1, 4)
    {
        stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
    }

    GaussSeidel(const GaussSeidel& c): stencil(-3, -3, -3, +4, +4, +4, false, 1, 4)
    {
        stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
    }
#endif

    template<typename TLab, typename TBlock>
    inline void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    Real p[15];
                    const int span = 5;

                    for(int dir=0;dir<3;dir++)
                        for(int i=0;i<span;i++)
                            p[i+dir*span] = lab(ix+ (dir==0? i-2 : 0 ), iy+ (dir==1? i-2 : 0), iz+ (dir==2? i-2 : 0)).energy;

                    Real pressure_new = 0;
                    for(int dir=0;dir<3;dir++)
                        pressure_new += p[dir*span+1]+p[dir*span+3];

                    const double G1 = 1.0/(Simulation_Environment::GAMMA1-1);
                    const double G2 = 1.0/(Simulation_Environment::GAMMA2-1);

                    o(ix,iy,iz).dummy = (lab(ix,iy,iz).G < 0.5*(G1+G2)) ? pressure_new/(Real)6.0 : lab(ix,iy,iz).energy;
                }
    }
};

#endif /* 0 */

#endif /* VECTOROPERATOR_H_UBJ0NRPP */
