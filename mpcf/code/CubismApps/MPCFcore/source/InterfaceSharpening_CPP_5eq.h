/*
 *  InterfaceSharpening_CPP_5eq.h
 *  MPCFcore
 *
 *  Created by Ursula Rasthofer
 *  Date: October 2015
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

//-------------------------------------------------
// this imlementation ussed corner gradients as applied for the diffusion tensor
//-------------------------------------------------

#ifndef INTERFACESHARPENING_CPP_5EQ_H
#define INTERFACESHARPENING_CPP_5EQ_H

#include "SOA2D.h"

typedef SOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1> InputSOA_IS;
typedef RingSOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1, 2> RingInputSOA_IS;
typedef RingSOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1, 3> RingInputSOA_IS_ex;
typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_+1> TempSOA_IS;
typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_> TempGrXSOA_IS;
typedef SOA2D<0,_BLOCKSIZE_, 0,_BLOCKSIZE_+1> TempGrYSOA_IS;
typedef SOA2D<0,_BLOCKSIZE_, 0,_BLOCKSIZE_> TempGrZSOA_IS;
typedef RingSOA2D<0, _BLOCKSIZE_, 0,_BLOCKSIZE_, 2> RingTempGrZSOA_IS;
typedef RingSOA2D<0, _BLOCKSIZE_+1, 0,_BLOCKSIZE_+1, 2> RingTempSOA_IS;

class InterfaceSharpening_CPP_5eq
{
public:

  const Real a; //factor that affects the rhs, related to low-storage rk3 (see _copyback method)
  const Real dtinvh; //time-step length
  const Real epsilon, U0; //user defined interface thickness and regularization parameter
  const Real _g1, _g2, _pc1, _pc2; // gamma and correction pressure
  Real G1, G2, P1, P2; //Gamma and Pi for both fluids
  Real diff_G, diff_P; //difference of Gamma and Pi

  InterfaceSharpening_CPP_5eq(const Real ina = 1, const Real indtinvh = 1, const Real inepsilon=1,
                              const Real inU0=1, const Real inG1=1.4, const Real inG2=1.4, const Real inP1=0, const Real inP2=0):
  a(ina), dtinvh(indtinvh), epsilon(inepsilon*static_cast<Real>(0.0625)), U0(inU0), _g1(inG1), _g2(inG2), _pc1(inP1), _pc2(inP2)
  {
    G1 = static_cast<Real>(1.0)/(_g1-static_cast<Real>(1.0));
    G2 = static_cast<Real>(1.0)/(_g2-static_cast<Real>(1.0));
    P1 = _g1*_pc1/(_g1-static_cast<Real>(1.0));
    P2 = _g2*_pc2/(_g2-static_cast<Real>(1.0));
    diff_G = G2-G1;
    diff_P = P2-P1;
    return; 
  }

  //main loop, processing the z-slices in sequence. all the logic is in here
  void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
               Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts);

  static void printflops(const float PEAKPERF_CORE, const float PEAKBAND,
                         const int NCORES, const int NTIMES, const int NBLOCKS,
                         const float MEASUREDTIME) {std::cout << "No printflops for InterfaceSharpening_CPP, yet" << std::endl; return;}

protected:

  struct AssumedType { Real a1r1, a2r2, ru, rv, rw, e, a2, dummy;};

  RingInputSOA_IS ringu, ringv, ringw, ringp; //slices for working variables > values for u, v, w, and p may be reduced to one slice
  RingInputSOA_IS_ex ringa2, ringa1r1, ringa2r2;
  OutputSOA rhsa1r1, rhsa2r2, rhsru, rhsrv, rhsrw, rhse, rhsa2; //output slices
  OutputSOA nx, ny, nz; //normal vector cell
  TempGrXSOA_IS vxxa2, vxya2, vxza2, vxxa1r1, vxya1r1, vxza1r1, vxxa2r2, vxya2r2, vxza2r2;
  TempGrYSOA_IS vyxa2, vyya2, vyza2, vyxa1r1, vyya1r1, vyza1r1, vyxa2r2, vyya2r2, vyza2r2;
  RingTempGrZSOA_IS ringvzxa2, ringvzya2, ringvzza2, ringvzxa1r1, ringvzya1r1, ringvzza1r1, ringvzxa2r2, ringvzya2r2, ringvzza2r2;
  RingTempSOA_IS ringgrada2x, ringgrada2y, ringgrada2z, ringgrada1r1x, ringgrada1r1y, ringgrada1r1z, ringgrada2r2x, ringgrada2r2y, ringgrada2r2z; // for corner values

  void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);

  void _input_next() { ringa2.next(); ringa1r1.next(); ringa2r2.next();
                       ringu.next(); ringv.next(); ringw.next(); ringp.next(); return; }

  void _gradient_next() { ringvzxa2.next(); ringvzya2.next(); ringvzza2.next(); 
                          ringvzxa1r1.next(); ringvzya1r1.next(); ringvzza1r1.next();
                          ringvzxa2r2.next(); ringvzya2r2.next(); ringvzza2r2.next(); return; }

  void _corner_gradient_next() { ringgrada2x.next(); ringgrada2y.next(); ringgrada2z.next();
                                 ringgrada1r1x.next(); ringgrada1r1y.next(); ringgrada1r1z.next();
                                 ringgrada2r2x.next(); ringgrada2r2y.next(); ringgrada2r2z.next(); return; }

  void _corners(const InputSOA_IS& vf00, const InputSOA_IS& vf01, TempSOA_IS& gf0x, TempSOA_IS& gf0y, TempSOA_IS& gf0z,
                const InputSOA_IS& vf10, const InputSOA_IS& vf11, TempSOA_IS& gf1x, TempSOA_IS& gf1y, TempSOA_IS& gf1z,
                const InputSOA_IS& vf20, const InputSOA_IS& vf21, TempSOA_IS& gf2x, TempSOA_IS& gf2y, TempSOA_IS& gf2z);

  void _gradient_xface(const TempSOA_IS& a2x0, const TempSOA_IS& a2y0, const TempSOA_IS& a2z0,const TempSOA_IS& a2x1, const TempSOA_IS& a2y1, const TempSOA_IS& a2z1,
                       const TempSOA_IS& a1r1x0, const TempSOA_IS& a1r1y0, const TempSOA_IS& a1r1z0,const TempSOA_IS& a1r1x1, const TempSOA_IS& a1r1y1, const TempSOA_IS& a1r1z1,
                       const TempSOA_IS& a2r2x0, const TempSOA_IS& a2r2y0, const TempSOA_IS& a2r2z0,const TempSOA_IS& a2r2x1, const TempSOA_IS& a2r2y1, const TempSOA_IS& a2r2z1);

  void _gradient_yface(const TempSOA_IS& a2x0, const TempSOA_IS& a2y0, const TempSOA_IS& a2z0,const TempSOA_IS& a2x1, const TempSOA_IS& a2y1, const TempSOA_IS& a2z1,
                       const TempSOA_IS& a1r1x0, const TempSOA_IS& a1r1y0, const TempSOA_IS& a1r1z0,const TempSOA_IS& a1r1x1, const TempSOA_IS& a1r1y1, const TempSOA_IS& a1r1z1,
                       const TempSOA_IS& a2r2x0, const TempSOA_IS& a2r2y0, const TempSOA_IS& a2r2z0,const TempSOA_IS& a2r2x1, const TempSOA_IS& a2r2y1, const TempSOA_IS& a2r2z1);

  void _gradient_zface(const TempSOA_IS& a2x, const TempSOA_IS& a2y, const TempSOA_IS& a2z,
                       const TempSOA_IS& a1r1x, const TempSOA_IS& a1r1y, const TempSOA_IS& a1r1z,
                       const TempSOA_IS& a2r2x, const TempSOA_IS& a2r2y, const TempSOA_IS& a2r2z,
                       TempGrZSOA_IS& vzxa2, TempGrZSOA_IS& vzya2, TempGrZSOA_IS& vzza2,
                       TempGrZSOA_IS& vzxa1r1, TempGrZSOA_IS& vzya1r1, TempGrZSOA_IS& vzza1r1,
                       TempGrZSOA_IS& vzxa2r2, TempGrZSOA_IS& vzya2r2, TempGrZSOA_IS& vzza2r2);

  void _normal_cell(const InputSOA_IS& a2m1, const InputSOA_IS& a2, const InputSOA_IS& a2p1);

  void _rhs(const InputSOA_IS& a2m1, const InputSOA_IS& a1r1m1, const InputSOA_IS& a2r2m1,
            const InputSOA_IS& a2, const InputSOA_IS& a1r1, const InputSOA_IS& a2r2,
            const InputSOA_IS& a2p1, const InputSOA_IS& a1r1p1, const InputSOA_IS& a2r2p1,
            const InputSOA_IS& u, const InputSOA_IS& v, const InputSOA_IS& w, const InputSOA_IS& p,
            const TempGrZSOA_IS& vzxa2_0, const TempGrZSOA_IS& vzya2_0, const TempGrZSOA_IS& vzza2_0,
            const TempGrZSOA_IS& vzxa1r1_0, const TempGrZSOA_IS& vzya1r1_0, const TempGrZSOA_IS& vzza1r1_0,
            const TempGrZSOA_IS& vzxa2r2_0, const TempGrZSOA_IS& vzya2r2_0, const TempGrZSOA_IS& vzza2r2_0,
            const TempGrZSOA_IS& vzxa2_1, const TempGrZSOA_IS& vzya2_1, const TempGrZSOA_IS& vzza2_1,
            const TempGrZSOA_IS& vzxa1r1_1, const TempGrZSOA_IS& vzya1r1_1, const TempGrZSOA_IS& vzza1r1_1,
            const TempGrZSOA_IS& vzxa2r2_1, const TempGrZSOA_IS& vzya2r2_1, const TempGrZSOA_IS& vzza2r2_1);

  void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);
};

#endif
