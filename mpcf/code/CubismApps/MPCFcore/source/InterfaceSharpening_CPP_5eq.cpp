/*
 *  InterfaceSharpening_CPP_5eq.cpp
 *  MPCFcore
 *
 *  Created by Ursula Rasthofer
 *  Date: October 2015
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#include "InterfaceSharpening_CPP_5eq.h"

void InterfaceSharpening_CPP_5eq::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
  InputSOA_IS& a2=ringa2.ref(), &a1r1=ringa1r1.ref(), &a2r2=ringa2r2.ref();
  InputSOA_IS& u=ringu.ref(), &v=ringv.ref(), &w=ringw.ref(), &p=ringp.ref();

  for(int sy=0; sy<_BLOCKSIZE_+2; sy++)
    for(int sx=0; sx<_BLOCKSIZE_+2; sx++)
    {
      AssumedType pt = *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));

      const int dx = sx-1;
      const int dy = sy-1;

#ifdef _CONVERTCLIP_

      const Real a1r1clip = max(static_cast<Real>(0.0),pt.a1r1);
      const Real a2r2clip = max(static_cast<Real>(0.0),pt.a2r2);
      const Real alpha2 = max(static_cast<Real>(ALPHAEPS), min(pt.a2, static_cast<Real>(1.0-ALPHAEPS)));

      a1r1.ref(dx, dy) = a1r1clip;
      a2r2.ref(dx, dy) = a2r2clip;
      a2.ref(dx, dy)   = alpha2;

      const Real invrho = static_cast<Real>(1.0) / (a1r1clip + a2r2clip);
      // const Real alpha1 = max(static_cast<Real>(0.0), min(static_cast<Real>(1.0)-alpha2, static_cast<Real>(1.0)));
      const Real alpha1 = static_cast<Real>(1.0) - alpha2;

      // simple conversion to primitive variables
      u.ref(dx, dy) = pt.ru * invrho;
      v.ref(dx, dy) = pt.rv * invrho;
      w.ref(dx, dy) = pt.rw * invrho;

      const Real ke = static_cast<Real>(0.5) * (pt.ru*pt.ru + pt.rv*pt.rv + pt.rw*pt.rw) * invrho;
      const Real GmixInv = static_cast<Real>(1.0) / (alpha1*G1 + alpha2*G2);
      const Real Pmix = max(static_cast<Real>(0.0),alpha1*P1 + alpha2*P2);

      const Real actpres = GmixInv*pt.e - GmixInv*ke - GmixInv*Pmix;

#ifdef _NOK_
      const Real pthresh = Pmix*GmixInv/(GmixInv+1.0);
#else
      // TODO: (fabianw@mavt.ethz.ch; Mon 02 May 2016 02:41:56 PM CEST) What to
      // do here if ALPHAEPS is used?
      const Real pthresh = alpha2 == static_cast<Real>(1.0) ? _pc2 : (alpha2 == static_cast<Real>(0.0) ? _pc1 : min(_pc1, _pc2));
#endif
      //const Real deltap = pthresh < (-static_cast<Real>(4.0)*actpres + static_cast<Real>(PRESEPS)) ? (-static_cast<Real>(6.0)*actpres - pthresh + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
      const Real deltap =  pthresh <= (-actpres + static_cast<Real>(PRESEPS)) ? (-actpres - pthresh + static_cast<Real>(PRESEPS)) : static_cast<Real>(0.0);
/*
#ifdef _CLIPVERBOSE_
      if (deltap > static_cast<Real>(0.0))
      {
        cout << "Pressure correction in SHARPENING convert " << setprecision(12) << actpres << " " << deltap << " " << pthresh << " " << alpha1 << " " << alpha2  << endl;
      }
#endif
*/
      const Real pressure = actpres + deltap;
/*
      const Real pressure = actpres;
*/
      p.ref(dx,dy) = pressure;

#else

    const Real rInv = static_cast<Real>(1.0) / (pt.a1r1 + pt.a2r2);
#ifdef _ALPHACLIP_
    const Real alpha2 = max(static_cast<Real>(ALPHAEPS), min(pt.a2, static_cast<Real>(1.0-ALPHAEPS)));
#else
    const Real alpha2 = pt.a2;
#endif /* _ALPHACLIP_ */
    const Real alpha1 = static_cast<Real>(1.0) - alpha2;
    const Real GmixInv = static_cast<Real>(1.0) / (alpha1*G1 + alpha2*G2);
    const Real Pmix  = alpha1*P1 + alpha2*P2;
    const Real pres = GmixInv*pt.e - GmixInv*static_cast<Real>(0.5)*rInv*(pt.ru*pt.ru + pt.rv*pt.rv + pt.rw*pt.rw) - GmixInv*Pmix;

    a1r1.ref(dx, dy) = pt.a1r1;
    a2r2.ref(dx, dy) = pt.a2r2;
    a2.ref(dx, dy) = pt.a2;
    u.ref(dx, dy) = pt.ru*rInv;
    v.ref(dx, dy) = pt.rv*rInv;
    w.ref(dx, dy) = pt.rw*rInv;
    p.ref(dx,dy) = pres;

#endif

    }

  return;
}

void InterfaceSharpening_CPP_5eq::compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
                                          Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
{
  _convert(srcfirst, srcfloats, rowsrcs);
  _input_next();
  _convert(srcfirst + srcfloats*slicesrcs, srcfloats, rowsrcs);

  _corners(ringa2(-1), ringa2(), ringgrada2x.ref(), ringgrada2y.ref(), ringgrada2z.ref(),
           ringa1r1(-1), ringa1r1(), ringgrada1r1x.ref(), ringgrada1r1y.ref(), ringgrada1r1z.ref(),
           ringa2r2(-1), ringa2r2(), ringgrada2r2x.ref(), ringgrada2r2y.ref(), ringgrada2r2z.ref());

  _gradient_zface(ringgrada2x(), ringgrada2y(), ringgrada2z(),
                  ringgrada1r1x(), ringgrada1r1y(), ringgrada1r1z(),
                  ringgrada2r2x(), ringgrada2r2y(), ringgrada2r2z(),
                  ringvzxa2.ref(), ringvzya2.ref(), ringvzza2.ref(),
                  ringvzxa1r1.ref(), ringvzya1r1.ref(), ringvzza1r1.ref(),
                  ringvzxa2r2.ref(), ringvzya2r2.ref(), ringvzza2r2.ref());

  for (int islice=0; islice<_BLOCKSIZE_; islice++)
  {
    _gradient_next();
    _input_next();
    _corner_gradient_next();

    _convert(srcfirst + (islice+2)*srcfloats*slicesrcs, srcfloats, rowsrcs);

    _corners(ringa2(-1), ringa2(0), ringgrada2x.ref(), ringgrada2y.ref(), ringgrada2z.ref(),
             ringa1r1(-1), ringa1r1(0), ringgrada1r1x.ref(), ringgrada1r1y.ref(), ringgrada1r1z.ref(),
             ringa2r2(-1), ringa2r2(0), ringgrada2r2x.ref(), ringgrada2r2y.ref(), ringgrada2r2z.ref());


    _gradient_xface(ringgrada2x(-1), ringgrada2y(-1), ringgrada2z(-1), ringgrada2x(), ringgrada2y(), ringgrada2z(),
                    ringgrada1r1x(-1), ringgrada1r1y(-1), ringgrada1r1z(-1), ringgrada1r1x(), ringgrada1r1y(), ringgrada1r1z(),
                    ringgrada2r2x(-1), ringgrada2r2y(-1), ringgrada2r2z(-1), ringgrada2r2x(), ringgrada2r2y(), ringgrada2r2z());
    _gradient_yface(ringgrada2x(-1), ringgrada2y(-1), ringgrada2z(-1), ringgrada2x(), ringgrada2y(), ringgrada2z(),
                    ringgrada1r1x(-1), ringgrada1r1y(-1), ringgrada1r1z(-1), ringgrada1r1x(), ringgrada1r1y(), ringgrada1r1z(),
                    ringgrada2r2x(-1), ringgrada2r2y(-1), ringgrada2r2z(-1), ringgrada2r2x(), ringgrada2r2y(), ringgrada2r2z());
    _gradient_zface(ringgrada2x(), ringgrada2y(), ringgrada2z(),
                    ringgrada1r1x(), ringgrada1r1y(), ringgrada1r1z(),
                    ringgrada2r2x(), ringgrada2r2y(), ringgrada2r2z(),
                    ringvzxa2.ref(), ringvzya2.ref(), ringvzza2.ref(),
                    ringvzxa1r1.ref(), ringvzya1r1.ref(), ringvzza1r1.ref(),
                    ringvzxa2r2.ref(), ringvzya2r2.ref(), ringvzza2r2.ref());

    _normal_cell(ringa2(-2), ringa2(-1), ringa2(0));

    _rhs(ringa2.ref(-2), ringa1r1.ref(-2), ringa2r2.ref(-2), ringa2.ref(-1), ringa1r1.ref(-1), ringa2r2.ref(-1),
         ringa2.ref(0), ringa1r1.ref(0), ringa2r2.ref(0), ringu.ref(-1), ringv.ref(-1), ringw.ref(-1), ringp.ref(-1),
         ringvzxa2.ref(-1), ringvzya2.ref(-1), ringvzza2.ref(-1), ringvzxa1r1.ref(-1), ringvzya1r1.ref(-1), ringvzza1r1.ref(-1),
         ringvzxa2r2.ref(-1), ringvzya2r2.ref(-1), ringvzza2r2.ref(-1), ringvzxa2.ref(0), ringvzya2.ref(0), ringvzza2.ref(0),
         ringvzxa1r1.ref(0), ringvzya1r1.ref(0), ringvzza1r1.ref(0), ringvzxa2r2.ref(0), ringvzya2r2.ref(0), ringvzza2r2.ref(0));

    _copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);
  }

  return;
}


void InterfaceSharpening_CPP_5eq::_corners(const InputSOA_IS& vf00, const InputSOA_IS& vf01, TempSOA_IS& gf0x, TempSOA_IS& gf0y, TempSOA_IS& gf0z,
                                           const InputSOA_IS& vf10, const InputSOA_IS& vf11, TempSOA_IS& gf1x, TempSOA_IS& gf1y, TempSOA_IS& gf1z,
                                           const InputSOA_IS& vf20, const InputSOA_IS& vf21, TempSOA_IS& gf2x, TempSOA_IS& gf2y, TempSOA_IS& gf2z)
{
  for(int iy=0; iy<TempSOA_IS::NY; iy++)
    for(int ix=0; ix<TempSOA_IS::NX; ix++)
    {
      Real c000    = vf01(ix,iy);
      Real cm100   = vf01(ix-1,iy);
      Real c0m10   = vf01(ix,iy-1);
      Real cm1m10  = vf01(ix-1,iy-1);
      Real c00m1   = vf00(ix,iy);
      Real cm10m1  = vf00(ix-1,iy);
      Real c0m1m1  = vf00(ix,iy-1);
      Real cm1m1m1 = vf00(ix-1,iy-1);

      gf0x.ref(ix, iy) = c000-cm100 + c0m10-cm1m10 + c00m1-cm10m1 + c0m1m1-cm1m1m1;
      gf0y.ref(ix, iy) = cm100-cm1m10 + c000-c0m10 + cm10m1-cm1m1m1  + c00m1-c0m1m1;
      gf0z.ref(ix, iy) = cm1m10-cm1m1m1 + cm100-cm10m1 + c0m10-c0m1m1  + c000-c00m1;

      c000    = vf11(ix,iy);
      cm100   = vf11(ix-1,iy);
      c0m10   = vf11(ix,iy-1);
      cm1m10  = vf11(ix-1,iy-1);
      c00m1   = vf10(ix,iy);
      cm10m1  = vf10(ix-1,iy);
      c0m1m1  = vf10(ix,iy-1);
      cm1m1m1 = vf10(ix-1,iy-1);

      gf1x.ref(ix, iy) = c000-cm100 + c0m10-cm1m10 + c00m1-cm10m1 + c0m1m1-cm1m1m1;
      gf1y.ref(ix, iy) = cm100-cm1m10 + c000-c0m10 + cm10m1-cm1m1m1  + c00m1-c0m1m1;
      gf1z.ref(ix, iy) = cm1m10-cm1m1m1 + cm100-cm10m1 + c0m10-c0m1m1  + c000-c00m1;

      c000    = vf21(ix,iy);
      cm100   = vf21(ix-1,iy);
      c0m10   = vf21(ix,iy-1);
      cm1m10  = vf21(ix-1,iy-1);
      c00m1   = vf20(ix,iy);
      cm10m1  = vf20(ix-1,iy);
      c0m1m1  = vf20(ix,iy-1);
      cm1m1m1 = vf20(ix-1,iy-1);

      gf2x.ref(ix, iy) = c000-cm100 + c0m10-cm1m10 + c00m1-cm10m1 + c0m1m1-cm1m1m1;
      gf2y.ref(ix, iy) = cm100-cm1m10 + c000-c0m10 + cm10m1-cm1m1m1  + c00m1-c0m1m1;
      gf2z.ref(ix, iy) = cm1m10-cm1m1m1 + cm100-cm10m1 + c0m10-c0m1m1  + c000-c00m1;
    }

  return;
}

void InterfaceSharpening_CPP_5eq::_gradient_xface(const TempSOA_IS& a2x0, const TempSOA_IS& a2y0, const TempSOA_IS& a2z0, const TempSOA_IS& a2x1, const TempSOA_IS& a2y1, const TempSOA_IS& a2z1,
                                                  const TempSOA_IS& a1r1x0, const TempSOA_IS& a1r1y0, const TempSOA_IS& a1r1z0, const TempSOA_IS& a1r1x1, const TempSOA_IS& a1r1y1, const TempSOA_IS& a1r1z1,
                                                  const TempSOA_IS& a2r2x0, const TempSOA_IS& a2r2y0, const TempSOA_IS& a2r2z0, const TempSOA_IS& a2r2x1, const TempSOA_IS& a2r2y1, const TempSOA_IS& a2r2z1)
{
  // missing factor 1/(16h) = 1/4 (averaging corner values) / h (missing in corner values) * 1/4 (present average)  is included in rhs function below

  for (int iy=0; iy<TempGrXSOA_IS::NY; iy++)
    for(int ix=0; ix<TempGrXSOA_IS::NX; ix++)
    {
      vxxa2.ref(ix,iy) = (a2x0(ix,iy)+a2x0(ix,iy+1)+a2x1(ix,iy)+a2x1(ix,iy+1));
      vxya2.ref(ix,iy) = (a2y0(ix,iy)+a2y0(ix,iy+1)+a2y1(ix,iy)+a2y1(ix,iy+1));
      vxza2.ref(ix,iy) = (a2z0(ix,iy)+a2z0(ix,iy+1)+a2z1(ix,iy)+a2z1(ix,iy+1));

      vxxa1r1.ref(ix,iy) = (a1r1x0(ix,iy)+a1r1x0(ix,iy+1)+a1r1x1(ix,iy)+a1r1x1(ix,iy+1));
      vxya1r1.ref(ix,iy) = (a1r1y0(ix,iy)+a1r1y0(ix,iy+1)+a1r1y1(ix,iy)+a1r1y1(ix,iy+1));
      vxza1r1.ref(ix,iy) = (a1r1z0(ix,iy)+a1r1z0(ix,iy+1)+a1r1z1(ix,iy)+a1r1z1(ix,iy+1));

      vxxa2r2.ref(ix,iy) = (a2r2x0(ix,iy)+a2r2x0(ix,iy+1)+a2r2x1(ix,iy)+a2r2x1(ix,iy+1));
      vxya2r2.ref(ix,iy) = (a2r2y0(ix,iy)+a2r2y0(ix,iy+1)+a2r2y1(ix,iy)+a2r2y1(ix,iy+1));
      vxza2r2.ref(ix,iy) = (a2r2z0(ix,iy)+a2r2z0(ix,iy+1)+a2r2z1(ix,iy)+a2r2z1(ix,iy+1));
   }

  return;
}


void InterfaceSharpening_CPP_5eq::_gradient_yface(const TempSOA_IS& a2x0, const TempSOA_IS& a2y0, const TempSOA_IS& a2z0,const TempSOA_IS& a2x1, const TempSOA_IS& a2y1, const TempSOA_IS& a2z1,
                                                  const TempSOA_IS& a1r1x0, const TempSOA_IS& a1r1y0, const TempSOA_IS& a1r1z0,const TempSOA_IS& a1r1x1, const TempSOA_IS& a1r1y1, const TempSOA_IS& a1r1z1,
                                                  const TempSOA_IS& a2r2x0, const TempSOA_IS& a2r2y0, const TempSOA_IS& a2r2z0,const TempSOA_IS& a2r2x1, const TempSOA_IS& a2r2y1, const TempSOA_IS& a2r2z1)
{
  // see xface for missing factor

  for (int iy=0; iy<TempGrYSOA_IS::NY; iy++)
    for (int ix=0; ix<TempGrYSOA_IS::NX; ix++)
    {
      vyxa2.ref(ix,iy) = (a2x0(ix,iy)+a2x0(ix+1,iy)+a2x1(ix,iy)+a2x1(ix+1,iy));
      vyya2.ref(ix,iy) = (a2y0(ix,iy)+a2y0(ix+1,iy)+a2y1(ix,iy)+a2y1(ix+1,iy));
      vyza2.ref(ix,iy) = (a2z0(ix,iy)+a2z0(ix+1,iy)+a2z1(ix,iy)+a2z1(ix+1,iy));

      vyxa1r1.ref(ix,iy) = (a1r1x0(ix,iy)+a1r1x0(ix+1,iy)+a1r1x1(ix,iy)+a1r1x1(ix+1,iy));
      vyya1r1.ref(ix,iy) = (a1r1y0(ix,iy)+a1r1y0(ix+1,iy)+a1r1y1(ix,iy)+a1r1y1(ix+1,iy));
      vyza1r1.ref(ix,iy) = (a1r1z0(ix,iy)+a1r1z0(ix+1,iy)+a1r1z1(ix,iy)+a1r1z1(ix+1,iy));

      vyxa2r2.ref(ix,iy) = (a2r2x0(ix,iy)+a2r2x0(ix+1,iy)+a2r2x1(ix,iy)+a2r2x1(ix+1,iy));
      vyya2r2.ref(ix,iy) = (a2r2y0(ix,iy)+a2r2y0(ix+1,iy)+a2r2y1(ix,iy)+a2r2y1(ix+1,iy));
      vyza2r2.ref(ix,iy) = (a2r2z0(ix,iy)+a2r2z0(ix+1,iy)+a2r2z1(ix,iy)+a2r2z1(ix+1,iy));
   }

  return;
}


void InterfaceSharpening_CPP_5eq::_gradient_zface(const TempSOA_IS& a2x, const TempSOA_IS& a2y, const TempSOA_IS& a2z,
                                                  const TempSOA_IS& a1r1x, const TempSOA_IS& a1r1y, const TempSOA_IS& a1r1z,
                                                  const TempSOA_IS& a2r2x, const TempSOA_IS& a2r2y, const TempSOA_IS& a2r2z,
                                                  TempGrZSOA_IS& vzxa2, TempGrZSOA_IS& vzya2, TempGrZSOA_IS& vzza2,
                                                  TempGrZSOA_IS& vzxa1r1, TempGrZSOA_IS& vzya1r1, TempGrZSOA_IS& vzza1r1,
                                                  TempGrZSOA_IS& vzxa2r2, TempGrZSOA_IS& vzya2r2, TempGrZSOA_IS& vzza2r2)
{
  // see xface for missing factor

  for (int iy=0; iy<TempGrZSOA_IS::NY; iy++)
    for (int ix=0; ix<TempGrZSOA_IS::NX; ix++)
    {
      vzxa2.ref(ix,iy) = (a2x(ix,iy)+a2x(ix+1,iy)+a2x(ix,iy+1)+a2x(ix+1,iy+1));
      vzya2.ref(ix,iy) = (a2y(ix,iy)+a2y(ix+1,iy)+a2y(ix,iy+1)+a2y(ix+1,iy+1));
      vzza2.ref(ix,iy) = (a2z(ix,iy)+a2z(ix+1,iy)+a2z(ix,iy+1)+a2z(ix+1,iy+1));

      vzxa1r1.ref(ix,iy) = (a1r1x(ix,iy)+a1r1x(ix+1,iy)+a1r1x(ix,iy+1)+a1r1x(ix+1,iy+1));
      vzya1r1.ref(ix,iy) = (a1r1y(ix,iy)+a1r1y(ix+1,iy)+a1r1y(ix,iy+1)+a1r1y(ix+1,iy+1));
      vzza1r1.ref(ix,iy) = (a1r1z(ix,iy)+a1r1z(ix+1,iy)+a1r1z(ix,iy+1)+a1r1z(ix+1,iy+1));

      vzxa2r2.ref(ix,iy) = (a2r2x(ix,iy)+a2r2x(ix+1,iy)+a2r2x(ix,iy+1)+a2r2x(ix+1,iy+1));
      vzya2r2.ref(ix,iy) = (a2r2y(ix,iy)+a2r2y(ix+1,iy)+a2r2y(ix,iy+1)+a2r2y(ix+1,iy+1));
      vzza2r2.ref(ix,iy) = (a2r2z(ix,iy)+a2r2z(ix+1,iy)+a2r2z(ix,iy+1)+a2r2z(ix+1,iy+1));
    }

  return;
}

void InterfaceSharpening_CPP_5eq::_normal_cell(const InputSOA_IS& a2m1, const InputSOA_IS& a2, const InputSOA_IS& a2p1)
{
  for(int iy=0; iy<OutputSOA::NY; iy++)
    for(int ix=0; ix<OutputSOA::NX; ix++)
    {
      nx.ref(ix,iy) = (a2(ix+1,iy) - a2(ix-1,iy));
      ny.ref(ix,iy) = (a2(ix,iy+1) - a2(ix,iy-1));
      nz.ref(ix,iy) = (a2p1(ix,iy) - a2m1(ix,iy));
    }

  //nomalization by mag(n) in rhs function below

  return;
}


void InterfaceSharpening_CPP_5eq::_rhs(const InputSOA_IS& a2m1, const InputSOA_IS& a1r1m1, const InputSOA_IS& a2r2m1,
                                       const InputSOA_IS& a2, const InputSOA_IS& a1r1, const InputSOA_IS& a2r2,
                                       const InputSOA_IS& a2p1, const InputSOA_IS& a1r1p1, const InputSOA_IS& a2r2p1,
                                       const InputSOA_IS& u, const InputSOA_IS& v, const InputSOA_IS& w, const InputSOA_IS& p,
                                       const TempGrZSOA_IS& vzxa2_0, const TempGrZSOA_IS& vzya2_0, const TempGrZSOA_IS& vzza2_0,
                                       const TempGrZSOA_IS& vzxa1r1_0, const TempGrZSOA_IS& vzya1r1_0, const TempGrZSOA_IS& vzza1r1_0,
                                       const TempGrZSOA_IS& vzxa2r2_0, const TempGrZSOA_IS& vzya2r2_0, const TempGrZSOA_IS& vzza2r2_0,
                                       const TempGrZSOA_IS& vzxa2_1, const TempGrZSOA_IS& vzya2_1, const TempGrZSOA_IS& vzza2_1,
                                       const TempGrZSOA_IS& vzxa1r1_1, const TempGrZSOA_IS& vzya1r1_1, const TempGrZSOA_IS& vzza1r1_1,
                                       const TempGrZSOA_IS& vzxa2r2_1, const TempGrZSOA_IS& vzya2r2_1, const TempGrZSOA_IS& vzza2r2_1)
{
  Real inv_mag_h;
  Real mag_xm, mag_xp, mag_ym, mag_yp, mag_zm, mag_zp;
  Real eps_mag_xm, eps_mag_xp, eps_mag_ym, eps_mag_yp, eps_mag_zm, eps_mag_zp;
  Real f_xm_a1r1, f_xp_a1r1, f_ym_a1r1, f_yp_a1r1, f_zm_a1r1, f_zp_a1r1;
  Real f_xm_a2r2, f_xp_a2r2, f_ym_a2r2, f_yp_a2r2, f_zm_a2r2, f_zp_a2r2;
  Real sum_R, kappa;

  // missing factor in calculation of gradients is incorporated here (1/(16h))
  // const Real fac = static_cast<Real>(0.0625)/h;
  // as eps=epsilon*h defined in Tiwari et al 2013 is multplied by this factor,
  // h cancels and epsilon is the only remaining factor
  // epsilon is multiplied by static_cast<Real>(0.0625) in the constructor
  // factor rhs including dt
  const Real facdt = U0*dtinvh;

  for(int iy=0; iy<OutputSOA::NY; iy++)
    for(int ix=0; ix<OutputSOA::NX; ix++)
    {
      // TODO: is there a smarter way?
      if ((a2(ix,iy) > static_cast<Real>(1.0e-6)) && (a2(ix,iy) < static_cast<Real>(1.0-1.0e-6)))
      {
        const Real mag_n = mysqrt(nx(ix,iy)*nx(ix,iy) + ny(ix,iy)*ny(ix,iy) + nz(ix,iy)*nz(ix,iy));

        if (mag_n > static_cast<Real>(1.0e-9))
          inv_mag_h = facdt/mag_n;
        else inv_mag_h = static_cast<Real>(0.0);

        assert(!isnan(inv_mag_h));
        assert(!isinf(inv_mag_h));

        mag_xm = mysqrt(vxxa2(ix,iy)*vxxa2(ix,iy) + vxya2(ix,iy)*vxya2(ix,iy) + vxza2(ix,iy)*vxza2(ix,iy));
        mag_xp = mysqrt(vxxa2(ix+1,iy)*vxxa2(ix+1,iy) + vxya2(ix+1,iy)*vxya2(ix+1,iy) + vxza2(ix+1,iy)*vxza2(ix+1,iy));
        mag_ym = mysqrt(vyxa2(ix,iy)*vyxa2(ix,iy) + vyya2(ix,iy)*vyya2(ix,iy) + vyza2(ix,iy)*vyza2(ix,iy));
        mag_yp = mysqrt(vyxa2(ix,iy+1)*vyxa2(ix,iy+1) + vyya2(ix,iy+1)*vyya2(ix,iy+1) + vyza2(ix,iy+1)*vyza2(ix,iy+1));
        mag_zm = mysqrt(vzxa2_0(ix,iy)*vzxa2_0(ix,iy) + vzya2_0(ix,iy)*vzya2_0(ix,iy) + vzza2_0(ix,iy)*vzza2_0(ix,iy));
        mag_zp = mysqrt(vzxa2_1(ix,iy)*vzxa2_1(ix,iy) + vzya2_1(ix,iy)*vzya2_1(ix,iy) + vzza2_1(ix,iy)*vzza2_1(ix,iy));


        if (mag_xm > static_cast<Real>(1.0e-9))
          eps_mag_xm = epsilon/mag_xm;
        else eps_mag_xm = static_cast<Real>(0.0);
        if (mag_xp > static_cast<Real>(1.0e-9))
          eps_mag_xp = epsilon/mag_xp;
        else eps_mag_xp = static_cast<Real>(0.0);
        if (mag_ym > static_cast<Real>(1.0e-9))
          eps_mag_ym = epsilon/mag_ym;
        else eps_mag_ym = static_cast<Real>(0.0);
        if (mag_yp > static_cast<Real>(1.0e-9))
          eps_mag_yp = epsilon/mag_yp;
        else eps_mag_yp = static_cast<Real>(0.0);
        if (mag_zm > static_cast<Real>(1.0e-9))
          eps_mag_zm = epsilon/mag_zm;
        else eps_mag_zm = static_cast<Real>(0.0);
        if (mag_zp > static_cast<Real>(1.0e-9))
          eps_mag_zp = epsilon/mag_zp;
        else eps_mag_zp = static_cast<Real>(0.0);


        assert(!isnan(eps_mag_xm));
        assert(!isnan(eps_mag_xp));
        assert(!isnan(eps_mag_ym));
        assert(!isnan(eps_mag_yp));
        assert(!isnan(eps_mag_zm));
        assert(!isnan(eps_mag_zp));

        assert(!isinf(eps_mag_xm));
        assert(!isinf(eps_mag_xp));
        assert(!isinf(eps_mag_ym));
        assert(!isinf(eps_mag_yp));
        assert(!isinf(eps_mag_zm));
        assert(!isinf(eps_mag_zp));

        f_xm_a1r1 = (vxxa2(ix,iy)*vxxa1r1(ix,iy) + vxya2(ix,iy)*vxya1r1(ix,iy) + vxza2(ix,iy)*vxza1r1(ix,iy)) * eps_mag_xm;
        f_xp_a1r1 = (vxxa2(ix+1,iy)*vxxa1r1(ix+1,iy) + vxya2(ix+1,iy)*vxya1r1(ix+1,iy) + vxza2(ix+1,iy)*vxza1r1(ix+1,iy)) * eps_mag_xp;
        f_ym_a1r1 = (vyxa2(ix,iy)*vyxa1r1(ix,iy) + vyya2(ix,iy)*vyya1r1(ix,iy) + vyza2(ix,iy)*vyza1r1(ix,iy)) * eps_mag_ym;
        f_yp_a1r1 = (vyxa2(ix,iy+1)*vyxa1r1(ix,iy+1) + vyya2(ix,iy+1)*vyya1r1(ix,iy+1) + vyza2(ix,iy+1)*vyza1r1(ix,iy+1)) * eps_mag_yp;
        f_zm_a1r1 = (vzxa2_0(ix,iy)*vzxa1r1_0(ix,iy) + vzya2_0(ix,iy)*vzya1r1_0(ix,iy) + vzza2_0(ix,iy)*vzza1r1_0(ix,iy)) * eps_mag_zm;
        f_zp_a1r1 = (vzxa2_1(ix,iy)*vzxa1r1_1(ix,iy) + vzya2_1(ix,iy)*vzya1r1_1(ix,iy) + vzza2_1(ix,iy)*vzza1r1_1(ix,iy)) * eps_mag_zp;

        f_xm_a2r2 = (vxxa2(ix,iy)*vxxa2r2(ix,iy) + vxya2(ix,iy)*vxya2r2(ix,iy) + vxza2(ix,iy)*vxza2r2(ix,iy)) * eps_mag_xm;
        f_xp_a2r2 = (vxxa2(ix+1,iy)*vxxa2r2(ix+1,iy) + vxya2(ix+1,iy)*vxya2r2(ix+1,iy) + vxza2(ix+1,iy)*vxza2r2(ix+1,iy)) * eps_mag_xp;
        f_ym_a2r2 = (vyxa2(ix,iy)*vyxa2r2(ix,iy) + vyya2(ix,iy)*vyya2r2(ix,iy) + vyza2(ix,iy)*vyza2r2(ix,iy)) * eps_mag_ym;
        f_yp_a2r2 = (vyxa2(ix,iy+1)*vyxa2r2(ix,iy+1) + vyya2(ix,iy+1)*vyya2r2(ix,iy+1) + vyza2(ix,iy+1)*vyza2r2(ix,iy+1)) * eps_mag_yp;
        f_zm_a2r2 = (vzxa2_0(ix,iy)*vzxa2r2_0(ix,iy) + vzya2_0(ix,iy)*vzya2r2_0(ix,iy) + vzza2_0(ix,iy)*vzza2r2_0(ix,iy)) * eps_mag_zm;
        f_zp_a2r2 = (vzxa2_1(ix,iy)*vzxa2r2_1(ix,iy) + vzya2_1(ix,iy)*vzya2r2_1(ix,iy) + vzza2_1(ix,iy)*vzza2r2_1(ix,iy)) * eps_mag_zp;


        rhsa1r1.ref(ix,iy) = inv_mag_h * ( nx(ix,iy) * ((f_xp_a1r1 - f_xm_a1r1) - (static_cast<Real>(0.5) - a2(ix,iy)) * (a1r1(ix+1,iy) - a1r1(ix-1,iy)))
                                         + ny(ix,iy) * ((f_yp_a1r1 - f_ym_a1r1) - (static_cast<Real>(0.5) - a2(ix,iy)) * (a1r1(ix,iy+1) - a1r1(ix,iy-1)))
                                         + nz(ix,iy) * ((f_zp_a1r1 - f_zm_a1r1) - (static_cast<Real>(0.5) - a2(ix,iy)) * (a1r1p1(ix,iy) - a1r1m1(ix,iy))) );

        rhsa2r2.ref(ix,iy) = inv_mag_h * ( nx(ix,iy) * ((f_xp_a2r2 - f_xm_a2r2) - (static_cast<Real>(0.5) - a2(ix,iy)) * (a2r2(ix+1,iy) - a2r2(ix-1,iy)))
                                         + ny(ix,iy) * ((f_yp_a2r2 - f_ym_a2r2) - (static_cast<Real>(0.5) - a2(ix,iy)) * (a2r2(ix,iy+1) - a2r2(ix,iy-1)))
                                         + nz(ix,iy) * ((f_zp_a2r2 - f_zm_a2r2) - (static_cast<Real>(0.5) - a2(ix,iy)) * (a2r2p1(ix,iy) - a2r2m1(ix,iy))) );

        sum_R = rhsa1r1(ix,iy) + rhsa2r2(ix,iy);

        rhsru.ref(ix,iy) = u(ix,iy) * sum_R;
        rhsrv.ref(ix,iy) = v(ix,iy) * sum_R;
        rhsrw.ref(ix,iy) = w(ix,iy) * sum_R;

        rhsa2.ref(ix,iy) = inv_mag_h * (nx(ix,iy) * (epsilon*(mag_xp - mag_xm) - static_cast<Real>(0.25) * ( (a2(ix+1,iy) + a2(ix,iy)) * (static_cast<Real>(2.0) - a2(ix+1,iy) - a2(ix,iy))
                                                                                                           - (a2(ix,iy) + a2(ix-1,iy)) * (static_cast<Real>(2.0) - a2(ix,iy) - a2(ix-1,iy)) ))
                                       +ny(ix,iy) * (epsilon*(mag_yp - mag_ym) - static_cast<Real>(0.25) * ( (a2(ix,iy+1) + a2(ix,iy)) * (static_cast<Real>(2.0) - a2(ix,iy+1) - a2(ix,iy))
                                                                                                           - (a2(ix,iy) + a2(ix,iy-1)) * (static_cast<Real>(2.0) - a2(ix,iy) - a2(ix,iy-1)) ))
                                       +nz(ix,iy) * (epsilon*(mag_zp - mag_zm) - static_cast<Real>(0.25) * ( (a2p1(ix,iy) + a2(ix,iy)) * (static_cast<Real>(2.0) - a2p1(ix,iy) - a2(ix,iy))
                                                                                                          - (a2(ix,iy) + a2m1(ix,iy)) * (static_cast<Real>(2.0) - a2(ix,iy) - a2m1(ix,iy)) )) );
        kappa = static_cast<Real>(0.5) * (u(ix,iy)*u(ix,iy) + v(ix,iy)*v(ix,iy) + w(ix,iy)*w(ix,iy));

        rhse.ref(ix,iy) = kappa * sum_R + (p(ix,iy)*diff_G + diff_P) * rhsa2(ix,iy);
      }
      else
      {
        rhsa1r1.ref(ix,iy) = static_cast<Real>(0.0);
        rhsa2r2.ref(ix,iy) = static_cast<Real>(0.0);
        rhsru.ref(ix,iy) = static_cast<Real>(0.0);
        rhsrv.ref(ix,iy) = static_cast<Real>(0.0);
        rhsrw.ref(ix,iy) = static_cast<Real>(0.0);
        rhse.ref(ix,iy) = static_cast<Real>(0.0);
        rhsa2.ref(ix,iy) = static_cast<Real>(0.0);
      }
    }

  return;
}

void InterfaceSharpening_CPP_5eq::_copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
{
  // dtinvh already incorporated in rhs function

  for(int iy=0; iy<OutputSOA::NY; iy++)
    for(int ix=0; ix<OutputSOA::NX; ix++)
    {
      AssumedType& rhs = *(AssumedType*)(gptfirst + gptfloats*(ix + iy*rowgpts));

      assert(!isnan(rhsa1r1(ix,iy)));
      assert(!isnan(rhsa2r2(ix,iy)));
      assert(!isnan(rhsru(ix,iy)));
      assert(!isnan(rhsrv(ix,iy)));
      assert(!isnan(rhsrw(ix,iy)));
      assert(!isnan(rhse(ix,iy)));
      assert(!isnan(rhsa2(ix,iy)));

      rhs.a1r1 = a*rhs.a1r1 + rhsa1r1(ix, iy);
      rhs.a2r2 = a*rhs.a2r2 + rhsa2r2(ix, iy);
      rhs.ru   = a*rhs.ru   + rhsru(ix,iy);
      rhs.rv   = a*rhs.rv   + rhsrv(ix,iy);
      rhs.rw   = a*rhs.rw   + rhsrw(ix,iy);
      rhs.e    = a*rhs.e    + rhse(ix,iy);
      rhs.a2   = a*rhs.a2   + rhsa2(ix,iy);

      // cout << rhs.a1r1 << " " << rhs.a2r2  << " " << rhs.ru << " " << rhs.rv << " " << rhs.rw << " " << rhs.e << " " << rhs.a2 << endl;
    }
  return;
}
