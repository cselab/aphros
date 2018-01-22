/*
 *  Test_Utils_5eq.h
 *  MPCFcore
 *
 *  Created by Ursula Rasthofer on 2016/03/02.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */

#ifndef TEST_UTILS_5EQ_H
#define TEST_UTILS_5EQ_H

#include <ArgumentParser.h>
#include "TestTypes.h"
#include <fstream>

Real get_distance_to_interface(string interfacetype, const Real p[], const Real rb=0.0)
{
  Real d = 0.0;
  if (interfacetype == "bubble-segment")
    d = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) -  rb;
  else
  {
    // inclined plane cutting unit cube at (0, 0 , 0.4), (0.5, 0, 1) and (0, 0.6, 1)
    // normal = (0.36, 0.3, -0.3)
    d = abs(static_cast<Real>(0.36) * p[0] + static_cast<Real>(0.3) * p[1] - static_cast<Real>(0.3) * p[2] + static_cast<Real>(0.12)) / static_cast<Real>(sqrt(0.36*0.36 + 2.0 * 0.3 *0.3));
  }

  return d;
}

void get_velocity(string velocitytype, const Real p[], Real u[], const Real rb=0.0)
{
  if (velocitytype == "bubble")
  {
    const Real norm = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    const Real val = (static_cast<Real>(0.5) * norm /(rb*rb) -  static_cast<Real>(6.0) * norm / rb);
    u[0] = val * p[0];
    u[1] = val * p[1];
    u[2] = val * p[2];
  }
  else if (velocitytype == "inclined-interface")
  {
    // inclined plane cutting unit cube at (0, 0 , 0.4), (0.5, 0, 1) and (0, 0.6, 1)
    const Real n[3] = {0.36, 0.3, -0.3};
    const Real norm = 0.25*static_cast<Real>(sqrt(n[0]*n[0] + n[1] * n[1] + n[2] * n[2])); 
    u[0] = norm * n[0];
    u[1] = norm * n[1];
    u[2] = norm * n[2];
  }
  else
  {
    u[0] = static_cast<Real>(1.5) * p[1]*p[1];
    u[1] = static_cast<Real>(0.3) * exp(-static_cast<Real>(0.1)*p[2]);
    u[2] = -static_cast<Real>(3.8) * p[0]*p[2] + static_cast<Real>(1.5);
  }

  return;
}

void _initialize_block(Block& block)
{
  for(int iz = 0; iz<_BLOCKSIZE_; iz++)
    for(int iy = 0; iy<_BLOCKSIZE_; iy++)
      for(int ix = 0; ix<_BLOCKSIZE_; ix++)
      {
        block(ix, iy, iz).clear();

        block(ix, iy, iz).s.a1r1 = 1.0;
        block(ix, iy, iz).s.a2r2 = 0.0;
        block(ix, iy, iz).s.ru = 1.0;
        block(ix, iy, iz).s.rv = 0.1;
        block(ix, iy, iz).s.rw = -0.4;
        block(ix, iy, iz).s.e = 100.0;
        block(ix, iy, iz).s.a2 = 0.0;

        block(ix, iy, iz).dsdt.a1r1 = 1.0;
        block(ix, iy, iz).dsdt.a2r2 = 0.0;
        block(ix, iy, iz).dsdt.ru = 1.0;
        block(ix, iy, iz).dsdt.rv = 0.1;
        block(ix, iy, iz).dsdt.rw = -0.4;
        block(ix, iy, iz).dsdt.e = 100.0;
        block(ix, iy, iz).dsdt.a2 = 0.0;
      }
  return;
}

void _initialize_block(Block& block, ArgumentParser& parser, bool fill_s=true)
{
  for(int iz = 0; iz<_BLOCKSIZE_; iz++)
    for(int iy = 0; iy<_BLOCKSIZE_; iy++)
      for(int ix = 0; ix<_BLOCKSIZE_; ix++)
      {
        block(ix, iy, iz).clear();
        // assume a block of unit length
        const Real h = static_cast<Real>(1.0)/static_cast<Real>(_BLOCKSIZE_);
        // get cell center position
        Real p[3];
        p[0] = static_cast<Real>(0.5)*h*(static_cast<Real>(1.0) + static_cast<Real>(2.0) * static_cast<Real>(ix));
        p[1] = static_cast<Real>(0.5)*h*(static_cast<Real>(1.0) + static_cast<Real>(2.0) * static_cast<Real>(iy));
        p[2] = static_cast<Real>(0.5)*h*(static_cast<Real>(1.0) + static_cast<Real>(2.0) * static_cast<Real>(iz));
        // get physical parameters
        const Real g1 = parser("-g1").asDouble(6.59);
        const Real g2 = parser("-g2").asDouble(1.4);
        const Real g1m1Inv = static_cast<Real>(1.0)/(g1-static_cast<Real>(1.0));
        const Real g2m1Inv = static_cast<Real>(1.0)/(g2-static_cast<Real>(1.0));
        const Real pc1 = parser("-pc1").asDouble(4.069e3);
        const Real pc2 = parser("-pc2").asDouble(0.01);
        const Real rho1 = parser("-rho1").asDouble(858.3);
        const Real rho2 = parser("-rho2").asDouble(1.16);
        const Real p1 = parser("-p1").asDouble(35.8);
        const Real p2 = parser("-p2").asDouble(5.0);
        const Real eps = parser("-mollfactor").asDouble(1.0);

        // compute distance to interface
        const string interfacetype = parser("-interface-shape").asString("bubble-segment");
        const string velocitytype = parser("-velocity").asString("bubble");
        const Real rb = 0.5;
        const Real d = get_distance_to_interface(interfacetype, p, rb);

        // compute conserved variables
        const Real phi = M_PI*min(static_cast<Real>(1.0), max(static_cast<Real>(0.0), (d + eps)/(static_cast<Real>(2.0) * eps)));
        const Real alpha2 = static_cast<Real>(0.5) + static_cast<Real>(0.5) * cos(phi);
        const Real alpha1 = static_cast<Real>(1.0) -  alpha2;
        const Real alpha1rho1 = alpha1 * rho1;
        const Real alpha2rho2 = alpha2 * rho2;
        const Real dens = alpha1rho1 + alpha2rho2;
        const Real GmixInv = alpha1*g1m1Inv + alpha2*g2m1Inv;
        const Real Pmix = max(static_cast<Real>(0.0),alpha1*g1m1Inv*g1*pc1 + alpha2*g2m1Inv*g2*pc2);
        Real pressure = 0.0;
        if (interfacetype == "bubble-segment")
          pressure = d<=static_cast<Real>(0.0) ? p2 : (p1-rb*(p1-p2)/(rb+d));
        else
        {
          if (velocitytype == "inclined-interface")
            pressure = p2;
          else
            pressure = static_cast<Real>(65.8)  * p[0] * p[1] + static_cast<Real>(20.3) * p[2] * p[0] * p[0] - static_cast<Real>(8.9) * p[2] + static_cast<Real>(42.7);
         }
        const Real inte = pressure * GmixInv + Pmix;
        // compute velocity
        Real u[3];
        get_velocity(velocitytype, p, u, rb);
        const Real kine = 0.5 * dens * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]); 

        // initialize with some default values
        block(ix, iy, iz).s.a1r1 = 1.0;
        block(ix, iy, iz).s.a2r2 = 0.0;
        block(ix, iy, iz).s.ru = 1.0;
        block(ix, iy, iz).s.rv = 0.1;
        block(ix, iy, iz).s.rw = -0.4;
        block(ix, iy, iz).s.e = 100.0;
        block(ix, iy, iz).s.a2 = 0.0;

        block(ix, iy, iz).dsdt.a1r1 = 1.0;
        block(ix, iy, iz).dsdt.a2r2 = 0.0;
        block(ix, iy, iz).dsdt.ru = 1.0;
        block(ix, iy, iz).dsdt.rv = 0.1;
        block(ix, iy, iz).dsdt.rw = -0.4;
        block(ix, iy, iz).dsdt.e = 100.0;
        block(ix, iy, iz).dsdt.a2 = 0.0;

       // initialize block accrding to computed values
       if (fill_s)
       {
         block(ix, iy, iz).s.a1r1 = alpha1rho1;
         block(ix, iy, iz).s.a2r2 = alpha2rho2;
         block(ix, iy, iz).s.ru = dens * u[0];
         block(ix, iy, iz).s.rv = dens * u[1];
         block(ix, iy, iz).s.rw = dens * u[2];
         block(ix, iy, iz).s.e = kine + inte;
         block(ix, iy, iz).s.a2 = alpha2;
       }
       else
       {
         block(ix, iy, iz).dsdt.a1r1 = alpha1rho1;
         block(ix, iy, iz).dsdt.a2r2 = alpha2rho2;
         block(ix, iy, iz).dsdt.ru = dens * u[0];
         block(ix, iy, iz).dsdt.rv = dens * u[1];
         block(ix, iy, iz).dsdt.rw = dens * u[2];
         block(ix, iy, iz).dsdt.e = kine + inte;
         block(ix, iy, iz).dsdt.a2 = alpha2;
       }
    }

  return;
}

void _initialize_block(Block& block, string fname, const bool fill_s)
{
  std::ifstream infile(fname.c_str());
  Real ia1r1, ia2r2, iru, irv, irw, iE, ia2;

  int c=0;
  while (infile >> ia1r1 >> ia2r2 >> iru >> irv >> irw >> iE >> ia2)
  {
    const int ix = (c % _BLOCKSIZE_);
    const int iy = (c/_BLOCKSIZE_ % _BLOCKSIZE_);
    const int iz = (c/_BLOCKSIZE_ / _BLOCKSIZE_);

    block(ix, iy, iz).clear();

   // initialize with some default values
    block(ix, iy, iz).s.a1r1 = 1.0;
    block(ix, iy, iz).s.a2r2 = 0.0;
    block(ix, iy, iz).s.ru = 1.0;
    block(ix, iy, iz).s.rv = 0.1;
    block(ix, iy, iz).s.rw = -0.4;
    block(ix, iy, iz).s.e = 100.0;
    block(ix, iy, iz).s.a2 = 0.0;

    block(ix, iy, iz).dsdt.a1r1 = 1.0;
    block(ix, iy, iz).dsdt.a2r2 = 0.0;
    block(ix, iy, iz).dsdt.ru = 1.0;
    block(ix, iy, iz).dsdt.rv = 0.1;
    block(ix, iy, iz).dsdt.rw = -0.4;
    block(ix, iy, iz).dsdt.e = 100.0;
    block(ix, iy, iz).dsdt.a2 = 0.0;

    // initialize block accrding to computed values
    if (fill_s)
    {
      block(ix, iy, iz).s.a1r1 = ia1r1;
      block(ix, iy, iz).s.a2r2 = ia2r2;
      block(ix, iy, iz).s.ru = iru;
      block(ix, iy, iz).s.rv = irv;
      block(ix, iy, iz).s.rw = irw;
      block(ix, iy, iz).s.e = iE;
      block(ix, iy, iz).s.a2 = ia2;
    }
    else
    {
      block(ix, iy, iz).dsdt.a1r1 = ia1r1;
      block(ix, iy, iz).dsdt.a2r2 = ia2r2;
      block(ix, iy, iz).dsdt.ru = iru;
      block(ix, iy, iz).dsdt.rv = irv;
      block(ix, iy, iz).dsdt.rw = irw;
      block(ix, iy, iz).dsdt.e = iE;
      block(ix, iy, iz).dsdt.a2 = ia2;
    }

    ++c;
  }

  infile.close();

  return;
}

void write_block_to_file(Block& block, string fname, const bool write_s)
{
  ofstream outfile (fname.c_str());

  for(int iz = 0; iz<_BLOCKSIZE_; iz++)
    for(int iy = 0; iy<_BLOCKSIZE_; iy++)
      for(int ix = 0; ix<_BLOCKSIZE_; ix++)
      {
        if (write_s)
          outfile << scientific << setprecision(16) << block(ix, iy, iz).s.a1r1 << " " << block(ix, iy, iz).s.a2r2 << " " << block(ix, iy, iz).s.ru << " " << block(ix, iy, iz).s.rv << " " << block(ix, iy, iz).s.rw << " " << block(ix, iy, iz).s.e << " " << block(ix, iy, iz).s.a2 << endl;
        else
          outfile << scientific << setprecision(16) << block(ix, iy, iz).dsdt.a1r1 << " " << block(ix, iy, iz).dsdt.a2r2 << " " << block(ix, iy, iz).dsdt.ru << " " << block(ix, iy, iz).dsdt.rv << " " << block(ix, iy, iz).dsdt.rw << " " << block(ix, iy, iz).dsdt.e << " " << block(ix, iy, iz).dsdt.a2 << endl;
      }

  outfile.close(); 

  return;
}


template<int WIDTH>
void _initialize_lab(TestLab<WIDTH>& lab)
{
  for(int iz = -WIDTH; iz<_BLOCKSIZE_+WIDTH; iz++)
    for(int iy = -WIDTH; iy<_BLOCKSIZE_+WIDTH; iy++)
      for(int ix = -WIDTH; ix<_BLOCKSIZE_+WIDTH; ix++)
      {
        lab(ix, iy, iz).clear();

        lab(ix, iy, iz).s.a1r1 = 1.0;
        lab(ix, iy, iz).s.a2r2 = 0.0;
        lab(ix, iy, iz).s.ru = 1.0;
        lab(ix, iy, iz).s.rv = 0.1;
        lab(ix, iy, iz).s.rw = -0.4;
        lab(ix, iy, iz).s.e = 100.0;
        lab(ix, iy, iz).s.a2 = 0.0;
      }
  return;
}

template<int WIDTH>
void _initialize_lab(TestLab<WIDTH>& lab, ArgumentParser& parser)
{
  for(int iz = -WIDTH; iz<_BLOCKSIZE_+WIDTH; iz++)
    for(int iy = -WIDTH; iy<_BLOCKSIZE_+WIDTH; iy++)
      for(int ix = -WIDTH; ix<_BLOCKSIZE_+WIDTH; ix++)
      {
        lab(ix, iy, iz).clear();
        // assume a block of unit length
        const Real h = static_cast<Real>(1.0)/static_cast<Real>(_BLOCKSIZE_);
        // get cell center position
        Real p[3];
        p[0] = static_cast<Real>(0.5)*h*(static_cast<Real>(1.0) + static_cast<Real>(2.0) * static_cast<Real>(ix));
        p[1] = static_cast<Real>(0.5)*h*(static_cast<Real>(1.0) + static_cast<Real>(2.0) * static_cast<Real>(iy));
        p[2] = static_cast<Real>(0.5)*h*(static_cast<Real>(1.0) + static_cast<Real>(2.0) * static_cast<Real>(iz));
        // get physical parameters
        const Real g1 = parser("-g1").asDouble(6.59);
        const Real g2 = parser("-g2").asDouble(1.4);
        const Real g1m1Inv = static_cast<Real>(1.0)/(g1-static_cast<Real>(1.0));
        const Real g2m1Inv = static_cast<Real>(1.0)/(g2-static_cast<Real>(1.0));
        const Real pc1 = parser("-pc1").asDouble(4.069e3);
        const Real pc2 = parser("-pc2").asDouble(0.01);
        const Real rho1 = parser("-rho1").asDouble(858.3);
        const Real rho2 = parser("-rho2").asDouble(1.16);
        const Real p1 = parser("-p1").asDouble(35.8);
        const Real p2 = parser("-p2").asDouble(5.0);
        const Real eps = parser("-mollfactor").asDouble(1.0);

        // compute distance to interface
        const string interfacetype = parser("-interface-shape").asString("bubble-segment");
        const string velocitytype = parser("-velocity").asString("bubble");
        const Real rb = 0.5;
        const Real d = get_distance_to_interface(interfacetype, p, rb);

        // compute conserved variables
        const Real phi = M_PI*min(static_cast<Real>(1.0), max(static_cast<Real>(0.0), (d + eps)/(static_cast<Real>(2.0) * eps)));
        const Real alpha2 = static_cast<Real>(0.5) + static_cast<Real>(0.5) * cos(phi);
        const Real alpha1 = static_cast<Real>(1.0) -  alpha2;
        const Real alpha1rho1 = alpha1 * rho1;
        const Real alpha2rho2 = alpha2 * rho2;
        const Real dens = alpha1rho1 + alpha2rho2;
        const Real GmixInv = alpha1*g1m1Inv + alpha2*g2m1Inv;
        const Real Pmix = max(static_cast<Real>(0.0),alpha1*g1m1Inv*g1*pc1 + alpha2*g2m1Inv*g2*pc2);
        Real pressure = 0.0;
        if (interfacetype == "bubble-segment")
          pressure = d<=static_cast<Real>(0.0) ? p2 : (p1-rb*(p1-p2)/(rb+d));
        else
        {
          if (velocitytype == "inclined-interface")
            pressure = p2;
          else
            pressure = static_cast<Real>(65.8)  * p[0] * p[1] + static_cast<Real>(20.3) * p[2] * p[0] * p[0] - static_cast<Real>(8.9) * p[2] + static_cast<Real>(42.7);
         }
        const Real inte = pressure * GmixInv + Pmix;
        // compute velocity
        Real u[3];
        get_velocity(velocitytype, p, u, rb);
        const Real kine = 0.5 * dens * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

       // initialize block accrding to computed values
       lab(ix, iy, iz).s.a1r1 = alpha1rho1;
       lab(ix, iy, iz).s.a2r2 = alpha2rho2;
       lab(ix, iy, iz).s.ru = dens * u[0];
       lab(ix, iy, iz).s.rv = dens * u[1];
       lab(ix, iy, iz).s.rw = dens * u[2];
       lab(ix, iy, iz).s.e = kine + inte;
       lab(ix, iy, iz).s.a2 = alpha2;
    }

  return;
}


#endif /* TEST_UTILS_%EQ_H */
