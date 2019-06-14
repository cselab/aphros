#pragma once

#include <functional>
#include <stdexcept>

#include "solver/fluid.h"
#include "parse/vars.h"


template <class M>
FieldCell<typename M::Scal> GetBcField(
    MapFace<std::shared_ptr<solver::CondFaceFluid>>& mf, const M& m) {
  FieldCell<typename M::Scal> fc(m, 0);
  for (auto it : mf) {
    IdxFace f = it.GetIdx();
    auto* b = it.GetValue().get();
    size_t nci = b->GetNci();
    IdxCell c = m.GetNeighbourCell(f, nci);
    if (dynamic_cast<solver::fluid_condition::NoSlipWall<M>*>(b)) {
      fc[c] = 1.;
    } else if (dynamic_cast<solver::fluid_condition::SlipWall<M>*>(b)) {
      fc[c] = 2.;
    } else if (dynamic_cast<solver::fluid_condition::Inlet<M>*>(b)) {
      fc[c] = 3.;
    } else if (dynamic_cast<solver::fluid_condition::Outlet<M>*>(b)) {
      fc[c] = 4.;
    } else {
      fc[c] = -1.;
    }
  }
  return fc;
}

template <class M>
FieldCell<typename M::Vect> GetVort(const FieldCell<typename M::Vect>& fcv,
                       const MapFace<std::shared_ptr<solver::CondFace>>& mf,
                       M& m) {
  auto ffv = solver::Interpolate(fcv, mf, m);

  auto d0 = solver::Gradient(GetComponent(ffv, 0), m);
  auto d1 = solver::Gradient(GetComponent(ffv, 1), m);
  auto d2 = solver::Gradient(GetComponent(ffv, 2), m);

  FieldCell<typename M::Vect> r(m);
  for (auto c : m.Cells()) {
    r[c][0] = d2[c][1] - d1[c][2];
    r[c][1] = d0[c][2] - d2[c][0];
    r[c][2] = d1[c][0] - d0[c][1];
  }

  return r;
}

template <class M>
void InitVel(FieldCell<typename M::Vect>& fcv, const Vars& var, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const std::string vi = var.String["vel_init"];
  if (vi == "taylor-green") {
    for (auto i : m.AllCells()) {
      auto& v = fcv[i];
      auto x = m.GetCenter(i);
      if (var.Int["dim"] == 2) {
        x[2] = 0.;
      }
      v[0] = std::sin(x[0]) * std::cos(x[1]) * std::cos(x[2]);
      v[1] = -std::cos(x[0]) * std::sin(x[1]) * std::cos(x[2]);
      v[2] = 0.;
    }
  } else if (vi == "vortex") {
    Scal pi = M_PI;
    Vect xc1(var.Vect["vort_x1"]);
    Vect xc2(var.Vect["vort_x2"]);
    Scal g1 = var.Double["vort_g1"];
    Scal g2 = var.Double["vort_g2"];
    Scal vmax = var.Double["vort_vmax"];
    std::vector<Vect> xxc = {xc1, xc2};
    std::vector<Scal> gg = {g1, g2};
    for (size_t i = 0; i < xxc.size(); ++i) {
      Vect xc = xxc[i];
      Scal g = gg[i];
      for (auto c : m.AllCells()) {
        Vect x = m.GetCenter(c);
        Scal r = xc.dist(x);
        Scal r2 = r * r;
        Vect v(0);
        v[0] = -(x[1] - xc[1]) * g / r2 * 2. * pi;
        v[1] = (x[0] - xc[0]) * g / r2 * 2. * pi;
        if (v.norm() > vmax) {
          v *= vmax / v.norm();
        }
        fcv[c][0] += v[0];
        fcv[c][1] += v[1];
      }
    }
  } else if (vi == "vortexring") { // XXX: use with initvort=1
    Vect xc(var.Vect["ring_c"]); // center of ring
    Vect n(var.Vect["ring_n"]);  // normal
    n /= n.norm();
    Scal om = var.Double["ring_om"];  // vorticity in crosssection
    Scal r0 = var.Double["ring_r0"];  // inner radius
    Scal r1 = var.Double["ring_r1"];  // outer radius
    Scal qr = (r1 - r0) * 0.5;   // radius
    Vect qc((r1 + r0) * 0.5, 0., 0.);   // center
    const Scal eps = 1e-10;
    for (auto c : m.AllCells()) {
      Vect x = m.GetCenter(c) - xc;
      // along axis
      Scal xn = n.dot(x);
      // along plane
      Scal xt = (x - n * xn).norm();
      // unit along plane
      Vect t = (x - n * xn) / std::max(eps, xt);
      // unit along circle 
      Vect s = n.cross(t);
      Vect q(xt, xn, 0.);
      fcv[c] = ((q - qc).sqrnorm() <= sqr(qr) ? s * om : Vect(0));
    }
  } else if (vi == "vortexgauss") {
    // XXX: Bergdorf 2007, Direct numerical simulations of vortex rings
    Scal pi = M_PI;
    Scal g(var.Double["ring_gamma"]); // circulation
    Scal sig(var.Double["ring_sigma"]); // sigma
    Scal rad(var.Double["ring_r"]); // radius
    Scal nfr(var.Double["ring_noise_freq"]); // noise angular frequency
    Scal namp(var.Double["ring_noise_amp"]); // noise amp relative to r
    Scal nfr2(var.Double["ring_noise2_freq"]); // noise angular frequency
    Scal namp2(var.Double["ring_noise2_amp"]); // noise amp relative to r
    Vect xc(var.Vect["ring_c"]); // center
    Vect n(var.Vect["ring_n"]);  // normal
    n /= n.norm();
    const Scal eps = 1e-10;
    for (auto c : m.AllCells()) {
      Vect x = m.GetCenter(c) - xc;
      // along axis
      Scal xn = n.dot(x);
      // select direction 
      Vect v(0);
      v[n.abs().argmin()] = 1.;
      // unit vectors in plane
      Vect vx = n.cross(v);
      Vect vy = n.cross(vx);
      // angle
      Scal a = std::atan2(vx.dot(x), vy.dot(x));
      // noise
      xn += rad * namp * std::sin(a * nfr);
      // along plane
      Scal xt = (x - n * xn).norm();
      // noise2
      xt += rad * namp2 * std::sin(a * nfr2);
      // unit radial along plane
      Vect et = (x - n * xn) / std::max(eps, xt);
      // unit along circle 
      Vect es = n.cross(et);
      Scal s2 = sqr(xn) + sqr(xt - rad);
      Scal sig2 = sqr(sig);
      Scal om = g / (pi * sig2) * std::exp(-s2 / sig2);
      fcv[c] = es * om;
    }
  } else if (vi == "rot") {
    Vect xc(var.Vect["rot_c"]);
    Scal om = var.Double["rot_om"];
    for (auto c : m.AllCells()) {
      Vect x = m.GetCenter(c);
      Vect v(0);
      v[0] = -(x[1] - xc[1]) * om * 0.5;
      v[1] = (x[0] - xc[0]) * om * 0.5;
      fcv[c] = v;
    }
  } else if (vi == "grad") {
    Vect xc(var.Vect["grad_c"]);
    Vect vx(var.Vect["grad_vx"]); // gradient of vel[0]
    Vect vy(var.Vect["grad_vy"]); // gradient of vel[1]
    Vect vz(var.Vect["grad_vz"]); // gradient of vel[2]
    for (auto c : m.AllCells()) {
      Vect x = m.GetCenter(c) - xc;
      Vect v;
      v[0] = vx.dot(x);
      v[1] = vy.dot(x);
      v[2] = vz.dot(x);
      fcv[c] = v;
    }
  } else if (vi == "pois" || vi == "poisy") {
    // Poiseuille with walls in y
    auto gs = m.GetGlobalSize(); // global mesh size
    Scal ext = var.Double["extent"]; // TODO: revise
    Vect gh = Vect(gs) * ext / gs.max(); // global domain length
    Scal pv = var.Double["poisvel"]; // centerline velocity

    for (auto i : m.AllCells()) {
      Scal y = m.GetCenter(i)[1] / gh[1];
      fcv[i][0] = y * (1. - y) * 4. * pv;
    }
  } else if (vi == "poisyz") {
    // Poiseuille with walls in y and z
    // Spiga 1994: Symmetric solution for velocity in rectangular ducts
    auto gs = m.GetGlobalSize(); // global mesh size
    Scal ext = var.Double["extent"]; // TODO: revise
    Vect gh = Vect(gs) * ext / gs.max(); // global domain length
    Scal mu = var.Double["poismu"]; // viscosity
    Scal pg = var.Double["poisgrad"]; // pressure gradient
    int im = var.Int["poisiter"]; // depth to evaluate series
    bool wym = var.Int["poiswym"]; // wallym
    bool wyp = var.Int["poiswyp"]; // wallyp
    bool wzm = var.Int["poiswzm"]; // wallzm
    bool wzp = var.Int["poiswzp"]; // wallzp
    Scal pi = M_PI;

    Scal ly = gh[1];
    Scal lz = gh[2];

    if ((!wym && !wyp) || (!wzm && !wzp)) {
      throw std::runtime_error("poisyz: can't remove both walls");
    }

    Scal oy = 0.;
    if (!wym) {
      oy = 0.5;
      ly *= 2;
    } else if (!wyp) {
      oy = 0;
      ly *= 2;
    }

    Scal oz = 0.;
    if (!wzm) {
      oz = 0.5;
      lz *= 2;
    } else if (!wzp) {
      oz = 0;
      lz *= 2;
    }


    Scal p = sqr(ly) * pg / mu;
    Scal b = lz / ly;
    Scal k = 16. * sqr(b) / std::pow(pi, 4);

    // TODO: tests for gh[1] != gh[2]
    for (auto i : m.AllCells()) {
      Scal y = m.GetCenter(i)[1] / ly;
      y = y + oy;
      Scal z = m.GetCenter(i)[2] / lz;
      z = z + oz;
      Scal s = 0.;
      for (int iy = 1; iy < im * 2; iy += 2) {
        for (int iz = 1; iz < im * 2; iz += 2) {
          s += std::sin(iy * pi * y) *  std::sin(iz * pi * z) / 
              (iy * iz * (sqr(b) * sqr(iy) + sqr(iz)));
        }
      }
      fcv[i][0] = p * s * k;
    }
  } else if (vi == "uniform") {
    Vect v(var.Vect["vel"]);
    for (auto i : m.AllCells()) {
      fcv[i] = v;
    }
  } else if (vi == "box") {
    Vect v(var.Vect["vel"]);
    Vect a(var.Vect["vel_box_a"]);
    Vect b(var.Vect["vel_box_b"]);
    Rect<Vect> r(a, b);
    for (auto c : m.AllCells()) {
      if (r.IsInside(m.GetCenter(c))) {
        fcv[c] = v;
      }
    }
  } else if (vi == "solitoncos") {
    Scal xc(var.Double["soliton_xc"]);
    Scal yc(var.Double["soliton_yc"]);
    Scal yh(var.Double["soliton_yh"]);
    Scal xw(var.Double["soliton_xw"]);
    auto sinh = [](Scal t) { return std::sinh(t); };
    auto cosh = [](Scal t) { return std::cosh(t); };
    auto tanh = [](Scal t) { return std::tanh(t); };
    auto cos = [](Scal t) { return std::cos(t); };
    auto sin = [](Scal t) { return std::sin(t); };
    Scal g(Vect(var.Vect["gravity"]).norm());
    Scal p2 = 2 * M_PI;
    Scal L = xw;
    Scal D = yc;
    Scal C = std::sqrt(g * L * tanh(p2 * D / L) / p2);
    Scal T = L / C;
    Scal H = yh;
    for (auto c : m.AllCells()) {
      Vect xx = m.GetCenter(c);
      Scal x = xx[0] - xc;
      Scal y = xx[1];
      Scal Ux = H * p2 / T * cosh((p2 / L) * y) * 
          cos(p2 * x / L) / sinh(p2 * D / L);
      Scal Uy = H * p2 / T * sinh((p2 / L) * y) * 
          sin(p2 * x / L) / sinh(p2 * D / L);
      fcv[c] = Vect(Ux, Uy, 0.);
    }
  } else if (vi == "zero") {
    // nop
  } else  {
    throw std::runtime_error("Init(): unknown vel_init=" + vi);
  }
}

template <class M>
void GetFluidFaceCond(
    const Vars& var, const M& m,
    MapFace<std::shared_ptr<solver::CondFaceFluid>>& mfvel,
    MapFace<std::shared_ptr<solver::CondFace>>& mfvf) {
  using Dir = typename M::Dir;
  using MIdx = typename M::MIdx;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  size_t edim = var.Int["dim"];
  auto& fi = m.GetIndexFaces();
  MIdx gs = m.GetGlobalSize();

  // boundary xm of global mesh
  auto gxm = [&fi](IdxFace i) -> bool {
    return fi.GetDir(i) == Dir::i && fi.GetMIdx(i)[0] == 0;
  };
  auto gxp = [&fi,gs](IdxFace i) -> bool {
    return fi.GetDir(i) == Dir::i && fi.GetMIdx(i)[0] == gs[0];
  };
  auto gym = [&fi](IdxFace i) -> bool {
    return fi.GetDir(i) == Dir::j && fi.GetMIdx(i)[1] == 0;
  };
  auto gyp = [&fi,gs](IdxFace i) -> bool {
    return fi.GetDir(i) == Dir::j && fi.GetMIdx(i)[1] == gs[1];
  };
  auto gzm = [&fi,edim](IdxFace i) -> bool {
    return edim >= 3 && fi.GetDir(i) == Dir::k && fi.GetMIdx(i)[2] == 0;
  };
  auto gzp = [&fi,gs,edim](IdxFace i) -> bool {
    return edim >= 3 && fi.GetDir(i) == Dir::k && fi.GetMIdx(i)[2] == gs[2];
  };

  // Set condition bc for face i on global box boundary
  // choosing proper neighbour cell id (nci)
  // Return true if on global boundary
  auto set_bc = [&](IdxFace i, std::string bc) -> bool {
    if (gxm(i) || gym(i) || gzm(i)) {
      mfvel[i] = solver::Parse(bc, i, 1, m);
      return true;
    } else if (gxp(i) || gyp(i) || gzp(i)) {
      mfvel[i] = solver::Parse(bc, i, 0, m);
      return true;
    }
    return false;
  };

  // Boundary conditions for fluid 
  auto ff = m.AllFaces();
  std::vector<std::pair<std::string, std::function<bool(IdxFace)>>> pp = 
      {{"bc_xm", gxm}, {"bc_xp", gxp},
       {"bc_ym", gym}, {"bc_yp", gyp},
       {"bc_zm", gzm}, {"bc_zp", gzp}};

  for (auto p : pp) {
    if (auto bc = var.String(p.first)) {
      for (auto i : ff) {
        p.second(i) && set_bc(i, *bc);
      }
    }
  }

  // boundary conditions for advection
  for (auto it : mfvel) {
    IdxFace i = it.GetIdx();
    solver::CondFaceFluid* cb = it.GetValue().get();
    if (dynamic_cast<solver::fluid_condition::SlipWall<M>*>(cb)) {
      mfvf[i] = std::make_shared<solver::
          CondFaceReflect>(it.GetValue()->GetNci());
    } else {
      mfvf[i] = std::make_shared<solver::
          CondFaceGradFixed<Scal>>(Scal(0), it.GetValue()->GetNci());
    }
  }
  // selection boxes
  // Parameters (N>=0):
  // string boxN -- bc description
  // vect boxN_a -- lower corner
  // vect boxN_b -- upper corner
  // double boxN_vf -- inlet volume fraction
  // Check at least first nmax indices and all contiguous
  {
    int n = 0;
    const int nmax = 100;
    while (true) {
      std::string k = "box" + std::to_string(n);
      if (auto p = var.String(k)) {
        Vect a(var.Vect[k + "_a"]);
        Vect b(var.Vect[k + "_b"]);
        Scal vf = var.Double[k + "_vf"];
        Rect<Vect> r(a, b);
        for (auto i : m.AllFaces()) {
          Vect x = m.GetCenter(i);
          if (r.IsInside(x)) {
            if (set_bc(i, *p)) {
              auto b = mfvel[i];
              mfvf[i] = std::make_shared
                  <solver::CondFaceValFixed<Scal>>(vf, b->GetNci());
            }
          }
        }
      } else if (n > nmax) { 
        break;
      }
      ++n;
    }
  }
  // selection faces
  // Parameters (N>=0):
  // string faceN -- bc description
  // vect faceN_a -- lower corner
  // vect faceN_b -- upper corner
  // double faceN_vf -- inlet volume fraction
  // Check at least first nmax indices and all contiguous
  {
    int n = -1;
    const int nmax = 100;
    while (true) {
      ++n;
      std::string k = "face" + std::to_string(n);
      if (auto p = var.String(k)) {
        // set boundary conditions on faces of box (a,b) 
        // normal to d with outer normal towards (b-a)[d]
        Vect a(var.Vect[k + "_a"]); 
        Vect b(var.Vect[k + "_b"]);
        int d(var.Int[k + "_dir"]); // direction: 0:x, 1:y, 2:z
        Scal vf = var.Double[k + "_vf"];
        Rect<Vect> r(a, b);
        Vect h = m.GetCellSize();
        auto& cb = m.GetInBlockCells();
        auto& fi = m.GetIndexFaces();
        // indices of [a,b), [begin,end)
        Vect xd(0);
        xd[d] = 1.;
        // round to faces
        a = Vect(MIdx((a + h * 0.5) / h)) * h;
        b = Vect(MIdx((b + h * 0.5) / h)) * h;
        // indices
        MIdx wa(a / h + xd * 0.5);
        MIdx wb(b / h + xd * 0.5);
        wb[d] = wa[d] + 1; // size 1 in direction d
        // direction
        MIdx wd(0);
        wd[d] = 1;
        // direction of neighbour cell
        int nci = ((b - a)[d] > 0. ? 0 : 1);
        // box of valid indices
        MIdx w0 = cb.GetBegin();
        MIdx w1 = cb.GetEnd() + wd;
        // clip (a,b) to valid indices
        wa = wa.clip(w0, w1);
        wb = wb.clip(w0, w1);
        // size of local block
        MIdx ws = wb - wa;
        if (ws.prod() == 0) {
          continue;
        }
        typename M::BlockCells bb(wa, ws);
        for (auto w : bb) {
          IdxFace f = fi.GetIdx(w, Dir(d));
          mfvel[f] = solver::Parse(*p, f, nci, m);
        }
        (void) vf;
      } else if (n > nmax) { 
        break;
      }
    }
  }
  // selection spheres
  // Parameters (N>=0):
  // string sphN -- bc description ("inlet")
  // vect sphN_c -- center
  // double sphN_r0 -- radius inner
  // double sphN_r1 -- radius outer
  // double sphN_vf -- inlet volume fraction
  // vect sphN_vel -- velocity
  // Check at least first nmax indices and all contiguous
  {
    int n = 0;
    const int nmax = 100;
    while (true) {
      std::string k = "sph" + std::to_string(n);
      if (auto p = var.String(k)) {
        if (*p == "inlet") {
          Vect xc(var.Vect[k + "_c"]);
          Scal r0 = var.Double[k + "_r0"];
          Scal r1 = var.Double[k + "_r1"];
          Scal vf = var.Double[k + "_vf"];
          Vect vel(var.Vect[k + "_vel"]);
          for (auto i : m.AllFaces()) {
            Vect x = m.GetCenter(i);
            if (xc.dist(x) < r1) {
              Scal r = xc.dist(x);
              Vect v = vel * std::min(1., (r1 - r) / (r1 - r0));
              std::string a = "inlet";
              a += " " + std::to_string(v[0]);
              a += " " + std::to_string(v[1]);
              a += " " + std::to_string(v[2]);
              if (set_bc(i, a)) {
                auto b = mfvel[i];
                mfvf[i] = std::make_shared
                    <solver::CondFaceValFixed<Scal>>(vf, b->GetNci());
              }
            }
          }
        } else {
          throw std::runtime_error("unknown selection sphere cond: " + *p);
        }
      } else if (n > nmax) { 
        break;
      }
      ++n;
    }
  }
}

template <class M>
void GetFluidCellCond(
    const Vars& var, M& m,
    MapCell<std::shared_ptr<solver::CondCellFluid>>& mcvel,
    std::pair<typename M::Scal, int>& pdist) {
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;

  auto sem = m.GetSem(__func__);
  if (sem("reduce")) {
    if (auto* p = var.Double("pfixed")) {
      Vect x(var.Vect["pfixed_x"]);
      // Find cell nearest to pfixed_x
      IdxCell c = m.FindNearestCell(x);
      pdist.first = m.GetCenter(c).dist(x);
      pdist.second = m.GetId();
      m.Reduce(std::make_shared<typename M::OpMinloc>(&pdist));
    }
  }

  if (sem("set")) {
    // Fix pressure at one cell
    if (auto* p = var.Double("pfixed")) {
      if (pdist.second == m.GetId()) {
        Vect x(var.Vect["pfixed_x"]);
        IdxCell c = m.FindNearestCell(x);
        mcvel[c] = std::make_shared
            <solver::fluid_condition::GivenPressureFixed<M>>(*p);
        std::cout 
            << "pfixed id=" << pdist.second 
            << " dist=" << pdist.first << std::endl;
      }
    }

    // exclude cells
    int n = -1;
    const int nmax = 100;
    while (true) {
      ++n;
      std::string k = "cellbox" + std::to_string(n);
      if (auto p = var.String(k)) {
        // set boundary conditions on faces of box (a,b) 
        // normal to d with outer normal towards (b-a)[d]
        Vect a(var.Vect[k + "_a"]); 
        Vect b(var.Vect[k + "_b"]);
        Rect<Vect> r(a, b);
        Vect h = m.GetCellSize();
        auto& cb = m.GetInBlockCells();
        auto& ci = m.GetIndexCells();
        // round to faces
        a = Vect(MIdx((a + h * 0.5) / h)) * h;
        b = Vect(MIdx((b + h * 0.5) / h)) * h;
        // indices
        MIdx wa(a / h);
        MIdx wb(b / h);
        // box of valid indices
        MIdx w0 = cb.GetBegin();
        MIdx w1 = cb.GetEnd();
        // clip (a,b) to valid indices
        wa = wa.clip(w0, w1);
        wb = wb.clip(w0, w1);
        // size of local block
        MIdx ws = wb - wa;
        if (ws.prod() == 0) {
          continue;
        }
        typename M::BlockCells bb(wa, ws);
        for (auto w : bb) {
          IdxCell c = ci.GetIdx(w);
          mcvel[c] = std::make_shared
              <solver::fluid_condition::
              GivenVelocityAndPressureFixed<M>>(Vect(0), 0.);
        }
      } else if (n > nmax) { 
        break;
      }
    }
  }
}
