// Created by Petr Karnakov on 15.07.2020
// Copyright 2020 ETH Zurich

#include "young/young.h"

class YoungMarangoni {
  YoungParam GetYoungPar() const {
    YoungParam q;
    q.rhov = var.Double["rho1"];
    q.rhou = var.Double["rho2"];
    q.muv = var.Double["mu1"];
    q.muu = var.Double["mu2"];
    q.hv = 1.;
    q.hu = 1.;
    q.gamma0 = var.Double["sigma"];
    q.gamma1 = var.Vect["sig_grad"][0];
    q.T0 = 0.;
    q.T1 = 1.;
    q.R = var.Double["youngbc_r"];
    return q;
  }
  void InitYoung() {
    young_ini(GetYoungPar());
  }
  Vect GetYoungVel(Vect x) const {
    x -= Vect(0.5);
    // 0: streamwise
    // 1,2: crossstream with symmetry around 0-axis
    Scal r = Vect(0, x[1], x[2]).norm();
    Vect v(0);
    Scal vc; // vel cross
    young_fields(r, x[0], &vc, &v[0], nullptr, nullptr);
    Scal e1 = x[1] / r;
    Scal e2 = x[2] / r;
    v[1] = vc * e1;
    v[2] = vc * e2;
    return v;
  }
  void Hydro<M>::InitFluid(const FieldCell<Vect>& fc_vel) {
    if (var.Int["youngbc"]) {
      fcyv_.Reinit(m);
      InitYoung();
      for (auto c : m.Cells()) {
        Vect x = m.GetCenter(c);
        fcyv_[c] = GetYoungVel(x);
      }
    }
  }
  void Hydro<M>::CalcStat() {
    if (sem("young")) {
      if (var.Int["youngbc"]) {
        InitYoung();
        for (auto& p : mebc_fluid_.GetMapFace()) {
          const IdxFace f = p.first;
          auto& bc = p.second;
          if (bc.type == BCondFluidType::wall) {
            const Vect x = m.GetCenter(f);
            bc.velocity = GetYoungVel(x);
          }
        }
      }
    }
  }
  void Dump() {
    if (var.Int["youngbc"]) {
      dumpv(fcyv_, 0, "yvx");
      dumpv(fcyv_, 1, "yvy");
      dumpv(fcyv_, 2, "yvz");
    }
  }

  FieldCell<Vect> fcyv_; // Young velocity
}

void Force() {
}

void Init() {
  if (var.Int["wavelamb_vort"] && sem("wavelamb")) {
    FieldCell<Vect> fc(m);
    Vars vr;
    vr.String.Set("vel_init", "wavelamb");
    vr.Double.Set("wavelamb_a0", var.Double["wavelamb_a0"]);
    vr.Double.Set("wavelamb_xc", var.Double["wavelamb_xc"]);
    vr.Double.Set("wavelamb_h", var.Double["wavelamb_h"]);
    vr.Double.Set("wavelamb_k", var.Double["wavelamb_k"]);
    vr.Vect.Set("gravity", var.Vect["gravity"]);
    InitVel(fc, vr, m);
    for (auto c : m.AllCells()) {
      if (fc[c].sqrnorm() != 0) {
        fcvel[c] = fc[c];
      }
    }
  }

  if (var.Int["vel_init_noise"] && sem("noise")) {
    Vect vel0(var.Vect["noise_vel0"]);
    Vect vel1(var.Vect["noise_vel1"]);
    Vect vel2(var.Vect["noise_vel2"]);
    Vect per0(var.Vect["noise_per0"]);
    Vect per1(var.Vect["noise_per1"]);
    Vect per2(var.Vect["noise_per2"]);
    Vect k0 = Vect(2 * M_PI) / (per0 * m.GetCellSize());
    Vect k1 = Vect(2 * M_PI) / (per1 * m.GetCellSize());
    Vect k2 = Vect(2 * M_PI) / (per2 * m.GetCellSize());
    for (auto c : m.AllCells()) {
      auto x = m.GetCenter(c);
      fcvel[c] += vel0 * std::sin(k0.dot(x));
      fcvel[c] += vel1 * std::sin(k1.dot(x));
      fcvel[c] += vel2 * std::sin(k2.dot(x));
    }
  }
}

void CalcStat() {
  if (sem("vfslip")) {
    auto kslip = var.Double["kslip"];
    if (kslip != 0) {
      Vect slipvel(var.Vect["slipvel"]);
      // XXX: adhoc, overwrite wall conditions
      auto& fa = as_->GetField();
      for (auto& p : mebc_fluid_.GetMapFace()) {
        const IdxFace f = p.first;
        auto& bc = p.second;
        if (bc.type == BCondFluidType::wall) {
          const IdxCell c = m.GetCell(f, bc.nci);
          bc.velocity = slipvel * kslip * fa[c];
        }
      }
    }

    // Slip velocity penalization.
    const Scal penalslip = var.Double["penalslip"];
    if (penalslip != 0) {
      const Scal dt = fs_->GetTimeStep();
      const Vect slipvel(var.Vect["slipvel"]);
      // const auto& fa = as_->GetField();
      const auto& fa = fc_smvf_;
      const auto& fv = fs_->GetVelocity();
      for (auto& p : mebc_fluid_.GetMapFace()) {
        const IdxFace f = p.first;
        const auto& bc = p.second;
        if (bc.type == BCondFluidType::wall) {
          const IdxCell c = m.GetCell(f, bc.nci);
          Scal sgn = (slipvel - fv[c]).dot(slipvel);
          if (sgn > 0) {
            fc_force_[c] +=
                (slipvel - fv[c]) * (fc_rho_[c] * penalslip * fa[c] / dt);
          }
        }
      }
    }
    // Repulsive force from walls.
    const Scal slipnormal = var.Double["slipnormal"];
    if (slipnormal != 0) {
      const Scal slipnormal_dist = var.Double["slipnormal_dist"];
      const Scal dt = fs_->GetTimeStep();
      const Scal h = m.GetCellSize()[0];
      const auto& fcvf = fc_smvf_;
      if (eb_) {
        auto& eb = *eb_;
        fc_dist_.Reinit(eb, 0);
        for (auto c : eb.SuCells()) {
          fc_dist_[c] = eb.GetSignedDistance(c);
        }
        fc_phi_.Reinit(eb, 0); // potential [length]
        for (auto c : eb.SuCells()) {
          if (!IsNan(fc_dist_[c])) {
            const Scal d0 = slipnormal_dist * h;
            const Scal d = std::max(0., d0 - fc_dist_[c]);
            fc_phi_[c] += slipnormal * d;
          }
        }
        const auto ffg = UEB::Gradient(fc_phi_, {}, eb); // potential grad [-]
        const auto ff_rho = UEB::InterpolateHarmonic(fc_rho_, {}, eb);
        if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
          for (auto f : eb.Faces()) {
            const IdxCell cm = m.GetCell(f, 0);
            const IdxCell cp = m.GetCell(f, 1);
            if (eb.GetType(cm) == EB::Type::excluded ||
                eb.GetType(cp) == EB::Type::excluded) {
              continue;
            }
            const auto& fccl = as->GetColor();
            const auto& fcu = as->GetFieldM();
            std::set<Scal> colors;
            for (auto l : layers) {
              const Scal clm = (*fccl[l])[cm];
              const Scal clp = (*fccl[l])[cp];
              if (clm != kClNone) colors.insert(clm);
              if (clp != kClNone) colors.insert(clp);
            }
            for (auto cl : colors) {
              Scal um = 0;
              Scal up = 0;
              for (auto l : layers) {
                if ((*fccl[l])[cm] == cl) um = (*fcu[l])[cm];
                if ((*fccl[l])[cp] == cl) up = (*fcu[l])[cp];
              }
              const Scal uf = (um + up) * 0.5;
              if (uf != 0) {
                febp_[f] -= ff_rho[f] * ffg[f] * uf * h / sqr(dt);
              }
            }
          }
        }
      } else {
        for (auto& p : mebc_fluid_.GetMapFace()) {
          const IdxFace f = p.first;
          const auto& bc = p.second;
          const auto nci = bc.nci;
          if (bc.type == BCondFluidType::wall) {
            const IdxCell c = m.GetCell(f, bc.nci);
            fc_force_[c] += m.GetNormal(f) * ((nci == 1 ? 1 : -1) * fc_rho_[c] *
                                              slipnormal * fcvf[c]);
          }
        }
      }
    }
  }
}

void CalcMixture() {
  if (sem.Nested("smoothnode")) {
    SmoothenNode(fc_smvf_, m, var.Int["vfsmoothnode"]);
  }

  if (sem("calc")) {
    // vortex force
    Scal force_vort = var.Double["force_vort_g"];
    if (force_vort != 0) {
      Scal r = var.Double["force_vort_r"];
      Vect xc(var.Vect["force_vort_c"]);
      for (auto c : m.Cells()) {
        Vect x = m.GetCenter(c);
        Vect dx = x - xc;
        Scal q = std::exp(-dx.sqrnorm() / sqr(r)) * force_vort;
        Scal fx = -dx[1];
        Scal fy = dx[0];
        fc_force_[c] = Vect(fx, fy, 0.) * q * fc_rho_[c];
      }
    }

    // moving force on the interface
    auto fmov = [this, &a](std::string pre) {
      Vect force_mov(var.Vect[pre]);
      if (force_mov.sqrnorm()) {
        Vect x0(var.Vect[pre + "_x0"]);
        Vect x1(var.Vect[pre + "_x1"]);
        Scal t0 = var.Double[pre + "_t0"];
        Scal t1 = var.Double[pre + "_t1"];
        Vect sig(var.Vect[pre + "_sig"]);
        Scal pi = M_PI;
        int edim = var.Int["dim"];

        Scal t = (st_.t - t0) / (t1 - t0);
        if (t >= 0 && t <= 1) {
          if (var.Int[pre + "_parab"]) {
            t = t * t;
          }
          for (auto c : m.Cells()) {
            Vect xt = x0 + (x1 - x0) * t;
            Vect r = (xt - m.GetCenter(c)) / sig;
            Scal k = std::exp(-r.sqrnorm() * 0.5) /
                     (std::pow(2 * pi, 1. / edim) * sig.prod());
            Scal vf = a[c];
            fc_force_[c] += force_mov * (k * fc_rho_[c] * vf * (1 - a[c]) * 2);
          }
        }
      }
    };

    fmov("force_mov");
    fmov("force_mov2");
    // Kolmogorov forcing
    if (var.Int["force_kolm"]) {
      for (auto f : m.AllFaces()) {
        Vect n = m.GetNormal(f);
        Vect x = m.GetCenter(f);
        febp_[f] += Vect(std::sin(x[1]), 0., 0.).dot(n);
      }
    }

    // Kolmogorov forcing as acceleration
    if (var.Int["force_kolm_accel"]) {
      for (auto f : m.AllFaces()) {
        Vect n = m.GetNormal(f);
        Vect x = m.GetCenter(f);
        febp_[f] += Vect(std::sin(x[1]), 0., 0.).dot(n) * ff_rho[f];
      }
    }
  }
}
