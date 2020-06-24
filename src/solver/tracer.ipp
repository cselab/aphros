// Created by Petr Karnakov on 24.06.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>

#include "approx.h"
#include "approx_eb.h"
#include "util/vof.h"

#include "tracer.h"

template <class EB_>
struct Tracer<EB_>::Imp {
  using Owner = Tracer<EB_>;
  using UEB = UEmbed<M>;

  Imp(Owner* owner, M& m, const EB& eb_,
      const Multi<const FieldCell<Scal>*>& vfcvf,
      const MapEmbed<BCond<Scal>>& mebc, Scal time, Conf conf_) {
      : owner_(owner)
      , conf(conf_)
      , m(m)
      , eb(eb_)
      , layers(conf.layers)
      , mebc_(mebc)
      , time_(time) {
        fcu_.time_curr = fcu;
      }
      void AdvAulisa(Sem & sem) {
        // directions, format: {dir EI, dir LE, ...}
        std::vector<size_t> dd;
        Scal vsc; // scaling factor for time step
        if (conf.dim == 3) { // 3d
          if (count_ % 3 == 0) {
            dd = {0, 1, 1, 2, 2, 0};
          } else if (count_ % 3 == 1) {
            dd = {1, 2, 2, 0, 0, 1};
          } else {
            dd = {2, 0, 0, 1, 1, 2};
          }
          vsc = 0.5;
        } else { // 2d
          if (count_ % 2 == 0) {
            dd = {0, 1};
          } else {
            dd = {1, 0};
          }
          vsc = 1.0;
        }
        for (size_t id = 0; id < dd.size(); ++id) {
          size_t d = dd[id]; // direction as index
          if (sem("copyface")) {
            if (id % 2 == 1) { // copy fluxes for Lagrange Explicit step
              auto& ffv =
                  owner_->fev_->GetFieldFace(); // [f]ield [f]ace [v]olume flux
              fcfm_.Reinit(m);
              fcfp_.Reinit(m);
              for (auto c : eb.Cells()) {
                fcfm_[c] = ffv[m.GetFace(c, 2 * d)];
                fcfp_[c] = ffv[m.GetFace(c, 2 * d + 1)];
              }
              m.Comm(&fcfm_);
              m.Comm(&fcfp_);
            }
          }
          auto& uc = fcu_.iter_curr;
          if (sem("sweep")) {
            Sweep(
                uc, d, owner_->fev_->GetFieldFace(), fccl_, fcim_, fcn_, fca_,
                &me_vf_, id % 2 == 0 ? SweepType::EI : SweepType::LE, &fcfm_,
                &fcfp_, nullptr, owner_->GetTimeStep() * vsc, conf.clipth, eb);
          }
          CommRec(sem, uc, fccl_, fcim_);
        }
      }
      static void BcMarchFill(FieldCell<Scal> & fcu, Scal fill, const M& m) {
        for (auto c : m.AllCells()) {
          auto x = m.GetCenter(c);
          if (!(Vect(0) <= x && x <= m.GetGlobalLength())) {
            fcu[c] = fill;
          }
        }
      }
      void AdvPlain(Sem & sem, SweepType type) {
        std::vector<size_t> dd; // sweep directions
        if (conf.dim == 3) { // 3d
          if (count_ % 3 == 0) {
            dd = {0, 1, 2};
          } else if (count_ % 3 == 1) {
            dd = {1, 2, 0};
          } else {
            dd = {2, 0, 1};
          }
        } else { // 2d
          if (count_ % 2 == 0) {
            dd = {0, 1};
          } else {
            dd = {1, 0};
          }
        }
        for (size_t id = 0; id < dd.size(); ++id) {
          auto& uc = fcu_.iter_curr;
          if (sem("sweep")) {
            Sweep(
                uc, dd[id], owner_->fev_->GetFieldFace(), fccl_, fcim_, fcn_,
                fca_, &me_vf_, type, nullptr, nullptr, &fcuu_,
                owner_->GetTimeStep(), conf.clipth, eb);
          }
          CommRec(sem, uc, fccl_, fcim_);
          if (conf.extrapolate_boundaries) {
            ExtrapolatePlic(sem, uc);
          }
        }
      }
      void Sharpen() {
        auto sem = m.GetSem("sharp");
        std::vector<size_t> dd; // sweep directions
        if (conf.dim == 3) { // 3d
          if (count_ % 3 == 0) {
            dd = {0, 0, 1, 1, 2, 2};
          } else if (count_ % 3 == 1) {
            dd = {1, 1, 2, 2, 0, 0};
          } else {
            dd = {2, 2, 0, 0, 1, 1};
          }
        } else { // 2d
          if (count_ % 2 == 0) {
            dd = {0, 0, 1, 1};
          } else {
            dd = {1, 1, 0, 0};
          }
        }
        for (size_t id = 0; id < dd.size(); ++id) {
          size_t d = dd[id]; // direction as index
          auto& uc = fcu_.iter_curr;
          if (sem("sweep")) {
            const Scal sgn = (id % 2 == count_ / conf.dim % 2 ? -1 : 1);
            FieldFace<Scal> ffv(m, 0);
            for (auto f : eb.Faces()) {
              const IdxCell cm = m.GetCell(f, 0);
              const IdxCell cp = m.GetCell(f, 1);
              ffv[f] = std::min(eb.GetVolume(cm), eb.GetVolume(cp)) * sgn *
                       conf.sharpen_cfl;
            }
            // zero flux on boundaries
            for (const auto& it : me_vf_.GetMapFace()) {
              ffv[it.first] = 0;
            }
            Sweep(
                uc, d, ffv, fccl_, fcim_, fcn_, fca_, &me_vf_,
                SweepType::weymouth, nullptr, nullptr, &fcuu_, 1., conf.clipth,
                eb);
          }
          CommRec(sem, uc, fccl_, fcim_);
        }
      }
      void Step(
          Scal dt, const FieldEmbed<Scal>& fe_flux,
          const Multi<FieldCell<Scal>*>& vfc_src) {
        auto sem = m.GetSem("iter");
        if (sem("init")) {
        }
        if (sem("stat")) {
          time_ += dt;
        }
      }

      Owner* owner_;
      Conf conf;
      M& m;
      const EB& eb;
      GRange<size_t> layers;

      StepData<FieldCell<Scal>> fcu_;
      FieldCell<Scal> fcuu_; // volume fraction for Weymouth div term

      // boundary conditions
      MapEmbed<BCond<Scal>> mebc_; // field
  };

  template <class EB_>
  Tracer<EB_>::Tracer(
      M& m, const EB& eb, const FieldCell<Scal>& fcu,
      const FieldCell<Scal>& fccl, const MapEmbed<BCondAdvection<Scal>>& mfc,
      const FieldEmbed<Scal>* fev, const FieldCell<Scal>* fcs, double t,
      double dt, Conf conf)
      : AdvectionSolver<M>(t, dt, m, fev, fcs)
      , imp(new Imp(this, m, eb, vfcu, mebc, t, conf)) {}

  template <class EB_>
  Tracer<EB_>::~Tracer() = default;

  template <class EB_>
  auto Tracer<EB_>::GetEmbed() const -> const EB& {
    return imp->eb;
  }

  template <class EB_>
  auto Tracer<EB_>::GetConf() const -> const Conf& {
    return imp->conf;
  }

  template <class EB_>
  void Tracer<EB_>::SetConf(Conf conf) {
    imp->conf = conf;
  }
};

template <class EB_>
auto Tracer<EB_>::GetField(Step l) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
}

template <class EB_>
auto Tracer<EB_>::GetView() const -> Plic {
  return {imp->layers, &GetField(), &GetAlpha(), &GetNormal(),
          &GetMask(),  nullptr,     imp->mfc_};
}
