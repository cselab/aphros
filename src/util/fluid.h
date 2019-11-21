#pragma once

#include "geom/mesh.h"
#include "solver/fluid.h"
#include "solver/cond.h"

template <class M_>
class UFluid {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  // TODO: Consider seperate channels in one domain
  // fcw: velocity
  // mfc: fluid face conditions
  // fcsv: volume source
  static void UpdateOutletBaseConditions(
      M& m, const FieldCell<Vect>& fcw,
      MapCondFaceFluid& mfc, const FieldCell<Scal>& fcsv) {
    using namespace solver;
    using namespace solver::fluid_condition;

    auto sem = m.GetSem("outlet");

    struct {
      Scal fi; // total inlet volume flux
      Scal fo; // total outlet volume flux
      Scal ao; // total outlet area
    }* ctx(sem);
    auto& fi = ctx->fi;
    auto& fo = ctx->fo;
    auto& ao = ctx->ao;

    if (sem("local")) {
      fi = 0.;
      fo = 0.;
      ao = 0.;

      // Extrapolate velocity to outlet from neighbour cells,
      // and compute total fluxes
      for (auto p : mfc) {
        IdxFace f = p.GetIdx();
        auto& cb = p.GetValue(); // cond base

        size_t id = cb->GetNci();
        IdxCell c = m.GetCell(f, id);
        if (m.IsInner(c)) {
          if (auto cd = cb.Get<Outlet<M>>()) {
            Scal w = (id == 0 ? 1. : -1.);
            Vect vc = fcw[c];
            Vect s = m.GetSurface(f);
            // clip normal component, let only positive
            vc -= s * (w * std::min(0., vc.dot(s) * w)  / s.dot(s));
            cd->SetVelocity(vc);
            fo += cd->GetVelocity().dot(s) * w;
            ao += m.GetArea(f);
          } else if (auto cd = cb.Get<Inlet<M>>()) {
            Scal w = (id == 0 ? -1. : 1.);
            fi += cd->GetVelocity().dot(m.GetSurface(f)) * w;
          }
        }
      }

      // Append volume source to inlet flux
      for (auto c : m.Cells()) {
        fi += fcsv[c] * m.GetVolume(c);
      }

      m.Reduce(&fi, "sum");
      m.Reduce(&fo, "sum");
      m.Reduce(&ao, "sum");
    }

    if (sem("corr")) {
      Scal velcor = (fi - fo) / ao; // Additive correction for velocity

      // Apply correction on outlet faces
      for (auto it : mfc) {
        IdxFace f = it.GetIdx();
        auto& cb = it.GetValue(); // cond base

        if (auto cd = cb.Get<Outlet<M>>()) {
          size_t id = cd->GetNci();
          Scal w = (id == 0 ? 1. : -1.);
          Vect n = m.GetNormal(f);
          cd->SetVelocity(cd->GetVelocity() + n * (velcor * w));
        }
      }
    }
  }
};
