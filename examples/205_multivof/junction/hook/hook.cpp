// Created by Petr Karnakov on 30.05.2020
// Copyright 2020 ETH Zurich

#include <fstream>
#include <iostream>
#include <limits>

#include <func/init_u.h>
#include <func/init_vel.h>
#include <kernel/hydro.h>
#include <solver/multi.h>
#include <util/format.h>
#include <util/posthook.h>

// listpath: path to list of primitives b.dat
template <class M>
void InitVolumeFractionFromList(
    std::string listpath, const Multi<FieldCell<typename M::Scal>*>& fcu,
    const Multi<FieldCell<typename M::Scal>*>& fccl,
    const GRange<size_t>& layers, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Primitive = typename UPrimList<Vect>::Primitive;
  const auto kClNone = Vofm<M>::kClNone;
  std::ifstream fin(listpath);
  const std::vector<Primitive> primitives =
      UPrimList<Vect>::GetPrimitives(fin, m.flags.edim);
  for (auto l : layers) {
    fcu[l]->Reinit(m, 0);
    fccl[l]->Reinit(m, kClNone);
  }
  const auto h = m.GetCellSize();
  for (auto c : m.AllCellsM()) {
    const size_t none = -1;
    size_t imax = none; // primitive with largest volume fraction in c
    size_t imax2 = none; // primitive with second largest volume fraction in c
    Scal lsmax = -std::numeric_limits<Scal>::max();
    Scal lsmax2 = -std::numeric_limits<Scal>::max();
    for (size_t i = 0; i < primitives.size(); ++i) {
      const Scal ls = primitives[i].ls(c.center);
      if (ls > lsmax) {
        imax = i;
        lsmax = ls;
      }
    }
    for (size_t i = 0; i < primitives.size(); ++i) {
      const Scal ls = primitives[i].ls(c.center);
      if (ls > lsmax2 && i != imax) {
        imax2 = i;
        lsmax2 = ls;
      }
    }
    const Scal umax =
        (imax == none ? 0
                      : GetLevelSetVolume<Scal, M::dim>(
                            primitives[imax].ls, c.center, h));
    const Scal umax2 =
        (imax2 == none ? 0
                       : GetLevelSetVolume<Scal, M::dim>(
                             primitives[imax2].ls, c.center, h));
    if (umax + umax2 >= 1) { // two or more positive level-sets
      auto f = [&](const Vect& x) -> Scal {
        return primitives[imax].ls(x) - primitives[imax2].ls(x);
      };
      auto f2 = [&](const Vect& x) -> Scal {
        return primitives[imax2].ls(x) - primitives[imax].ls(x);
      };
      const Scal u = GetLevelSetVolume<Scal, M::dim>(f, c.center, h);
      const Scal u2 = GetLevelSetVolume<Scal, M::dim>(f2, c.center, h);
      if (u > 0) {
        for (auto l : layers) {
          if ((*fccl[l])[c] == kClNone) {
            (*fccl[l])[c] = imax;
            (*fcu[l])[c] = u;
            break;
          }
        }
      }
      if (u2 > 0) {
        for (auto l : layers) {
          if ((*fccl[l])[c] == kClNone) {
            (*fccl[l])[c] = imax2;
            (*fcu[l])[c] = u2;
            break;
          }
        }
      }
    } else if (1 || umax > 0) {
      (*fccl[0])[c] = imax;
      (*fcu[0])[c] =
          GetLevelSetVolume<Scal, M::dim>(primitives[imax].ls, c.center, h);
    } else if (umax2 > 0) {
      (*fccl[0])[c] = imax2;
      (*fcu[0])[c] =
          GetLevelSetVolume<Scal, M::dim>(primitives[imax2].ls, c.center, h);
    }
  }
}

template <class M>
class Vortex : public ModuleInitVelocity<M> {
 public:
  using Vect = typename M::Vect;
  Vortex() : ModuleInitVelocity<M>("vortex") {}
  void operator()(FieldCell<Vect>& fcv, const Vars&, const M& m) override {
    using Scal = typename M::Scal;
    for (auto c : m.AllCellsM()) {
      const Scal x = c.center[0] * M_PI;
      const Scal y = c.center[1] * M_PI;
      Scal& u = fcv[c][0];
      Scal& v = fcv[c][1];
      u = std::sin(x) * std::cos(y);
      v = -std::cos(x) * std::sin(y);
    }
  }
};

template <class M>
void FluidDummyHook(
    FieldCell<typename M::Vect>& fcvel, FieldFace<typename M::Scal>& ffv,
    typename M::Scal t, typename M::Scal dt, const Vars& var, const M& m) {
  using Scal = typename M::Scal;
  const Scal rt = var.Double["reverse_t"];
  if (t >= rt && t < rt + dt) {
    for (auto f : m.Faces()) {
      ffv[f] *= -1;
    }
    for (auto c : m.AllCells()) {
      fcvel[c] *= -1;
    }
  }
}

using M = MeshCartesian<double, 3>;

bool reg[] = {
    RegisterModule<Vortex<M>>(),
};

template <class M>
void InitHook(Hydro<M>* hydro) {
  auto& var = hydro->var;
  auto mod = [&var](auto& fcu, auto& fccl, auto layers, auto& m) {
    InitVolumeFractionFromList(var.String["list_path"], fcu, fccl, layers, m);
  };
  if (auto* as = dynamic_cast<Vofm<M>*>(hydro->as_.get())) {
    as->AddModifier(mod);
  }
}

template void InitHook(Hydro<M>*);
template void FluidDummyHook(
    FieldCell<typename M::Vect>&, FieldFace<typename M::Scal>&,
    typename M::Scal t, typename M::Scal dt, const Vars&, const M&);
