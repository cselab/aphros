// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <map>
#include <memory>

#include "approx.h"
#include "debug/isnan.h"
#include "embed.h"
#include "geom/block.h"
#include "normal.h"
#include "parse/curv.h"
#include "reconst.h"
#include "trackerm.h"
#include "util/format.h"
#include "util/height.h"
#include "util/vof.h"
#include "vof.h"
#include "vofm.h"

#include "curv.h"

namespace curvature {

template <class M_>
struct Particles<M_>::Imp {
  Imp(M& m, const typename PartStrMeshM<M>::Par& par,
      const GRange<size_t>& layers)
      : partstrmeshm_(new PartStrMeshM<M>(m, par, layers)) {}

  template <class EB>
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const EB& eb) {
    const auto& layers = plic.layers;
    auto sem = m.GetSem();

    if (sem.Nested("part")) {
      partstrmeshm_->Part(plic, eb);
    }
    if (sem("copy")) {
      fck.assert_size(layers);
      for (auto l : layers) {
        (*fck[l]) = *partstrmeshm_->GetCurv()[l];
      }
    }
  }

  std::unique_ptr<PartStrMeshM<M>> partstrmeshm_;
};

template <class M_>
Particles<M_>::Particles(
    M& m, const typename PartStrMeshM<M>::Par& par,
    const GRange<size_t>& layers)
    : imp(new Imp(m, par, layers)) {}

template <class EB_>
Particles<EB_>::~Particles() = default;

template <class M_>
void Particles<M_>::CalcCurvature(
    const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m, const M& eb) {
  imp->CalcCurvature(fck, plic, m, eb);
}

template <class M_>
void Particles<M_>::CalcCurvature(
    const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
    const Embed<M>& eb) {
  imp->CalcCurvature(fck, plic, m, eb);
}

template <class M_>
std::unique_ptr<PartStrMeshM<M_>> Particles<M_>::ReleaseParticles() {
  return std::move(imp->partstrmeshm_);
}

template <class M_>
const PartStrMeshM<M_>* Particles<M_>::GetParticles() const {
  return imp->partstrmeshm_.get();
}

template <class M_>
void Particles<M_>::DumpAux(std::string, int, M&) {}

template <class M_>
struct Heights<M_>::Imp {
  Imp() {}
  static constexpr Scal kClNone{Vofm<M>::kClNone};

  // fcu: volume fraction [a]
  // Output:
  // fcdu2: difference of volume fractions with step 2 (u[xpp]-u[xmm], ...) [i]
  // fccl: updated with colors from neighbors
  static void CalcDiff2(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Vect>*>& fcdu2, const Multi<FieldCell<Scal>*>& fccl,
      const M& m) {
    for (auto l : layers) {
      fcdu2[l]->Reinit(m, Vect(0));
    }
    for (auto c : m.CellsM()) {
      for (auto dr : m.dirs) {
        const auto d = m.direction(dr);
        const auto cmm = c - d - d;
        const auto cpp = c + d + d;
        for (int sign : {-1, 1}) {
          const auto cn = sign > 0 ? cpp : cmm;
          for (auto ln : layers) {
            const auto cl = (*fccl[ln])[cn];
            if (cl != kClNone) {
              for (auto l : layers) {
                if ((*fccl[l])[c] == cl) {
                  (*fcdu2[l])[c][d] += sign * (*fcu[ln])[cn];
                }
              }
              for (auto l : layers) {
                if ((*fccl[l])[c] == kClNone) {
                  (*fcdu2[l])[c][d] += sign * (*fcu[ln])[cn];
                  (*fccl[l])[c] = cl;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }
  // fcu: volume fraction [a]
  // fcdu2: difference of volume fractions with step 2 (u[xpp]-u[xmm], ...) [a]
  // Output:
  // fcdu4: difference of volume fractions with step 4 (u[xp4]-u[xm4], ...) [i]
  // fccl: updated with colors from neighbors
  static void CalcDiff4(
      const GRange<size_t>& layers, const Multi<const FieldCell<Vect>*>& fcdu2,
      const Multi<FieldCell<Vect>*>& fcdu4, const Multi<FieldCell<Scal>*>& fccl,
      const M& m) {
    for (auto l : layers) {
      fcdu4[l]->Reinit(m, Vect(0));
    }
    for (auto c : m.CellsM()) {
      for (auto dr : m.dirs) {
        const auto d = m.direction(dr);
        const auto cmm = c - d - d;
        const auto cpp = c + d + d;
        for (int sign : {-1, 1}) {
          const auto cn = sign > 0 ? cpp : cmm;
          for (auto ln : layers) {
            const auto cl = (*fccl[ln])[cn];
            if (cl != kClNone) {
              for (auto l : layers) {
                if ((*fccl[l])[c] == cl) {
                  (*fcdu4[l])[c][d] += sign * (*fcdu2[ln])[cn][d];
                }
              }
              for (auto l : layers) {
                if ((*fccl[l])[c] == kClNone) {
                  (*fcdu4[l])[c][d] += sign * (*fcdu2[ln])[cn][d];
                  (*fccl[l])[c] = cl;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  // Computes height functions in interfacial cells.
  // fcu: volume fraction [a]
  // fcdu2: difference of volume fractions with step 2 [a]
  // fcdu4: difference of volume fractions with step 4 [a]
  // fccl: colors
  // Output:
  // fch: fch[c][d] is absolute position of the interface
  //      from a column in direction d starting from an interfacial cell c
  //      NaN if undefined
  static void CalcHeight(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Vect>*>& fcdu2,
      const Multi<const FieldCell<Vect>*>& fcdu4,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<FieldCell<Vect>*>& fch, M& m) {
    const size_t S = 2;
    auto I = [](Scal a) { return a > 0 && a < 1; }; // interface

    fch.assert_size(layers);
    for (auto l : layers) {
      fch[l]->Reinit(m, GetNan<Vect>());
    }

    auto get = [&](auto& fc, Scal cl, IdxCell c) {
      for (auto l : layers) {
        if ((*fccl[l])[c] == cl) {
          return (*fc[l])[c];
        }
      }
      using T = decltype((*fc[0])[c]);
      return T{0};
    };

    for (auto l : layers) {
      for (auto c : m.CellsM()) {
        const Scal cl = (*fccl[l])[c];
        if (cl == kClNone || !I((*fcu[l])[c])) {
          continue;
        }
        for (size_t dr = 0; dr < m.GetEdim(); ++dr) {
          const auto d = m.direction(dr);
          const auto cm = c - d;
          const auto cp = c + d;
          const auto cmm = c - d - d;
          const auto cpp = c + d + d;

          const size_t si = (S + 1) * 4 + 1;
          const size_t sih = (S + 1) * 2;

          std::array<Scal, si> uu;

          uu[sih] = get(fcu, cl, c);
          uu[sih - 1] = get(fcu, cl, cm);
          uu[sih + 1] = get(fcu, cl, cp);
          uu[sih - 2] = get(fcu, cl, cmm);
          uu[sih + 2] = get(fcu, cl, cpp);

          for (size_t s = 0; s < S; ++s) {
            const size_t q = (s + 1) * 2;
            const size_t ia = q + 1;
            const size_t ib = q + 2;
            auto& fcdu = s == 0 ? fcdu2 : fcdu4;

            uu[sih - ia] = uu[sih - ia + 2 * q] - get(fcdu, cl, cm)[d];
            uu[sih + ia] = uu[sih + ia - 2 * q] + get(fcdu, cl, cp)[d];
            uu[sih - ib] = uu[sih - ib + 2 * q] - get(fcdu, cl, cmm)[d];
            uu[sih + ib] = uu[sih + ib - 2 * q] + get(fcdu, cl, cpp)[d];
          }

          // |cm6|cm5|cm4|cm3|cmm| cm| c |cp |cpp|cp3|cp4|cp5|cp6|
          // |   |   |   |   | * |   | c |   | * |   |   |   |   |
          // |   |   | * |   |   |   | c |   |   |   | * |   |   |

          const Scal offset = UHeight<Scal>::Good(uu);
          (*fch[l])[c][d] = c.center[d] + offset * m.GetCellSize()[d];
        }
      }
    }
  }

  // Computes curvature from height functions.
  // fcu: volume fraction
  // fcn: normal, antigradient of fcu
  // Output: modified in cells with fci=1, resized to m
  // fck: curvature [i]
  static void CalcCurvHeight(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Vect>*>& fch,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<FieldCell<Scal>*>& fck, M& m) {
    auto I = [](Scal a) { return a > 0 && a < 1; }; // interface

    auto shift = [](auto c, auto d, int i) {
      return i == 0 ? c : //
                 i == 1 ? c + d : //
                     i == -1 ? c - d : //
                         i == 2 ? c + d + d : c - d - d;
    };

    fck.assert_size(layers);
    for (auto l : layers) {
      fck[l]->Reinit(m, GetNan<Scal>());
    }

    for (auto l : layers) {
      for (auto c : m.CellsM()) {
        const Scal cl = (*fccl[l])[c];
        if (cl == kClNone || !I((*fcu[l])[c])) {
          continue;
        }

        // best direction
        const auto dn = m.direction((*fcn[l])[c].abs().argmax());
        // directions of tangents
        const auto dx = dn.next(1);
        const auto dy = dn.next(2);

        // mesh step
        const Scal lx = m.GetCellSize()[dx];
        const Scal ly = m.GetCellSize()[dy];

        // Evaluates height function from nearest interface
        // o: offset from w
        auto hh = [&](auto cb) -> Scal {
          const int si = 5; // stencil size
          const int sih = si / 2;
          int i = sih; // closest interface to center
          while (i < si) {
            const auto cn = shift(cb, dn, i - sih);
            for (auto ln : layers) {
              if ((*fccl[ln])[cn] == cl && I((*fcu[ln])[cn])) {
                return (*fch[ln])[cn][dn];
              }
            }
            if (i > sih) {
              i = si - i - 1;
            } else {
              i = si - i;
            }
          }
          return GetNan<Scal>();
        };

        // height function
        const Scal hcc = hh(c);
        const Scal hmc = hh(c - dx);
        const Scal hpc = hh(c + dx);
        const Scal hcm = hh(c - dy);
        const Scal hcp = hh(c + dy);
        // corners: hxy
        const Scal hmm = hh(c - dx - dy);
        const Scal hmp = hh(c - dx + dy);
        const Scal hpm = hh(c + dx - dy);
        const Scal hpp = hh(c + dx + dy);

        // first derivative (slope)
        const Scal hx = (hpc - hmc) / (2 * lx);
        const Scal hy = (hcp - hcm) / (2 * ly);
        // second derivative
        const Scal fl = 0.2; // filter factor (Basilisk: fl=0.2)
        const Scal hxx = ((hpm - 2 * hcm + hmm) * fl + (hpc - 2 * hcc + hmc) +
                          (hpp - 2 * hcp + hmp) * fl) /
                         ((1 + 2 * fl) * lx * lx);
        const Scal hyy = ((hmp - 2 * hmc + hmm) * fl + (hcp - 2 * hcc + hcm) +
                          (hpp - 2 * hpc + hpm) * fl) /
                         ((1 + 2 * fl) * ly * ly);
        const Scal hxy = ((hpp - hmp) - (hpm - hmm)) / (4 * lx * ly);
        // curvature
        Scal k =
            (2 * hx * hy * hxy - (sqr(hy) + 1) * hxx - (sqr(hx) + 1) * hyy) /
            std::pow(sqr(hx) + sqr(hy) + 1, 3. / 2.);

        if ((*fcn[l])[c][dn] < 0) {
          k *= -1;
        }

        (*fck[l])[c] = k;
      }
    }
  }

  template <class EB>
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m, const EB&) {
    auto sem = m.GetSem();
    struct {
      Multi<FieldCell<Scal>> fccl; // extended colors, contain colors
                                   // from current cell and neighbors
      Multi<FieldCell<Vect>> fcdu2; // difference of volume
                                    // fractions with step 2
      Multi<FieldCell<Vect>> fcdu4; // difference of volume
                                    // fractions with step 4
      Multi<FieldCell<Vect>> fch; // height function
    } * ctx(sem);
    auto& t = *ctx;
    auto& layers = plic.layers;

    if (sem("init")) {
      t.fccl.resize(layers);
      if (!plic.vfccl[0]) { // no colors defined
        fassert_equal(layers.size(), 1);
        t.fccl[0].Reinit(m, 0);
      } else {
        for (auto l : layers) {
          t.fccl[l] = *plic.vfccl[l];
        }
      }
    }
    if (sem("diff2")) {
      t.fcdu2.resize(layers);
      CalcDiff2(layers, plic.vfcu, t.fcdu2, t.fccl, m);
      for (auto l : layers) {
        m.Comm(&t.fcdu2[l]);
        m.Comm(&t.fccl[l]);
      }
    }
    if (sem("diff4")) {
      t.fcdu4.resize(layers);
      CalcDiff4(layers, t.fcdu2, t.fcdu4, t.fccl, m);
      for (auto l : layers) {
        m.Comm(&t.fcdu4[l]);
        m.Comm(&t.fccl[l]);
      }
    }
    if (sem("height")) {
      t.fch.resize(layers);
      CalcHeight(layers, plic.vfcu, t.fcdu2, t.fcdu4, t.fccl, t.fch, m);
      for (auto l : layers) {
        m.Comm(&t.fch[l]);
        m.Comm(&t.fccl[l]);
      }
    }
    if (sem("curvcomm")) {
      CalcCurvHeight(layers, plic.vfcu, t.fch, plic.vfcn, t.fccl, fck, m);
      for (auto l : layers) {
        m.Comm(fck[l]);
      }
    }
    if (sem("savelast")) {
      fc_height_ = t.fch;
      fccl_ = t.fccl;
    }
  }

  Multi<FieldCell<Vect>> fc_height_; // last height to be used in DumpAux
  Multi<FieldCell<Scal>> fccl_; // last colors to be used in DumpAux
};

template <class M_>
Heights<M_>::Heights() : imp(new Imp()) {}

template <class EB_>
Heights<EB_>::~Heights() = default;

template <class M_>
void Heights<M_>::CalcCurvature(
    const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m, const M& eb) {
  imp->CalcCurvature(fck, plic, m, eb);
}

template <class M_>
void Heights<M_>::CalcCurvature(
    const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
    const Embed<M>& eb) {
  imp->CalcCurvature(fck, plic, m, eb);
}

template <class M_>
void Heights<M_>::DumpAux(std::string, int frame, M& m) {
  auto sem = m.GetSem();
  struct {
    std::map<std::string, std::vector<Scal>> data;
  } * ctx(sem);
  auto& t = *ctx;
  constexpr Scal kClNone{Vofm<M>::kClNone};
  if (sem()) {
    auto& d_x = t.data["x"];
    auto& d_y = t.data["y"];
    auto& d_z = t.data["z"];
    auto& d_hx = t.data["hx"];
    auto& d_hy = t.data["hy"];
    auto& d_hz = t.data["hz"];
    auto& d_dir = t.data["dir"];
    auto& d_cl = t.data["cl"];
    const GRange<size_t> layers(imp->fc_height_.size());
    using Vect3 = generic::Vect<Scal, 3>;
    for (auto l : layers) {
      for (auto c : m.CellsM()) {
        auto cl = imp->fccl_[l][c];
        const Vect3 h(imp->fc_height_[l][c]);
        const Vect3 center(c.center());
        if (cl != kClNone) {
          for (auto d : m.dirs) {
            if (!IsNan(h[d])) {
              d_x.push_back(center[0]);
              d_y.push_back(center[1]);
              d_z.push_back(center[2]);
              d_dir.push_back(d);
              auto hh = center;
              hh[d] = h[d];
              d_hx.push_back(hh[0]);
              d_hy.push_back(hh[1]);
              d_hz.push_back(hh[2]);
              d_cl.push_back(cl);
            }
          }
        }
      }
    }
    for (auto& p : t.data) {
      m.Reduce(&p.second, Reduction::concat);
    }
  }
  if (sem()) {
    if (m.IsRoot()) {
      std::ofstream csv(util::Format("h_{:04d}.csv", frame));
      size_t size = 0;
      { // write header
        bool first = true;
        for (auto& p : t.data) {
          if (!first) {
            csv << ',';
          } else {
            first = false;
          }
          csv << p.first;
          size = p.second.size();
        }
        csv << '\n';
      }
      { // write rows
        for (size_t i = 0; i < size; ++i) {
          bool first = true;
          for (auto& p : t.data) {
            if (!first) {
              csv << ',';
            } else {
              first = false;
            }
            csv << p.second[i];
          }
          csv << '\n';
        }
      }
    }
  }
}

template <class M_>
struct Hybrid<M_>::Imp {
  Imp(M& m, const typename PartStrMeshM<M>::Par& par,
      const GRange<size_t>& layers)
      : heights_(new Heights<M>())
      , particles_(new Particles<M>(m, par, layers)) {}

  template <class EB>
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const EB& eb) {
    const auto& layers = plic.layers;
    auto sem = m.GetSem();
    struct {
      Multi<FieldCell<Scal>> fck;
    } * ctx(sem);
    auto& t = *ctx;
    if (sem.Nested("part")) {
      heights_->CalcCurvature(fck, plic, m, eb);
    }
    // TODO: only seed particles in cells with undefined curvature from heights
    if (sem("resize")) {
      t.fck.resize(layers);
    }
    if (sem.Nested("part")) {
      particles_->CalcCurvature(t.fck, plic, m, eb);
    }
    if (sem("update")) {
      for (auto l : layers) {
        for (auto c : m.AllCells()) {
          if (IsNan((*fck[l])[c])) {
            (*fck[l])[c] = t.fck[l][c];
          }
        }
      }
    }
  }

  std::unique_ptr<Heights<M>> heights_;
  std::unique_ptr<Particles<M>> particles_;
};

template <class M_>
Hybrid<M_>::Hybrid(
    M& m, const typename PartStrMeshM<M>::Par& par,
    const GRange<size_t>& layers)
    : imp(new Imp(m, par, layers)) {}

template <class EB_>
Hybrid<EB_>::~Hybrid() = default;

template <class M_>
void Hybrid<M_>::CalcCurvature(
    const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m, const M& eb) {
  imp->CalcCurvature(fck, plic, m, eb);
}

template <class M_>
void Hybrid<M_>::CalcCurvature(
    const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
    const Embed<M>& eb) {
  imp->CalcCurvature(fck, plic, m, eb);
}

template <class M_>
std::unique_ptr<PartStrMeshM<M_>> Hybrid<M_>::ReleaseParticles() {
  return std::move(imp->particles_->ReleaseParticles());
}

template <class M_>
const PartStrMeshM<M_>* Hybrid<M_>::GetParticles() const {
  return imp->particles_->GetParticles();
}

template <class M_>
void Hybrid<M_>::DumpAux(std::string request, int frame, M& m) {
  imp->heights_->DumpAux(request, frame, m);
}

template <class M>
std::unique_ptr<Estimator<M>> MakeEstimator(
    const Vars& var, M& m, const GRange<size_t>& layers) {
  using Scal = typename M::Scal;
  const auto name = var.String("curvature", "particles");
  if (name == "particles") {
    const auto ps = ParsePar<PartStr<Scal>>()(m.GetCellSize()[0], var);
    const auto psm = ParsePar<PartStrMeshM<M>>()(ps, var);
    return std::make_unique<curvature::Particles<M>>(m, psm, layers);
  } else if (name == "heights") {
    return std::make_unique<curvature::Heights<M>>();
  } else if (name == "hybrid") {
    const auto ps = ParsePar<PartStr<Scal>>()(m.GetCellSize()[0], var);
    const auto psm = ParsePar<PartStrMeshM<M>>()(ps, var);
    return std::make_unique<curvature::Hybrid<M>>(m, psm, layers);
  }
  fassert(false, util::Format("Unknown curvature estimator '{}'", name));
}

} // namespace curvature
