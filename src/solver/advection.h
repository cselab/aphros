// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include "geom/mesh.h"
#include "multi.h"
#include "solver/embed.h"
#include "solver/solver.h"

template <class Scal>
struct BCondAdvection {
  enum class Halo { fill, reflect };

  BCondAdvection() = default;
  BCondAdvection(size_t nci_) : nci(nci_) {}
  size_t GetNci() const {
    return nci;
  }

  size_t nci;
  Scal clear0 = 0; // snap to 0 if vf<clear0
  Scal clear1 = 1; // snap to 1 if vf>clear1
  Halo halo = Halo::reflect;
  Scal fill_vf; // volume fraction in halo cells if halo=fill
  Scal fill_cl; // color in halo cells if halo=fill
  Scal contang = -1; // contact angle [0..pi], negative to disable
};

namespace generic {

// Piecewise Linear Interface Characterization.
// Describes the multilayer interface.
template <class Vect_>
struct Plic {
  using Vect = Vect_;
  using Scal = typename Vect::Scal;
  static constexpr size_t dim = Vect::dim;
  using MIdx = generic::MIdx<dim>;
  GRange<size_t> layers;
  Multi<const FieldCell<Scal>*> vfcu; // volume fraction
  Multi<const FieldCell<Scal>*> vfca; // plane constant
  Multi<const FieldCell<Vect>*> vfcn; // normal
  Multi<const FieldCell<bool>*> vfci; // interface mask
  Multi<const FieldCell<Scal>*> vfccl; // color
  // XXX Field `vfccl_stat` may have value kClNone even in cells with u>0
  // close to the interface and should be used for statitics only,
  // and not for filtering the cells (e.g. in surface tension).
  // Required for single-field VOF solver where the criterion for merging the
  // cells is more strict to distinguish closely located bubbles.
  Multi<const FieldCell<Scal>*> vfccl_stat; // color for statistics only
  Multi<const FieldCell<MIdx>*> vfcim; // image vector
  const MapEmbed<BCondAdvection<Scal>>& me_adv; // boundary conditions
};

} // namespace generic

template <class M_>
class AdvectionSolver : public UnsteadyIterativeSolver {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Plic = generic::Plic<Vect>;
  // fev: volume flux
  // fcs: source
  // mebc: boundary conditions
  AdvectionSolver(
      double t, double dt, M& m_, const FieldEmbed<Scal>* fev,
      const FieldCell<Scal>* fcs, const MapEmbed<BCondAdvection<Scal>>& mebc)
      : UnsteadyIterativeSolver(t, dt)
      , m(m_)
      , fev_(fev)
      , fcs_(fcs)
      , mebc_(mebc) {}
  // Postprocessing after time step (dumps)
  virtual void PostStep() {}
  // Volume fraction
  virtual const FieldCell<Scal>& GetField(Step) const = 0;
  // Volume fraction at last time step
  virtual const FieldCell<Scal>& GetField() const {
    return GetField(Step::time_curr);
  }
  virtual Plic GetPlic() const = 0;
  virtual void DumpInterface(
      std::string /*filename*/,
      std::vector<Multi<const FieldCell<Scal>*>> /*extra_fields*/,
      std::vector<std::string> /*extra_names*/) const {
    fassert(false, "not implemented");
  };
  virtual void DumpInterfaceMarch(std::string /*filename*/) const {
    fassert(false, "not implemented");
  };

 protected:
  M& m;
  const FieldEmbed<Scal>* fev_; // volume flux [velocity*area]
  const FieldCell<Scal>* fcs_; // source [value/time]
  const MapEmbed<BCondAdvection<Scal>>& mebc_; // boundary conditions
};
