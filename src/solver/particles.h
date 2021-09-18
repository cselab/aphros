// Created by Petr Karnakov on 26.06.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <functional>
#include <memory>
#include <unordered_set>

#include "cond.h"
#include "geom/mesh.h"
#include "multi.h"
#include "solver/embed.h"

enum class ParticlesMode {
  tracer, // particle velocity equals mixture velocity
  stokes, // linear relaxation with mixture velocity,
          // relaxation time from Stokes' law
          // based on gravity, mixture density, mixture viscosity in Conf,
          // and particle density `rho`, particle radius `r`
  termvel, // linear relaxation with mixture velocity,
           // relaxation time based on given terminal velocity `termvel`
};

namespace generic {

template <class Vect_>
struct ParticlesView {
  using Scal = typename Vect_::Scal;
  using Vect = Vect_;
  std::vector<Vect>& x; // positions
  std::vector<bool>& inner; // true if particle is owned by block
  std::vector<Vect>& v; // velocity
  std::vector<Scal>& id; // identifier
  std::vector<Scal>& r; // radius
  std::vector<Scal>& source; // volume source
  std::vector<Scal>& rho; // density
  std::vector<Scal>& termvel; // terminal settling velocity
  std::vector<Scal>& removed; // 1 if particle is removed, else 0
                              // TODO: use type bool, needs support by Comm()
  std::vector<std::vector<Scal>*> GetAttrScal() const {
    return {&id, &r, &source, &rho, &termvel, &removed};
  }
  std::vector<std::vector<Vect>*> GetAttrVect() const {
    return {&v};
  }
};

} // namespace generic

template <class M>
typename M::CommPartRequest GetCommPartRequest(
    const generic::ParticlesView<typename M::Vect>& view) {
  typename M::CommPartRequest res;
  res.x = &view.x;
  res.inner = &view.inner;
  res.attr_scal = view.GetAttrScal();
  res.attr_vect = view.GetAttrVect();
  return res;
}

template <class M_>
class ParticlesInterface {
 public:
  using M = M_;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using ParticlesView = generic::ParticlesView<Vect>;

  struct Conf {
    Scal mixture_density{0};
    Scal mixture_viscosity{0};
    Vect gravity{0};
    ParticlesMode mode{ParticlesMode::tracer};
  };

  virtual ~ParticlesInterface() {}
  virtual const Conf& GetConf() const = 0;
  virtual void SetConf(Conf) = 0;
  // Makes one time step.
  // dt: time step
  // fe_flux: mixture volume flux
  // velocity_hook: function called after velocity is projected to particles
  virtual void Step(
      Scal dt, const FieldEmbed<Scal>& fe_flux,
      const MapEmbed<BCond<Vect>>& mebc_velocity,
      std::function<void(const ParticlesView&)> velocity_hook) = 0;
  // Appends current state with particles from `view`.
  virtual void Append(const ParticlesView& view) = 0;
  // Updates `id` with a unique contigous index for particles appended since
  // last call of SetUniqueIdForAppended().
  virtual void SetUniqueIdForAppended(M& m) = 0;
  // Returns view with pointers to fields.
  virtual ParticlesView GetView() const = 0;
  virtual Scal GetTime() const = 0;
  // Returns the number of particles received at the last communication.
  virtual size_t GetNumRecv() const = 0;
  // Dumps current particles to CSV file.
  // path: path to output file
  // sel: set of fields to dump
  virtual void DumpCsv(
      const std::string& path,
      const std::unordered_set<std::string>& sel) const = 0;
};

template <class EB_>
class Particles : public ParticlesInterface<typename EB_::M> {
 public:
  using M = typename EB_::M;
  static constexpr size_t dim = M::dim;
  using Base = ParticlesInterface<M>;
  using EB = EB_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using ParticlesView = generic::ParticlesView<Vect>;
  using Conf = typename Base::Conf;
  using UEB = UEmbed<M>;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;

  // Constructor
  // vfcu: initial volume fraction
  Particles(M& m, const EB& eb, const ParticlesView& init, Scal t, Conf conf);
  ~Particles();
  const Conf& GetConf() const override;
  void SetConf(Conf) override;
  void Step(
      Scal dt, const FieldEmbed<Scal>& fe_flux,
      const MapEmbed<BCond<Vect>>& mebc_velocity,
      std::function<void(const ParticlesView&)> velocity_hook) override;
  void Append(const ParticlesView&) override;
  void SetUniqueIdForAppended(M& m) override;
  ParticlesView GetView() const override;
  Scal GetTime() const override;
  size_t GetNumRecv() const override;
  void DumpCsv(
      const std::string& path,
      const std::unordered_set<std::string>& sel) const override;
  struct ReadCsvStatus {
    // Columns read from the CSV file.
    bool x = false;
    bool inner = false;
    bool v = false;
    bool id = false;
    bool r = false;
    bool source = false;
    bool rho = false;
    bool termvel = false;
    bool removed = false;
  };
  // Reads recognized columns (position, velocity, radius, ...) from CSV file.
  static ReadCsvStatus ReadCsv(
      std::istream& fin, const ParticlesView& view, char deliimiter = ',');
  static ReadCsvStatus ReadCsv(
      const std::string& path, const ParticlesView& view,
      char deliimiter = ',');

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
