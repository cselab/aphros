// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "advection.h"
#include "geom/mesh.h"
#include "parse/vars.h"
#include "partstrmeshm.h"
#include "util/format.h"

namespace curvature {

template <class M_>
class Estimator {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Plic = generic::Plic<Vect>;

  Estimator() = default;
  Estimator(const Estimator&) = delete;
  virtual ~Estimator() {}
  // Computes curvature from volume fractions
  // and a piecewise linear reconstruction
  virtual void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const M& eb) = 0;
  virtual void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const Embed<M>& eb) = 0;
};

// Curvature estimator using particles fitted to piecewise linear interface
template <class M_>
class Particles : public Estimator<M_> {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Plic = generic::Plic<Vect>;

  Particles(
      M& m, const typename PartStrMeshM<M>::Par& par,
      const GRange<size_t>& layers);
  Particles(const Particles&) = delete;
  ~Particles();
  // Computes curvature from volume fractions
  // and a piecewise linear reconstruction
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const M& eb) override;
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const Embed<M>& eb) override;
  std::unique_ptr<PartStrMeshM<M>> ReleaseParticles();
  const PartStrMeshM<M>* GetParticles() const;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

// Curvature estimator using height functions
template <class M_>
class Heights : public Estimator<M_> {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Plic = generic::Plic<Vect>;

  Heights();
  ~Heights();
  // Computes curvature from volume fractions
  // and a piecewise linear reconstruction
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const M& eb) override;
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const Embed<M>& eb) override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

// Hybrid curvature estimator using height functions where defined
// and particles otherwise.
template <class M_>
class Hybrid : public Estimator<M_> {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Plic = generic::Plic<Vect>;

  Hybrid(
      M& m, const typename PartStrMeshM<M>::Par& par,
      const GRange<size_t>& layers);
  Hybrid(const Hybrid&) = delete;
  ~Hybrid();
  // Computes curvature from volume fractions
  // and a piecewise linear reconstruction
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const M& eb) override;
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const Embed<M>& eb) override;
  std::unique_ptr<PartStrMeshM<M>> ReleaseParticles();
  const PartStrMeshM<M>* GetParticles() const;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

template <class M>
std::unique_ptr<Estimator<M>> MakeEstimator(
    const Vars& var, M& m, const GRange<size_t>& layers);

} // namespace curvature
