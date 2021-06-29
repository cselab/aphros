// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "advection.h"
#include "geom/mesh.h"
#include "partstrmeshm.h"

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

  Particles(const typename PartStrMeshM<M>::Par& par);
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

  Heights() = default;
  // Computes curvature from volume fractions
  // and a piecewise linear reconstruction
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const M& eb) override;
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const Embed<M>& eb) override;
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

  Hybrid(const typename PartStrMeshM<M>::Par& par);
  // Computes curvature from volume fractions
  // and a piecewise linear reconstruction
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const M& eb) override;
  void CalcCurvature(
      const Multi<FieldCell<Scal>*>& fck, const Plic& plic, M& m,
      const Embed<M>& eb) override;
};

} // namespace curvature

template <class M_>
struct UCurv {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Plic = generic::Plic<Vect>;

  // Computes curvature with height functions.
  // fcu: volume fraction [a]
  // fcn: normal
  // edim: effective dimension (2 or 3)
  // Output:
  // fck: curvature
  static void CalcCurvHeight(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn, size_t edim,
      FieldCell<Scal>& fck, M& m);

 private:
  struct Imp;
};
