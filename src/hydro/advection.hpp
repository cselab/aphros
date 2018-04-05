#pragma once

#include <exception>
#include <fstream>

#include "mesh.hpp"
#include "linear.hpp"
#include "solver.hpp"

namespace solver {

template <class M>
geom::FieldCell<typename M::Vect> GetDeformingVelocity(const M& m) {
  using Vect = typename M::Vect;
  geom::FieldCell<Vect> r(m, 0);
  for (auto c : m.Cells()) {
    auto x = m.GetCenter(c);
    r[c][0] = -std::cos(x[1]) * std::sin(x[0]);
    r[c][1] = std::cos(x[0]) * std::sin(x[1]);
  }
  return r;
}

template <class Mesh, class Vel>
class AdvectionSolver : public UnsteadyIterativeSolver {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;

 protected:
  const Vel* p_f_velocity_;
  const geom::FieldCell<Scal>* p_fc_source_;

 public:
  AdvectionSolver(double time, double time_step,
                  const Vel* p_f_velocity,
                  const geom::FieldCell<Scal>* p_fc_source)
      : UnsteadyIterativeSolver(time, time_step)
      , p_f_velocity_(p_f_velocity)
      , p_fc_source_(p_fc_source)
  {}
  virtual void AssignVelocity(const Vel* p_f_velocity) {
    p_f_velocity_ = p_f_velocity;
  }
  virtual const geom::FieldCell<Scal>& GetField() = 0;
};

template <class Mesh, class Vel>
class AdvectionSolverExplicit :
    public AdvectionSolver<Mesh, Vel> {
  Mesh& mesh;
  static constexpr size_t dim = Mesh::dim;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  LayersData<geom::FieldCell<Scal>> fc_u_;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  geom::MapFace<std::shared_ptr<ConditionFace>> mf_u_cond_;
  // Common buffers:
  geom::FieldFace<Vect> ff_velocity_;
  geom::FieldFace<Scal> ff_volume_flux_;
  geom::FieldFace<Scal> ff_flux_;
  geom::FieldFace<Scal> ff_u_;
  geom::FieldFace<Scal> sharp_af_;
  geom::FieldCell<Vect> sharp_gc_;
  geom::FieldFace<Vect> sharp_gf_;
  Scal sharp_;
  Scal sharpo_;
  Scal sharp_max_;

 public:
  AdvectionSolverExplicit(
      Mesh& mesh,
      const geom::FieldCell<Scal>& fc_u_initial,
      const geom::MapFace<std::shared_ptr<ConditionFace>>& mf_u_cond_,
      const Vel* p_fn_velocity,
      const geom::FieldCell<Scal>* p_fc_source,
      double time, double time_step,
      Scal sharp, Scal sharpo, Scal sharp_max)
      : AdvectionSolver<Mesh, Vel>(
          time, time_step, p_fn_velocity, p_fc_source)
      , mesh(mesh)
      , mf_u_cond_(mf_u_cond_)
      , ff_volume_flux_(mesh)
      , ff_flux_(mesh)
      , sharp_(sharp)
      , sharpo_(sharpo)
      , sharp_max_(sharp_max)
  {
    fc_u_.time_curr = fc_u_initial;
  }
  void StartStep() override {
    this->ClearIter();
    fc_u_.time_prev = fc_u_.time_curr;
    fc_u_.iter_curr = fc_u_.time_prev;
  }
  // Correct the inconsistency with arguments: velocity vs volume flux
  geom::FieldFace<Scal> ConvertVolumeFlux(
      const geom::FieldFace<Scal>* p_f_velocity) {
    return *p_f_velocity;
  }
  geom::FieldFace<Scal> ConvertVolumeFlux(
      const geom::FieldNode<Vect>* p_f_velocity)  {
    auto ff_velocity = Interpolate(*p_f_velocity, mesh);

    geom::FieldFace<Scal> ff_volume_flux(mesh);
    for (IdxFace idxface : mesh.AllFaces()) {
      ff_volume_flux[idxface] =
          ff_velocity[idxface].dot(mesh.GetSurface(idxface));
    }
    return ff_volume_flux;
  }
  void MakeIteration() override {
    auto sem = mesh.GetSem();
    if (sem()) {
      auto& prev = fc_u_.iter_prev;
      auto& curr = fc_u_.iter_curr;
      prev = curr;

      ff_volume_flux_ = ConvertVolumeFlux(this->p_f_velocity_);

      ff_u_ = InterpolateSuperbee(
          prev,
          Gradient(Interpolate(prev, mf_u_cond_, mesh), mesh),
          mf_u_cond_, ff_volume_flux_, mesh);

      for (auto idxface : mesh.Faces()) {
        ff_flux_[idxface] = ff_u_[idxface] * ff_volume_flux_[idxface];
      }

      // Interface sharpening
      // append to ff_flux_
      if (std::abs(sharp_) != 0.) {
        // zero-derivative bc for Vect
        geom::MapFace<std::shared_ptr<ConditionFace>> mfvz;
        for (auto it : mf_u_cond_) {
          IdxFace i = it.GetIdx();
          mfvz[i] = std::make_shared<
              ConditionFaceDerivativeFixed<Vect>>(
                  Vect(0), it.GetValue()->GetNci());
        }
        auto& af = sharp_af_;
        auto& gc = sharp_gc_;
        auto& gf = sharp_gf_;
        af = Interpolate(curr, mf_u_cond_, mesh);
        gc = Gradient(af, mesh);
        gf = Interpolate(gc, mfvz, mesh);
        for (auto i : mesh.Faces()) {
          // normal to interface
          const Vect ni = gf[i] / (gf[i].norm() + 1e-6); 
          const Vect n = mesh.GetNormal(i);
          const Vect s = mesh.GetSurface(i);
          const Scal a = mesh.GetArea(i);
          //const Scal nf = n.dot(mesh.GetNormal(i));
          //const Scal uf = ff_volume_flux_[i];
          //const Scal uf = 1.;
          //const Scal am = sharp_max_;
          //const Scal eh = sharp_ * mesh.GetArea(i);
          //ff_flux_[i] -= sharpo_ * std::abs(uf * nf) * nf * 
          //    (eh * gf[i].norm() - af[i] * (1. - af[i] / am));
          ff_flux_[i] += 
            sharp_ * af[i] * (1. - af[i]) * ni.dot(s)
            - sharpo_ * gf[i].dot(s);
        }
      }

      for (auto idxcell : mesh.Cells()) {
        Scal flux_sum = 0.;
        for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
          IdxFace idxface = mesh.GetNeighbourFace(idxcell, i);
          flux_sum += ff_flux_[idxface] * mesh.GetOutwardFactor(idxcell, i);
        }

        curr[idxcell] = fc_u_.time_prev[idxcell] -
            this->GetTimeStep() / mesh.GetVolume(idxcell) * flux_sum +
            this->GetTimeStep() * (*this->p_fc_source_)[idxcell];
      }
      mesh.Comm(&curr);

      this->IncIter();
    }
  }
  void FinishStep() override {
    fc_u_.time_curr = fc_u_.iter_curr;
    this->IncTime();
  }
  double GetError() const override {
    if (this->GetError() == 0) {
      return 1.;
    }
    return CalcDiff<geom::FieldCell<Scal>, Mesh>(
        fc_u_.iter_curr, fc_u_.iter_prev, mesh);
  }
  const geom::FieldCell<Scal>& GetField() override {
    return fc_u_.time_curr;
  }
};

} // namespace solver
