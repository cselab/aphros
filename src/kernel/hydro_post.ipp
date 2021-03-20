// Created by Petr Karnakov on 20.03.2021
// Copyright 2021 ETH Zurich

#include "hydro_post.h"

template <class M>
struct HydroPost<M>::Imp {
  template <class MEB>
  static FieldCell<Scal> GetVortScal(
      const FieldCell<Vect>& fcvel, const MapEmbed<BCond<Vect>>& me_vel,
      MEB& eb) {
    constexpr auto dim = M::dim;
    auto& m = eb.GetMesh();
    using UEB = UEmbed<M>;

    std::array<FieldCell<Vect>, dim> grad;
    for (size_t d = 0; d < dim; ++d) {
      grad[d].Reinit(m, Vect(0));
      const auto mebc = GetScalarCond(me_vel, d, m);
      const FieldCell<Scal> fcu = GetComponent(fcvel, d);
      const FieldFace<Scal> ffg = UEB::Gradient(fcu, mebc, m);
      grad[d] = UEB::AverageGradient(ffg, m);
    }

    FieldCell<Scal> res(m, 0);
    for (auto c : m.Cells()) {
      res[c] = grad[1][c][0] - grad[0][c][1];
    }
    return res;
  }

  static FieldCell<Scal> GetField(
      const Hydro<M>* hydro, std::string name, const M& m) {
    if (name == "p" || name == "pressure") {
      return hydro->fs_->GetPressure();
    }
    if (name == "omz" || name == "vorticity") {
      if (hydro->eb_) {
        return GetVortScal(
            hydro->fs_->GetVelocity(), hydro->fs_->GetVelocityCond(),
            *hydro->eb_);
      }
      return GetVortScal(
          hydro->fs_->GetVelocity(), hydro->fs_->GetVelocityCond(), m);
    }
    if (name == "ebvf" || name == "embed fraction") {
      FieldCell<Scal> fc(m, 1);
      if (hydro->eb_) {
        auto& eb = *hydro->eb_;
        for (auto c : m.SuCells()) {
          fc[c] = eb.GetVolumeFraction(c);
        }
      }
      return fc;
    }
    if (name == "vf" || name == "volume fraction") {
      return hydro->as_->GetField();
    }
    if (name == "vmagn" || name == "velocity magnitude") {
      FieldCell<Scal> fc(m, 0);
      for (auto c : m.SuCells()) {
        fc[c] = hydro->fs_->GetVelocity()[c].norm();
      }
      return fc;
    }
    if (name == "vx" || name == "velocity x") {
      FieldCell<Scal> fc(m, 0);
      for (auto c : m.SuCells()) {
        fc[c] = hydro->fs_->GetVelocity()[c][0];
      }
      return fc;
    }
    if (name == "vy" || name == "velocity y") {
      FieldCell<Scal> fc(m, 0);
      for (auto c : m.SuCells()) {
        fc[c] = hydro->fs_->GetVelocity()[c][1];
      }
      return fc;
    }
    if (m.IsRoot()) {
      std::cerr << "Unknown field '" + name + "'\n";
    }
    return FieldCell<Scal>(m, 0);
  };
};

template <class M>
auto HydroPost<M>::GetField(const Hydro<M>* hydro, std::string name, const M& m)
    -> FieldCell<Scal> {
  return Imp::GetField(hydro, name, m);
}
