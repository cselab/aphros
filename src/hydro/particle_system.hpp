/*
 * particle_system.hpp
 *
 *  Created on: Feb 11, 2016
 *      Author: Petr Karnakov (petr.karnakov@gmail.com)
 */

#pragma once

#include <vector>
#include <array>
#include "mesh.hpp"
#include "solver.hpp"

namespace solver {

template<class Vect, class Value>
Value GetBilinearInterpolation(Vect x, Vect xlb, Vect xrt,
                               Value alb, Value arb, Value art, Value alt)
{
  Vect lower = x - xlb, upper = xrt - x, size = xrt - xlb;
  Value blb = alb * upper[0] * upper[1];
  Value brb = arb * lower[0] * upper[1];
  Value brt = art * lower[0] * lower[1];
  Value blt = alt * upper[0] * lower[1];
  return (blb + brb + brt + blt) / (size[0] * size[1]);
}

template <class Mesh>
class ParticleSystem : public UnsteadySolver {
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using IdxCell = geom::IdxCell;
  using IdxNode = geom::IdxNode;
  using IdxParticle = geom::IdxGeneric<20160211>;
  template <class T>
  using FieldCell = geom::FieldCell<T>;
  template <class T>
  using FieldNode = geom::FieldNode<T>;
  template <class T>
  using FieldParticle = geom::FieldGeneric<T, IdxParticle>;

  const Mesh& mesh;

  geom::SearchMesh<Mesh> search_mesh_;

  const Scal kSpawningGapRelative;
  const Scal kParticleRadiusFactor;
  const size_t kMinNumParticlesInCell;
  const size_t kMaxNumParticlesInCell;
  const Scal back_relaxation_factor_;

  std::vector<FieldCell<Scal>> fc_fields_;
  std::vector<FieldParticle<Scal>> fp_fields_;

  geom::Range<IdxParticle> particles_;

  FieldParticle<Vect> fp_position_;
  FieldParticle<IdxCell> fp_cell_;
  FieldParticle<Scal> fp_cell_center_distance_;
  FieldParticle<bool> fp_is_removed_;

  FieldNode<Vect>* p_fn_velocity_;
  std::vector<const geom::FieldCell<Scal>*> v_p_fc_source_;
  FieldCell<IdxParticle> fc_nearest_particle_;
  FieldCell<Scal> fc_nearest_particle_distance_;
  FieldCell<size_t> fc_num_particles_;
  FieldCell<Scal> fc_particle_weight_sum_;
  FieldCell<Scal> fc_base_radius_;

  // Temporary buffers:
  FieldParticle<Vect> fp_velocity_;

  bool is_consistent_;
  bool is_fields_relevant_;

  void MarkInconsistent() {
    is_consistent_ = false;
  }
  void MarkConsistent() {
    is_consistent_ = true;
  }
  bool IsConsistent() const {
    return is_consistent_;
  }
  void MarkFieldsIrrelevant() {
    is_fields_relevant_ = false;
  }
  void MarkFieldsRelevant() {
    is_fields_relevant_ = true;
  }
  bool IsFieldsRelevant() const {
    return is_fields_relevant_;
  }
  bool MayAddParticle(IdxCell idxcell) const {
    return fc_num_particles_[idxcell] == 0 ||
        (fc_num_particles_[idxcell] < kMinNumParticlesInCell &&
        fc_nearest_particle_distance_[idxcell] >=
        kSpawningGapRelative * fc_base_radius_[idxcell]);
  }
  bool MayRemoveParticle(IdxCell idxcell) const {
    return fc_num_particles_[idxcell] > kMaxNumParticlesInCell;
  }
  IdxParticle AddParticle(Vect position, IdxCell idxcell) {
    IdxParticle idxpart(fp_position_.size());
    fp_position_.push_back(position);
    particles_ = geom::Range<IdxParticle>(0, fp_position_.size());
    fp_velocity_.push_back(Vect());
    fp_cell_.push_back(idxcell);
    fp_cell_center_distance_.push_back(position.dist(mesh.GetCenter(idxcell)));
    fp_is_removed_.push_back(false);
    for (auto& fp_f : fp_fields_) {
      fp_f.push_back(Scal());
    }
    MarkInconsistent();
    MarkFieldsIrrelevant();
    return idxpart;
  }
  void MarkRemoved(IdxParticle idx) {
    fp_is_removed_[idx] = true;
    MarkInconsistent();
    MarkFieldsIrrelevant();
  }
  size_t GetNumLiveParticles() const {
    size_t res = 0;
    for (auto idx : particles_) {
      if (!fp_is_removed_[idx]) {
        ++res;
      }
    }
    return res;
  }
  void RestoreConsistency() {
    if (fp_position_.empty()) {
      return;
    }
    if (IsConsistent()) {
      return;
    }

    const Scal infinity = 1e16;

    // Find a cell for each particle
    for (auto idxpart : particles_) {
      if (!fp_is_removed_[idxpart]) {
        const Vect particle_position = fp_position_[idxpart];
        auto oldcell = fp_cell_[idxpart];
        if (!mesh.IsInside(oldcell, particle_position)) {
          IdxCell newcell = search_mesh_.FindCell(particle_position);
          if (newcell == IdxCell::None()) {
            fp_is_removed_[idxpart] = true;
          } else {
            fp_cell_[idxpart] = newcell;
            fp_cell_center_distance_[idxpart] =
                fp_position_[idxpart].dist(mesh.GetCenter(newcell));
          }
        }
      }
    }

    // Find the nearest particle for each cell
    for (auto idxcell : mesh.Cells()) {
      fc_num_particles_[idxcell] = 0;
      fc_nearest_particle_[idxcell] = IdxParticle::None();
      fc_nearest_particle_distance_[idxcell] = infinity;
    }
    for (auto idxpart : particles_) {
      if (!fp_is_removed_[idxpart]) {
        const auto idxcell = fp_cell_[idxpart];
        const Scal distance = fp_cell_center_distance_[idxpart];
        ++fc_num_particles_[idxcell];
        if (fc_nearest_particle_distance_[idxcell] < distance) {
          fc_nearest_particle_[idxcell] = idxpart;
          fc_nearest_particle_distance_[idxcell] = distance;
        }
      }
    }
  }
  void AddParticles() {
    RestoreConsistency();
    for (auto idxcell : mesh.Cells()) {
      if (MayAddParticle(idxcell)) {
        IdxParticle idxpart = AddParticle(mesh.GetCenter(idxcell), idxcell);
        for (size_t n = 0; n < fp_fields_.size(); ++n) {
          fp_fields_[n][idxpart] = fc_fields_[n][idxcell];
        }
      }
    }
    MarkFieldsIrrelevant();
    MarkInconsistent();
  }
  void RemoveParticles() {
    RestoreConsistency();
    for (auto it = particles_.end(); it != particles_.begin(); ) {
      --it;
      auto idxpart = *it;
      if (MayRemoveParticle(fp_cell_[idxpart])) {
        MarkRemoved(idxpart);
      }
    }
    MarkFieldsIrrelevant();
    MarkInconsistent();
  }
  // TODO: Think of which trajectories the particles should travel
  // in a cell to ensure volume conservation
  // for a given face volumetric fluxes
  void CalcParticleVelocity() {
    RestoreConsistency();
    for (auto idxpart : particles_) {
      IdxCell idxcell = fp_cell_[idxpart];
      std::array<IdxNode, 4> nodes;
      for (size_t i = 0; i < mesh.GetNumNeighbourNodes(idxcell); ++i) {
        nodes[i] = mesh.GetNeighbourNode(idxcell, i);
      }
      fp_velocity_[idxpart] = GetBilinearInterpolation(
          fp_position_[idxpart],
          mesh.GetNode(nodes[0]), mesh.GetNode(nodes[2]),
          (*p_fn_velocity_)[nodes[0]],
          (*p_fn_velocity_)[nodes[1]],
          (*p_fn_velocity_)[nodes[2]],
          (*p_fn_velocity_)[nodes[3]]);
    }
  }
  void AdvanceParticles(double time_step_factor = 1.,
                        bool apply_source = true) {
    const Scal time_step = this->GetTimeStep() * time_step_factor;
    for (auto idxpart : particles_) {
      fp_position_[idxpart] += fp_velocity_[idxpart] * time_step;
    }
    if (apply_source) {
      for (auto idxpart : particles_) {
        IdxCell idxcell = fp_cell_[idxpart];
        for (size_t n = 0; n < fc_fields_.size(); ++n) {
          fp_fields_[n][idxpart] += (*v_p_fc_source_[n])[idxcell] * time_step;
          //fp_fields_[n][idxpart] =
          //    std::min(1., std::max(0., fp_fields_[n][idxpart]));
        }
      }
    }
    MarkInconsistent();
    MarkFieldsIrrelevant();
  }
  void MakeBackRelaxation() {
    for (auto idxpart : particles_) {
      IdxCell idxcell = fp_cell_[idxpart];
      for (size_t n = 0; n < fc_fields_.size(); ++n) {
        fp_fields_[n][idxpart] =
            back_relaxation_factor_ * fp_fields_[n][idxpart] +
            (1. - back_relaxation_factor_) * fc_fields_[n][idxcell];
      }
    }
  }
  template <class T>
  void EraseRemovedParticlesField(FieldParticle<T>& fp_field,
                             geom::Range<IdxParticle> particles_new) {
    FieldParticle<T> res(particles_new);
    size_t it = 0;
    for (auto idxpart : particles_) {
      if (!fp_is_removed_[idxpart]) {
        res[IdxParticle(it++)] = fp_field[idxpart];
      }
    }
    assert(it == res.size());
    fp_field = res;
  }
  void EraseRemovedParticles() {
    geom::Range<IdxParticle> particles_new(0, GetNumLiveParticles());
    EraseRemovedParticlesField(fp_position_, particles_new);
    EraseRemovedParticlesField(fp_cell_, particles_new);
    EraseRemovedParticlesField(fp_cell_center_distance_, particles_new);
    for (size_t n = 0; n < fp_fields_.size(); ++n) {
      EraseRemovedParticlesField(fp_fields_[n], particles_new);
    }
    EraseRemovedParticlesField(fp_is_removed_, particles_new);
    particles_ = particles_new;
    fp_velocity_.resize(particles_.size());
    RestoreConsistency();
  }
  Scal CalcBaseRadius(IdxCell idxcell) const {
    Scal res = 1e16;
    for (size_t i = 0; i < mesh.GetNumNeighbourFaces(idxcell); ++i) {
      res = std::min(
          res,
          mesh.GetCenter(mesh.GetNeighbourFace(idxcell, i)).dist(
              mesh.GetCenter(idxcell)));
    }
    return res;
  }
  void CalcFields() {
    RestoreConsistency();

    for (auto idxcell : mesh.Cells()) {
      fc_particle_weight_sum_[idxcell] = 0;
    }
    for (size_t n = 0; n < fc_fields_.size(); ++n) {
      for (auto idxcell : mesh.Cells()) {
        fc_fields_[n][idxcell] = 0.;
      }
    }

    const Scal particle_radius =
        fc_base_radius_[*mesh.Cells().begin()] * kParticleRadiusFactor;

    auto weight = [particle_radius](Vect target, Vect particle_position) {
      const Scal relative_distance =
          particle_position.dist(target) / particle_radius;
      return std::max(0., 1. - relative_distance);
    };

    for (auto idxpart : particles_) {
      if (!fp_is_removed_[idxpart]) {
        auto nearby_cells = search_mesh_.GetNearbyCells(
            fp_position_[idxpart], particle_radius);
        for (auto idxcell : nearby_cells) {
          fc_particle_weight_sum_[idxcell] +=
              weight(mesh.GetCenter(idxcell), fp_position_[idxpart]);

          for (size_t n = 0; n < fc_fields_.size(); ++n) {
            fc_fields_[n][idxcell] += fp_fields_[n][idxpart] *
                weight(mesh.GetCenter(idxcell), fp_position_[idxpart]);
          }
        }

      }
    }

/*
    // Nearest neighbour rule:
    for (auto idxpart : particles_) {
      fc_particle_weight_sum_[fp_cell_[idxpart]] += 1.0;
    }
    for (size_t n = 0; n < fc_fields_.size(); ++n) {
      for (auto idxpart : particles_) {
        fc_fields_[n][fp_cell_[idxpart]] += fp_fields_[n][idxpart];
      }
    }
*/

    for (size_t n = 0; n < fc_fields_.size(); ++n) {
      for (auto idxcell : mesh.Cells()) {
        fc_fields_[n][idxcell] /= fc_particle_weight_sum_[idxcell];
      }
    }

    MarkFieldsRelevant();
  }

  // TODO sort particles by fp_cell_.GetRaw() (cache usage optimization)
  // TODO eliminate unnecessary RestoreConsistency() calls
  //      and measure earned performance gains

 public:
  ParticleSystem(const Mesh& mesh,
                 const std::vector<FieldCell<Scal>>& fc_initial,
                 const std::vector<const geom::FieldCell<Scal>*>& v_p_fc_source,
                 FieldNode<Vect>* p_fn_velocity,
                 double time, double time_step,
                 Scal spawning_gap_relative = 0.2,
                 Scal particle_radius_factor = 2.,
                 size_t min_num_particle_in_cell = 3,
                 size_t max_num_particle_in_cell = 5,
                 Scal back_relaxation_factor = 1.)
      : UnsteadySolver(time, time_step)
      , mesh(mesh)
      , search_mesh_(mesh, 4, mesh.GetNumCells() * 16, 1)
      , kSpawningGapRelative(spawning_gap_relative)
      , kParticleRadiusFactor(particle_radius_factor)
      , kMinNumParticlesInCell(min_num_particle_in_cell)
      , kMaxNumParticlesInCell(max_num_particle_in_cell)
      , back_relaxation_factor_(back_relaxation_factor)
      , fc_fields_(fc_initial)
      , fp_fields_(fc_fields_.size())
      , particles_(0, 0)
      , p_fn_velocity_(p_fn_velocity)
      , v_p_fc_source_(v_p_fc_source)
      , fc_nearest_particle_(mesh)
      , fc_nearest_particle_distance_(mesh)
      , fc_num_particles_(mesh)
      , fc_particle_weight_sum_(mesh)
      , fc_base_radius_(mesh)
      , is_consistent_(false)
      , is_fields_relevant_(false)
  {
    for (auto idxcell : mesh.Cells()) {
      fc_base_radius_[idxcell] = CalcBaseRadius(idxcell);
    }
    AddParticles();
    RestoreConsistency();
  }
  void CalcStep() override {
    auto fp_tmp = fp_position_;
    CalcParticleVelocity();
    AdvanceParticles(0.5, false);
    CalcParticleVelocity();
    fp_position_ = fp_tmp;
    AdvanceParticles(1.);
    RemoveParticles();
    AddParticles();
    EraseRemovedParticles();
    if (!IsFieldsRelevant()) {
      CalcFields();
    }
    MakeBackRelaxation();
    std::cout << "\nnum_live_particles = " << GetNumLiveParticles();
    std::cout << "\nnum_all_particles = " << fp_position_.size();
    std::cout << std::endl;
  }
  const FieldCell<Scal>& GetField(size_t n) {
    return fc_fields_[n];
  }
};

} // namespace solver
