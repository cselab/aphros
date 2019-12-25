#pragma once

#include <memory>

#include "advection.h"
#include "dump/dumper.h"
#include "partstrmeshm.h"

namespace solver {

template <class M_>
class Vof final : public AdvectionSolver<M_> {
 public:
  using M = M_;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;

  struct Par {
    size_t dim = 3; // dimension (dim=2 assumes zero velocity in z)
    bool part = false; // particles
    Scal part_relax = 0.5;
    Scal part_h = 4.; // particle string equilibrium length
    int part_verb = 0; // debug output
    bool part_n = false; // normal from particles
    // curvature from particles
    // if true, GetCurv returns fckp_
    bool part_k = false; // curvature from particles
    size_t part_dump_fr = 100; // num frames dump
    size_t part_report_fr = 100; // num frames report
    Scal clipth = 1e-6; // vf clipping threshold
    std::unique_ptr<Dumper> dmp; // dumper for particles
    bool dumppoly = false; // dump reconstructed interface (cut polygons)
    bool dumppolymarch = false; // dump reconstructed interface (marching cube)
    bool dumppart = false; // dump particles
    bool dumppartinter = false; // dump interface for particles
    bool bcc_reflectpoly = true; // reflection for DumpPolyMarch
    Scal dumppolymarch_fill = -1; // fill cells outside

    Scal part_segcirc = 1.; // factor for circular segment
    size_t part_np = 11; // number of particles per string
    size_t part_ns = 2; // number of strings per cell
    size_t part_itermax = 100; // particles itermax
    Scal part_tol = 0.01; // tolerance
    bool part_dn = false;
    Scal part_maxr = 0;
    int verb = 0;
    bool recolor_unionfind = true;
    bool recolor_reduce = true;
    bool recolor_grid = true;
    bool vtkbin = false;
    bool vtkmerge = true;
    Scal vtkiso = 0.5;
    bool sharpen = false;
    Scal sharpen_cfl = 0.5;
    enum class Scheme { plain, aulisa, weymouth };
    Scheme scheme;
    Scal avgnorm0 = 1; // original normal with sum(u)<avgnorm0
    Scal avgnorm1 = 1; // overriden normal with sum(u)>=acgnorm1
    Scal clfixed = -1; // if >= 0, value for color at point clfixed_x
    Vect clfixed_x = Vect(1e10);
    bool cloverride = false; // XXX adhoc if clear1<1, override color with 0
    size_t layers = 4;
    Scal coalth = 1.5;
  };

  // Constructor
  Vof(M& m, const FieldCell<Scal>& fcu, const FieldCell<Scal>& fccl,
      const MapCondFaceAdvection<Scal>& mfc, const FieldFace<Scal>* ffv,
      const FieldCell<Scal>* fcs, double t, double dt,
      std::shared_ptr<Par> par);
  ~Vof();
  // Parameters
  Par* GetPar();
  // ...
  void StartStep() override;
  // ...
  void MakeIteration() override;
  // ...
  void FinishStep() override;
  // ...
  void PostStep() override;
  // Volume fraction
  const FieldCell<Scal>& GetField(Layers l) const override;
  // ...
  using P::GetField;
  // Plane constant
  const FieldCell<Scal>& GetAlpha() const;
  // Normal to interface
  const FieldCell<Vect>& GetNormal() const;
  // Default curvature
  const FieldCell<Scal>& GetCurv() const override;
  // Curvature from height functions
  const FieldCell<Scal>& GetCurvH() const;
  // Curvature from particles
  const FieldCell<Scal>& GetCurvP() const;
  // Color
  const FieldCell<Scal>& GetColor() const;
  static constexpr Scal kClNone = -1; // no color
  // Image vector, number of passes through periodic boundaries
  MIdx GetImage(IdxCell c) const;

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};

} // namespace solver
