#pragma once

#include <memory>

#include "advection.h"
#include "dump/dumper.h"
#include "partstrmeshm.h"

namespace solver {

template <class M_>
class Vof : public AdvectionSolver<M_> {
 public:
  using M = M_;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  struct Par {
    size_t dim = 3; // dimension (dim=2 assumes zero velocity in z)
    bool curvgrad = false; // compute curvature using gradient
    bool part = false; // particles
    Scal part_relax = 0.5; 
    Scal part_h = 4.;  // particle string equilibrium length
    int part_verb = 0; // debug output
    Scal part_kstr = 1.; // stretching
    Scal part_kattr = 1.; // attraction to reconstructed interface
    Scal part_kbend = 1.; // bending
    bool part_bendmean = true; // bending to mean angle (fit circle)
    bool part_n = false; // normal from particles
    // curvature from particles
    // if true, GetCurv returns fckp_
    bool part_k = false; // curvature from particles
    size_t part_dump_fr = 100; // num frames dump
    size_t part_report_fr = 100; // num frames report
    Scal part_intth = 1e-5; // interface threshold for particle seed
    Scal poly_intth = 0.; // interface threshold for DumpPoly
    Scal clipth = 1e-6; // vf clipping threshold
    std::unique_ptr<Dumper> dmp; // dumper for particles
    bool dumppoly = false; // dump reconstructed interface (cut polygons)
    bool dumppolymarch = false; // dump reconstructed interface (marching cube)
    bool dumppart = false; // dump particles
    bool dumppartinter = false; // dump interface for particles
    bool bcc_reflect = false; // reflection at boundaries
    Scal bcc_fill = 0;        // fill value for halo cells

    int part_constr = 0; // 0: no constraints
                         // 1: fixed distance, constant angle
                         // 2: fixed distance, linear angle
    Scal part_segcirc = 1.; // factor for circular segment
    size_t part_np = 11; // number of particles per string
    size_t part_ns = 2; // number of strings per cell
    size_t part_itermax = 100; // particles itermax
    Scal part_tol = 0.01; // tolerance
    Scal part_tmax = 180.; 
    Scal part_dtmax = 10.; 
    Scal part_anglim = 90.; 
    bool part_dn = false;
    Scal part_maxr = 0;
    using AF = typename solver::PartStrMeshM<M>::AF;
    using AR = typename solver::PartStrMeshM<M>::AR;
    AF part_attrforce = AF::line;
    AR part_attrreconst = AR::line;
    int verb = 0;
    bool vtkbin = false;
  };

  // Constructor
  // fcu: initial volume fraction
  // fccl: initial color
  // mfc: boundary conditions for volume fraction
  // ffv: pointer to mixture flux
  // fcs: poitner to volume sources
  // t,dt: initial time and timestep
  // par: parameters
  Vof(M& m, const FieldCell<Scal>& fcu, const FieldCell<Scal>& fccl,
      const MapFace<std::shared_ptr<CondFace>>& mfc,
      const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs,
      double t, double dt, std::shared_ptr<Par> par);
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
  const FieldCell<Scal>& GetField(Layers l, size_t i) const;
  const FieldCell<Scal>& GetField(size_t i) const;
  // ...
  using P::GetField;
  // Plane constant
  const FieldCell<Scal>& GetAlpha(size_t i) const;
  // Normal to interface
  const FieldCell<Vect>& GetNormal(size_t i) const;
  // Number of layers
  size_t GetNumLayers() const;
  // Layer mask
  const FieldCell<bool>& GetMask(size_t i) const;
  // Dependent cells
  const FieldCell<bool>& GetDepend(size_t i) const;
  // Color
  const FieldCell<Scal>& GetColor(size_t i) const;
  const FieldCell<Scal>& GetColor() const;
  // Height function
  const FieldCell<Vect>& GetHeight() const;
  // Default curvature 
  const FieldCell<Scal>& GetCurv() const override;
  // Default curvature 
  const FieldCell<Scal>& GetCurv(size_t i) const;
  // Curvature from height functions
  const FieldCell<Scal>& GetCurvH() const;
  // Curvature from particles
  const FieldCell<Scal>& GetCurvP() const;

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};

} // namespace solver
