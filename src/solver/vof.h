#pragma once

#include <memory>

#include "advection.h"
#include "dump/dumper.h"

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
    bool part_verb = false; // debug output
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
    // XXX: adhoc
    Scal bcc_k0 = 1.;   // mul corrections to bc 
    Scal bcc_k1 = 1.;   
    Scal bcc_t0 = -1.;   // duration of phases (one negative to disable)
    Scal bcc_t1 = -1.;   
    Scal bcc_y0 = -1e10; // overwrite u=0 if y<y0 or y>y1
    Scal bcc_y1 = 1e10;  // (to remove periodic conditions)
    int part_constr = 0; // 0: no constraints
                         // 1: fixed distance, constant angle
                         // 2: fixed distance, linear angle
    Scal part_segcirc = 1.; // factor for circular segment
    size_t part_np = 11; // number of particles per string
    size_t part_ns = 4; // number of strings per cell
    size_t part_itermax = 100; // particles itermax
    Scal part_tol = 0.01; // tolerance
    Scal part_tmax = 180.; 
    Scal part_dtmax = 10.; 
    Scal part_anglim = 90.; 
    enum class Attr { // attraction force type
        interface  // to reconstructed interface
      , volume // to fluid volume
    }; 
    Attr part_attr = Attr::interface;
  };

  // Constructor
  Vof(M& m, const FieldCell<Scal>& fcu,
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

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};

} // namespace solver
