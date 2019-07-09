#pragma once

#include <memory>

#include "geom/mesh.h"
#include "partstr.h"
#include "dump/dumper.h"
#include "cond.h"

namespace solver {

template <class M_>
class PartStrMesh {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using PS = PartStr<Scal>;

  enum class AF { // attraction force type
      line   // nearest line
    , center // nearest line center
    , volume // fluid volume
  };

  enum class AR { // attraction reconstruction type
      line   // interface line
    , volume // fluid volume
  };

  struct Par {
    std::shared_ptr<typename PS::Par> ps;
    Scal intth = 0; // interface threshold for particle seed
    size_t dump_fr = -1; // num frames dump
    size_t ns = 2; // number of strings per cell
    Scal tol = 0.01; // tolerance
    size_t itermax = 20;
    int verb = 0; // debug output
    size_t dim = 3;
    bool bcc_reflect = false;
    AR attrreconst = AR::line;
  };

  PartStrMesh(M& m, std::shared_ptr<Par> par);
  ~PartStrMesh();

  // Computes curvature with particles.
  // fcu: volume fraction
  // fca: plane constant
  // fcn: normal
  // fci: interface mask (1: contains interface)
  void Part(const FieldCell<Scal>& uc, 
      FieldCell<Scal>& fca, FieldCell<Vect>& fcn, FieldCell<bool>& fci,
      FieldCell<Scal>& fck,
      const MapFace<std::shared_ptr<CondFace>>& mfc);
  // Dump particles to csv.
  // fca: plane constant
  // fcn: normal
  // n: frame index
  // t: time
  void DumpParticles(FieldCell<Scal>& fca, FieldCell<Vect>& fcn,
                     size_t id, Scal t);
  void DumpPartInter(FieldCell<Scal>& fca, FieldCell<Vect>& fcn,
                     size_t id, Scal t);
  // Returns curvature field from last call of Part()
  const FieldCell<Scal>& GetCurv();

 private:
  struct Imp; // implementation
  //std::unique_ptr<Imp> imp;
  std::unique_ptr<Imp> imp;
};

} // namespace solver
