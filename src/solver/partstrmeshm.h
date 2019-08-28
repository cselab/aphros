#pragma once

#include <memory>

#include "geom/mesh.h"
#include "partstr.h"
#include "dump/dumper.h"
#include "cond.h"

namespace solver {

template <class M_>
class PartStrMeshM {
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

  PartStrMeshM(M& m, std::shared_ptr<Par> par);
  ~PartStrMeshM();

  // Computes curvature with particles.
  // vfcu: volume fraction
  // vfca: plane constant
  // vfcn: normal
  // vfci: interface mask (1: contains interface)
  void Part(
      std::vector<const FieldCell<Scal>*> vfcu,
      std::vector<const FieldCell<Scal>*> vfca,
      std::vector<const FieldCell<Vect>*> vfcn,
      std::vector<const FieldCell<bool>*> vfci,
      const MapFace<std::shared_ptr<CondFace>>& mfc);
  // Dump particles to csv.
  // vfca: plane constant
  // vfcn: normal
  // n: frame index
  // t: time
  void DumpParticles(std::vector<const FieldCell<Scal>*> vfca,
                     std::vector<const FieldCell<Vect>*> vfcn,
                     size_t id, Scal t);
  void DumpPartInter(std::vector<const FieldCell<Scal>*> vfca,
                     std::vector<const FieldCell<Vect>*> vfcn,
                     size_t id, Scal t);
  // Returns curvature field from last call of Part()
  const FieldCell<Scal>& GetCurv(size_t l);

 private:
  struct Imp; // implementation
  //std::unique_ptr<Imp> imp;
  std::unique_ptr<Imp> imp;
};

} // namespace solver
