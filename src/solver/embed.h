#pragma once

#include <limits>

#include "solver.h"
#include "reconst.h"
#include "dump/vtk.h"

namespace solver {

// Embedded boundaries.
template <class M_>
class Embed {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using R = Reconst<Scal>;
  static constexpr size_t dim = M::dim;

 public:
  // fnf: level-set function on nodes, interface at fnf=0
  Embed(M& m, const FieldNode<Scal>& fnf)
      : m(m), fnf_(fnf) {
    InitInterface(fnf_);
  }
  const FieldCell<char>& GetCellType() const { return fct_; }
  const FieldCell<Scal>& GetNormal() const { return fcn_; }
  const FieldCell<Vect>& GetPlane() const { return fca_; }

  // Dump cut polygons
  void DumpPoly() {
    auto sem = m.GetSem("dumppoly");
    if (sem("local")) {
      dl_.clear();
      dlc_.clear();
      auto h = m.GetCellSize();

      auto Q = [](Vect x0, Vect x1, Scal f0, Scal f1) {
        return (x0 * f1 - x1 * f0) / (f1 - f0);
      };
      auto QB = [this](IdxNode n0, IdxNode n1) -> bool {
        return (fnf_[n0] > 0) != (fnf_[n1] > 0);
      };
      auto QV = [this,Q](IdxNode n0, IdxNode n1) -> Vect {
        return Q(m.GetNode(n0), m.GetNode(n1), fnf_[n0], fnf_[n1]);
      };
      auto MI = [this](IdxFace f) -> size_t {
        return GVect<Scal, 4>(
            fnf_[m.GetNeighbourNode(f, 0)], 
            fnf_[m.GetNeighbourNode(f, 1)], 
            fnf_[m.GetNeighbourNode(f, 2)], 
            fnf_[m.GetNeighbourNode(f, 3)]).argmin();
      };
      (void) MI;

      for (auto f : m.Faces()) {
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        if (fct_[cm] == 1 || fct_[cp] == 1) {
          IdxNode n0 = m.GetNeighbourNode(f, 0);
          IdxNode n1 = m.GetNeighbourNode(f, 1);
          IdxNode n2 = m.GetNeighbourNode(f, 2);
          IdxNode n3 = m.GetNeighbourNode(f, 3);
          std::vector<Vect> xx;
          if (QB(n0,n1)) { xx.push_back(QV(n0,n1)); }
          if (QB(n1,n2)) { xx.push_back(QV(n1,n2)); }
          if (QB(n2,n3)) { xx.push_back(QV(n2,n3)); }
          if (QB(n3,n0)) { xx.push_back(QV(n3,n0)); }
          if (xx.size() == 2) {
            dl_.push_back(xx);
          }
        }
      }

      if(0)
      for (auto c : m.Cells()) {
        if (fct_[c] == 1) {
          dl_.push_back(R::GetCutPoly(m.GetCenter(c), fcn_[c], fca_[c], h));
          //dlc_.push_back(m.GetHash(c));
        }
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl_));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dlc_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::string fn = GetDumpName("eb", ".vtk", 0);
        std::cout << std::fixed << std::setprecision(8)
            << "dump" 
            << " to " << fn << std::endl;
        //WriteVtkPoly(fn, dl_, {&dlc_}, {"c"}, "Embedded boundary", true);
        WriteVtkPoly(fn, dl_, {}, {}, "Embedded boundary", true);
      }
    }
  }

 private:
  void InitInterface(const FieldNode<Scal>& fnf) {
    fcn_.Reinit(m, Vect(1));
    fca_.Reinit(m, 0);
    fct_.Reinit(m);
    for (auto c : m.Cells()) {
      size_t q = 0; // number of nodes with f > 0
      const size_t mi = m.GetNumNeighbourNodes(c);
      for (size_t i = 0; i < mi; ++i) {
        IdxNode n = m.GetNeighbourNode(c, i);
        if (fnf[n] > 0) {
          ++q;
        }
      }
      fct_[c] = (q == 0 ? 0 : q < mi ? 1 : 2);
    }
  }

  M& m;
  FieldCell<char> fct_; // cell type [a] (0: regular, 1: cut, 2: excluded)
  FieldCell<Vect> fcn_; // normal [a]
  FieldCell<Scal> fca_; // plane constant [a]
  FieldNode<Scal> fnf_; // level-set [a]
  std::vector<std::vector<Vect>> dl_; // dump poly
  std::vector<Scal> dlc_; // dump poly
};


} // namespace solver
