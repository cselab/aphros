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
  enum class Type {regular, cut, excluded};
  // fnf: level-set function on nodes, interface at fnf=0
  Embed(M& m, const FieldNode<Scal>& fnf)
      : m(m), fnf_(fnf) {
    InitInterface(fnf_);
  }
  const FieldCell<Type>& GetCellType() const { return fct_; }
  const FieldCell<Scal>& GetNormal() const { return fcn_; }
  const FieldCell<Vect>& GetPlane() const { return fca_; }

  // Dump cut polygons
  void DumpPoly() {
    auto sem = m.GetSem("dumppoly");
    if (sem("local")) {
      dl_.clear();
      dlc_.clear();
      ffs_.Reinit(m, 0.);


      for (auto f : m.Faces()) {
        const size_t em = m.GetNumNeighbourNodes(f);
        std::vector<Vect> xx;
        for (size_t e = 0; e < em; ++e) {
          size_t ep = (e + 1) % em;
          IdxNode n = m.GetNeighbourNode(f, e);
          IdxNode np = m.GetNeighbourNode(f, ep);
          Scal f = fnf_[n];
          Scal fp = fnf_[np];
          Vect x = m.GetNode(n);
          Vect xp = m.GetNode(np);
          if (f < iso_) {
            xx.push_back(x);
          }
          if ((f < iso_) != (fp < iso_)) {
            xx.push_back(Q(x, xp, f, fp));
          }
        }
        ffs_[f] = std::abs(R::GetArea(xx, m.GetNormal(f))) / m.GetArea(f);

        if (fct_[m.GetNeighbourCell(f, 0)] == 1 ||
            fct_[m.GetNeighbourCell(f, 1)] == 1) {
          dl_.push_back(xx);
          dlc_.push_back(size_t(m.GetIndexFaces().GetDir(f)));
        }
      }

      for (auto c : m.Cells()) {
        if (fct_[c] == 1) {
          // calc normal
          Vect n(0);
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetNeighbourFace(c, q);
            n += m.GetOutwardSurface(c, q) * ffs_[f];
          }
          n /= n.norm1();

          // calc plane constant
          Scal a = 0;
          Scal aw = 0;
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetNeighbourFace(c, q);
            auto V = [this,Q,f](size_t i, size_t ip) {
              IdxNode n = m.GetNeighbourNode(f, i);
              IdxNode np = m.GetNeighbourNode(f, ip);
              return Q(m.GetNode(n), m.GetNode(np), fnf_[n], fnf_[np]);
            };
            auto F = [this,f](size_t i) -> bool {
              IdxNode n = m.GetNeighbourNode(f, i);
              return fnf_[n] > 0;
            };
            for (size_t i = 0; i < 4; ++i) {
              size_t ip = (i + 1) % 4;
              if (F(i) != F(ip)) {
                a += n.dot(V(i, ip) - m.GetCenter(c));
                aw += 1.;
              }
            }
          }
          if (aw != 0) {
            a /= aw;
          }

          dl_.push_back(R::GetCutPoly(m.GetCenter(c), n, a, m.GetCellSize()));
          dlc_.push_back(3);
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
        WriteVtkPoly(fn, dl_, {&dlc_}, {"dir"}, "Embedded boundary", true);
      }
    }
  }

 private:
  // Returns point at which interpolant has value 0.
  // x0,x1: points
  // f0,f1: values
  static Vect GetIso(Vect x0, Vect x1, Scal f0, Scal f1) {
    return (x0 * f1 - x1 * f0) / (f1 - f0);
  }
  // Determines the type of faces, constructs polygons and computes fractions.
  // fnf: f on nodes
  // fft: type of faces
  // ffpoly: if fft=1, polygon representing f < 0; otherwise empty
  // ffs: fraction of face for which f > 0
  static void InitFaces(const FieldNode<Scal>& fnf, 
                        FieldFace<Type>& fft,
                        FieldFace<std::vector<Vect>>& ffpoly,
                        const M& m) {
    fft.Reinit(m);
    ffpoly.Reinit(m);
    ffs.Reinit(m);
    for (auto f : m.Faces()) {
      const size_t em = m.GetNumNeighbourNodes(f);
      std::vector<Vect> xx;
      bool cut = false;
      for (size_t e = 0; e < em; ++e) {
        size_t ep = (e + 1) % em;
        IdxNode n = m.GetNeighbourNode(f, e);
        IdxNode np = m.GetNeighbourNode(f, ep);
        Scal f = fnf[n];
        Scal fp = fnf[np];
        Vect x = m.GetNode(n);
        Vect xp = m.GetNode(np);
        if (f < 0) {
          xx.push_back(x);
        }
        if ((f < 0) != (fp < 0)) {
          xx.push_back(GetIso(x, xp, f, fp, iso_));
          cut = true;
        }
      }
      fft[f] = (cut ? Type::cut : xx.empty() ? Type::excluded : Type::regular);
      if (fft[f] == Type::cut) {
        ffpoly[f] = xx;
      }
    }
  }
  void InitInterface(const FieldNode<Scal>& fnf) {
    fcn_.Reinit(m, Vect(1));
    fca_.Reinit(m, 0);
    fct_.Reinit(m);
    for (auto c : m.Cells()) {
      size_t q = 0; // number of nodes with f > 0
      const size_t mi = m.GetNumNeighbourNodes(c);
      for (size_t i = 0; i < mi; ++i) {
        IdxNode n = m.GetNeighbourNode(c, i);
        if (fnf[n] > iso_) {
          ++q;
        }
      }
      fct_[c] = (q == 0 ? 0 : q < mi ? 1 : 2);
    }
  }

  M& m;
  Scal iso_ = 0.;        // isovalue
  FieldCell<Type> fct_;  // cell type [a] (0: regular, 1: cut, 2: excluded)
  FieldFace<Type> fft_;  // face type [a] (0: regular, 1: cut, 2: excluded)
  FieldCell<Vect> fcn_;  // normal [a]
  FieldCell<Scal> fca_;  // plane constant [a]
  FieldNode<Scal> fnf_;  // level-set [a]
  FieldFace<Scal> ffs_;  // inner face area fraction [a]
  FieldFace<std::vector<Vect>> ffpoly_;  // polygon representing f < 0
  std::vector<std::vector<Vect>> dl_;  // tmp: dump poly
  std::vector<Scal> dlc_;              // tmp: dump poly
};


} // namespace solver
