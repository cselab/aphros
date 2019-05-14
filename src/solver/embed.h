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
    InitFaces(fnf_, fft_, ffpoly_, ffs_, m);
    InitCells(fnf_, ffs_, fct_, fcn_, fca_, fcs_, m);
  }
  const FieldCell<Type>& GetCellType() const { return fct_; }
  const FieldCell<Vect>& GetNormal() const { return fcn_; }
  const FieldCell<Scal>& GetPlane() const { return fca_; }
  const FieldFace<Scal>& GetFaceArea() const { return ffs_; }
  const FieldCell<Scal>& GetCellArea() const { return fcs_; }

  // Dump cut polygons
  void DumpPoly() {
    auto sem = m.GetSem("dumppoly");
    if (sem("local")) {
      dl_.clear();
      dld_.clear();
      dls_.clear();

      for (auto f : m.Faces()) {
        if (fct_[m.GetNeighbourCell(f, 0)] == Type::cut ||
            fct_[m.GetNeighbourCell(f, 1)] == Type::cut) {
          size_t d = m.GetIndexFaces().GetDir(f);
          if (fft_[f] == Type::cut || fft_[f] == Type::regular) {
            if (fft_[f] == Type::cut) {
              dl_.push_back(ffpoly_[f]);
            } else {
              dl_.push_back(GetPoly(f, m));
            }
            dld_.push_back(d);
            dls_.push_back(ffs_[f]);
          }
        }
      }

      for (auto c : m.Cells()) {
        if (fct_[c] == Type::cut) {
          auto xx = R::GetCutPoly(
              m.GetCenter(c), fcn_[c], fca_[c], m.GetCellSize());
          dl_.push_back(xx);
          dld_.push_back(3);
          dls_.push_back(fcs_[c]);
        }
      }

      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl_));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dld_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::string fn = GetDumpName("eb", ".vtk", 0);
        std::cout << std::fixed << std::setprecision(8)
            << "dump" 
            << " to " << fn << std::endl;
        WriteVtkPoly(fn, dl_, {&dld_, &dls_}, {"dir", "area"}, 
                     "Embedded boundary", true);
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
  // Returns face polygon.
  // x0,x1: points
  // f0,f1: values
  static std::vector<Vect> GetPoly(IdxFace f, const M& m) {
    std::vector<Vect> xx;
    for (size_t e = 0; e < m.GetNumNeighbourNodes(f); ++e) {
      auto n = m.GetNeighbourNode(f, e);
      xx.push_back(m.GetNode(n));
    }
    return xx;
  }
  // Determines the face types, constructs polygons and computes fractions.
  // fnf: f on nodes
  // fft: type of faces
  // ffpoly: if fft=1, polygon representing f < 0; otherwise empty
  // ffs: face area for which f > 0
  static void InitFaces(const FieldNode<Scal>& fnf, 
                        FieldFace<Type>& fft,
                        FieldFace<std::vector<Vect>>& ffpoly,
                        FieldFace<Scal>& ffs,
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
          xx.push_back(GetIso(x, xp, f, fp));
          cut = true;
        }
      }
      fft[f] = (cut ? Type::cut : xx.empty() ? Type::excluded : Type::regular);
      switch (fft[f]) {
        case Type::regular:
          ffs[f] = m.GetArea(f);
          break;
        case Type::cut:
          ffs[f] = std::abs(R::GetArea(xx, m.GetNormal(f)));
          ffpoly[f] = xx;
          break;
        case Type::excluded:
          ffs[f] = 0;
          break;
      }
    }
  }
  // Determines the cell types, normals and plane constants.
  // fnf: f on nodes
  // ffs: face area for which f > 0
  // Output:
  // fct: cell types
  // fcn: normals
  // fca: plane constants
  // fcs: polygon area
  static void InitCells(const FieldNode<Scal>& fnf, const FieldFace<Scal>& ffs,
                        FieldCell<Type>& fct, FieldCell<Vect>& fcn,
                        FieldCell<Scal>& fca, FieldCell<Scal>& fcs, 
                        const M& m) {
    fct.Reinit(m);
    fcn.Reinit(m);
    fca.Reinit(m);
    fcs.Reinit(m);
    for (auto c : m.Cells()) {
      size_t q = 0; // number of nodes with f > 0
      const size_t mi = m.GetNumNeighbourNodes(c);
      for (size_t i = 0; i < mi; ++i) {
        IdxNode n = m.GetNeighbourNode(c, i);
        if (fnf[n] > 0) {
          ++q;
        }
      }
      fct[c] = (q == 0 ? Type::regular : q < mi ? Type::cut : Type::excluded);

      if (fct[c] == Type::cut) {
        // calc normal
        {
          Vect n(0);
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetNeighbourFace(c, q);
            n += m.GetOutwardSurface(c, q) * ffs[f];
          }
          fcn[c] = n / n.norm1();
        }

        // calc plane constant
        {
          Scal a = 0;
          Scal aw = 0;
          for (auto q : m.Nci(c)) {
            // FIXME: edges traversed twice
            IdxFace f = m.GetNeighbourFace(c, q);
            const size_t em = m.GetNumNeighbourNodes(f);
            for (size_t e = 0; e < em; ++e) {
              size_t ep = (e + 1) % em;
              IdxNode n = m.GetNeighbourNode(f, e);
              IdxNode np = m.GetNeighbourNode(f, ep);
              Scal f = fnf[n];
              Scal fp = fnf[np];
              Vect x = m.GetNode(n);
              Vect xp = m.GetNode(np);

              if ((f < 0) != (fp < 0)) {
                a += fcn[c].dot(GetIso(x, xp, f, fp) - m.GetCenter(c));
                aw += 1.;
              }
            }
          }
          fca[c] = a / aw;
        }

        auto xx = R::GetCutPoly(
            m.GetCenter(c), fcn[c], fca[c], m.GetCellSize());
        fcs[c] = std::abs(R::GetArea(xx, fcn[c]));
      }
    }
  }

  M& m;
  // nodes
  FieldNode<Scal> fnf_;  // level-set
  // faces
  FieldFace<Type> fft_;  // face type (0: regular, 1: cut, 2: excluded)
  FieldFace<std::vector<Vect>> ffpoly_;  // polygon representing f < 0
  FieldFace<Scal> ffs_;  // area for which f > 0
  // cells
  FieldCell<Type> fct_;  // cell type (0: regular, 1: cut, 2: excluded)
  FieldCell<Vect> fcn_;  // normal
  FieldCell<Scal> fca_;  // plane constant
  FieldCell<Scal> fcs_;  // area of polygon
  // tmp
  std::vector<std::vector<Vect>> dl_;  // dump poly, polygon
  std::vector<Scal> dld_;              // dump poly, direction
  std::vector<Scal> dls_;              // dump poly, area
};


} // namespace solver
