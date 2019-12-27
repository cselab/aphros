#pragma once

#include <limits>

#include "dump/vtk.h"
#include "reconst.h"
#include "solver.h"

// Embedded boundaries.
template <class M_>
class Embed {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using R = Reconst<Scal>;
  static constexpr size_t dim = M::dim;

 public:
  enum class Type { regular, cut, excluded };
  // fnf: level-set function on nodes, interface at fnf=0
  Embed(M& m, const FieldNode<Scal>& fnf) : m(m), fnf_(fnf) {
    InitFaces(fnf_, fft_, ffpoly_, ffs_, m);
    InitCells(fnf_, ffs_, fct_, fcn_, fca_, fcs_, fcd_, fcv_, m);
  }
  const FieldCell<Type>& GetCellType() const {
    return fct_;
  }
  const FieldCell<Vect>& GetNormal() const {
    return fcn_;
  }
  const FieldCell<Scal>& GetPlane() const {
    return fca_;
  }
  const FieldFace<Type>& GetFaceType() const {
    return fft_;
  }
  const FieldFace<Scal>& GetFaceArea() const {
    return ffs_;
  }
  const FieldCell<Scal>& GetCellArea() const {
    return fcs_;
  }
  const FieldCell<Scal>& GetCellOffset() const {
    return fcd_;
  }
  const FieldCell<Scal>& GetCellVolume() const {
    return fcv_;
  }
  const FieldFace<std::vector<Vect>>& GetPoly() const {
    return ffpoly_;
  }
  Scal GetCellVolume(IdxCell c) const {
    return fcv_[c];
  }
  Vect GetFaceCenter(IdxCell c) const {
    if (fct_[c] == Type::cut) {
      return R::GetCenter(
          R::GetCutPoly(m.GetCenter(c), fcn_[c], fca_[c], m.GetCellSize()));
    }
    return GetNan<Vect>();
  }
  Vect GetFaceCenter(IdxFace f) const {
    switch (fft_[f]) {
      case Type::regular:
        return m.GetCenter(f);
      case Type::cut:
        return R::GetCenter(ffpoly_[f]);
      default:
        return GetNan<Vect>();
    }
  }
  Vect GetCellCenter(IdxCell c) const {
    switch (fct_[c]) {
      case Type::regular:
        return m.GetCenter(c);
      case Type::cut: {
        const Scal w = GetCellArea()[c];
        Vect sum = GetFaceCenter(c) * w;
        Scal sumw = w;
        for (auto q : m.Nci(c)) {
          auto f = m.GetFace(c, q);
          if (fft_[f] != Type::excluded) {
            const Scal w = GetFaceArea()[f];
            sum += GetFaceCenter(f) * w;
            sumw += w;
          }
        }
        return sum / sumw;
      }
      default:
        return GetNan<Vect>();
    }
  }

  // Dump cut polygons
  // ffs: face area for which f > 0
  // fft: type of faces
  // fcs: polygon area
  // fct: cell types
  // fcn: normals
  // fca: plane constants
  // ffpoly: polygon representing f < 0
  static void DumpPoly(
      const FieldFace<Scal>& ffs, const FieldFace<Type>& fft,
      const FieldCell<Scal>& fcs, const FieldCell<Type>& fct,
      const FieldCell<Vect>& fcn, const FieldCell<Scal>& fca,
      const FieldFace<std::vector<Vect>>& ffpoly, M& m) {
    auto sem = m.GetSem("dumppoly");
    struct {
      std::vector<std::vector<Vect>> dl; // polygon
      std::vector<Scal> dld; // direction
      std::vector<Scal> dls; // area
    } * ctx(sem);
    auto& dl = ctx->dl;
    auto& dld = ctx->dld;
    auto& dls = ctx->dls;
    if (sem("local")) {
      for (auto f : m.Faces()) {
        if (fct[m.GetCell(f, 0)] == Type::cut ||
            fct[m.GetCell(f, 1)] == Type::cut) {
          size_t d(m.GetIndexFaces().GetDir(f));
          if (fft[f] == Type::cut || fft[f] == Type::regular) {
            if (fft[f] == Type::cut) {
              dl.push_back(ffpoly[f]);
            } else {
              dl.push_back(GetPoly(f, m));
            }
            dld.push_back(d);
            dls.push_back(ffs[f]);
          }
        }
      }

      for (auto c : m.Cells()) {
        if (fct[c] == Type::cut) {
          auto xx =
              R::GetCutPoly(m.GetCenter(c), fcn[c], fca[c], m.GetCellSize());
          dl.push_back(xx);
          dld.push_back(3);
          dls.push_back(fcs[c]);
        }
      }

      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dld));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::string fn = GetDumpName("eb", ".vtk", 0);
        std::cout << std::fixed << std::setprecision(8) << "dump"
                  << " to " << fn << std::endl;
        WriteVtkPoly<Vect>(
            fn, dl, nullptr, {&dld, &dls}, {"dir", "area"}, "Embedded boundary",
            true, true, true);
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
    for (size_t e = 0; e < m.GetNumNodes(f); ++e) {
      auto n = m.GetNode(f, e);
      xx.push_back(m.GetNode(n));
    }
    return xx;
  }
  // Determines the face types, constructs polygons and computes fractions.
  // fnf: f on nodes
  // fft: type of faces
  // ffpoly: if fft=1, polygon representing f < 0; otherwise empty
  // ffs: face area for which f > 0
  static void InitFaces(
      const FieldNode<Scal>& fnf, FieldFace<Type>& fft,
      FieldFace<std::vector<Vect>>& ffpoly, FieldFace<Scal>& ffs, const M& m) {
    fft.Reinit(m);
    ffpoly.Reinit(m);
    ffs.Reinit(m, 0);
    for (auto f : m.Faces()) {
      const size_t em = m.GetNumNodes(f);
      std::vector<Vect> xx;
      bool cut = false;
      for (size_t e = 0; e < em; ++e) {
        size_t ep = (e + 1) % em;
        IdxNode n = m.GetNode(f, e);
        IdxNode np = m.GetNode(f, ep);
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
  // ffpoly: polygon representing f < 0
  // Output:
  // fct: cell types
  // fcn: normals
  // fca: plane constants
  // fcs: polygon area
  static void InitCells(
      const FieldNode<Scal>& fnf, const FieldFace<Scal>& ffs,
      FieldCell<Type>& fct, FieldCell<Vect>& fcn, FieldCell<Scal>& fca,
      FieldCell<Scal>& fcs, FieldCell<Scal>& fcd, FieldCell<Scal>& fcv,
      const M& m) {
    fct.Reinit(m);
    fcn.Reinit(m);
    fca.Reinit(m);
    fcs.Reinit(m, 0);
    fcd.Reinit(m);
    fcv.Reinit(m, 0);
    for (auto c : m.Cells()) {
      size_t q = 0; // number of nodes with f > 0
      const size_t mi = m.GetNumNodes(c);
      for (size_t i = 0; i < mi; ++i) {
        IdxNode n = m.GetNode(c, i);
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
            IdxFace f = m.GetFace(c, q);
            n += m.GetNormal(f) * ffs[f] * m.GetOutwardFactor(c, q);
          }
          fcn[c] = -n / n.norm();
        }

        // calc plane constant
        {
          Scal a = 0;
          Scal aw = 0;
          for (auto q : m.Nci(c)) {
            // FIXME: edges traversed twice
            IdxFace f = m.GetFace(c, q);
            const size_t em = m.GetNumNodes(f);
            for (size_t e = 0; e < em; ++e) {
              size_t ep = (e + 1) % em;
              IdxNode n = m.GetNode(f, e);
              IdxNode np = m.GetNode(f, ep);
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

        const auto h = m.GetCellSize();

        auto xx = R::GetCutPoly(m.GetCenter(c), fcn[c], fca[c], h);
        fcs[c] = std::abs(R::GetArea(xx, fcn[c]));

        // normal distance from center to face
        auto xc = R::GetCenter(xx);
        fcd[c] = (xc - m.GetCenter(c)).dot(fcn[c]);
        // cut cell volume
        fcv[c] = R::GetLineU(fcn[c], fca[c], h) * m.GetVolume(c);
      } else if (fct[c] == Type::regular) {
        fcv[c] = m.GetVolume(c);
      }
    }
  }

  M& m;
  // nodes
  FieldNode<Scal> fnf_; // level-set
  // faces
  FieldFace<Type> fft_; // face type (0: regular, 1: cut, 2: excluded)
  FieldFace<std::vector<Vect>> ffpoly_; // polygon representing f < 0
  FieldFace<Scal> ffs_; // area for which f > 0
  // cells
  FieldCell<Type> fct_; // cell type (0: regular, 1: cut, 2: excluded)
  FieldCell<Vect> fcn_; // unit outer normal
  FieldCell<Scal> fca_; // plane constant
  FieldCell<Scal> fcs_; // area of polygon
  FieldCell<Scal> fcd_; // normal component of displacement from cell center
  FieldCell<Scal> fcv_; // volume of cut cell
};
