// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#include "embed.h"
#include "dump/dumper.h"
#include "dump/vtk.h"
#include "func/init_u.h"

template <class M>
void Embed<M>::Init(const FieldNode<Scal>& fnl) {
  auto sem = m.GetSem("init");
  if (sem()) {
    fnl_ = fnl;
    InitFaces(fnl_, fft_, ffpoly_, ffs_, m);
    InitCells(fnl_, ffs_, fct_, fcn_, fca_, fcs_, fcv_, m);
    InitRedistr(fct_, fcv_, fcs_, fc_redistr_, mc_redistr_, m);
    m.Comm(&fc_redistr_);

    // volume of neighbor cells
    fcvst3_.Reinit(m, 0);
    for (auto c : eb.Cells()) {
      for (auto cn : eb.Stencil(c)) {
        fcvst3_[c] += eb.GetVolume(cn);
      }
    }
    m.Comm(&fcvst3_);
  }
}

template <class M>
void Embed<M>::DumpPoly(
    const std::string fn, const FieldFace<Scal>& ffs,
    const FieldFace<Type>& fft, const FieldCell<Scal>& fcs,
    const FieldCell<Type>& fct, const FieldCell<Vect>& fcn,
    const FieldCell<Scal>& fca, const FieldFace<std::vector<Vect>>& ffpoly,
    bool vtkbin, bool vtkmerge, M& m) {
  auto sem = m.GetSem("dumppoly");
  struct {
    std::vector<std::vector<Vect>> dl; // polygon
    std::vector<Scal> dld; // direction
    std::vector<Scal> dls; // area
    std::vector<Scal> dlf; // 0: cell, 1: face
  } * ctx(sem);
  auto& dl = ctx->dl;
  auto& dld = ctx->dld;
  auto& dls = ctx->dls;
  auto& dlf = ctx->dlf;
  if (sem("local")) {
    for (auto f : m.Faces()) {
      if (fct[m.GetCell(f, 0)] == Type::cut ||
          fct[m.GetCell(f, 1)] == Type::cut) {
        size_t d(m.GetIndexFaces().GetDir(f));
        if (fft[f] == Type::cut || fft[f] == Type::regular) {
          if (fft[f] == Type::cut) {
            dl.push_back(ffpoly[f]);
          } else {
            dl.push_back(GetRegularFacePoly(f, m));
          }
          dld.push_back(d);
          dls.push_back(ffs[f]);
          dlf.push_back(1);
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
        dlf.push_back(0);
      }
    }

    using TV = typename M::template OpCatVT<Vect>;
    m.Reduce(std::make_shared<TV>(&dl));
    using TS = typename M::template OpCatT<Scal>;
    m.Reduce(std::make_shared<TS>(&dld));
    m.Reduce(std::make_shared<TS>(&dls));
    m.Reduce(std::make_shared<TS>(&dlf));
  }
  if (sem("write")) {
    if (m.IsRoot()) {
      std::cout << std::fixed << std::setprecision(8) << "dump"
                << " to " << fn << std::endl;
      WriteVtkPoly<Vect>(
          fn, dl, nullptr, {&dld, &dls, &dlf}, {"dir", "area", "face"},
          "Embedded boundary", true, vtkbin, vtkmerge);
    }
  }
}

template <class M>
void Embed<M>::DumpPlaneSection(
    const std::string fn, const FieldCell<Type>& fct,
    const FieldCell<Vect>& fcn, const FieldCell<Scal>& fca, Vect plane_x,
    Vect plane_n, M& m) {
  auto sem = m.GetSem("dumpplanesection");
  struct {
    std::vector<std::vector<Vect>> dl;
  } * ctx(sem);
  auto& dl = ctx->dl;
  if (sem("local")) {
    for (auto c : m.Cells()) {
      if (fct[c] == Type::cut) {
        auto xx =
            R::GetCutPoly(m.GetCenter(c), fcn[c], fca[c], m.GetCellSize());
        std::array<Vect, 2> e; // ends of intersection
        if (R::GetInterPoly(xx, plane_x, plane_n, e)) {
          dl.push_back({e[0], e[1]});
        }
      }
    }
    using TV = typename M::template OpCatVT<Vect>;
    m.Reduce(std::make_shared<TV>(&dl));
  }
  if (sem("write")) {
    if (m.IsRoot()) {
      std::cout << std::fixed << std::setprecision(8) << "dump"
                << " to " << fn << std::endl;
      std::ofstream f(fn);
      for (auto e : dl) {
        f << e[0][0] << ' ' << e[0][1] << ' ' << e[0][2] << ' ';
        f << e[1][0] << ' ' << e[1][1] << ' ' << e[1][2] << '\n';
      }
    }
  }
}

template <class M>
void Embed<M>::InitFaces(
    const FieldNode<Scal>& fnl, FieldFace<Type>& fft,
    FieldFace<std::vector<Vect>>& ffpoly, FieldFace<Scal>& ffs, const M& m) {
  fft.Reinit(m);
  ffpoly.Reinit(m);
  ffs.Reinit(m, 0);
  for (auto f : m.AllFaces()) {
    const size_t em = m.GetNumNodes(f);
    std::vector<Vect> xx;
    bool cut = false;
    for (size_t e = 0; e < em; ++e) {
      const size_t ep = (e + 1) % em;
      const IdxNode n = m.GetNode(f, e);
      const IdxNode np = m.GetNode(f, ep);
      const Scal l = fnl[n];
      const Scal lp = fnl[np];
      const Vect x = m.GetNode(n);
      const Vect xp = m.GetNode(np);
      if (l >= 0) {
        xx.push_back(x);
      }
      if (l * lp < 0) {
        xx.push_back(GetIso(x, xp, l, lp));
        cut = true;
      }
    }
    fft[f] =
        (cut ? Type::cut : xx.size() <= 2 ? Type::excluded : Type::regular);
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

template <class M>
void Embed<M>::InitCells(
    const FieldNode<Scal>& fnl, const FieldFace<Scal>& ffs,
    FieldCell<Type>& fct, FieldCell<Vect>& fcn, FieldCell<Scal>& fca,
    FieldCell<Scal>& fcs, FieldCell<Scal>& fcv, const M& m) {
  fct.Reinit(m), GetNan<Scal>();
  fcn.Reinit(m), GetNan<Scal>();
  fca.Reinit(m, GetNan<Scal>());
  fcs.Reinit(m, GetNan<Scal>());
  fcv.Reinit(m, GetNan<Scal>());
  for (auto c : m.AllCells()) {
    size_t q = 0; // number of nodes with fnl > 0
    const size_t mi = m.GetNumNodes(c);
    for (size_t i = 0; i < mi; ++i) {
      IdxNode n = m.GetNode(c, i);
      if (fnl[n] >= 0) {
        ++q;
      }
    }
    fct[c] = (q == mi ? Type::regular : q > 0 ? Type::cut : Type::excluded);

    if (fct[c] == Type::cut) {
      // calc normal
      {
        Vect n(0);
        for (auto q : m.Nci(c)) {
          const IdxFace f = m.GetFace(c, q);
          n += m.GetNormal(f) * ffs[f] * m.GetOutwardFactor(c, q);
        }
        fcn[c] = -n / n.norm();
      }

      // calc plane constant
      {
        Scal a = 0;
        Scal aw = 0;
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetFace(c, q);
          const size_t em = m.GetNumNodes(f);
          for (size_t e = 0; e < em; ++e) {
            size_t ep = (e + 1) % em;
            IdxNode n = m.GetNode(f, e);
            IdxNode np = m.GetNode(f, ep);
            Scal l = fnl[n];
            Scal lp = fnl[np];
            Vect x = m.GetNode(n);
            Vect xp = m.GetNode(np);

            if (l * lp < 0) {
              a += fcn[c].dot(GetIso(x, xp, l, lp) - m.GetCenter(c));
              aw += 1.;
            }
          }
        }
        fca[c] = a / aw;
      }

      const auto h = m.GetCellSize();

      Vect sum(0);
      for (auto q : m.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        sum += m.GetNormal(f) * ffs[f] * m.GetOutwardFactor(c, q);
      }
      fcs[c] = -sum.dot(fcn[c]);

      fcv[c] = R::GetLineU(fcn[c], fca[c], h) * m.GetVolume(c);
    } else if (fct[c] == Type::regular) {
      fcv[c] = m.GetVolume(c);
    }
  }
}

template <class M>
void Embed<M>::InitRedistr(
    const FieldCell<Type>& fct, const FieldCell<Scal>& fcv,
    const FieldCell<Scal>& fcs, FieldCell<Scal>& fc_redistr,
    MapCell<std::vector<std::pair<IdxCell, Scal>>>& mc_redistr, const M& m) {
  fc_redistr.Reinit(m, 1);
  const Scal hreg = m.GetCellSize()[0]; // XXX adhoc cubic cells
  for (auto c : m.Cells()) {
    if (fct[c] == Type::cut) {
      const Scal a = std::min(1., fcv[c] / (hreg * fcs[c]));
      std::vector<std::pair<IdxCell, Scal>> vp;
      for (auto cn : m.Stencil(c)) {
        if (fct[cn] != Type::excluded && cn != c) {
          vp.push_back({cn, 1});
        }
      }
      Scal sum = 0;
      for (auto& p : vp) {
        sum += p.second * fcv[p.first];
      }
      for (auto& p : vp) {
        p.second *= fcv[c] * (1 - a) / sum;
      }
      fc_redistr[c] = a;
      mc_redistr[c] = vp;
    }
  }
}