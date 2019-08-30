#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>
#include <set>

#include <march.h>

#include "vof.h"
#include "geom/block.h"
#include "dump/vtk.h"
#include "reconst.h"
#include "normal.h"
#include "debug/isnan.h"
#include "partstr.h"
#include "multi.h"

namespace solver {

// Copies values from layer 0 to layers [1:] in cells with mask=0
// in stencil around each cell with mask=1.
// mfc: fields to propagate
// mask: if true, copy values from layer 0 in stencil 1.
// sw: stencil halfwidth
template <class M>
void Propagate(Multi<FieldCell<typename M::Scal>*> mfc,
               const Multi<const FieldCell<bool>*>& mask, size_t sw, M& m) {

  if (!mfc.size()) {
    return;
  }

  using MIdx = typename M::MIdx;
  auto& bc = m.GetIndexCells();
  const int sn = sw * 2 + 1;
  // block of offsets
  GBlock<IdxCell, M::dim> bo(MIdx(-sw), MIdx(sn)); 

  auto& fc0 = *mfc[0];
  for (size_t i = 1; i < mfc.size(); ++i) {
    auto& fc = *mfc[i];
    auto& ms = *mask[i];
    for (auto c : m.Cells()) {
      if (ms[c]) {
        MIdx w = bc.GetMIdx(c);
        for (MIdx wo : bo) {
          IdxCell cn = bc.GetIdx(w + wo); 
          if (!ms[cn]) {
            fc[cn] = fc0[cn];
          }
        }
      }
    }
  }
}

// Same as Propagate() but copies only from cells of the same color,
// writes fill otherwise.
// fccl: color
template <class M>
void Propagate(Multi<FieldCell<typename M::Scal>*> mfc,
               const Multi<const FieldCell<bool>*>& mfcm, 
               const Multi<const FieldCell<typename M::Scal>*>& mfccl, 
               size_t sw, M& m) {

  if (!mfc.size()) {
    return;
  }

  using MIdx = typename M::MIdx;
  auto& bc = m.GetIndexCells();
  const int sn = sw * 2 + 1;
  // block of offsets
  GBlock<IdxCell, M::dim> bo(MIdx(-sw), MIdx(sn)); 

  auto& fc0 = *mfc[0];
  auto& fccl0 = *mfccl[0];
  for (size_t i = 1; i < mfc.size(); ++i) {
    auto& fc = *mfc[i];
    auto& fcm = *mfcm[i];
    auto& fccl = *mfccl[i];
    for (auto c : m.Cells()) {
      if (fcm[c]) {
        MIdx w = bc.GetMIdx(c);
        for (MIdx wo : bo) {
          IdxCell cn = bc.GetIdx(w + wo); 
          if (!fcm[cn]) {
            if (fccl[c] == fccl0[cn]) {
              fc[cn] = fc0[cn];
            }
          }
        }
      }
    }
  }
}

// sw: stencil halfwidth
template <class M>
void Propagate(FieldCell<bool>& mask, size_t sw, M& m) {
  using MIdx = typename M::MIdx;
  auto& bc = m.GetIndexCells();
  const int sn = sw * 2 + 1;
  // block of offsets
  GBlock<IdxCell, M::dim> bo(MIdx(-sw), MIdx(sn)); 

  auto ms = mask;
  for (auto c : m.Cells()) {
    if (ms[c]) {
      MIdx w = bc.GetMIdx(c);
      for (MIdx wo : bo) {
        IdxCell cn = bc.GetIdx(w + wo);
        mask[cn] = true;
      }
    }
  }
}

template <class M>
void Propagate(const Multi<FieldCell<bool>*>& mfc, size_t sw, M& m) {
  for (auto p : mfc.data()) {
    Propagate(*p, sw, m);
  }
}

template <class T>
Multi<FieldCell<T>*> GetLayer(Multi<LayersData<FieldCell<T>>>& u, Layers l) {
  Multi<FieldCell<T>*> r(u.size());
  for (size_t i = 0; i < u.size(); ++i) {
    r[i] = &u[i].Get(l);
  }
  return r;
}

template <class M>
FieldCell<bool> And(const FieldCell<bool>& u, 
                    const FieldCell<bool>& v, const M& m) {
  FieldCell<bool> r = u;
  for (auto c : m.AllCells()) {
    r[c] = r[c] && v[c];
  }
  return r;
}

template <class M>
FieldCell<bool> Or(const FieldCell<bool>& u, 
                   const FieldCell<bool>& v, const M& m) {
  FieldCell<bool> r = u;
  for (auto c : m.AllCells()) {
    r[c] = r[c] || v[c];
  }
  return r;
}

template <class M_>
struct Vof<M_>::Imp {
  using Owner = Vof<M_>;
  using R = Reconst<Scal>;
  using PS = PartStr<Scal>;
  using PSM = PartStrMeshM<M>;
  static constexpr Scal kClNone = -1;
  static constexpr size_t dim = M::dim;
  using Vect2 = GVect<Scal, 2>;

  Imp(Owner* owner, const FieldCell<Scal>& fcu0, const FieldCell<Scal>& fccl0,
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m)
      , mfc_(mfc), fch_(m, Vect(0)), layers(0, 4)
  {
    fcu_.resize(layers.size());
    fcn_.resize(layers.size());
    fca_.resize(layers.size());
    fci_.resize(layers.size());
    fcdp_.resize(layers.size());
    fcm_.resize(layers.size());
    fck_.resize(layers.size());
    fccl_.resize(layers.size());
    fcclt_.resize(layers.size());
    ffvu_.resize(layers.size());
    ffcl_.resize(layers.size());
    ffi_.resize(layers.size());


    fcn_.InitAll(FieldCell<Vect>(m, Vect(0)));
    fca_.InitAll(FieldCell<Scal>(m, 0));
    fck_.InitAll(FieldCell<Scal>(m, 0));
    fcdp_.InitAll(FieldCell<bool>(m, false));
    fcm_.InitAll(FieldCell<bool>(m, false));
    fci_.InitAll(FieldCell<bool>(m, false));

    fcus_.time_curr = fcu0;

    fcu_[0].time_curr.Reinit(m, 0);
    fcu_.InitAll(fcu_[0]);

    fccl_.InitAll(FieldCell<Scal>(m, kClNone));
    ffvu_.InitAll(FieldFace<Scal>(m, 0));
    ffcl_.InitAll(FieldFace<Scal>(m, kClNone));
    ffi_.InitAll(FieldFace<bool>(m, false));


    fcu_[0].time_curr = fcu0;
    fccl_[0] = fccl0;
    fccls_ = fccl0;

    for (auto it : mfc_) {
      IdxFace f = it.GetIdx();
      mfvz_[f] = std::make_shared<
          CondFaceGradFixed<Vect>>(Vect(0), it.GetValue()->GetNci());
    }

    // particles
    auto ps = std::make_shared<typename PS::Par>();
    Update(ps.get());
    auto psm = std::make_shared<typename PSM::Par>();
    psm->ps = ps;
    psm->intth = par->part_intth;
    psm->ns = par->part_ns;
    psm->tol = par->part_tol;
    psm->itermax = par->part_itermax;
    psm->verb = par->part_verb;
    psm->dim = par->dim;
    psm->bcc_reflect = par->bcc_reflect;
    psm->dump_fr = par->part_dump_fr;
    psm->maxr = par->part_maxr;
    psm_ = std::unique_ptr<PSM>(new PSM(m, psm));
  }
  void Update(typename PS::Par* p) const {
    Scal hc = m.GetCellSize().norminf(); // cell size

    p->leq = par->part_h;
    p->relax = par->part_relax;
    p->constr = par->part_constr;
    p->npmax = par->part_np;
    p->segcirc = par->part_segcirc;
    p->hc = hc;

    p->kstr = par->part_kstr;
    p->kbend = par->part_kbend;
    p->kattr = par->part_kattr;
    p->bendmean = par->part_bendmean;
    p->ka = par->part_kattr;
    p->kt = par->part_kbend;
    p->kx = par->part_kstr;
    p->tmax = par->part_tmax;
    p->dtmax = par->part_dtmax;
    p->anglim = par->part_anglim;
    p->dn = par->part_dn;

    {
      using FT = typename PS::FT;
      switch (par->part_attrforce) { 
        case Par::AF::line:
          p->forcetype = FT::line;
          break;
        case Par::AF::center:
          p->forcetype = FT::center;
          break;
        case Par::AF::volume:
          p->forcetype = FT::volume;
          break;
        default:
          throw std::runtime_error("Update(): Unknown part_attrforce");
      }
    }
  }

  // Dump cut polygons
  void DumpPoly() {
    auto sem = m.GetSem("dumppoly");
    if (sem("local")) {
      dl_.clear();
      dlc_.clear();
      dll_.clear();
      dlcl_.clear();
      auto h = m.GetCellSize();
      for (auto i : layers) {
        auto& fcn = fcn_[i];
        auto& fca = fca_[i];
        auto& fci = fci_[i];
        auto& fccl = fccl_[i];
        auto& fcu = fcu_[i].iter_curr;
        for (auto c : m.Cells()) {
          Scal u = fcu[c];
          if (IsNan(u) || IsNan(fcn[c]) || IsNan(fca[c])) {
            continue;
          }
          const Scal th = par->poly_intth;
          if (fci[c] && u > th && u < 1. - th) {
            dl_.push_back(R::GetCutPoly(m.GetCenter(c), fcn[c], fca[c], h));
            dlc_.push_back(m.GetHash(c));
            dll_.push_back(i);
            dlcl_.push_back(fccl[c]);
          }
        }
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl_));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dlc_));
      m.Reduce(std::make_shared<TS>(&dll_));
      m.Reduce(std::make_shared<TS>(&dlcl_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::string fn = GetDumpName("s", ".vtk", par->dmp->GetN());
        std::cout << std::fixed << std::setprecision(8)
            << "dump" 
            << " t=" << owner_->GetTime() + owner_->GetTimeStep()
            << " to " << fn << std::endl;
        WriteVtkPoly(fn, dl_, {&dlc_, &dll_, &dlcl_}, {"c", "l", "cl"}, 
            "Reconstructed linear interface", true, par->vtkbin);
      }
    }
  }
  // uu: volume fraction in nodes
  // xc: cell center
  // h: cell size
  std::vector<std::vector<Vect>> GetMarchTriangles(
      const std::array<Scal, 8>& uu, const Vect& xc, const Vect& h) {
    std::array<double, 8> uuz = uu;
    for (auto& u : uuz) {
      u -= 0.5;
    }
    int n;
    std::array<double, 45> tri;
    march_cube(uuz.data(), &n, tri.data());
    assert(n * 3 * 3 <= tri.size());
    std::vector<std::vector<Vect>> vv;

    vv.resize(n);
    size_t i = 0;
    for (auto& v : vv) {
      v.resize(3);
      for (auto& x : v) {
        x[0] = tri[i++] - 0.5;
        x[1] = tri[i++] - 0.5;
        x[2] = tri[i++] - 0.5;
        x = xc + h * x;
      }
    }
    return vv;
  }
  void DumpPolyMarch() {
    auto sem = m.GetSem("dumppolymarch");
    if (sem("local")) {
      dl_.clear();
      dlc_.clear();
      dll_.clear();
      dlcl_.clear();
      auto h = m.GetCellSize();
      for (auto i : layers) {
        auto& fci = fci_[i];
        auto& fccl = fccl_[i];
        auto& fcu = fcu_[i].iter_curr;
        for (auto c : m.Cells()) {
          Scal u = fcu[c];
          const Scal th = par->poly_intth;
          if (fci[c] && u > th && u < 1. - th) {
            auto uu = GetStencil<1>(
                GetLayer(fcu_, Layers::iter_curr), c, fccl[c]);
            auto uun = ToNodes<1>(uu);
            auto vv = GetMarchTriangles(uun, m.GetCenter(c), h);
            for (auto& v : vv) {
              dl_.push_back(v);
              dlc_.push_back(m.GetHash(c));
              dll_.push_back(i);
              dlcl_.push_back(fccl[c]);
            }
          }
        }
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl_));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dlc_));
      m.Reduce(std::make_shared<TS>(&dll_));
      m.Reduce(std::make_shared<TS>(&dlcl_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::string fn = GetDumpName("sm", ".vtk", par->dmp->GetN());
        std::cout << std::fixed << std::setprecision(8)
            << "dump" 
            << " t=" << owner_->GetTime() + owner_->GetTimeStep()
            << " to " << fn << std::endl;
        WriteVtkPoly(fn, dl_, {&dlc_, &dll_, &dlcl_}, {"c", "l", "cl"}, 
            "Reconstructed linear interface", true, par->vtkbin);
      }
    }
  }
  // fcu: volume fraction [a]
  // fcud2: volume fraction difference double (xpp-xmm, ypp-ymm, zpp-zmm) [i]
  void CalcDiff2(const FieldCell<Scal>& fcu, FieldCell<Vect>& fcud2) {
    fcud2.Reinit(m);
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetNeighbourCell(c, 2 * d);
        auto cmm = m.GetNeighbourCell(cm, 2 * d);
        auto cp = m.GetNeighbourCell(c, 2 * d + 1);
        auto cpp = m.GetNeighbourCell(cp, 2 * d + 1);
        fcud2[c][d] = fcu[cpp] - fcu[cmm];
      }
    }
  }
  // fcu: volume fraction [a]
  // fcud2: volume fraction difference double (xpp-xmm, ...) [a]
  // Output:
  // fcud4: volume fraction difference quad (xp4-xm4, ...) [i]
  void CalcDiff4(const FieldCell<Scal>& fcu, 
      const FieldCell<Vect>& fcud2, FieldCell<Vect>& fcud4) {
    fcud4.Reinit(m);
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetNeighbourCell(c, 2 * d);
        auto cmm = m.GetNeighbourCell(cm, 2 * d);
        auto cp = m.GetNeighbourCell(c, 2 * d + 1);
        auto cpp = m.GetNeighbourCell(cp, 2 * d + 1);
        Scal um4 = fcu[c] - fcud2[cmm][d];
        Scal up4 = fcud2[cpp][d] + fcu[c];
        fcud4[c][d] = up4 - um4;
      }
    }
  }
  // fcu: volume fraction [a]
  // fcud4: volume fraction difference double (xp4-xm4, ...) [a]
  // Output:
  // fcud6: volume fraction difference quad (xp6-xm6, ...) [i]
  void CalcDiff6(const FieldCell<Scal>& fcu, 
      const FieldCell<Vect>& fcud4, FieldCell<Vect>& fcud6) {
    fcud6.Reinit(m);
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetNeighbourCell(c, 2 * d);
        auto cmm = m.GetNeighbourCell(cm, 2 * d);
        auto cp = m.GetNeighbourCell(c, 2 * d + 1);
        auto cpp = m.GetNeighbourCell(cp, 2 * d + 1);
        Scal um6 = fcu[cpp] - fcud4[cmm][d];
        Scal up6 = fcud4[cpp][d] + fcu[cmm];
        fcud6[c][d] = up6 - um6;
      }
    }
  }
  // u: volume fraction, array of size 3x3x3
  static Vect GetNormalYoungs(const std::array<Scal, 27>& u) {
    auto q = [&u](int dx, int dy, int dz) {
      const int w = 3;  // stencil width
      int i = (dx + 1) + (dy + 1) * w + (dz + 1) * w * w;
      return u[i];
    };

    Vect n;
    // generated by gen/normal.py
    n[0] = (1.0/32.0)*(-q(-1,-1,-1) - 2*q(-1,-1,0) -
    q(-1,-1,1) - 2*q(-1,0,-1) - 4*q(-1,0,0) - 2*q(-1,0,1) -
    q(-1,1,-1) - 2*q(-1,1,0) - q(-1,1,1) + q(1,-1,-1) +
    2*q(1,-1,0) + q(1,-1,1) + 2*q(1,0,-1) + 4*q(1,0,0) +
    2*q(1,0,1) + q(1,1,-1) + 2*q(1,1,0) + q(1,1,1));
    n[1] = (1.0/32.0)*(-q(-1,-1,-1) - 2*q(-1,-1,0) -
    q(-1,-1,1) + q(-1,1,-1) + 2*q(-1,1,0) + q(-1,1,1) -
    2*q(0,-1,-1) - 4*q(0,-1,0) - 2*q(0,-1,1) + 2*q(0,1,-1) +
    4*q(0,1,0) + 2*q(0,1,1) - q(1,-1,-1) - 2*q(1,-1,0) -
    q(1,-1,1) + q(1,1,-1) + 2*q(1,1,0) + q(1,1,1));
    n[2] = (1.0/32.0)*(-q(-1,-1,-1) + q(-1,-1,1) -
    2*q(-1,0,-1) + 2*q(-1,0,1) - q(-1,1,-1) + q(-1,1,1) -
    2*q(0,-1,-1) + 2*q(0,-1,1) - 4*q(0,0,-1) + 4*q(0,0,1) -
    2*q(0,1,-1) + 2*q(0,1,1) - q(1,-1,-1) + q(1,-1,1) -
    2*q(1,0,-1) + 2*q(1,0,1) - q(1,1,-1) + q(1,1,1));

    n /= -n.norm1();
    return n;
  }
  // Returns values over stencil centered at cell c with color cl.
  // Values for neighbors without color cl are filled with 0.
  // sw: stencil half-width
  template <size_t sw, size_t sn=sw*2+1>
  std::array<Scal, sn*sn*sn> GetStencil(
      const Multi<const FieldCell<Scal>*>& fc, IdxCell c, Scal cl) {
    using MIdx = typename M::MIdx;
    auto& bc = m.GetIndexCells();
    GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sn));
    MIdx w = bc.GetMIdx(c);
    std::array<Scal, sn*sn*sn> uu;
    size_t k = 0;
    for (MIdx wo : bo) {
      IdxCell cn = bc.GetIdx(w + wo);
      Scal u = 0;
      for (auto j : layers) {
        if (fccl_[j][cn] == cl) {
          u = (*fc[j])[cn];
          break;
        }
      }
      uu[k++] = u;
    }
    return uu;
  }
  // Interpolates from cells to nodes.
  // stencil half-width
  template <int sw, int sn=sw*2+1, int snn=sw*2>
  std::array<Scal, snn*snn*snn> ToNodes(const std::array<Scal, sn*sn*sn>& uu) {
    std::array<Scal, snn*snn*snn> uun;
    size_t i = 0;
    for (int z = 0; z < snn; ++z) {
      for (int y = 0; y < snn; ++y) {
        for (int x = 0; x < snn; ++x) {
          auto u = [&uu,x,y,z](int dx, int dy, int dz) {
            return uu[(z+dz)*sn*sn + (y+dy)*sn + (x+dx)];
          };
          uun[i++] = (1. / 8.) * (
              u(0,0,0) + u(1,0,0) + u(0,1,0) + u(1,1,0) +
              u(0,0,1) + u(1,0,1) + u(0,1,1) + u(1,1,1));
        }
      }
    }
    return uun;
  }
  // reconstruct interface
  void Rec(const Multi<FieldCell<Scal>*>& uc) {
    auto sem = m.GetSem("rec");
    if (sem("detect")) {
      DetectInterface(uc);
    }
    if (sem("local")) {
      for (auto i : layers) {
        auto& fcn = fcn_[i];
        auto& fci = fci_[i];
        auto& fcu = *uc[i];
        auto& fccl = fccl_[i];
        if (par->bcc_reflect) {
          BcReflect(fcu, mfc_, par->bcc_fill, m);
        }

        const int sw = 1; // stencil halfwidth
        using MIdx = typename M::MIdx;
        auto& bc = m.GetIndexCells();
        GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sw * 2 + 1));
        for (auto c : m.Cells()) {
          if (fci[c]) {
            MIdx w = bc.GetMIdx(c);
            std::array<Scal, 27> uu;
            size_t k = 0;
            for (MIdx wo : bo) {
              IdxCell cn = bc.GetIdx(w + wo);
              Scal u = 0;
              for (auto j : layers) {
                if (fccl_[j][cn] == fccl[c]) {
                  u = (*uc[j])[cn];
                  break;
                }
              }
              uu[k++] = u;
            }
            fcn[c] = GetNormalYoungs(uu);
          }
        }
      }
      for (auto i : layers) {
        if (par->bcc_reflect) {
          auto& fcn = fcn_[i];
          BcReflect(fcn, mfc_, Vect(0), m);
        }
      }
      for (auto i : layers) {
        auto& fcn = fcn_[i];
        auto& fca = fca_[i];
        auto& fci = fci_[i];
        auto& fcu = *uc[i];
        for (auto c : m.SuCells()) {
          if (fci[c]) {
            fca[c] = R::GetLineA(fcn[c], fcu[c], m.GetCellSize());
          }
        }
        if (par->bcc_reflect) {
          BcReflect(fca, mfc_, Scal(0), m);
        }
      }
    }
    for (auto i : layers) {
      if (sem("comm")) {
        auto& fcn = fcn_[i];
        auto& fca = fca_[i];
        m.Comm(&fca);
        m.Comm(&fcn);
      }
    }
  }
  void DetectInterface(const Multi<const FieldCell<Scal>*>& uc) {
    for (auto i : layers) {
      auto& fci = fci_[i];
      auto& fcu = *uc[i];
      fci.Reinit(m, false);
      for (auto c : m.AllCells()) {
        Scal u = fcu[c];
        if (u > 0. && u < 1.) {
          fci[c] = true;
        }
      }
      /*
      // neighbour cell has different value but both are 0 or 1
      for (auto f : m.SuFaces()) {
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        if (fcm[cm] && fcm[cp]) {
          Scal um = fcu[cm];
          Scal up = fcu[cp];
          if ((um == 0. || um == 1.) && (up == 0. || up == 1.) && (um != up)) {
            fci[cm] = true;
            fci[cp] = true;
          }
        }
      }
      */
    }
  }
  void StartStep() {
    auto sem = m.GetSem("start");
    if (owner_->GetTime() == 0.) {
      if (sem.Nested("reconst")) {
        Rec(GetLayer(fcu_, Layers::time_curr));
      }
    }
    if (sem("rotate")) {
      owner_->ClearIter();
      for (auto& u : fcu_.data()) {
        u.time_prev = u.time_curr;
        u.iter_curr = u.time_prev;
      }
      fcus_.time_prev = fcus_.time_curr;
      fcus_.iter_curr = fcus_.time_prev;
    }
  }
  void Dump() {
    auto sem = m.GetSem("iter");
    bool dm = par->dmp->Try(owner_->GetTime(),
                            owner_->GetTimeStep());
    if (par->dumppoly && dm && sem.Nested()) {
      DumpPoly();
    }
    if (par->dumppolymarch && dm && sem.Nested()) {
      DumpPolyMarch();
    }
    if (par->dumppart && dm && sem.Nested("part-dump")) {
      psm_->DumpParticles(fca_, fcn_, par->dmp->GetN(), owner_->GetTime());
    }
    if (par->dumppartinter && dm && sem.Nested("partinter-dump")) {
      psm_->DumpPartInter(fca_, fcn_, par->dmp->GetN(), owner_->GetTime());
    }
  }
  void Recolor() {
    auto sem = m.GetSem("recolor");
    if (sem("init")) {
      fcclt_.InitAll(FieldCell<Scal>(m, kClNone));
      // initial unique color
      Scal q = m.GetId() * m.GetInBlockCells().size() * layers.size();
      for (auto i : layers) {
        for (auto c : m.Cells()) {
          if (fccl_[i][c] != kClNone) {
            fcclt_[i][c] = (q += 1);
          }
        }
        m.Comm(&fcclt_[i]);
      }
    }
    sem.LoopBegin();
    if (sem("min")) {
      const int sw = 1; // stencil halfwidth
      using MIdx = typename M::MIdx;
      auto& bc = m.GetIndexCells();
      GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sw * 2 + 1));
      size_t tries = 0;
      size_t cells = 0;
      while (true) {
        bool chg = false;
        for (auto i : layers) {
          for (auto c : m.Cells()) {
            if (fccl_[i][c] != kClNone) {
              MIdx w = bc.GetMIdx(c);
              for (MIdx wo : bo) {
                IdxCell cn = bc.GetIdx(w + wo);
                for (auto j : layers) {
                  if (fccl_[j][cn] == fccl_[i][c]) {
                    if (fcclt_[j][cn] < fcclt_[i][c]) {
                      chg = true;
                      ++cells;
                      fcclt_[i][c] = fcclt_[j][cn];
                    }
                  }
                }
              }
            }
          }
        }
        if (!chg) {
          break;
        }
        ++tries;
      }
      for (auto i : layers) {
        m.Comm(&fcclt_[i]);
      }
      recolor_tries_ = tries;
      m.Reduce(&recolor_tries_, "max");
    }
    if (sem("check")) {
      if (par->verb && m.IsRoot()) {
        std::cerr << "recolor:"
          << " max tries: " << recolor_tries_ 
          << std::endl;
      }
      if (!recolor_tries_) {
        sem.LoopBreak();
      }
    }
    sem.LoopEnd();
    if (sem("free")) {
      for (auto i : layers) {
        for (auto c : m.AllCells()) {
          fccl_[i][c] = fcclt_[i][c];
        }
        fcclt_[i].Free();
      }
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("iter");
    if (sem("init")) {
      for (auto i : layers) {
        auto& fcu = fcu_[i];
        auto& uc = fcu.iter_curr;
        const Scal dt = owner_->GetTimeStep();
        auto& fcs = *owner_->fcs_;
        for (auto c : m.Cells()) {
          uc[c] =
              fcu.time_prev[c] +  // previous time step
              dt * fcs[c]; // source
        }
      }
    }

    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    // directions, format: {dir LE, dir EI, ...}
    std::vector<size_t> dd; 
    Scal vsc; // scaling factor for ffv, used for splitting
    if (par->dim == 3) { // 3d
      if (count_ % 3 == 0) {
        dd = {0, 1, 1, 2, 2, 0};
      } else if (count_ % 3 == 1) {
        dd = {1, 2, 2, 0, 0, 1};
      } else {
        dd = {2, 0, 0, 1, 1, 2};
      }
      vsc = 0.5;
    } else {
      if (count_ % 2 == 0) {
        dd = {0, 1};
      } else {
        dd = {1, 0};
      } 
      vsc = 1.0;
    }
    for (size_t id = 0; id < dd.size(); ++id) {
      size_t d = dd[id]; // direction as index
      if (sem("copyface")) {
        if (id % 2 == 1) { // copy fluxes for Lagrange Explicit step
          auto& ffv = *owner_->ffv_; // [f]ield [f]ace [v]olume flux
          fcfm_.Reinit(m);
          fcfp_.Reinit(m);
          for (auto c : m.Cells()) {
            fcfm_[c] = ffv[m.GetNeighbourFace(c, 2 * d)];
            fcfp_[c] = ffv[m.GetNeighbourFace(c, 2 * d + 1)];
          }
          m.Comm(&fcfm_);
          m.Comm(&fcfp_);
        }
        for (auto i : layers) {
          ffvu_[i].Reinit(m, 0);
          ffcl_[i].Reinit(m, kClNone);
          ffi_[i].Reinit(m, false);
        }
      }
      for (auto i : layers) {
        if (sem("flux")) {
          Dir md(d); // direction as Dir
          MIdx wd(md); // offset in direction d
          auto& bf = m.GetIndexFaces();
          auto h = m.GetCellSize();
          const Scal dt = owner_->GetTimeStep();
          auto& ffv = *owner_->ffv_; // [f]ield [f]ace [v]olume flux

          // Returns phase 2 flux
          auto F = [this,vsc,id,h,dt,d,&ffv](
              IdxFace f, IdxCell cu, Scal u, Vect n, Scal a, bool in) {
            const Scal v = ffv[f] * vsc; // mixture flux
            if (in) {
              if (id % 2 == 0) { // Euler Implicit
                return R::GetLineFlux(n, a, h, v, dt, d);
              } else {           // Lagrange Explicit
                // upwind mixture flux
                Scal vu = (v > 0. ? fcfm_[cu] : fcfp_[cu]) * vsc;
                return R::GetLineFluxStr(n, a, h, v, vu, dt, d);
              }
            }
            return v * u;
          };

          auto& fcu = fcu_[i].iter_curr;
          auto& fcn = fcn_[i];
          auto& fca = fca_[i];
          auto& fci = fci_[i];
          auto& fccl = fccl_[i];
          auto& ffvu = ffvu_[i];
          auto& ffcl = ffcl_[i];
          auto& ffi = ffi_[i];

          for (auto f : m.Faces()) {
            if (bf.GetDir(f) != md) {
              continue;
            }

            // upwind cell
            IdxCell cu = m.GetNeighbourCell(f, ffv[f] > 0. ? 0 : 1);
            if (fccl[cu] != kClNone) {
              ffvu[f] = F(f, cu, fcu[cu], fcn[cu], fca[cu], fci[cu]);
              ffi[f] = fci[cu];
              ffcl[f] = fccl[cu];
            }
          }

          FieldFace<Scal> ffu(m);
          // interpolate field value to boundaries
          InterpolateB(fcu, mfc_, ffu, m);
          // override boundary upwind flux
          for (const auto& it : mfc_) {
            IdxFace f = it.GetIdx();
            CondFace* cb = it.GetValue().get(); 
            Scal v = ffv[f];
            if ((cb->GetNci() == 0) != (v > 0.)) {
              ffvu[f] = v * ffu[f];
              ffi[f] = true;
            }
          }
        }
      }
      if (sem("cell")) {
        for (auto i : layers) {
          Dir md(d); // direction as Dir
          MIdx wd(md); // offset in direction d
          const Scal dt = owner_->GetTimeStep();
          auto& bc = m.GetIndexCells();
          auto& bf = m.GetIndexFaces();
          auto& ffv = *owner_->ffv_; // mixture flux
          auto& fcu = fcu_[i].iter_curr;
          auto& fccl = fccl_[i];

          for (auto c : m.Cells()) {
            auto w = bc.GetMIdx(c);
            const Scal lc = m.GetVolume(c);
            Scal cl = fccl[c];
            if (cl == kClNone) {
              continue;
            }
            // faces
            IdxFace fm = bf.GetIdx(w, md);
            IdxFace fp = bf.GetIdx(w + wd, md);
            // mixture fluxes
            const Scal vm = ffv[fm] * vsc;
            const Scal vp = ffv[fp] * vsc;
            // mixture cfl
            const Scal sm = vm * dt / lc;
            const Scal sp = vp * dt / lc;
            const Scal ds = sp - sm;
            // phase 2 fluxes
            Scal qm = 0;
            Scal qp = 0;
            // interface
            bool im = false;
            bool ip = false;
            for (auto j : layers) {
              if (ffcl_[j][fm] == cl) {
                qm = ffvu_[j][fm];
                im = ffi_[j][fm];
                break;
              }
            }
            for (auto j : layers) {
              if (ffcl_[j][fp] == cl) {
                qp = ffvu_[j][fp];
                ip = ffi_[j][fp];
                break;
              }
            }
            if (!im && !ip) {
              //continue; // XXX
            }
            // phase 2 cfl
            const Scal lm = qm * dt / lc;
            const Scal lp = qp * dt / lc;
            const Scal dl = lp - lm;
            auto& u = fcu[c];
            if (id % 2 == 0) { // Euler Implicit
              u = (u - dl) / (1. - ds);
            } else { // Lagrange Explicit
              u = u * (1. + ds) - dl;
            }
            // clip
            const Scal th = par->clipth;
            if (u < th) {
              u = 0.;
            } else if (u > 1. - th) {
              u = 1.;
            } else if (IsNan(u)) {
              u = 0.;
            }
            // update color
            if (u == 0) {
              fccl[c] = kClNone;
            }
          }
        }

        for (auto i : layers) {
          Dir md(d); // direction as Dir
          MIdx wd(md); // offset in direction d
          const Scal dt = owner_->GetTimeStep();
          auto& bc = m.GetIndexCells();
          auto& bf = m.GetIndexFaces();
          auto& ffv = *owner_->ffv_; // mixture flux
          auto& ffcl = ffcl_[i];
          for (auto f : m.Faces()) {
            auto p = bf.GetMIdxDir(f);
            MIdx w = p.first;
            Dir df = p.second;
            if (df != md) {
              continue;
            }
            Scal cl = ffcl[f];
            if (cl == kClNone) {
              continue;
            }

            const Scal v = ffv[f] * vsc;

            IdxFace fm, fp;
            IdxCell c;
            if (v > 0) {
              fm = f;
              fp = bf.GetIdx(w + wd, md);
              c = bc.GetIdx(w);
            } else {
              fm = bf.GetIdx(w - wd, md);
              fp = f;
              c = bc.GetIdx(w - wd);
            }

            const Scal lc = m.GetVolume(c);

            // mixture fluxes
            const Scal vm = ffv[fm] * vsc;
            const Scal vp = ffv[fp] * vsc;
            // mixture cfl
            const Scal sm = vm * dt / lc;
            const Scal sp = vp * dt / lc;
            const Scal ds = sp - sm;
            // phase 2 fluxes
            Scal qm = 0;
            Scal qp = 0;
            // interface
            bool im = false;
            bool ip = false;
            // qm,im from layer with same color
            for (auto j : layers) {
              if (ffcl_[j][fm] == cl) {
                qm = ffvu_[j][fm];
                im = ffi_[j][fm];
                break;
              }
            }
            // qp,ip from layer with same color
            for (auto j : layers) {
              if (ffcl_[j][fp] == cl) {
                qp = ffvu_[j][fp];
                ip = ffi_[j][fp];
                break;
              }
            }
            if (!im && !ip) {
              //continue; // XXX
            }
            // phase 2 cfl
            const Scal lm = qm * dt / lc;
            const Scal lp = qp * dt / lc;
            const Scal dl = lp - lm;

            // find cell with same color
            const size_t jnone = -1;
            size_t j = jnone;
            for (auto jj : layers) {
              if (fccl_[jj][c] == cl) {
                j = jj;
                break;
              }
            }
            if (j == jnone) { // if not found, find first empty
              for (auto jj : layers) {
                if (fccl_[jj][c] == kClNone) {
                  j = jj;
                  break;
                }
              }
            }
            if (j == jnone) {
              //std::cerr << "warn: no layer i=" << i << " w=" << w << "\n";
              j = 0;
            }

            if (fccl_[j][c] != kClNone) {
              continue;
            }

            auto& u = fcu_[j].iter_curr[c];
            if (id % 2 == 0) { // Euler Implicit
              u = (u - dl) / (1. - ds);
            } else { // Lagrange Explicit
              u = u * (1. + ds) - dl;
            }

            // clip
            const Scal th = par->clipth;
            if (u < th) {
              u = 0.;
            } else if (u > 1. - th) {
              u = 1.;
            } else if (IsNan(u)) {
              u = 0.;
            }
            // update color
            auto& clc = fccl_[j][c];
            if (u == 0) {
              clc = kClNone;
            } else {
              clc = cl;
            }
          }
        }
        for (auto i : layers) {
          m.Comm(&fcu_[i].iter_curr);
          m.Comm(&fccl_[i]);
        }
      }
      if (sem.Nested("reconst")) {
        Rec(GetLayer(fcu_, Layers::iter_curr));
      }
    }
    if (sem.Nested("recolor")) {
      Recolor();
    }
    if (sem("stat")) {
      owner_->IncIter();
      ++count_;
      auto l = Layers::iter_curr;
      auto& fcus = fcus_.Get(l);
      fcus.Reinit(m, 0);
      fccls_.Reinit(m, kClNone);
      for (auto i : layers) {
        auto& fccl = fccl_[i];
        auto& fcu = fcu_[i].Get(l);
        for (auto c : m.AllCells()) {
          if (fccl[c] != kClNone) {
            fcus[c] += fcu[c];
            fccls_[c] = fccl[c];
          }
        }
      }
      for (auto c : m.AllCells()) {
        auto& u = fcus[c];
        if (!(u >= 0)) {
          u = 0;
        } else if (!(u <= 1.)) {
          u = 1;
        } 
      }
    }
  }
  void FinishStep() {
    for (auto& u : fcu_.data()) {
      u.time_curr = u.iter_curr;
    }
    fcus_.time_curr = fcus_.iter_curr;
    owner_->IncTime();
  }
  void PostStep() {
    auto sem = m.GetSem("iter");
    // Curvature from gradient of volume fraction
    if (par->curvgrad && sem("curv")) {
      for (auto i : layers) {
        auto ffu = Interpolate(fcu_[i].iter_curr, mfc_, m); // [s]
        auto fcg = Gradient(ffu, m); // [s]
        auto ffg = Interpolate(fcg, mfvz_, m); // [i]

        auto& fck = fck_[i];
        fck.Reinit(m, GetNan<Scal>()); // curvature [i]

        auto& fci = fci_[i];
        for (auto c : m.Cells()) {
          if (!fci[c]) {
            continue;
          }
          Scal s = 0.;
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetNeighbourFace(c, q);
            auto& g = ffg[f];
            auto n = g / g.norm();  // inner normal
            s += -n.dot(m.GetOutwardSurface(c, q));
          }
          fck[c] = s / m.GetVolume(c);
        }
      }
    }
    if (!par->curvgrad) {
      for (auto i : layers) {
        auto& uc = fcu_[i].iter_curr;
        if (sem("diff2")) {
          CalcDiff2(uc, fcud2_);
          m.Comm(&fcud2_);
        }
        if (sem("diff4")) {
          CalcDiff4(uc, fcud2_, fcud4_);
          m.Comm(&fcud4_);
        }
        if (sem("height")) {
          UNormal<M>::CalcHeight(m, uc, fcud2_, fcud4_, par->dim, fch_);
          m.Comm(&fch_);
        }
        if (sem("curvcomm")) {
          UNormal<M>::CalcCurvHeight(m, uc, fch_, fcn_[i], par->dim, fck_[i]);
          m.Comm(&fck_[i]);
        }
      }
    }
    if (par->part && sem.Nested("part")) {
      psm_->Part(GetLayer(fcu_, Layers::iter_curr),
                 fca_, fcn_, fci_, fccl_, mfc_);
    }
    if (sem.Nested("dump")) {
      Dump();
    }
  }


  Owner* owner_;
  std::shared_ptr<Par> par;
  M& m;

  Multi<LayersData<FieldCell<Scal>>> fcu_;
  LayersData<FieldCell<Scal>> fcus_;
  FieldCell<Scal> fccls_;
  MapFace<std::shared_ptr<CondFace>> mfc_;
  MapFace<std::shared_ptr<CondFace>> mfvz_; // zero-derivative bc for Vect

  Multi<FieldCell<Scal>> fca_; // alpha (plane constant)
  Multi<FieldCell<Vect>> fcn_; // n (normal to plane)
  Multi<FieldCell<Scal>> fck_; // curvature from height functions
  Multi<FieldCell<bool>> fci_; // interface mask (1: contains interface)
  Multi<FieldCell<bool>> fcm_; // layer mask
  Multi<FieldCell<bool>> fcdp_; // dependent cells with mask=0
  Multi<FieldCell<Scal>> fccl_;  // color
  Multi<FieldCell<Scal>> fcclt_;  // tmp color
  Multi<FieldFace<Scal>> ffvu_;  // flux: volume flux * field
  Multi<FieldFace<Scal>> ffcl_;  // flux color (from upwind cell)
  Multi<FieldFace<bool>> ffi_;   // interface mask (1: upwind cell contains interface)
  FieldCell<Vect> fcud2_; // volume fraction difference double
  FieldCell<Vect> fcud4_; // volume fraction difference quad
  FieldCell<Vect> fcud6_; // volume fraction difference 6th
  FieldCell<Vect> fch_; // curvature from height functions
  size_t count_ = 0; // number of MakeIter() calls, used for splitting

  std::vector<std::vector<Vect>> dl_; // dump poly
  std::vector<Scal> dlc_; // dump poly cell
  std::vector<Scal> dll_; // dump poly layer
  std::vector<Scal> dlcl_; // dump poly color

  std::unique_ptr<PSM> psm_;
  // tmp for MakeIteration, volume flux copied to cells
  FieldCell<Scal> fcfm_, fcfp_;
  GRange<size_t> layers;
  Scal recolor_tries_;
};

template <class M>
constexpr typename M::Scal Vof<M>::Imp::kClNone;

template <class M_>
Vof<M_>::Vof(
    M& m, const FieldCell<Scal>& fcu, const FieldCell<Scal>& fccl,
    const MapFace<std::shared_ptr<CondFace>>& mfc,
    const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs,
    double t, double dt, std::shared_ptr<Par> par)
    : AdvectionSolver<M>(t, dt, m, ffv, fcs)
    , imp(new Imp(this, fcu, fccl, mfc, par))
{}

template <class M_>
Vof<M_>::~Vof() = default;

template <class M_>
auto Vof<M_>::GetPar() -> Par* {
  return imp->par.get();
}

template <class M_>
void Vof<M_>::StartStep() {
  imp->StartStep();
}

template <class M_>
void Vof<M_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_>
void Vof<M_>::FinishStep() {
  imp->FinishStep();
}

template <class M_>
auto Vof<M_>::GetField(Layers l) const -> const FieldCell<Scal>& {
  return imp->fcus_.Get(l);
}

template <class M_>
auto Vof<M_>::GetField(Layers l, size_t i) const -> const FieldCell<Scal>& {
  return imp->fcu_[i].Get(l);
}

template <class M_>
auto Vof<M_>::GetField(size_t i) const -> const FieldCell<Scal>& {
  return imp->fcu_[i].Get(Layers::time_curr);
}

template <class M_>
auto Vof<M_>::GetAlpha(size_t i) const -> const FieldCell<Scal>& {
  (void) i;
  return imp->fca_[i];
}

template <class M_>
auto Vof<M_>::GetMask(size_t i) const -> const FieldCell<bool>& {
  return imp->fcm_[i];
}

template <class M_>
auto Vof<M_>::GetColor(size_t i) const -> const FieldCell<Scal>& {
  return imp->fccl_[i];
}

template <class M_>
auto Vof<M_>::GetColor() const -> const FieldCell<Scal>& {
  return imp->fccls_;
}

template <class M_>
auto Vof<M_>::GetDepend(size_t i) const -> const FieldCell<bool>& {
  return imp->fcdp_[i];
}

template <class M_>
size_t Vof<M_>::GetNumLayers() const {
  return imp->layers.size();
}

template <class M_>
auto Vof<M_>::GetNormal(size_t i) const -> const FieldCell<Vect>& {
  return imp->fcn_[i];
}

template <class M_>
auto Vof<M_>::GetHeight() const -> const FieldCell<Vect>& {
  return imp->fch_;
}

template <class M_>
auto Vof<M_>::GetCurv() const -> const FieldCell<Scal>& {
  return imp->par->part_k ? imp->psm_->GetCurv(0) : imp->fck_[0];
}

template <class M_>
auto Vof<M_>::GetCurv(size_t i) const -> const FieldCell<Scal>& {
  (void) i;
  return imp->par->part_k ? imp->psm_->GetCurv(i) : imp->fck_[i];
}

template <class M_>
void Vof<M_>::PostStep() {
  return imp->PostStep();
}


// curvature from height function
template <class M_>
auto Vof<M_>::GetCurvH() const -> const FieldCell<Scal>& {
  return imp->fck_[0];
}

// curvature from particles
template <class M_>
auto Vof<M_>::GetCurvP() const -> const FieldCell<Scal>& {
  return imp->psm_->GetCurv(0);
}

} // namespace solver
