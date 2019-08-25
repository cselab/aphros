#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>

#include "vof.h"
#include "geom/block.h"
#include "dump/vtk.h"
#include "reconst.h"
#include "normal.h"
#include "debug/isnan.h"
#include "partstr.h"

namespace solver {

template <class M_>
struct MultiMask<M_>::Imp {
  using Owner = MultiMask<M_>;
  static constexpr size_t dim = M::dim;

  Imp(Owner* owner, M& m, const FieldCell<Scal>* fccl, size_t edim)
      : owner_(owner)
      , m(m), fccl_(fccl), fcmask_(m, 0), fcmask2_(m, 0), edim_(edim)
  {}
  void Update(const FieldCell<Scal>& fcu) {
    using MIdx = typename M::MIdx;
    auto& bc = m.GetIndexCells();

    auto sem = m.GetSem("upd");

    if (sem("")) {
      const int sw = 2; // stencil halfwidth, [sw,sw]
      const int sn = sw * 2 + 1; // stencil size
      // block of offsets
      GBlock<IdxCell, dim> bo(MIdx(-sw, -sw, edim_ == 2 ? 0 : -sw), 
                              MIdx(sn, sn, edim_ == 2 ? 1 : sn)); 
      (void) bo;
      auto& cl = *fccl_;
      fcmask_.Reinit(m, kNone);
      fcmask2_.Reinit(m, kNone);
      auto I = [&fcu](IdxCell c) {
        return fcu[c] > 0 and fcu[c] < 1;
      };
      for (auto c : m.Cells()) {
        MIdx w = bc.GetMIdx(c);
        Scal max = kNone;
        Scal min = kNone;
        for (MIdx wo : bo) {
          IdxCell cn = bc.GetIdx(w + wo); 
          if (cl[cn] != kNone && I(cn)) {
            if (max == kNone || cl[cn] > max) {
              max = cl[cn];
            }
            if (min == kNone || cl[cn] < min) {
              min = cl[cn];
            }
          }
        }
        if (min != kNone) {
          fcmask_[c] = min;
        }
        if (max != min) {
          fcmask2_[c] = max;
        }
      }
      if(0)
      for (auto c : m.Cells()) {
        if (I(c) && cl[c] == fcmask2_[c] && fcmask_[c] != kNone &&
            fcmask_[c] != cl[c]) {
          MIdx w = bc.GetMIdx(c);
          for (MIdx wo : bo) {
            IdxCell cn = bc.GetIdx(w + wo); 
            fcmask2_[cn] = fcmask2_[c];
          }
        }
      }
    }
  }

  Owner* owner_;
  M& m;
  const FieldCell<Scal>* fccl_; // color
  FieldCell<Scal> fcmask_; // mask
  FieldCell<Scal> fcmask2_; // mask
  size_t edim_;
};

template <class M_>
MultiMask<M_>::MultiMask(M& m, const FieldCell<Scal>* fccl, size_t edim)
    : imp(new Imp(this, m, fccl, edim))
{}

template <class M_>
MultiMask<M_>::~MultiMask() = default;

template <class M_>
void MultiMask<M_>::Update(const FieldCell<Scal>& fcu) {
  imp->Update(fcu);
}

template <class M_>
auto MultiMask<M_>::GetMask() const -> const FieldCell<Scal>& { 
  return imp->fcmask_; 
}

template <class M_>
auto MultiMask<M_>::GetMask2() const -> const FieldCell<Scal>& { 
  return imp->fcmask2_; 
}

} // namespace solver
