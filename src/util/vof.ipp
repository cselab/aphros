#pragma once

#include <vector>
#include <memory>
#include <limits>
#include <iomanip>

#include "vof.h"
#include "geom/mesh.h"
#include "dump/vtk.h"
#include "solver/reconst.h"
#include "debug/isnan.h"

namespace solver {

template <class M_>
struct UVof<M_>::Imp {
  using R = Reconst<Scal>;
  static constexpr size_t dim = M::dim;

  Imp() = default;
  ~Imp() = default;

  void DumpPoly(
      const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fca,
      const Multi<const FieldCell<bool>*>& fci,
      std::string fn, Scal t, Scal th, bool bin, bool merge, M& m) {
    auto sem = m.GetSem("dumppoly");
    if (sem("local")) {
      dl_.clear();
      dlc_.clear();
      dll_.clear();
      dlcl_.clear();
      auto h = m.GetCellSize();
      for (auto i : layers) {
        for (auto c : m.Cells()) {
          Scal u = (*fcu[i])[c];
          if (IsNan(u) || IsNan((*fcn[i])[c]) || IsNan((*fca[i])[c])) {
            continue;
          }
          if ((*fci[i])[c] && u > th && u < 1 - th) {
            dl_.push_back(R::GetCutPoly(
                  m.GetCenter(c), (*fcn[i])[c], (*fca[i])[c], h));
            dlc_.push_back(m.GetHash(c));
            dll_.push_back(i);
            dlcl_.push_back(fccl[i] ? (*fccl[i])[c] : 0);
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
        std::cout << std::fixed << std::setprecision(8)
            << "dump" << " t=" << t << " to " << fn << std::endl;
        WriteVtkPoly(fn, dl_, {&dlc_, &dll_, &dlcl_}, {"c", "l", "cl"},
            "Interface from PLIC", true,
            bin, merge);
      }
    }
  }

  std::vector<std::vector<Vect>> dl_; // dump poly
  std::vector<Scal> dlc_; // dump poly cell
  std::vector<Scal> dll_; // dump poly layer
  std::vector<Scal> dlcl_; // dump poly color
};

template <class M_>
UVof<M_>::UVof() : imp(new Imp()) {}

template <class M_>
UVof<M_>::~UVof() = default;

template <class M_>
void UVof<M_>::DumpPoly(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<Scal>*>& fcu,
    const Multi<const FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<Vect>*>& fcn,
    const Multi<const FieldCell<Scal>*>& fca,
    const Multi<const FieldCell<bool>*>& fci,
    std::string fn, Scal t, Scal th, bool bin, bool merge, M& m) {
  imp->DumpPoly(layers, fcu, fccl, fcn, fca, fci, fn, t, th, bin, merge, m);
}

template <class M_>
void UVof<M_>::DumpPoly(
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn,
    const FieldCell<Scal>& fca, const FieldCell<bool>& fci,
    std::string fn, Scal t, Scal th, bool bin, bool merge, M& m) {
  GRange<size_t> layers(0, 1);
  const FieldCell<Scal>* fccl(nullptr);
  imp->DumpPoly(layers, &fcu, fccl, &fcn, &fca, &fci,
                fn, t, th, bin, merge, m);
}

} // namespace solver
