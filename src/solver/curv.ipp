#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>

#include "approx.h"
#include "debug/isnan.h"
#include "geom/block.h"
#include "normal.h"
#include "partstrmeshm.h"
#include "reconst.h"
#include "trackerm.h"
#include "util/vof.h"
#include "vof.h"

#include "curv.h"

using namespace solver;

template <class M_>
struct UCurv<M_>::Imp {
  using R = Reconst<Scal>;
  //using PS = PartStr<Scal>;
  //using PSM = PartStrMeshM<M>;
  static constexpr size_t dim = M::dim;
  using Vect2 = GVect<Scal, 2>;

  // fcu: volume fraction [a]
  // fcud2: volume fraction difference double (xpp-xmm, ypp-ymm, zpp-zmm) [i]
  static void CalcDiff2(
      const FieldCell<Scal>& fcu, FieldCell<Vect>& fcud2, const M& m) {
    fcud2.Reinit(m);
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetCell(c, 2 * d);
        auto cmm = m.GetCell(cm, 2 * d);
        auto cp = m.GetCell(c, 2 * d + 1);
        auto cpp = m.GetCell(cp, 2 * d + 1);
        fcud2[c][d] = fcu[cpp] - fcu[cmm];
      }
    }
  }
  // fcu: volume fraction [a]
  // fcud2: volume fraction difference double (xpp-xmm, ...) [a]
  // Output:
  // fcud4: volume fraction difference quad (xp4-xm4, ...) [i]
  static void CalcDiff4(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcud2,
      FieldCell<Vect>& fcud4, const M& m) {
    fcud4.Reinit(m);
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetCell(c, 2 * d);
        auto cmm = m.GetCell(cm, 2 * d);
        auto cp = m.GetCell(c, 2 * d + 1);
        auto cpp = m.GetCell(cp, 2 * d + 1);
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
  static void CalcDiff6(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcud4,
      FieldCell<Vect>& fcud6, const M& m) {
    fcud6.Reinit(m);
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetCell(c, 2 * d);
        auto cmm = m.GetCell(cm, 2 * d);
        auto cp = m.GetCell(c, 2 * d + 1);
        auto cpp = m.GetCell(cp, 2 * d + 1);
        Scal um6 = fcu[cpp] - fcud4[cmm][d];
        Scal up6 = fcud4[cpp][d] + fcu[cmm];
        fcud6[c][d] = up6 - um6;
      }
    }
  }
  static void CalcCurvHeight(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn, size_t edim,
      FieldCell<Scal>& fck, M& m) {
    auto sem = m.GetSem();
    struct {
      FieldCell<Vect> fcud2; // volume fraction difference double
      FieldCell<Vect> fcud4; // volume fraction difference quad
      FieldCell<Vect> fch; // height functions
    } * ctx(sem);
    auto& fcud2 = ctx->fcud2;
    auto& fcud4 = ctx->fcud4;
    auto& fch = ctx->fch;

    if (sem("diff2")) {
      CalcDiff2(fcu, fcud2, m);
      m.Comm(&fcud2);
    }
    if (sem("diff4")) {
      CalcDiff4(fcu, fcud2, fcud4, m);
      m.Comm(&fcud4);
    }
    if (sem("height")) {
      UNormal<M>::CalcHeight(m, fcu, fcud2, fcud4, edim, fch);
      m.Comm(&fch);
    }
    if (sem("curvcomm")) {
      UNormal<M>::CalcCurvHeight(m, fcu, fch, fcn, edim, fck);
      m.Comm(&fck);
    }
  }
};

template <class M>
void UCurv<M>::CalcCurvHeight(
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn, size_t edim,
    FieldCell<Scal>& fck, M& m) {
  Imp::CalcCurvHeight(fcu, fcn, edim, fck, m);
}
