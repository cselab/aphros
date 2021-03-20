// Created by Petr Karnakov on 28.02.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <iostream>
#include <tuple>
#include <vector>

#include "geom/field.h"
#include "geom/mesh.h"

template <class M_>
struct UDebug {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Expr = typename M::Expr;

  static std::vector<generic::Vect<Scal, 3>> GetNorms(
      const std::vector<const FieldCell<Scal>*>& vfc, M& m) {
    auto sem = m.GetSem(__func__);
    using V = generic::Vect<Scal, 3>;
    struct {
      std::vector<V> vr;
      Scal nc;
    } * ctx(sem);
    auto& vr = ctx->vr;
    auto& nc = ctx->nc;
    if (sem("local")) {
      vr.resize(vfc.size(), V(0));
      for (size_t i = 0; i < vfc.size(); ++i) {
        auto& fc = *vfc[i];
        auto& r = vr[i];
        auto& norm1 = r[0];
        auto& norm2 = r[1];
        auto& norminf = r[2];
        for (auto c : m.Cells()) {
          auto u = std::abs(fc[c]);
          norm1 += u;
          norm2 += u * u;
          norminf = std::max(norminf, u);
        }
        m.Reduce(&norm1, "sum");
        m.Reduce(&norm2, "sum");
        m.Reduce(&norminf, "max");
      }
      nc = 0;
      for (auto c : m.Cells()) {
        (void)c;
        nc += 1;
      }
      m.Reduce(&nc, "sum");
    }
    if (sem("reduce")) {
      for (auto& r : vr) {
        auto& norm1 = r[0];
        auto& norm2 = r[1];
        norm1 = norm1 / nc;
        norm2 = std::sqrt(norm2 / nc);
      }
      return vr;
    }
    return {};
  }

  static FieldCell<typename M::Scal> GetAsymmetryField(
      const FieldCell<typename M::Expr>& fce, M& m) {
    constexpr size_t vdim = Expr::dim - 2;
    static_assert(vdim == M::kCellNumNeighborFaces, "");
    auto sem = m.GetSem(__func__);
    struct {
      std::array<FieldCell<Scal>, vdim> vfck; // coefficient of neighbors
    } * ctx(sem);
    auto& vfck = ctx->vfck;
    if (sem("comm")) {
      for (size_t i = 0; i < vdim; ++i) {
        vfck[i] = GetComponent(fce, i + 1);
        m.Comm(&vfck[i]);
      }
    }
    if (sem()) {
      // FIXME: empty stage to finish communication in inner blocks
    }
    if (sem("local")) {
      FieldCell<Scal> fcr(m, 0);
      for (auto c : m.Cells()) {
        Scal sum = 0;
        for (auto q : m.Nci(c)) {
          auto cn = m.GetCell(c, q);
          sum += std::abs(vfck[q.raw()][c] - vfck[m.GetOpposite(q).raw()][cn]);
        }
        fcr[c] = sum;
      }
      return fcr;
    }
    return {};
  }

  static void CheckSymmetry(const FieldCell<Expr>& fce, M& m) {
    auto sem = m.GetSem(__func__);
    struct {
      FieldCell<Scal> fc_asymm;
      FieldCell<Scal> fc_diag;
      FieldCell<Scal> fc_nondiag;
      FieldCell<Scal> fc_const;
      std::vector<generic::Vect<Scal, 3>> vnorms;
      Scal sum_const;
    } * ctx(sem);
    auto& fc_asymm = ctx->fc_asymm;
    auto& fc_diag = ctx->fc_diag;
    auto& fc_nondiag = ctx->fc_nondiag;
    auto& fc_const = ctx->fc_const;
    auto& vnorms = ctx->vnorms;
    auto& sum_const = ctx->sum_const;
    if (sem.Nested()) {
      fc_asymm = UDebug<M>::GetAsymmetryField(fce, m);
    }
    if (sem("local")) {
      fc_diag.Reinit(m, 0);
      fc_nondiag.Reinit(m, 0);
      fc_const.Reinit(m, 0);
      for (auto c : m.Cells()) {
        fc_diag[c] = std::abs(fce[c][0]);
        fc_const[c] = std::abs(fce[c].back());
        for (auto q : m.Nci(c)) {
          fc_nondiag[c] += std::abs(fce[c][1 + q.raw()]);
        }
      }
    }
    if (sem.Nested()) {
      vnorms =
          UDebug<M>::GetNorms({&fc_asymm, &fc_diag, &fc_nondiag, &fc_const}, m);
    }
    if (sem("sum_const")) {
      sum_const = 0;
      for (auto c : m.Cells()) {
        sum_const += fce[c].back();
      }
      m.Reduce(&sum_const, "sum");
    }
    if (sem("print")) {
      if (m.IsRoot()) {
        auto flags = std::cerr.flags();
        std::cerr << std::scientific;
        std::cerr.precision(10);
        std::cerr << "CheckSymmetry() of system '" + fce.GetName() + "':"
                  << std::endl;
        std::cerr << "asymm: " << vnorms[0] << std::endl;
        std::cerr << "diag: " << vnorms[1] << std::endl;
        std::cerr << "nondiag: " << vnorms[2] << std::endl;
        std::cerr << "const: " << vnorms[3] << std::endl;
        std::cerr << "sum_const: " << sum_const << std::endl;
        std::cerr.flags(flags);
      }
      if (vnorms[0][2] > m.flags.check_symmetry_dump_threshold) {
        m.Dump(&fc_asymm, "asymm_" + fce.GetName());
      }
    }
    if (sem()) {
    }
  }
};
