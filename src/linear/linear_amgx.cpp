// Created by Petr Karnakov on 15.11.2020
// Copyright 2020 ETH Zurich

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "amgx.h"
#include "distr/commmap.h"
#include "distr/distrsolver.h"
#include "linear_amgx.h"

DECLARE_FORCE_LINK_TARGET(linear_amgx);

namespace linear {

template <class M>
struct SolverAmgx<M>::Imp {
  using Owner = SolverAmgx<M>;

  Imp(Owner* owner, const Extra& extra_)
      : owner_(owner), conf(owner_->conf), extra(extra_) {}
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) {
    auto sem = m.GetSem(__func__);
    struct {
      CommMap<M> comm_map;
      std::ofstream logfile;
      Info info;
    } * ctx(sem);
    static std::vector<Scal> buf_sol;
    static std::vector<Scal> buf_rhs;
    auto& t = *ctx;
    if (sem.Nested()) {
      t.comm_map.Init(m);
    }
    if (sem.Nested()) {
      t.comm_map.ConvertSystem(fc_system, m);
    }
    if (sem("solve")) {
      FieldCell<Scal> fc_rhs(m, 0);
      for (auto c : m.Cells()) {
        fc_rhs[c] = -fc_system[c].back();
      }
      if (fc_init) {
        fc_sol = *fc_init;
      } else {
        fc_sol.Reinit(m, 0);
      }
      t.comm_map.FieldToArray(fc_sol, buf_sol);
      t.comm_map.FieldToArray(fc_rhs, buf_rhs);
    }
    if (sem() && m.IsLead()) {
      const auto comm = m.GetMpiComm();
      auto& system = t.comm_map.GetSystem();

      const int gpu_id = 0;

      t.logfile.open(extra.log_path);
      Amgx::Library amgx(&t.logfile);

      const Amgx::Mode mode(extra.mode);
      Amgx::Config config(
          extra.config_path,
          "communicator=MPI, exception_handling=1" + extra.config_extra);
      Amgx::Resources resources(config, comm, gpu_id);

      Amgx::Matrix matrix(resources, mode);
      Amgx::Vector sol(resources, mode);
      Amgx::Vector rhs(resources, mode);

      const int n = system.row_ptrs.size() - 1;
      const int nnz = system.cols.size();
      fassert_equal(system.cols.size(), system.data.size());

      auto send_maps = const_cast<const int**>(system.send_maps.data());
      auto recv_maps = const_cast<const int**>(system.recv_maps.data());

      AMGXCALL(AMGX_matrix_comm_from_maps_one_ring(
          matrix, 1, system.neighbors.size(), system.neighbors.data(),
          system.send_sizes.data(), send_maps, system.recv_sizes.data(),
          recv_maps));

      AMGXCALL(AMGX_matrix_upload_all(
          matrix, n, nnz, 1, 1, system.row_ptrs.data(), system.cols.data(),
          system.data.data(), nullptr));

      sol.Bind(matrix);
      rhs.Bind(matrix);

      sol.Upload(buf_sol.data(), {(int)buf_sol.size(), 1});
      rhs.Upload(buf_rhs.data(), {(int)buf_rhs.size(), 1});

      Amgx::Solver solver(resources, mode, config);
      AMGXCALL(AMGX_solver_setup(solver, matrix));
      AMGXCALL(AMGX_solver_solve(solver, rhs, sol));

      t.info.residual = solver.GetResidual();
      t.info.iter = solver.GetNumIters();
      buf_sol = sol.Download<Scal>();
    }
    if (sem()) {
      t.comm_map.ArrayToField(buf_sol, fc_sol, m);
      m.Comm(&fc_sol);
      if (m.flags.linreport && m.IsRoot()) {
        std::cout << std::scientific;
        std::cout << "linear(amgx) '" + fc_system.GetName() + "':"
                  << " res=" << t.info.residual << " iter=" << t.info.iter
                  << std::endl;
      }
    }
    if (sem()) {
    }
    return t.info;
  }

 private:
  Owner* owner_;
  Conf& conf;
  Extra extra;
};

template <class M>
SolverAmgx<M>::SolverAmgx(const Conf& conf_, const Extra& extra)
    : Base(conf_), imp(new Imp(this, extra)) {}

template <class M>
SolverAmgx<M>::~SolverAmgx() = default;

template <class M>
auto SolverAmgx<M>::Solve(
    const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
    FieldCell<Scal>& fc_sol, M& m) -> Info {
  return imp->Solve(fc_system, fc_init, fc_sol, m);
}

template <class M>
class ModuleLinearAmgx : public ModuleLinear<M> {
 public:
  ModuleLinearAmgx() : ModuleLinear<M>("amgx") {}
  std::unique_ptr<Solver<M>> Make(
      const Vars& var, std::string prefix) override {
    auto addprefix = [prefix](std::string name) {
      return "amgx_" + prefix + "_" + name;
    };
    typename SolverAmgx<M>::Extra extra;
    extra.config_path = var.String[addprefix("config_path")];
    extra.config_extra = var.String[addprefix("config_extra")];
    extra.mode = var.String[addprefix("mode")];
    return std::make_unique<linear::SolverAmgx<M>>(
        this->GetConf(var, prefix), extra);
  }
};

using M = MeshStructured<double, 3>;

bool kReg_amgx[] = {
    RegisterModule<ModuleLinearAmgx<M>>(),
};

} // namespace linear
