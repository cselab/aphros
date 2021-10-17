// Created by Petr Karnakov on 15.11.2020
// Copyright 2020 ETH Zurich

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "distr/commmap.h"
#include "distr/distrsolver.h"
#include "linear_amgx.h"
#include "amgx.h"

DECLARE_FORCE_LINK_TARGET(linear_amgx);

namespace linear {

template <class M>
struct SolverAmgx<M>::Imp {
  using Owner = SolverAmgx<M>;

  Imp(Owner* owner, const Extra& extra_, const M& m)
      : owner_(owner), conf(owner_->conf), extra(extra_) {
    if (m.IsLead()) {
      if (extra.log_path.length()) {
        logfile_ = std::make_unique<std::ofstream>(extra.log_path);
      }
      amgx_ = std::make_unique<Amgx::Library>(logfile_.get());
    }
  }
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) {
    auto sem = m.GetSem(__func__);
    struct {
      Info info;
    } * ctx(sem);
    auto& t = *ctx;
    if (sem.Nested()) {
      if (!comm_map_init_called_) {
        comm_map_.Init(m);
      }
    }
    if (sem.Nested()) {
      comm_map_.ConvertSystem(fc_system, m);
    }
    if (sem("copyrhs")) {
      comm_map_init_called_ = true;
      FieldCell<Scal> fc_rhs(m, 0);
      for (auto c : m.Cells()) {
        fc_rhs[c] = -fc_system[c].back();
      }
      if (fc_init) {
        fc_sol = *fc_init;
      } else {
        fc_sol.Reinit(m, 0);
      }
      auto buf = comm_map_.GetBuffers();
      comm_map_.FieldToArray(fc_rhs, buf.rhs);
      comm_map_.FieldToArray(fc_sol, buf.sol);
    }
    if (sem("call") && m.IsLead()) {
      const auto comm = m.GetMpiComm();
      auto& system = comm_map_.GetSystem();

      const int gpu_id = 0;

      const Amgx::Mode mode(extra.mode);
      Amgx::Config config(
          extra.config_path, "communicator=MPI" + extra.config_extra);
      Amgx::Resources resources(config, comm, gpu_id);

      Amgx::Matrix matrix(resources, mode);
      Amgx::Vector sol(resources, mode);
      Amgx::Vector rhs(resources, mode);

      auto send_maps = const_cast<const int**>(system.send_maps.data());
      auto recv_maps = const_cast<const int**>(system.recv_maps.data());

      amgx_->matrix_comm_from_maps_one_ring(
          matrix, 1, system.neighbors.size(), system.neighbors.data(),
          system.send_sizes.data(), send_maps, system.recv_sizes.data(),
          recv_maps);

      amgx_->matrix_upload_all(
          matrix, system.n, system.nnz, 1, 1, system.row_ptrs.data(),
          system.cols.data(), system.data.data(), nullptr);

      sol.Bind(matrix);
      rhs.Bind(matrix);

      auto buf = comm_map_.GetBuffers();
      sol.Upload(buf.sol, {(int)buf.size, 1});
      rhs.Upload(buf.rhs, {(int)buf.size, 1});

      Amgx::Solver solver(resources, mode, config);
      amgx_->solver_setup(solver, matrix);
      amgx_->solver_solve(solver, rhs, sol);

      t.info.residual = solver.GetResidual();
      t.info.iter = solver.GetNumIters();
      sol.Download<Scal>(buf.sol);
    }
    if (sem("copysol")) {
      auto buf = comm_map_.GetBuffers();
      comm_map_.ArrayToField(buf.sol, fc_sol, m);
      m.Comm(&fc_sol);
      if (m.flags.linreport && m.IsRoot()) {
        std::cerr << std::scientific;
        std::cerr << "linear(amgx) '" + fc_system.GetName() + "':"
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
  CommMap<M> comm_map_;
  bool comm_map_init_called_ = false;
  std::unique_ptr<Amgx::Library> amgx_;
  std::unique_ptr<std::ofstream> logfile_;
};

template <class M>
SolverAmgx<M>::SolverAmgx(const Conf& conf_, const Extra& extra, const M& m)
    : Base(conf_), imp(new Imp(this, extra, m)) {}

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
      const Vars& var, std::string prefix, const M& m) override {
    auto addprefix = [prefix](std::string name) {
      return "amgx_" + prefix + "_" + name;
    };
    typename SolverAmgx<M>::Extra extra;
    extra.config_path = var.String[addprefix("config_path")];
    extra.config_extra = var.String[addprefix("config_extra")];
    extra.log_path = var.String[addprefix("log_path")];
    extra.mode = var.String[addprefix("mode")];
    return std::make_unique<linear::SolverAmgx<M>>(
        this->GetConf(var, prefix), extra, m);
  }
};

using M = MeshCartesian<double, 3>;

bool kReg_amgx[] = {
    RegisterModule<ModuleLinearAmgx<M>>(),
};

} // namespace linear
