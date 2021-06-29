// Created by Petr Karnakov on 05.03.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <array>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>

#include "debug/isnan.h"
#include "dump/dumper.h"
#include "dump/hdf.h"
#include "func/init.h"
#include "func/init_contang.h"
#include "geom/mesh.h"
#include "kernelmeshpar.h"
#include "linear/linear.h"
#include "parse/curv.h"
#include "parse/parser.h"
#include "parse/proj.h"
#include "parse/simple.h"
#include "parse/util.h"
#include "parse/vars.h"
#include "parse/vof.h"
#include "solver/advection.h"
#include "solver/approx.h"
#include "solver/approx2.h"
#include "solver/approx_eb.h"
#include "solver/curv.h"
#include "solver/electro.h"
#include "solver/embed.h"
#include "solver/fluid_dummy.h"
#include "solver/multi.h"
#include "solver/normal.h"
#include "solver/particles.h"
#include "solver/proj.h"
#include "solver/reconst.h"
#include "solver/simple.h"
#include "solver/solver.h"
#include "solver/tracer.h"
#include "solver/vof.h"
#include "solver/vofm.h"
#include "util/convdiff.h"
#include "util/events.h"
#include "util/filesystem.h"
#include "util/format.h"
#include "util/hydro.h"
#include "util/linear.h"
#include "util/metrics.h"
#include "util/posthook.h"
#include "util/stat.h"
#include "util/sysinfo.h"
#include "util/timer.h"

template <class M>
class Hydro;

template <class M>
struct GPar {
  void* ptr = nullptr; // opaque pointer to data
  std::function<void(void* ptr, Hydro<M>* hydro)> step_callback;
};

template <class M>
class ModulePostStep : public Module<ModulePostStep<M>> {
 public:
  using Vect = typename M::Vect;
  using Module<ModulePostStep>::Module;
  virtual ~ModulePostStep() = default;
  virtual void operator()(Hydro<M>*, M& m) = 0;
};

template <class M_>
class Hydro : public KernelMeshPar<M_, GPar<M_>> {
 public:
  using M = M_;
  using Par = GPar<M>;
  using P = KernelMeshPar<M_, Par>; // parent
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Sem = typename M::Sem;
  template <class T>
  using Multi = Multi<T>;
  using UEB = UEmbed<M>;
  using ParticlesView = typename ParticlesInterface<M>::ParticlesView;
  static constexpr size_t dim = M::dim;
  using BlockInfoProxy = generic::BlockInfoProxy<dim>;
  friend void StepHook<>(Hydro*);

  // TODO: issue warning if variable in Vars was not used
  // but differs from default (like in CMake)
  Hydro(Vars&, const BlockInfoProxy&, Par&);
  void Run() override;
  M& GetMesh() {
    return m;
  }

  // FIXME: make private and introduce public interface
 public:
  using P::m;
  using P::par_;
  using P::var;
  using P::var_mutable;

  // FIXME: make private and introduce public interface
  // Needed for ModulePostStep
 public:
  void Init();
  void InitEmbed();
  void InitStepwiseBody(FieldCell<bool>& fc_innermask);
  void InitTracer(Multi<FieldCell<Scal>>& vfcu);
  void InitTracerFields(Multi<FieldCell<Scal>>& vfcu);
  void InitElectro();
  void SpawnTracer();
  void InitParticles();
  void SpawnParticles(ParticlesView& view);
  void OverwriteBc();
  void InitFluid(const FieldCell<Vect>& fc_vel);
  void InitAdvection(const FieldCell<Scal>& fcvf, const FieldCell<Scal>& fccl);
  template <class MEB>
  void InitStat(const MEB& eb);
  Vect CalcPressureDrag(const FieldCell<Scal>& fcp, const Embed<M>& eb);
  Vect CalcViscousDrag(
      const FieldCell<Vect>& fcvel, const FieldCell<Scal>& fcmu,
      const Embed<M>& eb);
  void Dump(bool force);
  // Calc rho, mu and force based on volume fraction
  void CalcMixture(const FieldCell<Scal>& vf);
  // Clips v to given range, uses const_cast
  void Clip(const FieldCell<Scal>& v, Scal min, Scal max);
  void CalcStat();
  void CalcDt();
  void ReportStep();
  void ReportStepAdv();
  void ReportStepTracer();
  void ReportStepParticles();
  void ReportStepElectro();
  void ReportSysinfo(std::ostream& out);
  void ReportIter();
  // Issue sem.LoopBreak if abort conditions met
  void CheckAbort(Sem& sem, Scal& nabort);
  void StepFluid();
  void StepAdvection();
  void StepTracer();
  void StepParticles();
  void StepElectro();
  void StepBubgen();
  void StepEraseVolumeFraction(std::string prefix, Scal& last_t);
  void StepEraseColor(std::string prefix);

  using ASV = Vof<M>; // advection VOF
  using ASVM = Vofm<M>; // advection multi VOF
  using ASVEB = Vof<Embed<M>>; // advection VOF embed
  using ASVMEB = Vofm<Embed<M>>; // advection multi VOF embed
  using EB = Embed<M>;
  static constexpr Scal kClNone = ASVM::kClNone;

  void UpdateAdvectionPar();
  Scal GetSurfaceTensionDt() const;
  Scal GetViscosityDt() const;
  void CalcVort();
  FieldCell<Scal> CalcStrain(const FieldCell<Vect> fcvel) const;
  FieldCell<Scal> GetDiv() const;

  bool initialized_ = false;
  bool finished_ = false;
  bool silent_ = false;
  bool dumpstat_ = true; // write statistics to stat.dat, stat_summary

  GRange<size_t> layers;
  FieldCell<Scal> fc_mu_; // viscosity
  FieldCell<Scal> fc_rho_; // density
  FieldCell<Scal> fc_src_; // source of mixture volume
  FieldCell<Scal> fc_src2_; // source of second phase volume
  FieldCell<Scal> fc_srcm_; // mass source
  FieldCell<Vect> fc_force_; // force
  FieldCell<Scal> fc_phi_; // distance from eb
  FieldEmbed<Scal> febp_; // balanced force projections

  MapEmbed<BCondAdvection<Scal>> mebc_adv_;
  MapEmbed<BCond<Scal>> mebc_vfsm_;
  MapEmbed<BCondFluid<Vect>> mebc_fluid_;
  MapEmbed<BCondFluid<Vect>> mebc_fluid_orig_;
  MapEmbed<BCond<Scal>> mebc_electro_;
  MapEmbed<size_t> me_group_;
  MapEmbed<Scal> me_contang_;
  std::vector<std::string> bc_group_desc_;
  std::vector<std::map<std::string, Scal>> bc_group_custom_;
  MapCell<std::shared_ptr<CondCell>> mc_cond_;
  MapCell<std::shared_ptr<CondCellFluid>> mc_velcond_;

  std::unique_ptr<Embed<M>> eb_;
  std::unique_ptr<AdvectionSolver<M>> as_;
  std::unique_ptr<FluidSolver<M>> fs_;

  FieldCell<Scal> fc_smvf_; // smoothed volume fraction used by CalcMixture()
  FieldFace<Scal> ff_smvf_; // interpolated fc_smvf_

  FieldCell<Scal> fc_sig_; // surface tension sigma
  FieldCell<Scal> fc_contang_; // contact angle

  FieldCell<Vect> fcvm_; // velocity field time_prev // TODO: revise

  FieldCell<Vect> fcom_; // vorticity
  FieldCell<Scal> fcomm_; // vorticity magnitude
  FieldCell<Scal> fc_strain_; // double inner product of strain rate tensor
  Multi<FieldCell<Scal>> fck_; // curvature
  std::unique_ptr<curvature::Estimator<M>> curv_estimator_;

  Scal bgt_ = -1.; // bubgen last time
  Scal erasevf_last_t_ = -std::numeric_limits<Scal>::max();
  Scal erasevf2_last_t_ = -std::numeric_limits<Scal>::max();

  struct StatHydro {
    Scal dt = 0; // dt fluid
    Scal dta = 0; // dt advection
    Scal t = 0;
    size_t step = 0;
    size_t iter = 0;
    Vect meshvel = {};
  };
  StatHydro st_;
  std::ofstream fstat_;
  std::unique_ptr<Stat<M>> stat_;
  Dumper dumper_;
  Dumper dmptraj_; // dumper for traj
  Dumper dmptrep_; // dumper for timer report
  Dumper bubgen_;
  std::unique_ptr<Events> events_; // events from var
  SingleTimer timer_;
  std::shared_ptr<linear::Solver<M>> linsolver_symm_;

  std::unique_ptr<TracerInterface<M>> tracer_;
  Multi<FieldCell<Scal>> fc_tracer_source;

  std::set<IdxCell> nucl_cells_;

  std::unique_ptr<ParticlesInterface<M>> particles_;
  std::mt19937 randgen_;
  Scal tracer_dt_;
  Scal particles_dt_;
  std::string vf_save_state_path_;

  // electro
  std::unique_ptr<ElectroInterface<M>> electro_;

  ModulePostStep<M>* module_post_step_ = nullptr;
};
