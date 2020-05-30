// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <mpi.h>
#include <omp.h>
#include <sched.h>
#include <array>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "dump/dump.h"
#include "dump/dumper.h"
#include "geom/mesh.h"
#include "kernel/kernelmesh.h"
#include "linear/hypre.h"
#include "linear/hypresub.h"
#include "parse/vars.h"
#include "report.h"
#include "util/histogram.h"
#include "util/metrics.h"
#include "util/suspender.h"
#include "util/sysinfo.h"

// Abstract block processor aware of Mesh.
template <class M_>
class DistrMesh {
 public:
  using M = M_;
  static constexpr size_t dim = M::dim;
  using MIdx = typename M::MIdx;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  virtual void Run();
  virtual void Report();
  virtual void ReportOpenmp();
  virtual ~DistrMesh();
  virtual typename M::BlockCells GetGlobalBlock() const;
  virtual typename M::IndexCells GetGlobalIndex() const;
  // Returns data field i from buffer defined on global mesh
  virtual FieldCell<Scal> GetGlobalField(size_t i);

 protected:
  // TODO: remove comm, needed only by Hypre
  MPI_Comm comm_; // XXX: overwritten by Cubism<M>
  const Vars& var;
  Vars& var_mutable;
  const KernelMeshFactory<M>& kf_; // kernel factory
  Sampler samp_; // sampler accessible to derived classes

  int hl_; // number of halo cells (same in all directions)
  MIdx bs_; // block size
  MIdx p_; // number of ranks
  MIdx b_; // number of blocks
  Scal ext_; // extent (maximum over all directions)

  int stage_ = 0;
  size_t frame_ = 0; // current dump frame

  bool isroot_; // XXX: overwritten by Local<M> and Cubism<M>

  std::map<MIdx, std::unique_ptr<KernelMesh<M>>, typename MIdx::LexLess> mk;

  DistrMesh(MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var);
  // Performs communication and returns indices of blocks with updated halos.
  virtual std::vector<MIdx> GetBlocks(bool inner) = 0;
  virtual std::vector<MIdx> GetBlocks() {
    auto bbi = GetBlocks(true);
    auto bbh = GetBlocks(false);
    bbi.insert(bbi.end(), bbh.begin(), bbh.end());
    return bbi;
  }
  // Fill selected halo cells with garbage
  void ApplyNanFaces(const std::vector<MIdx>& bb);
  // Copy from communication buffer to fields
  virtual void ReadBuffer(const std::vector<MIdx>& bb) = 0;
  // Call kernels for current stage
  virtual void RunKernels(const std::vector<MIdx>& bb);
  // Copy from fields to communication buffer
  virtual void WriteBuffer(const std::vector<MIdx>& bb) = 0;
  // Reduce TODO: extend doc
  virtual void Reduce(const std::vector<MIdx>& bb) = 0;
  virtual void Scatter(const std::vector<MIdx>& bb) = 0;
  virtual void Bcast(const std::vector<MIdx>& bb) = 0;
  virtual void Solve(const std::vector<MIdx>& bb);
  // Writes dumps.
  virtual void DumpWrite(const std::vector<MIdx>& bb);
  virtual void ClearComm(const std::vector<MIdx>& bb);
  virtual void ClearDump(const std::vector<MIdx>& bb);
  // TODO: make Pending const
  virtual bool Pending(const std::vector<MIdx>& bb);
  // Create a kernel for each block and put into mk
  // Requires initialized isroot_;
  virtual void MakeKernels(const std::vector<MyBlockInfo>&);
  virtual void TimerReport(const std::vector<MIdx>& bb);
  virtual void ClearTimerReport(const std::vector<MIdx>& bb);

 private:
  Histogram hist_; // histogram sample collector
  MultiTimer<std::string> mt_; // timer all
  MultiTimer<std::string> mtp_; // timer partial
  using LS = typename M::LS;
  std::map<typename LS::T, std::unique_ptr<HypreSub>> mhp_; // hypre instances
};
