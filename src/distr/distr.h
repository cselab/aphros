// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

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
#include "parse/vars.h"
#include "util/metrics.h"
#include "util/module.h"
#include "util/mpi.h"
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
  using RedOp = typename M::Op;
  using BlockInfoProxy = generic::BlockInfoProxy<dim>;

  virtual void Run();
  virtual void Report();
  virtual void ReportOpenmp();
  virtual ~DistrMesh();

 protected:
  // TODO: remove comm, needed only by Hypre
  MPI_Comm comm_; // XXX: overwritten by Cubism<M>
  const Vars& var;
  Vars& var_mutable;
  const KernelMeshFactory<M>& kernelfactory_; // kernel factory

  int halos_; // number of halo cells (same in all directions)
  MIdx blocksize_; // block size
  MIdx nprocs_; // number of ranks
  MIdx nblocks_; // number of blocks
  Scal extent_; // extent (maximum over all directions)

  int stage_ = 0;
  size_t frame_ = 0; // current dump frame

  bool isroot_; // XXX: overwritten by Local<M> and Cubism<M>

  struct Hash {
    size_t operator()(const MIdx& x) const noexcept {
      const size_t h0 = std::hash<IntIdx>{}(x[0]);
      const size_t h1 = std::hash<IntIdx>{}(x[1]);
      const size_t h2 = std::hash<IntIdx>{}(x[2]);
      return h0 ^ (h1 << 1) ^ (h2 << 2);
    }
  };

  std::vector<std::unique_ptr<KernelMesh<M>>> kernels_;

  DistrMesh(MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var);
  // Performs communication and returns indices of blocks with updated halos.
  virtual std::vector<size_t> TransferHalos(bool inner) = 0;
  virtual std::vector<size_t> TransferHalos() {
    auto bbi = TransferHalos(true);
    auto bbh = TransferHalos(false);
    bbi.insert(bbi.end(), bbh.begin(), bbh.end());
    return bbi;
  }
  // Fill selected halo cells with garbage
  void ApplyNanFaces(const std::vector<size_t>& bb);
  // Call kernels for current stage
  virtual void RunKernels(const std::vector<size_t>& bb);
  // Performs reduction with a single request over all blocks.
  // block_request: request for each block, same dimension as `kernels_`
  virtual void ReduceSingleRequest(const std::vector<RedOp*>& blocks) = 0;
  void ReduceSingleRequest(const std::vector<std::unique_ptr<RedOp>>& blocks);
  // Performs reduction with all requests collected in m.GetReduce()
  // for all elements in `kernels_`
  virtual void Reduce(const std::vector<size_t>& bb);
  virtual void ReduceToLead(const std::vector<size_t>& bb);
  virtual void Scatter(const std::vector<size_t>& bb) = 0;
  virtual void Bcast(const std::vector<size_t>& bb) = 0;
  virtual void BcastFromLead(const std::vector<size_t>& bb);
  virtual void DumpWrite(const std::vector<size_t>& bb);
  virtual void ClearComm(const std::vector<size_t>& bb);
  virtual void ClearDump(const std::vector<size_t>& bb);
  // TODO: make Pending const
  virtual bool Pending(const std::vector<size_t>& bb);
  // Create a kernel for each block and put into kernels_
  // Requires initialized isroot_;
  virtual void MakeKernels(const std::vector<BlockInfoProxy>&);
  virtual void TimerReport(const std::vector<size_t>& bb);
  virtual void ClearTimerReport(const std::vector<size_t>& bb);

 private:
  MultiTimer<std::string> multitimer_all_;
  MultiTimer<std::string> multitimer_report_;
};

template <class M>
class ModuleDistr : public Module<ModuleDistr<M>> {
 public:
  using Module<ModuleDistr>::Module;
  virtual ~ModuleDistr() = default;
  virtual std::unique_ptr<DistrMesh<M>> Make(
      MPI_Comm, const KernelMeshFactory<M>&, Vars&) = 0;
};
