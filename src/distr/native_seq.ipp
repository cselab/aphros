// Created by Petr Karnakov on 02.01.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <cassert>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>

#include "comm_manager.h"
#include "distr.h"
#include "dump/dumper.h"
#include "util/format.h"
#include "util/mpi.h"

template <class M_>
class Native : public DistrMesh<M_> {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using RedOp = typename M::Op;
  using BlockInfoProxy = generic::BlockInfoProxy<M::dim>;

  Native(MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var);
  ~Native();

 private:
  using P = DistrMesh<M>; // parent
  using MIdx = typename M::MIdx;

  std::vector<size_t> TransferHalos(bool inner) override;
  void ReduceSingleRequest(const std::vector<RedOp*>& blocks) override;
  void Bcast(const std::vector<size_t>& bb) override;
  void Scatter(const std::vector<size_t>& bb) override;
  void DumpWrite(const std::vector<size_t>& bb) override;

  using P::blocksize_;
  using P::comm_;
  using P::dim;
  using P::extent_;
  using P::frame_;
  using P::halos_;
  using P::isroot_;
  using P::kernelfactory_;
  using P::kernels_;
  using P::nblocks_;
  using P::nprocs_;
  using P::stage_;
  using P::var;

  int commsize_;
  int commrank_;
  typename CommManager<dim>::Tasks tasks_;
  struct ReqTmp {
    size_t cnt = 0;
    std::vector<Scal> buf;
  };
  std::map<int, ReqTmp> tmp_send_; // rank to request
  std::map<int, ReqTmp> tmp_recv_;
};

template <class M>
Native<M>::Native(MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var_)
    : DistrMesh<M>(comm, kf, var_) {
  MpiWrapper mpi(comm_);
  commsize_ = mpi.GetCommSize();
  commrank_ = mpi.GetCommRank();
  isroot_ = (0 == commrank_); // XXX: overwrite isroot_

  if (commsize_ != nprocs_.prod()) {
    throw std::runtime_error(util::Format(
        "Number of MPI tasks {} does not match the number of subdomains {}",
        commsize_, nprocs_));
  }

  const MIdx globalsize = nprocs_ * nblocks_ * blocksize_;
  std::vector<BlockInfoProxy> proxies;
  GIndex<size_t, dim> procs(nprocs_);
  GIndex<size_t, dim> blocks(nblocks_);
  for (auto ib : blocks.Range()) {
    BlockInfoProxy p;
    p.index = nblocks_ * procs.GetMIdx(commrank_) + blocks.GetMIdx(ib);
    p.globalsize = globalsize;
    p.cellsize = Vect(extent_ / p.globalsize.max());
    p.blocksize = blocksize_;
    p.halos = halos_;
    p.isroot = (ib == 0 && isroot_);
    p.islead = (ib == 0);
    proxies.push_back(p);
  }

  this->MakeKernels(proxies);

  {
    std::vector<typename CommManager<dim>::Block> cm_blocks;
    for (auto& k : kernels_) {
      auto& m = k->GetMesh();
      cm_blocks.push_back({&m.GetInBlockCells(), &m.GetIndexCells()});
    }
    auto cell_to_rank = [&procs, this](MIdx w) -> int {
      return procs.GetIdx(w / blocksize_ / nblocks_);
    };
    const generic::Vect<bool, dim> is_periodic(true);
    tasks_ = CommManager<dim>::GetTasks(
        cm_blocks, cell_to_rank, globalsize, is_periodic, mpi);
  }
}

template <class M>
Native<M>::~Native() = default;

template <class M>
auto Native<M>::TransferHalos(bool inner) -> std::vector<size_t> {
  if (!inner) {
    return {};
  }
  std::vector<size_t> bb(kernels_.size());
  std::iota(bb.begin(), bb.end(), 0);
  auto& vcr = kernels_.front()->GetMesh().GetComm();
  if (vcr.empty()) {
    return bb;
  }

  using Task = typename CommManager<dim>::Task;
  using CommStencil = typename M::CommStencil;
  std::array<std::pair<const Task*, CommStencil>, 4> taskstencils{
      std::make_pair(&tasks_.full_two, CommStencil::full_two),
      std::make_pair(&tasks_.full_one, CommStencil::full_one),
      std::make_pair(&tasks_.direct_two, CommStencil::direct_two),
      std::make_pair(&tasks_.direct_one, CommStencil::direct_one),
  };

  for (auto& pair : taskstencils) {
    auto& task = *pair.first;
    auto stencil = pair.second;

    // indices of communication requests that have selected `stencil`
    std::vector<size_t> vcr_indices;
    for (size_t i = 0; i < vcr.size(); ++i) {
      if (vcr[i]->stencil == stencil) {
        vcr_indices.push_back(i);
      }
    }

    // Exchange data between blocks.
    for (auto i : vcr_indices) {
      if (dynamic_cast<typename M::CommRequestScal*>(vcr[i].get())) {
        std::vector<FieldCell<Scal>*> fields;
        for (auto& k : kernels_) {
          const auto* cr = static_cast<typename M::CommRequestScal*>(
              k->GetMesh().GetComm()[i].get());
          fields.push_back(cr->field);
        }
        auto& send = task.send.at(0);
        auto& recv = task.recv.at(0);
        fassert_equal(send.size(), recv.size());
        for (size_t ic = 0; ic < send.size(); ++ic) {
          (*fields[recv[ic].block])[recv[ic].cell] =
              (*fields[send[ic].block])[send[ic].cell];
        }
      }
      if (auto crd = dynamic_cast<typename M::CommRequestVect*>(vcr[i].get())) {
        std::vector<FieldCell<Vect>*> fields;
        for (auto& k : kernels_) {
          const auto* cr = static_cast<typename M::CommRequestVect*>(
              k->GetMesh().GetComm()[i].get());
          fields.push_back(cr->field);
        }
        auto& send = task.send.at(0);
        auto& recv = task.recv.at(0);
        fassert_equal(send.size(), recv.size());
        for (size_t ic = 0; ic < send.size(); ++ic) {
          const auto d = crd->d;
          if (d == -1) {
            (*fields[recv[ic].block])[recv[ic].cell] =
                (*fields[send[ic].block])[send[ic].cell];
          } else {
            fassert(0 <= d && d < int(M::dim));
            (*fields[recv[ic].block])[recv[ic].cell][d] =
                (*fields[send[ic].block])[send[ic].cell][d];
          }
        }
      }
    }
  }

  return bb;
}

template <class M>
void Native<M>::Bcast(const std::vector<size_t>& bb) {
  using OpConcat = typename UReduce<Scal>::OpCat;
  auto& vfirst = kernels_.front()->GetMesh().GetBcast();

  if (!vfirst.size()) {
    return;
  }

  for (auto b : bb) {
    fassert_equal(kernels_[b]->GetMesh().GetBcast().size(), vfirst.size());
  }

  for (size_t i = 0; i < vfirst.size(); ++i) {
    if (OpConcat* o = dynamic_cast<OpConcat*>(vfirst[i].get())) {
      std::vector<char> buf = o->Neutral();

      if (isroot_) {
        // read from root block
        for (auto b : bb) {
          auto& m = kernels_[b]->GetMesh();
          if (m.IsRoot()) {
            auto& v = m.GetBcast();
            OpConcat* ob = dynamic_cast<OpConcat*>(v[i].get());
            ob->Append(buf);
          }
        }
      }

      // write to all blocks
      for (auto b : bb) {
        auto& m = kernels_[b]->GetMesh();
        auto& v = m.GetBcast();
        OpConcat* ob = dynamic_cast<OpConcat*>(v[i].get());
        ob->Set(buf);
      }
    } else {
      throw std::runtime_error("Bcast: Unknown M::Op instance");
    }
  }

  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearBcast();
  }
}
template <class M>
void Native<M>::Scatter(const std::vector<size_t>& bb) {
  auto& vfirst = kernels_.front()->GetMesh().GetScatter();

  if (!vfirst.size()) {
    return;
  }

  for (auto b : bb) {
    fassert_equal(kernels_[b]->GetMesh().GetScatter().size(), vfirst.size());
  }

  for (size_t q = 0; q < vfirst.size(); ++q) {
    const size_t b_root = [&]() {
      for (auto b : bb) {
        if (kernels_[b]->GetMesh().IsRoot()) {
          return b;
        }
      }
      fassert(false);
    }();
    auto& req_root = kernels_.at(b_root)->GetMesh().GetScatter()[q];
    for (auto b : bb) {
      auto& req = kernels_[b]->GetMesh().GetScatter()[q];
      (*req.second) = (*req_root.first)[b];
    }
  }

  // Clear requests
  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearScatter();
  }
}

template <class M>
void Native<M>::ReduceSingleRequest(const std::vector<RedOp*>& blocks) {
  using OpScal = typename UReduce<Scal>::OpS;
  using OpScalInt = typename UReduce<Scal>::OpSI;
  using OpConcat = typename UReduce<Scal>::OpCat;

  auto* firstbase = blocks.front();

  if (auto* first = dynamic_cast<OpScal*>(firstbase)) {
    auto buf = first->Neutral(); // result
    // Reduce over blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScal*>(otherbase);
      other->Append(buf);
    }
    // Write results to all blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScal*>(otherbase);
      other->Set(buf);
    }
    return;
  }
  if (auto* first = dynamic_cast<OpScalInt*>(firstbase)) {
    auto buf = first->Neutral(); // result
    // Reduce over blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScalInt*>(otherbase);
      other->Append(buf);
    }
    // Write results to all blocks on current rank
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpScalInt*>(otherbase);
      other->Set(buf);
    }
    return;
  }
  if (auto* first = dynamic_cast<OpConcat*>(firstbase)) {
    auto buf = first->Neutral();
    // Reduce over local blocks
    for (auto otherbase : blocks) {
      auto* other = dynamic_cast<OpConcat*>(otherbase);
      other->Append(buf);
    }
    first->Set(buf);
    return;
  }
  fassert(false, "Unknown M::Op implementation");
}

template <class M>
void Native<M>::DumpWrite(const std::vector<size_t>& bb) {
  auto& mfirst = kernels_.front()->GetMesh();
  if (mfirst.GetDump().size()) {
    std::string dumpformat = var.String["dumpformat"];
    if (dumpformat == "default") {
      dumpformat = "hdf";
    }

    if (dumpformat == "hdf") {
      // nop
    } else {
      P::DumpWrite(bb);
    }
  }
}

template <class M>
std::unique_ptr<DistrMesh<M>> CreateNative(
    MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var) {
  return std::unique_ptr<DistrMesh<M>>(new Native<M>(comm, kf, var));
}
