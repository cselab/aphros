// Created by Petr Karnakov on 02.01.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <mpi.h>
#include <cassert>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>

#include "distr.h"
#include "dump/dumper.h"
#include "util/format.h"

#include "native.h"

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

  std::map<MIdx, size_t> midx_to_kernel_;
  int commsize_;
};

template <class M>
Native<M>::Native(MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var_)
    : DistrMesh<M>(comm, kf, var_) {
  MPI_Comm_size(comm, &commsize_);
  int rank;
  MPI_Comm_rank(comm, &rank);
  isroot_ = (0 == rank); // XXX: overwrite isroot_

  if (commsize_ != nprocs_.prod()) {
    throw std::runtime_error(util::Format(
        "Number of MPI tasks {} does not match the number of subdomains {}",
        commsize_, nprocs_));
  }

  std::vector<BlockInfoProxy> proxies;
  GIndex<size_t, dim> procs(nprocs_);
  GIndex<size_t, dim> blocks(nblocks_);
  const MIdx wproc = procs.GetMIdx(rank);
  for (auto ib : blocks.Range()) {
    BlockInfoProxy p;
    p.index = nblocks_ * wproc + blocks.GetMIdx(ib);
    p.globalsize = nprocs_ * nblocks_ * blocksize_;
    p.cellsize = Vect(extent_ / p.globalsize.max());
    p.blocksize = blocksize_;
    p.halos = halos_;
    p.isroot = (p.index == MIdx(0) && isroot_);
    p.islead = (p.index == MIdx(0));
    proxies.push_back(p);
  }
  comm_ = comm; // XXX: overwrite comm_
  this->MakeKernels(proxies);
}

template <class M>
Native<M>::~Native() = default;

template <class M>
auto Native<M>::TransferHalos(bool inner) -> std::vector<size_t> {
  std::vector<size_t> bb;
  std::iota(bb.begin(), bb.end(), 0);
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
      std::vector<char> r = o->Neutral(); // buffer

      if (isroot_) {
        // read from root block
        for (auto b : bb) {
          auto& m = kernels_[b]->GetMesh();
          if (m.IsRoot()) {
            auto& v = m.GetBcast();
            OpConcat* ob = dynamic_cast<OpConcat*>(v[i].get());
            ob->Append(r);
          }
        }
      }

      int s = r.size(); // size

      // broadcast size
      MPI_Bcast(&s, 1, MPI_INT, 0, comm_);

      // resize
      r.resize(s);

      // broadcast data
      MPI_Bcast(r.data(), r.size(), MPI_CHAR, 0, comm_);

      // write to all blocks
      for (auto b : bb) {
        auto& m = kernels_[b]->GetMesh();
        auto& v = m.GetBcast();
        OpConcat* ob = dynamic_cast<OpConcat*>(v[i].get());
        ob->Set(r);
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
    int recvcount;
    int sizes_recvcount;
    std::vector<Scal> rbuf;
    std::vector<int> sizes_rbuf;

    MPI_Datatype mscal = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);

    if (isroot_) {
      std::vector<Scal> buf;
      std::vector<int> dis(commsize_, 0);
      std::vector<int> cnt(commsize_, 0);
      std::vector<int> sizes_buf;
      std::vector<int> sizes_dis(commsize_, 0);
      std::vector<int> sizes_cnt(commsize_, 0);
      // find root block
      for (auto b : bb) {
        auto& m = kernels_[b]->GetMesh();
        if (m.IsRoot()) {
          auto& req = m.GetScatter()[q];
          size_t i = 0;
          // concatenate data for all blocks in buf
          for (int rank = 0; rank < commsize_; ++rank) {
            dis[rank] = buf.size();
            sizes_dis[rank] = sizes_buf.size();
            // XXX assuming the same number of blocks on all ranks
            for (size_t k = 0; k < bb.size(); ++k) {
              auto& v = (*req.first)[i];
              buf.insert(buf.end(), v.begin(), v.end());
              sizes_buf.push_back(v.size());
              ++i;
            }
            sizes_cnt[rank] = sizes_buf.size() - sizes_dis[rank];
            cnt[rank] = buf.size() - dis[rank];
          }
        }
      }
      // data recvcount
      MPI_Scatter(cnt.data(), 1, MPI_INT, &recvcount, 1, MPI_INT, 0, comm_);
      rbuf.resize(recvcount);
      // data
      MPI_Scatterv(
          buf.data(), cnt.data(), dis.data(), mscal, rbuf.data(), recvcount,
          mscal, 0, comm_);
      // sizes recvcount
      MPI_Scatter(
          sizes_cnt.data(), 1, MPI_INT, &sizes_recvcount, 1, MPI_INT, 0, comm_);
      sizes_rbuf.resize(sizes_recvcount);
      // sizes
      MPI_Scatterv(
          sizes_buf.data(), sizes_cnt.data(), sizes_dis.data(), MPI_INT,
          sizes_rbuf.data(), sizes_recvcount, MPI_INT, 0, comm_);
    } else {
      // data recvcount
      MPI_Scatter(nullptr, 0, MPI_INT, &recvcount, 1, MPI_INT, 0, comm_);
      rbuf.resize(recvcount);
      // data
      MPI_Scatterv(
          nullptr, nullptr, nullptr, mscal, rbuf.data(), recvcount, mscal, 0,
          comm_);
      // sizes recvcount
      MPI_Scatter(nullptr, 0, MPI_INT, &sizes_recvcount, 1, MPI_INT, 0, comm_);
      sizes_rbuf.resize(sizes_recvcount);
      // sizes
      MPI_Scatterv(
          nullptr, nullptr, nullptr, MPI_INT, sizes_rbuf.data(),
          sizes_recvcount, MPI_INT, 0, comm_);
    }

    // write to blocks on current rank
    size_t off = 0;
    for (size_t k = 0; k < bb.size(); ++k) {
      auto& v = *kernels_[bb[k]]->GetMesh().GetScatter()[q].second;
      v = std::vector<Scal>(
          rbuf.data() + off, rbuf.data() + off + sizes_rbuf[k]);
      off += sizes_rbuf[k];
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

    MPI_Op mpiop;
    if (dynamic_cast<typename UReduce<Scal>::OpSum*>(first)) {
      mpiop = MPI_SUM;
    } else if (dynamic_cast<typename UReduce<Scal>::OpProd*>(first)) {
      mpiop = MPI_PROD;
    } else if (dynamic_cast<typename UReduce<Scal>::OpMax*>(first)) {
      mpiop = MPI_MAX;
    } else if (dynamic_cast<typename UReduce<Scal>::OpMin*>(first)) {
      mpiop = MPI_MIN;
    } else {
      fassert(false, "Unknown reduction");
    }
    MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);

    // Reduce over ranks
    MPI_Allreduce(MPI_IN_PLACE, &buf, 1, mt, mpiop, comm_);

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

    MPI_Op mpiop;
    if (dynamic_cast<typename UReduce<Scal>::OpMinloc*>(first)) {
      mpiop = MPI_MINLOC;
    } else if (dynamic_cast<typename UReduce<Scal>::OpMaxloc*>(first)) {
      mpiop = MPI_MAXLOC;
    } else {
      fassert(false, "Unknown reduction");
    }

    MPI_Datatype mt = (sizeof(Scal) == 8 ? MPI_DOUBLE_INT : MPI_FLOAT_INT);

    // Reduce over all ranks
    MPI_Allreduce(MPI_IN_PLACE, &buf, 1, mt, mpiop, comm_);

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

    int bufsize = buf.size();

    if (isroot_) {
      std::vector<int> sizes(commsize_);

      MPI_Gather(&bufsize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, comm_);

      int size_all = 0;
      std::vector<int> offsets = {0};
      for (auto& q : sizes) {
        size_all += q;
        offsets.push_back(offsets.back() + q);
      }
      offsets.pop_back();
      assert(sizes.size() == offsets.size());

      std::vector<char> buf_all(size_all);

      MPI_Gatherv(
          buf.data(), buf.size(), MPI_CHAR, buf_all.data(), sizes.data(),
          offsets.data(), MPI_CHAR, 0, comm_);

      // Write results to root block only (assume first)
      // FIXME get IsRoot() from kernel_ after using std::vector for kernels
      first->Set(buf_all);
    } else {
      MPI_Gather(&bufsize, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm_);

      MPI_Gatherv(
          buf.data(), buf.size(), MPI_CHAR, nullptr, nullptr, nullptr, MPI_CHAR,
          0, comm_);
    }
    return;
  }
  fassert(false, "Unknown M::Op implementation");
}

template <class M>
void Native<M>::DumpWrite(const std::vector<size_t>& bb) {
  // nop
}

template <class M>
std::unique_ptr<DistrMesh<M>> CreateNative(
    MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var) {
  return std::unique_ptr<DistrMesh<M>>(new Native<M>(comm, kf, var));
}

