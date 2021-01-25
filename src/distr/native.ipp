// Created by Petr Karnakov on 02.01.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <mpi.h>
#include <cassert>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>

#include "comm_manager.h"
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

  int commsize_;
  int commrank_;
  typename CommManager<dim>::Tasks tasks_;
  struct HaloTmp {
    MPI_Request req;
    size_t cnt = 0;
    std::vector<Scal> buf;
  };
  std::map<int, HaloTmp> tmp_send_;
  std::map<int, HaloTmp> tmp_recv_;
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
    p.isroot = (p.index == MIdx(0) && isroot_);
    p.islead = (p.index == MIdx(0));
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
    const std::array<bool, dim> is_periodic{true, true, true};
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

  auto& task = tasks_.full_two;

  size_t nscal = 0; // number of scalar fields to transfer
  for (size_t i = 0; i < vcr.size(); ++i) {
    if (dynamic_cast<typename M::CommRequestScal*>(vcr[i].get())) {
      nscal += 1;
    }
    if (auto crd = dynamic_cast<typename M::CommRequestVect*>(vcr[i].get())) {
      nscal += (crd->d == -1 ? dim : 1);
    }
  }

  // Clear buffers
  for (auto& p : tmp_recv_) {
    auto& tmp = p.second;
    tmp.buf.resize(0);
    tmp.cnt = 0;
  }
  for (auto& p : tmp_send_) {
    auto& tmp = p.second;
    tmp.buf.resize(0);
    tmp.cnt = 0;
  }

  // Determine the size of receive buffer
  for (auto& p : task.recv) {
    auto& rank = p.first;
    auto& cells = p.second;
    tmp_recv_[rank].cnt += cells.size() * nscal;
  }

  auto type = sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT;
  const int tag = 0;

  // Receive
  for (auto& p : tmp_recv_) {
    auto& rank = p.first;
    auto& tmp = p.second;
    tmp.buf.resize(tmp.cnt);
    MPI_Irecv(tmp.buf.data(), tmp.cnt, type, rank, tag, comm_, &tmp.req);
  }
  // Determine the size of send buffer
  for (auto& p : task.send) {
    auto& rank = p.first;
    auto& cells = p.second;
    tmp_send_[rank].cnt += cells.size() * nscal;
  }
  // Reserve send buffer
  for (auto& p : tmp_send_) {
    auto& tmp = p.second;
    tmp.buf.reserve(tmp.cnt);
  }
  // Collect data to send
  for (size_t i = 0; i < vcr.size(); ++i) {
    if (dynamic_cast<typename M::CommRequestScal*>(vcr[i].get())) {
      std::vector<FieldCell<Scal>*> fields;
      for (auto& k : kernels_) {
        const auto* cr = static_cast<typename M::CommRequestScal*>(
            k->GetMesh().GetComm()[i].get());
        fields.push_back(cr->field);
      }
      for (auto& p : task.send) {
        auto& rank = p.first;
        auto& cells = p.second;
        auto& tmp = tmp_send_[rank];
        for (auto bc : cells) {
          tmp.buf.push_back((*fields[bc.block])[bc.cell]);
        }
      }
    }
    if (auto crd = dynamic_cast<typename M::CommRequestVect*>(vcr[i].get())) {
      std::vector<FieldCell<Vect>*> fields;
      for (auto& k : kernels_) {
        const auto* cr = static_cast<typename M::CommRequestVect*>(
            k->GetMesh().GetComm()[i].get());
        fields.push_back(cr->field);
      }
      for (auto& p : task.send) {
        auto& rank = p.first;
        auto& cells = p.second;
        auto& tmp = tmp_send_[rank];
        for (auto bc : cells) {
          if (crd->d == -1) {
            for (size_t d = 0; d < dim; ++d) {
              tmp.buf.push_back((*fields[bc.block])[bc.cell][d]);
            }
          } else {
            tmp.buf.push_back((*fields[bc.block])[bc.cell][crd->d]);
          }
        }
      }
    }
  }
  // Send
  for (auto& p : tmp_send_) {
    auto& rank = p.first;
    auto& tmp = p.second;
    tmp.buf.resize(tmp.cnt);
    MPI_Isend(tmp.buf.data(), tmp.cnt, type, rank, tag, comm_, &tmp.req);
  }
  // Wait for send and receive
  for (auto& p : tmp_send_) {
    MPI_Wait(&p.second.req, MPI_STATUS_IGNORE);
  }
  for (auto& p : tmp_recv_) {
    MPI_Wait(&p.second.req, MPI_STATUS_IGNORE);
  }
  // Copy received data to fields, use `cnt` for position
  for (auto& p : tmp_recv_) {
    p.second.cnt = 0;
  }
  for (size_t i = 0; i < vcr.size(); ++i) {
    if (dynamic_cast<typename M::CommRequestScal*>(vcr[i].get())) {
      std::vector<FieldCell<Scal>*> fields;
      for (auto& k : kernels_) {
        const auto* cr = static_cast<typename M::CommRequestScal*>(
            k->GetMesh().GetComm()[i].get());
        fields.push_back(cr->field);
      }
      for (auto& p : task.recv) {
        auto& rank = p.first;
        auto& cells = p.second;
        auto& tmp = tmp_recv_[rank];
        for (auto bc : cells) {
          (*fields[bc.block])[bc.cell] = tmp.buf[tmp.cnt++];
        }
      }
    }
    if (auto crd = dynamic_cast<typename M::CommRequestVect*>(vcr[i].get())) {
      std::vector<FieldCell<Vect>*> fields;
      for (auto& k : kernels_) {
        const auto* cr = static_cast<typename M::CommRequestVect*>(
            k->GetMesh().GetComm()[i].get());
        fields.push_back(cr->field);
      }
      for (auto& p : task.recv) {
        auto& rank = p.first;
        auto& cells = p.second;
        auto& tmp = tmp_recv_[rank];
        for (auto bc : cells) {
          if (crd->d == -1) {
            for (size_t d = 0; d < dim; ++d) {
              (*fields[bc.block])[bc.cell][d] = tmp.buf[tmp.cnt++];
            }
          } else {
            (*fields[bc.block])[bc.cell][crd->d] = tmp.buf[tmp.cnt++];
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

