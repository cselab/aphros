// Created by Petr Karnakov on 13.12.2019
// Copyright 2019 ETH Zurich

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <memory>
#include <stdexcept>

#include "distr.h"
#include "dump/raw.h"
#include "dump/xmf.h"
#include "reduce.h"
#include "report.h"
#include "util/filesystem.h"
#include "util/format.h"

template <class M>
M DistrMesh<M>::CreateSharedMesh(
    MIdx block, Vect cellsize, int halos, bool isroot,
    const DomainInfo& domain) {
  const MIdx bs = domain.blocksize;
  const MIdx global_size = bs * domain.nblocks * domain.nprocs;
  const MIdx global_blocks = global_size / bs;
  const MIdx shared_size = bs * domain.nblocks;
  const Vect h(cellsize);
  const MIdx begin = block * bs;
  const Rect<Vect> rect(Vect(begin) * h, Vect(begin + shared_size) * h);
  const int id = M::Flags::GetIdFromBlock(block, global_blocks);
  const bool islead = true;
  M m(begin, shared_size, rect, halos, isroot, islead, global_size, id);
  m.flags.global_origin = Vect(0);
  m.flags.global_blocks = global_blocks;
  m.flags.block_length = h * Vect(bs);
  return m;
}

template <class M>
void DistrMesh<M>::MakeKernels(const std::vector<BlockInfoProxy>& proxies) {
  {
    fassert(!proxies.empty());
    auto& p = proxies.front();
    mshared_ = std::make_unique<M>(
        CreateSharedMesh(p.index, p.cellsize, p.halos, isroot_, domain_));
  }

  for (auto proxy : proxies) {
    std::unique_ptr<KernelMesh<M>> kernel(
        kernelfactory_.Make(var_mutable, proxy));
    auto& m = kernel->GetMesh();
    m.flags.comm = comm_;
    m.flags.is_periodic[0] = var.Int["hypre_periodic_x"];
    m.flags.is_periodic[1] = var.Int["hypre_periodic_y"];
    if (M::dim > 2) {
      m.flags.is_periodic[2] = var.Int["hypre_periodic_z"];
    }
    m.SetShared(mshared_.get());
    kernels_.emplace_back(kernel.release());
  }
}

template <size_t dim, class Map>
auto GetMIdx(const Map& map, std::string prefix) {
  generic::MIdx<dim> res;
  for (size_t d = 0; d < dim; ++d) {
    res[d] = map[prefix + GDir<dim>(d).letter()];
  }
  return res;
}

template <class M>
DistrMesh<M>::DistrMesh(
    MPI_Comm comm, const KernelMeshFactory<M>& kf, Vars& var0)
    : comm_(comm)
    , var(var0)
    , var_mutable(var0)
    , kernelfactory_(kf)
    , domain_(
          var.Int["hl"], GetMIdx<dim>(var.Int, "bs"),
          GetMIdx<dim>(var.Int, "p"), GetMIdx<dim>(var.Int, "b"),
          var.Double["extent"]) {}

template <class M>
DistrMesh<M>::~DistrMesh() {}

template <class M>
void DistrMesh<M>::RunKernels(const std::vector<size_t>& bb) {
#pragma omp parallel for schedule(dynamic, 1)
  for (size_t i = 0; i < bb.size(); ++i) {
    kernels_[bb[i]]->Run();
  }
}

template <class M>
void DistrMesh<M>::ClearComm(const std::vector<size_t>& bb) {
  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearComm();
  }
}

template <class M>
void DistrMesh<M>::ClearDump(const std::vector<size_t>& bb) {
  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearDump();
  }
}

template <class M>
void DistrMesh<M>::TimerReport(const std::vector<size_t>& bb) {
  auto& m = kernels_.front()->GetMesh();
  std::string fn = m.GetTimerReport();
  if (fn.length()) {
    std::ofstream out;
    out.open(fn);
    out << "mem=" << (sysinfo::GetMem() / double(1 << 20)) << " MB"
        << std::endl;
    ParseReport(multitimer_report_.GetMap(), out);
    multitimer_report_.Reset();
  }
  ClearTimerReport(bb);
}

template <class M>
void DistrMesh<M>::ClearTimerReport(const std::vector<size_t>& bb) {
  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearTimerReport();
  }
}

template <class M>
void DistrMesh<M>::Reduce(const std::vector<size_t>& bb) {
  const size_t nreqs = kernels_.front()->GetMesh().GetReduce().size();
  if (!nreqs) {
    return;
  }

  for (auto b : bb) {
    fassert_equal(kernels_[b]->GetMesh().GetReduce().size(), nreqs);
  }

  for (size_t i = 0; i < nreqs; ++i) {
    std::vector<RedOp*> blocks;
    for (auto b : bb) {
      blocks.push_back(kernels_[b]->GetMesh().GetReduce()[i].get());
    }
    ReduceSingleRequest(blocks);
  }

  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearReduce();
  }
}

template <class M>
void DistrMesh<M>::ReduceShared(const std::vector<size_t>&) {
  const size_t nreqs = mshared_->GetReduce().size();
  if (!nreqs) {
    return;
  }
  for (size_t i = 0; i < nreqs; ++i) {
    ReduceSingleRequest({mshared_->GetReduce()[i].get()});
  }
  mshared_->ClearReduce();
}

template <class M>
void DistrMesh<M>::ReduceToLead(const std::vector<size_t>& bb) {
  using ReductionScal = typename UReduce<Scal>::OpS;
  using ReductionScalInt = typename UReduce<Scal>::OpSI;
  using ReductionConcat = typename UReduce<Scal>::OpCat;
  auto& vfirst = kernels_.front()->GetMesh().GetReduceToLead();

  if (!vfirst.size()) {
    return;
  }

  for (auto b : bb) {
    fassert_equal(
        kernels_[b]->GetMesh().GetReduceToLead().size(), vfirst.size());
  }

  auto append = [&bb, this](auto* derived, auto& buf, size_t i) {
    for (auto b : bb) {
      const auto& v = kernels_[b]->GetMesh().GetReduceToLead();
      dynamic_cast<decltype(derived)>(v[i].get())->Append(buf);
    }
  };
  auto set = [&bb, this](auto* derived, auto& buf, size_t i) {
    for (auto b : bb) {
      const auto& m = kernels_[b]->GetMesh();
      if (m.IsLead()) {
        const auto& v = m.GetReduceToLead();
        dynamic_cast<decltype(derived)>(v[i].get())->Set(buf);
      }
    }
  };
  auto apply = [append, set](auto* derived, size_t i) {
    if (derived) {
      auto buf = derived->Neutral();
      append(derived, buf, i);
      set(derived, buf, i);
      return true;
    }
    return false;
  };

  for (size_t i = 0; i < vfirst.size(); ++i) {
    auto* request = vfirst[i].get();
    if (apply(dynamic_cast<ReductionScal*>(request), i) ||
        apply(dynamic_cast<ReductionScalInt*>(request), i) ||
        apply(dynamic_cast<ReductionConcat*>(request), i)) {
    } else {
      fassert(false, "ReduceToLead: Unknown M::Op implementation");
    }
  }

  // Clear reduce requests
  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearReduceToLead();
  }
}

template <class M>
void DistrMesh<M>::BcastFromLead(const std::vector<size_t>& bb) {
  using OpConcat = typename UReduce<Scal>::OpCat;
  auto& vfirst = kernels_.front()->GetMesh().GetBcastFromLead();

  if (!vfirst.size()) {
    return;
  }

  for (auto b : bb) {
    fassert_equal(
        kernels_[b]->GetMesh().GetBcastFromLead().size(), vfirst.size());
  }

  for (size_t i = 0; i < vfirst.size(); ++i) {
    if (OpConcat* o = dynamic_cast<OpConcat*>(vfirst[i].get())) {
      std::vector<char> buf = o->Neutral();

      // Read from lead block
      for (auto b : bb) {
        const auto& m = kernels_[b]->GetMesh();
        if (m.IsLead()) {
          const OpConcat* ob =
              static_cast<OpConcat*>(m.GetBcastFromLead()[i].get());
          ob->Append(buf);
        }
      }

      // Write to all blocks
      for (auto b : bb) {
        const auto& m = kernels_[b]->GetMesh();
        OpConcat* ob = static_cast<OpConcat*>(m.GetBcastFromLead()[i].get());
        ob->Set(buf);
      }
    } else {
      fassert(false, "BcastFromLead: Unknown M::Op instance");
    }
  }

  for (auto b : bb) {
    kernels_[b]->GetMesh().ClearBcastFromLead();
  }
}

template <class M>
void DistrMesh<M>::DumpWrite(const std::vector<size_t>& bb) {
  auto& mfirst = kernels_.front()->GetMesh();
  if (mfirst.GetDump().size()) {
    std::string dumpformat = var.String["dumpformat"];
    if (dumpformat == "default") {
      dumpformat = "raw";
    }
    if (dumpformat == "plain") {
      const auto& dumpfirst = mfirst.GetDump();
      for (size_t idump = 0; idump < dumpfirst.size(); ++idump) {
        const size_t nblocks = kernels_.size();
        std::vector<std::vector<Scal>> b_data(nblocks);
        std::vector<std::vector<MIdx>> b_origin(nblocks);
        std::vector<std::vector<MIdx>> b_size(nblocks);
        std::vector<std::unique_ptr<RedOp>> reduce_data;
        std::vector<std::unique_ptr<RedOp>> reduce_origin;
        std::vector<std::unique_ptr<RedOp>> reduce_size;
        size_t iblock = 0;
        for (auto b : bb) {
          auto& m = kernels_[b]->GetMesh();
          auto& data = b_data[iblock];
          auto& origin = b_origin[iblock];
          auto& size = b_size[iblock];
          data.reserve(m.GetInBlockCells().size() * kernels_.size());
          const auto* req = m.GetDump()[idump].first.get();
          if (auto* req_scal =
                  dynamic_cast<const typename M::CommRequestScal*>(req)) {
            for (auto c : m.Cells()) {
              data.push_back((*req_scal->field)[c]);
            }
          } else if (
              auto* req_vect =
                  dynamic_cast<const typename M::CommRequestVect*>(req)) {
            fassert(
                req_vect->d != -1,
                "Dump supports vector fields with only selected one component");
            for (auto c : m.Cells()) {
              data.push_back((*req_vect->field)[c][req_vect->d]);
            }
          } else {
            fassert(false);
          }
          origin.push_back(m.GetInBlockCells().GetBegin());
          size.push_back(m.GetInBlockCells().GetSize());
          using R = UReduce<Scal>;
          reduce_data.push_back(R::Make(&data, Reduction::concat));
          reduce_origin.push_back(R::Make(&origin, Reduction::concat));
          reduce_size.push_back(R::Make(&size, Reduction::concat));
          ++iblock;
        }
        ReduceSingleRequest(reduce_data);
        ReduceSingleRequest(reduce_origin);
        ReduceSingleRequest(reduce_size);
        if (isroot_) {
          const std::string path =
              GetDumpName(dumpfirst[idump].second, ".dat", frame_);
          typename M::IndexCells indexc(MIdx(0), mfirst.GetGlobalSize());
          typename M::BlockCells blockc(MIdx(0), mfirst.GetGlobalSize());
          auto& data = b_data[0]; // XXX assume root is first block
          auto& origin = b_origin[0];
          auto& size = b_size[0];
          const size_t nblocks_global = origin.size();
          fassert_equal(data.size(), indexc.size());
          fassert_equal(
              nblocks_global * mfirst.GetInBlockCells().size(), indexc.size());
          fassert_equal(size.size(), nblocks_global);
          FieldCell<Scal> fc_global(indexc);
          size_t i = 0;
          for (size_t b = 0; b < nblocks_global; ++b) {
            for (auto w : GBlock<IntIdx, dim>(origin[b], size[b])) {
              fc_global[indexc.GetIdx(w)] = data[i++];
            }
          }
          Dump(fc_global, indexc, blockc, path);
        }
      }
      if (isroot_) {
        std::cerr << "Dump " << frame_ << ": format=" << dumpformat
                  << std::endl;
      }
      ++frame_;
    } else if (dumpformat == "raw") {
      std::vector<MIdx> starts;
      std::vector<MIdx> sizes;
      for (auto b : bb) {
        const auto& m = kernels_[b]->GetMesh();
        const auto bc = m.GetInBlockCells();
        starts.push_back(bc.GetBegin());
        sizes.push_back(bc.GetSize());
      }

      const auto& dumpfirst = mfirst.GetDump();
      for (size_t idump = 0; idump < dumpfirst.size(); ++idump) {
        std::vector<std::vector<Scal>> data;
        const auto* req = dumpfirst[idump].first.get();
        if (auto* req_scal =
                dynamic_cast<const typename M::CommRequestScal*>(req)) {
          for (auto b : bb) {
            const auto& m = kernels_[b]->GetMesh();
            const auto bc = m.GetInBlockCells();
            data.emplace_back();
            auto& dat = data.back();
            dat.reserve(bc.size());
            auto& field = *static_cast<const typename M::CommRequestScal*>(
                               m.GetDump()[idump].first.get())
                               ->field;
            for (auto c : m.Cells()) {
              dat.push_back(field[c]);
            }
          }

        } else if (
            auto* req_vect =
                dynamic_cast<const typename M::CommRequestVect*>(req)) {
          fassert(
              req_vect->d != -1,
              "Dump only supports vector fields with selected one component");
          for (auto b : bb) {
            const auto& m = kernels_[b]->GetMesh();
            const auto bc = m.GetInBlockCells();
            data.emplace_back();
            auto& dat = data.back();
            dat.reserve(bc.size());
            auto& field = *static_cast<const typename M::CommRequestVect*>(
                               m.GetDump()[idump].first.get())
                               ->field;
            for (auto c : m.Cells()) {
              dat.push_back(field[c][req_vect->d]);
            }
          }
        } else {
          fassert(false);
        }

        const std::string path =
            GetDumpName(dumpfirst[idump].second, ".raw", frame_);

        using Xmf = dump::Xmf<Vect>;
        using Vect3 = generic::Vect<Scal, 3>;
        using MIdx3 = generic::MIdx<3>;
        using Xmf3 = dump::Xmf<Vect3>;

        typename Xmf3::Meta meta3;
        {
          auto meta = Xmf::GetMeta(mfirst);
          meta3.binpath = path;
          meta3.name = dumpfirst[idump].second;
          meta3.type = dump::Type::Float64;
          meta3.dimensions = MIdx3(1).max(MIdx3(meta.dimensions));
          meta3.count = meta3.dimensions;
          meta3.spacing = Vect3(meta.spacing.min());
        }

        dump::Raw<M>::Write(
            path, starts, sizes, data, mfirst.GetGlobalSize(), meta3.type,
            MpiWrapper(comm_));

        if (isroot_) {
          Xmf3::WriteXmf(util::SplitExt(path)[0] + ".xmf", meta3);
        }
      }
      ++frame_;
    } else {
      fassert(false, "Unknown dumpformat=" + dumpformat);
    }
  }
}

template <class M>
bool DistrMesh<M>::Pending(const std::vector<size_t>& bb) const {
  size_t pending = 0;
  for (auto b : bb) {
    if (kernels_[b]->GetMesh().Pending()) {
      ++pending;
    }
  }
  fassert(pending == 0 || pending == bb.size());
  return pending;
}

template <class M>
void DistrMesh<M>::ApplyNanFaces(const std::vector<size_t>& bb) {
  for (auto b : bb) {
    auto& m = kernels_[b]->GetMesh();
    for (auto& o : m.GetComm()) {
      if (auto* os = dynamic_cast<typename M::CommRequestScal*>(o.get())) {
        m.ApplyNanFaces(*os->field);
      } else if (
          auto* ov = dynamic_cast<typename M::CommRequestVect*>(o.get())) {
        m.ApplyNanFaces(*ov->field);
      } else {
        fassert(false, "Distr::Run(): unknown field type for nan");
      }
    }
  }
}

template <class M>
void DistrMesh<M>::Run() {
  if (var.Int["verbose_openmp"]) {
    ReportOpenmp();
  }
  while (true) {
    multitimer_all_.Push();
    multitimer_report_.Push();

    std::vector<size_t> bb;
    if (kernels_.front()->GetMesh().GetDump().size() > 0) {
      bb = TransferHalos(); // all blocks, sync communication
      // ApplyNanFaces(bb);
      DumpWrite(bb);
      ClearDump(bb);
      ClearComm(bb);
      mshared_->ClearComm();
      RunKernels(bb);
    } else {
      auto bbi = TransferHalos(true); // inner blocks, async communication
      ClearComm(bbi);
      mshared_->ClearComm();
      RunKernels(bbi);

      auto bbh = TransferHalos(false); // halo blocks, wait for communication
      ClearComm(bbh);
      RunKernels(bbh);

      bb = bbi;
      bb.insert(bb.end(), bbh.begin(), bbh.end());
    }

    stage_ += 1;

    // Print current stage name
    if (isroot_ && var.Int["verbose"]) {
      const auto& mf = kernels_.front()->GetMesh();
      std::cerr << "*** STAGE"
                << " #" << stage_ << " depth=" << mf.GetSuspender().GetDepth()
                << " " << mf.GetSuspender().GetNameSequence() << " ***"
                << std::endl;
      // Check stage name is the same for all kernels
      for (auto b : bb) {
        auto& m = kernels_[b]->GetMesh();
        fassert(
            m.GetSuspender().GetNameSequence() ==
                mf.GetSuspender().GetNameSequence(),
            util::Format(
                "Blocks {} and {} diverged to different stages {} and {}",
                m.GetId(), mf.GetId(), m.GetSuspender().GetNameSequence(),
                mf.GetSuspender().GetNameSequence()));
      }
    }

    Reduce(bb);
    ReduceToLead(bb);
    ReduceShared(bb);
    Scatter(bb);
    Bcast(bb);
    BcastFromLead(bb);

    const std::string nameseq =
        kernels_.front()->GetMesh().GetSuspender().GetNameSequence();
    multitimer_all_.Pop(nameseq);
    multitimer_report_.Pop(nameseq);
    TimerReport(bb);

    if (!Pending(bb)) {
      break;
    }
  }

  if (var.Int["verbose_time"]) {
    Report();
  }
}

template <class M>
void DistrMesh<M>::Report() {
  if (isroot_) {
    double total = 0.;
    for (auto e : multitimer_all_.GetMap()) {
      total += e.second;
    }

    if (var.Int["verbose_stages"]) {
      std::cerr << std::fixed;
      auto& map = multitimer_all_.GetMap();
      std::cerr << "mem=" << (sysinfo::GetMem() / double(1 << 20)) << " MB"
                << std::endl;
      ParseReport(map, std::cerr);
    }

    const auto& m = kernels_.front()->GetMesh();
    const size_t nc = m.GetGlobalSize().prod();
    const size_t nt = var.Int["max_step"];
    const size_t ni = var.Int("iter", 1);

    // Returns: hour, minute, second, millisecond
    auto get_hmsm = [](double t) {
      std::array<int, 4> r;
      r[0] = int(t / 3600);
      t -= r[0] * 3600;
      r[1] = int(t / 60);
      t -= r[1] * 60;
      r[2] = int(t);
      r[3] = int((t - int(t)) * 1000) % 1000;
      return r;
    };

    const auto hmsm = get_hmsm(total);
    std::cerr << util::Format(
                     "cells = {}\n"
                     "steps = {}\n"
                     "iters = {}\n"
                     "total = {:.3f} s = {:02d}:{:02d}:{:02d}.{:03d}\n"
                     "time/cell/iter = {:e} s\n",
                     nc, nt, ni, total, hmsm[0], hmsm[1], hmsm[2], hmsm[3],
                     total / (nc * ni))
              << std::endl;
  }
}

template <class M>
void DistrMesh<M>::ReportOpenmp() {
#ifdef _OPENMP
  if (isroot_) {
#pragma omp parallel
    {
#pragma omp single
      std::cerr << "OpenMP threads" << std::endl;
#pragma omp for ordered
      for (int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp ordered
        {
          std::cerr << "thread=" << std::setw(2) << omp_get_thread_num();
          std::cerr << std::setw(8) << " cpu=" << std::setw(2)
                    << sched_getcpu();
          std::cerr << std::endl;
        }
      }
    }
  }
#endif
}

template <class M>
void DistrMesh<M>::ReduceSingleRequest(
    const std::vector<std::unique_ptr<RedOp>>& blocks) {
  std::vector<RedOp*> raw;
  for (auto& p : blocks) {
    raw.push_back(p.get());
  }
  ReduceSingleRequest(raw);
}
