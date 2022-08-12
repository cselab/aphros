// Created by Petr Karnakov on 28.08.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <cassert>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

#include "distr.h"
#include "util/format.h"
#include "util/mpi.h"

template <class M>
void TransferParticlesLocal(
    const std::vector<std::unique_ptr<KernelMesh<M>>>& kernels) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;

  const size_t nreq = kernels.front()->GetMesh().GetCommPart().size();
  for (auto& kernel : kernels) {
    fassert_equal(kernel->GetMesh().GetCommPart().size(), nreq);
  }

  for (size_t q = 0; q < nreq; ++q) {
    auto& mroot = kernels.front()->GetMesh();
    auto& reqroot = mroot.GetCommPart()[q];
    const size_t nattr_scal = reqroot.attr_scal.size();
    const size_t nattr_vect = reqroot.attr_vect.size();
    // Check that the number of attributes is the same in all blocks.
    // Check that the size of attributes equals the number of particles
    for (auto& kernel : kernels) {
      auto& req = kernel->GetMesh().GetCommPart()[q];
      fassert(req.x);
      fassert(req.inner);
      fassert_equal(req.inner->size(), req.x->size());
      fassert_equal(req.attr_scal.size(), nattr_scal);
      fassert_equal(req.attr_vect.size(), nattr_vect);
      for (size_t a = 0; a < nattr_scal; ++a) {
        fassert(req.attr_scal[a]);
        fassert_equal(req.attr_scal[a]->size(), req.x->size());
      }
      for (size_t a = 0; a < nattr_vect; ++a) {
        fassert(req.attr_vect[a]);
        fassert_equal(req.attr_vect[a]->size(), req.x->size());
      }
    }
    // Temporary buffers for positions and attributes
    std::vector<std::vector<Vect>> tmp_x(kernels.size());
    std::vector<std::vector<std::vector<Scal>>> tmp_attr_scal(
        kernels.size(), std::vector<std::vector<Scal>>(nattr_scal));
    std::vector<std::vector<std::vector<Vect>>> tmp_attr_vect(
        kernels.size(), std::vector<std::vector<Vect>>(nattr_vect));
    const auto halorad =
        mroot.flags.particles_halo_radius * mroot.GetCellSize();
    // Traverse all inner particles and add them to blocks within `halorad`
    for (auto& kernel : kernels) {
      auto& m = kernel->GetMesh();
      auto& req = m.GetCommPart()[q];
      for (size_t i = 0; i < req.x->size(); ++i) {
        const Vect x = (*req.x)[i];
        if (!(*req.inner)[i]) {
          continue;
        }
        const Vect halobox_xm = x - halorad;
        const Vect halobox_xp = x + halorad;
        const MIdx block_min(
            ((halobox_xm - m.flags.global_origin) / m.flags.block_length)
                .floor());
        const MIdx block_max(
            ((halobox_xp - m.flags.global_origin) / m.flags.block_length)
                .floor());
        const GBlock<int, M::dim> blocks(
            block_min, block_max - block_min + MIdx(1));
        for (MIdx block : blocks) {
          Vect xtrans = x;
          for (size_t d : m.dirs) {
            if (m.flags.is_periodic[d]) {
              // canonical index of block from periodic conditions
              const int canon =
                  mod_positive(block[d], m.flags.global_blocks[d]);
              if (canon != block[d]) {
                const int image = (canon - block[d]) / m.flags.global_blocks[d];
                xtrans[d] += image * m.GetGlobalLength()[d];
                block[d] = canon;
              }
            }
          }
          if (MIdx(0) <= block && block < m.flags.global_blocks) {
            const size_t bdest = m.flags.GetIdFromBlock(block);
            fassert(bdest < kernels.size());
            tmp_x[bdest].push_back(xtrans);
            for (size_t a = 0; a < nattr_scal; ++a) {
              tmp_attr_scal[bdest][a].push_back((*req.attr_scal[a])[i]);
            }
            for (size_t a = 0; a < nattr_vect; ++a) {
              tmp_attr_vect[bdest][a].push_back((*req.attr_vect[a])[i]);
            }
          }
        }
      }
    }
    // Copy particles from tmp_* to the request buffers
    for (auto& kernel : kernels) {
      auto& m = kernel->GetMesh();
      auto& req = m.GetCommPart()[q];
      const size_t id = m.GetId();
      fassert(id < kernels.size());
      (*req.x) = tmp_x[id];
      req.inner->resize(req.x->size());
      for (size_t i = 0; i < req.x->size(); ++i) {
        (*req.inner)[i] = m.IsInnerPoint((*req.x)[i]);
      }
      for (size_t a = 0; a < nattr_scal; ++a) {
        (*req.attr_scal[a]) = tmp_attr_scal[id][a];
      }
      for (size_t a = 0; a < nattr_vect; ++a) {
        (*req.attr_vect[a]) = tmp_attr_vect[id][a];
      }
    }
  }

  for (auto& kernel : kernels) {
    kernel->GetMesh().ClearCommPart();
  }
}

#if USEFLAG(MPI)
template <class M>
void TransferParticlesMpi(
    const std::vector<std::unique_ptr<KernelMesh<M>>>& kernels, MPI_Comm comm,
    const std::function<int(int)>& rank_from_id) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;

  const size_t nreq = kernels.front()->GetMesh().GetCommPart().size();
  for (auto& kernel : kernels) {
    fassert_equal(kernel->GetMesh().GetCommPart().size(), nreq);
  }

  for (size_t q = 0; q < nreq; ++q) {
    auto& mroot = kernels.front()->GetMesh();
    auto& reqroot = mroot.GetCommPart()[q];
    const size_t nattr_scal = reqroot.attr_scal.size();
    const size_t nattr_vect = reqroot.attr_vect.size();
    // Check that the number of attributes is the same in all blocks.
    // Check that the size of attributes equals the number of particles
    for (auto& kernel : kernels) {
      auto& req = kernel->GetMesh().GetCommPart()[q];
      fassert(req.x);
      fassert(req.inner);
      fassert_equal(req.inner->size(), req.x->size());
      fassert_equal(req.attr_scal.size(), nattr_scal);
      fassert_equal(req.attr_vect.size(), nattr_vect);
      for (size_t a = 0; a < nattr_scal; ++a) {
        fassert(req.attr_scal[a]);
        fassert_equal(req.attr_scal[a]->size(), req.x->size());
      }
      for (size_t a = 0; a < nattr_vect; ++a) {
        fassert(req.attr_vect[a]);
        fassert_equal(req.attr_vect[a]->size(), req.x->size());
      }
    }
    // Mappings from block id to temporary buffers.
    std::unordered_map<int, std::vector<Vect>> tmp_x;
    std::unordered_map<int, std::vector<std::vector<Scal>>> tmp_attr_scal;
    std::unordered_map<int, std::vector<std::vector<Vect>>> tmp_attr_vect;
    const auto halorad =
        mroot.flags.particles_halo_radius * mroot.GetCellSize();
    // Traverse all local inner particles and add to blocks within `halorad`
    for (auto& kernel : kernels) {
      auto& m = kernel->GetMesh();
      auto& req = m.GetCommPart()[q]; // Request on block `b`
      for (size_t i = 0; i < req.x->size(); ++i) {
        const Vect x = (*req.x)[i];
        if (!(*req.inner)[i]) {
          continue;
        }
        const Vect halobox_xm = x - halorad;
        const Vect halobox_xp = x + halorad;
        const MIdx block_min(
            ((halobox_xm - m.flags.global_origin) / m.flags.block_length)
                .floor());
        const MIdx block_max(
            ((halobox_xp - m.flags.global_origin) / m.flags.block_length)
                .floor());
        const GBlock<int, M::dim> blocks(
            block_min, block_max - block_min + MIdx(1));
        for (MIdx block : blocks) {
          Vect xtrans = x;
          for (size_t d : m.dirs) {
            if (m.flags.is_periodic[d]) {
              // Canonical index of block from periodic conditions
              const int canon =
                  mod_positive(block[d], m.flags.global_blocks[d]);
              if (canon != block[d]) {
                const int image = (canon - block[d]) / m.flags.global_blocks[d];
                xtrans[d] += image * m.GetGlobalLength()[d];
                block[d] = canon;
              }
            }
          }
          if (MIdx(0) <= block && block < m.flags.global_blocks) {
            const int dest = m.flags.GetIdFromBlock(block);
            tmp_x[dest].push_back(xtrans);
            tmp_attr_scal[dest].resize(nattr_scal);
            for (size_t a = 0; a < nattr_scal; ++a) {
              tmp_attr_scal[dest][a].push_back((*req.attr_scal[a])[i]);
            }
            tmp_attr_vect[dest].resize(nattr_vect);
            for (size_t a = 0; a < nattr_vect; ++a) {
              tmp_attr_vect[dest][a].push_back((*req.attr_vect[a])[i]);
            }
          }
        }
      }
    }
    MpiWrapper mpi(comm);

    const int myrank = mpi.GetCommRank();

    // Message for remote rank
    struct Msg {
      std::vector<char> buf; // Serial buffer
      MPI_Request req;
    };
    std::unordered_map<int, Msg> rank_to_msg;

    auto serialize = [](std::vector<char>& buf, const auto elem) {
      const char* ptr = reinterpret_cast<const char*>(&elem);
      buf.insert(buf.end(), ptr, ptr + sizeof(elem));
    };

    // Serialize send buffers for remote blocks to `rank_to_buf`
    for (auto& p : tmp_x) {
      const int id = p.first;
      const int rank = rank_from_id(id);
      if (rank == myrank) { // skip local blocks
        continue;
      }
      auto& buf = rank_to_msg[rank].buf;
      // Remote block id
      serialize(buf, id);
      // Number of particles
      serialize(buf, size_t(tmp_x[id].size()));
      // Number of scalar attributes
      serialize(buf, size_t(tmp_attr_scal[id].size()));
      // Number of vector attributes
      serialize(buf, size_t(tmp_attr_vect[id].size()));
      // Positions
      for (const auto& x : tmp_x[id]) {
        serialize(buf, x);
      }
      // Scalar attributes
      for (const auto& v : tmp_attr_scal[id]) {
        for (const auto& a : v) {
          serialize(buf, a);
        }
      }
      // Vector attributes
      for (const auto& v : tmp_attr_vect[id]) {
        for (const auto& a : v) {
          serialize(buf, a);
        }
      }
    }

    std::vector<int> rank_to_count(mpi.GetCommSize());

    // Send messages
    const int tag = 276;
    for (auto& p : rank_to_msg) {
      const int rank = p.first;
      Msg& msg = p.second;
      rank_to_count[rank] = 1;
      MPICALL(MPI_Isend(
          msg.buf.data(), msg.buf.size(), MPI_CHAR, rank, tag, comm, &msg.req));
    }

    // Compute the number of messages to receive
    // from other ranks in `rank_to_count[myrank]`
    MPICALL(MPI_Allreduce(
        MPI_IN_PLACE, rank_to_count.data(), rank_to_count.size(), MPI_INT,
        MPI_SUM, comm));

    // Initialize the request buffers with local particles
    for (auto& kernel : kernels) {
      auto& m = kernel->GetMesh();
      auto& req = m.GetCommPart()[q];
      const int id = m.GetId();
      if (tmp_x.count(id)) {
        (*req.x) = tmp_x.at(id);
        req.inner->resize(req.x->size());
        for (size_t i = 0; i < req.x->size(); ++i) {
          (*req.inner)[i] = m.IsInnerPoint((*req.x)[i]);
        }
        for (size_t a = 0; a < nattr_scal; ++a) {
          (*req.attr_scal[a]) = tmp_attr_scal.at(id)[a];
        }
        for (size_t a = 0; a < nattr_vect; ++a) {
          (*req.attr_vect[a]) = tmp_attr_vect.at(id)[a];
        }
      } else {
        req.x->clear();
        req.inner->clear();
        for (size_t a = 0; a < nattr_scal; ++a) {
          req.attr_scal[a]->clear();
        }
        for (size_t a = 0; a < nattr_vect; ++a) {
          req.attr_vect[a]->clear();
        }
      }
    }

    // Receive messages and append received particles to the request buffers
    for (int cnt = 0; cnt < rank_to_count[myrank]; ++cnt) {
      MPI_Status status;
      MPICALL(MPI_Probe(MPI_ANY_SOURCE, tag, comm, &status));
      int count; // Message length
      MPICALL(MPI_Get_count(&status, MPI_CHAR, &count));
      const int rank = status.MPI_SOURCE;
      std::vector<char> buf(count); // Serial buffer
      MPICALL(MPI_Recv(
          buf.data(), buf.size(), MPI_CHAR, rank, tag, comm,
          MPI_STATUS_IGNORE));
      // Current position in serial buffer
      auto pos = buf.begin();
      // Reads element from serial buffer to `elem`
      auto deserialize = [&buf, &pos](auto& elem) {
        char* ptr = reinterpret_cast<char*>(&elem);
        std::copy(pos, pos + sizeof(elem), ptr);
        pos += sizeof(elem);
      };

      std::unordered_map<int, KernelMesh<M>*> id_to_kernel;
      for (auto& kernel : kernels) {
        id_to_kernel[kernel->GetMesh().GetId()] = kernel.get();
      }

      // The serial buffer is a concatenation of particle data for all blocks
      // Traverse all of them and append to the request buffers.
      while (pos < buf.end()) {
        int id; // Block id
        deserialize(id);
        size_t npart; // Number of particles
        deserialize(npart);
        size_t nscal; // Number of scalar attributes
        deserialize(nscal);
        size_t nvect; // Number of vector attributes
        deserialize(nvect);

        fassert_equal(rank_from_id(id), myrank);
        fassert_equal(nscal, nattr_scal);
        fassert_equal(nvect, nattr_vect);
        auto& m = id_to_kernel.at(id)->GetMesh();
        auto& req = m.GetCommPart()[q];

        // Positions
        auto& x = *req.x;
        const size_t oldsize = x.size();
        const size_t newsize = x.size() + npart;
        x.resize(newsize);
        for (size_t i = oldsize; i < newsize; ++i) {
          deserialize(x[i]);
        }
        auto& inner = *req.inner;
        inner.resize(newsize);
        for (size_t i = oldsize; i < newsize; ++i) {
          inner[i] = m.IsInnerPoint(x[i]);
        }
        for (size_t a = 0; a < nattr_scal; ++a) {
          auto& v = *req.attr_scal[a];
          v.resize(newsize);
          for (size_t i = oldsize; i < newsize; ++i) {
            deserialize(v[i]);
          }
        }
        for (size_t a = 0; a < nattr_vect; ++a) {
          auto& v = *req.attr_vect[a];
          v.resize(newsize);
          for (size_t i = oldsize; i < newsize; ++i) {
            deserialize(v[i]);
          }
        }
      }
      // Check that particle data is not cut in the middle
      fassert(pos == buf.end());
    }

    // Wait for sends
    for (auto& p : rank_to_msg) {
      MPI_Wait(&p.second.req, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(comm);
  }

  for (auto& kernel : kernels) {
    kernel->GetMesh().ClearCommPart();
  }
}
#endif

template <class M>
void TransferParticles(
    const std::vector<std::unique_ptr<KernelMesh<M>>>& kernels, MPI_Comm comm,
    const std::function<int(int)>& rank_from_id) {
#if USEFLAG(MPI)
  TransferParticlesMpi(kernels, comm, rank_from_id);
#else
  (void)comm;
  (void)rank_from_id;
  TransferParticlesLocal(kernels);
#endif
}
