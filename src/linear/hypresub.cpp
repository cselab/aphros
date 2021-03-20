// Created by Petr Karnakov on 27.11.2019
// Copyright 2019 ETH Zurich

#include <cassert>
#include <chrono>
#include <iostream>
#include <map>
#include <stdexcept>
#include <thread>
#include <utility>

#include "hypresub.h"
#include "util/histogram.h"
#include "util/subcomm.h"

//#define DEB(x) x
#define DEB(x)

#define MSI MPI_STATUS_IGNORE

DEB(static std::ostream& operator<<(
    std::ostream& out, const HypreSub::MIdx& v) {
  out << "(";
  for (auto& a : v) {
    out << a << ",";
  }
  out << ")";
  return out;
})

DEB(template <class T> static std::ostream& operator<<(
    std::ostream& out, const std::vector<T>& v) {
  out << "[";
  for (auto& a : v) {
    out << a << " ";
  }
  out << "]";
  return out;
})

DEB(static std::ostream& operator<<(
    std::ostream& out, const HypreSub::Block& b) {
  out << "b.l=" << b.l << " b.u=" << b.u
      << " b.stencil.size()=" << b.stencil.size()
      << " b.a.size()=" << b.a->size() << " b.r.size()=" << b.r->size()
      << " b.x.size()=" << b.x->size();
  return out;
})

#define EV(x) (#x) << "=" << (x) << " "

struct HypreSub::Imp {
  enum class Cmd {
    construct,
    destruct,
    update,
    solve,
    get_residual,
    get_iter,
    exit
  };

  // Block with data, only default-constructable once
  class BlockBuffer {
   public:
    BlockBuffer() {
      b.a = &a;
      b.r = &r;
      b.x = &x;
    }
    BlockBuffer(const BlockBuffer&) = delete;
    BlockBuffer(BlockBuffer&&) = delete;
    BlockBuffer& operator=(const BlockBuffer&) = delete;
    BlockBuffer& operator=(BlockBuffer&&) = delete;
    // Copy data from
    void Update(const Block& b0) {
      b = b0;
      a = *b0.a;
      r = *b0.r;
      x = *b0.x;
      b.a = &a;
      b.r = &r;
      b.x = &x;
    }
    // Copy data from
    void Update(const BlockBuffer& buf) {
      a = buf.a;
      r = buf.r;
      x = buf.x;
    }
    const Block& Get() const {
      return b;
    }

   private:
    Block b; // b.a,b.r,b.x always point to a,r,x
    std::vector<Scal> a;
    std::vector<Scal> r;
    std::vector<Scal> x;
  };

  static std::vector<Block> GetBlocks(const std::vector<BlockBuffer>& vbuf) {
    std::vector<Block> bb;
    for (auto& buf : vbuf) {
      bb.push_back(buf.Get());
    }
    return bb;
  }
  static std::vector<BlockBuffer> GetBlockBuffers(
      const std::vector<Block>& bb) {
    std::vector<BlockBuffer> vbuf(bb.size());
    for (size_t i = 0; i < bb.size(); ++i) {
      vbuf[i].Update(bb[i]);
    }
    return vbuf;
  }

  class Instance {
   public:
    Instance() = delete;
    Instance(const Instance&) = delete;
    Instance(Instance&&) = delete;
    Instance& operator=(const Instance&) = delete;
    Instance& operator=(Instance&&) = delete;
    Instance(std::vector<BlockBuffer>&& vbuf0, MIdx gs, MIdx per)
        : vbuf(std::move(vbuf0))
        , hypre(state.comm, Imp::GetBlocks(vbuf), gs, per) {}
    ~Instance() {}
    void Update(const std::vector<Block>& src) {
      assert(src.size() == vbuf.size());
      for (size_t i = 0; i < vbuf.size(); ++i) {
        vbuf[i].Update(src[i]);
      }
    }
    void Update(const std::vector<BlockBuffer>& src) {
      assert(src.size() == vbuf.size());
      for (size_t i = 0; i < vbuf.size(); ++i) {
        vbuf[i].Update(src[i]);
      }
    }
    Hypre& GetHypre() {
      return hypre;
    }
    std::vector<Block> GetBlocks() const {
      return Imp::GetBlocks(vbuf);
    }

   private:
    std::vector<BlockBuffer> vbuf;
    Hypre hypre;
  };
  struct ServerState {
    std::vector<int> v;
    MPI_Comm comm;
    MPI_Comm commsub;
    int rank;
    int ranksub;
    int size;
    int sizesub;
    std::map<int, Instance> minst;
  };

  static ServerState state;
  static auto GetMpiScal() {
    return sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT;
  }
  static constexpr int tag = 1;

  // Returns partition of blocks for one rank
  std::vector<Block> GetPart(int ranksub) {
    std::vector<Block> bb;
    for (auto i : partidx_[ranksub]) {
      bb.push_back(bbarg_[i]);
    }
    return bb;
  }
  static Cmd GetCmd(std::string s) {
#define GETCMD(x) \
  if (s == #x) return Cmd::x;
    GETCMD(construct);
    GETCMD(destruct);
    GETCMD(update);
    GETCMD(solve);
    GETCMD(get_residual);
    GETCMD(get_iter);
    GETCMD(exit);
    fassert(false, "GetCmd: unknown name '" + s + "'");
  }
  static std::string GetString(Cmd cmd) {
#define GETSTR(x) \
  if (cmd == Cmd::x) return #x;
    GETSTR(construct);
    GETSTR(destruct);
    GETSTR(update);
    GETSTR(solve);
    GETSTR(get_residual);
    GETSTR(get_iter);
    GETSTR(exit);
    fassert(false, "GetString: unknown cmd");
  }
  static void InitServer(MPI_Comm comm, MPI_Comm commsub) {
    state.comm = comm;
    state.commsub = commsub;
    MPI_Comm_rank(comm, &state.rank);
    MPI_Comm_rank(commsub, &state.ranksub);
    MPI_Comm_size(comm, &state.size);
    MPI_Comm_size(commsub, &state.sizesub);
  }
  static void RunServer(Sampler& samp) {
    samp.SeedSample();
    while (true) {
      auto comm = state.commsub;
      Cmd cmd;
      int id;
      int root = 0;

      Notify(comm);
      Recv(cmd, root, comm);
      Recv(id, root, comm);

      if (cmd == Cmd::construct) {
        samp.SeedSample();
        MIdx gs;
        MIdx per;
        auto vbuf = RecvBlocks(root, comm);
        Recv(gs, root, comm);
        Recv(per, root, comm);
        DEB(std::cerr << "recv construct " << EV(id) << EV(state.rank)
                      << GetBlocks(vbuf) << std::endl;)
        state.minst.emplace(
            std::piecewise_construct, std::forward_as_tuple(id),
            std::forward_as_tuple(std::move(vbuf), gs, per));
        samp.CollectSample("RunServer::Construct");
      } else if (cmd == Cmd::update) {
        samp.SeedSample();
        DEB(std::cerr << "recv update " << EV(id) << std::endl;)
        auto vbuf = RecvBlocks(root, comm);
        auto& inst = state.minst.at(id);
        inst.Update(vbuf);
        inst.GetHypre().Update();
        samp.CollectSample("RunServer::Update");
      } else if (cmd == Cmd::solve) {
        samp.SeedSample();
        Scal tol;
        int print;
        std::string solver;
        int maxiter;
        Recv(tol, root, comm);
        Recv(print, root, comm);
        Recv(solver, root, comm);
        Recv(maxiter, root, comm);
        auto& inst = state.minst.at(id);
        DEB(std::cerr << "recv solve " << EV(id) << EV(solver) << std::endl;)
        inst.GetHypre().Solve(tol, print, solver, maxiter);
        SendBlocks(inst.GetBlocks(), root, comm);
        samp.CollectSample("RunServer::Solve");
      } else if (cmd == Cmd::destruct) {
        samp.SeedSample();
        DEB(std::cerr << "recv destruct " << EV(id) << EV(state.rank)
                      << std::endl;)
        state.minst.erase(id);
        samp.CollectSample("RunServer::Destruct");
      } else if (cmd == Cmd::exit) {
        DEB(std::cerr << "recv exit " << EV(id) << EV(state.rank) << std::endl;)
        break;
      }
    }
    samp.CollectSample("RunServer::Loop");

    fassert(
        state.minst.size() <= 0,
        "RunServer: Cmd::exit received on rank " + std::to_string(state.rank) +
            ", but state.minst still containts " +
            std::to_string(state.minst.size()) + " instances");
  }
  static void Notify(MPI_Comm comm) {
    MPI_Request req;
    MPI_Ibarrier(comm, &req);
    if (1) {
      while (true) {
        int flag = 0;
        MPI_Test(&req, &flag, MSI);
        if (flag) {
          break;
        }
        // sched_yield();
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
      }
    } else {
      MPI_Wait(&req, MSI);
    }
  }
  static void StopServer() {
    Notify(state.commsub);
    for (int rank = 1; rank < state.sizesub; ++rank) {
      Send(Cmd::exit, rank);
      Send(-1, rank, state.commsub); // dummy id, expected by RunServer
    }
  }
  static void Send(const MIdx& w, int rank, MPI_Comm comm) {
    MPI_Send(w.data(), dim, MPI_INT, rank, tag, comm);
  }
  static void Recv(MIdx& w, int rank, MPI_Comm comm) {
    MPI_Recv(w.data(), dim, MPI_INT, rank, tag, comm, MSI);
  }
  static void Send(const int& a, int rank, MPI_Comm comm) {
    MPI_Send(&a, 1, MPI_INT, rank, tag, comm);
  }
  static void Recv(int& a, int rank, MPI_Comm comm) {
    MPI_Recv(&a, 1, MPI_INT, rank, tag, comm, MSI);
  }
  static void Send(const Scal& a, int rank, MPI_Comm comm) {
    MPI_Send(&a, 1, GetMpiScal(), rank, tag, comm);
  }
  static void Recv(Scal& a, int rank, MPI_Comm comm) {
    MPI_Recv(&a, 1, GetMpiScal(), rank, tag, comm, MSI);
  }
  static void Send(const std::string& s, int rank, MPI_Comm comm) {
    int size = s.size();
    MPI_Send(&size, 1, MPI_INT, rank, tag, comm);
    MPI_Send(s.data(), size, MPI_CHAR, rank, tag, comm);
  }
  static void Recv(std::string& s, int rank, MPI_Comm comm) {
    int size;
    MPI_Recv(&size, 1, MPI_INT, rank, tag, comm, MSI);
    std::vector<char> v(size);
    MPI_Recv(v.data(), size, MPI_CHAR, rank, tag, comm, MSI);
    s = std::string(v.begin(), v.end());
  }
  static void Send(const std::vector<Scal>& v, int rank, MPI_Comm comm) {
    int size = v.size();
    MPI_Send(&size, 1, MPI_INT, rank, tag, comm);
    MPI_Send(v.data(), size, GetMpiScal(), rank, tag, comm);
  }
  static void Recv(std::vector<Scal>& v, int rank, MPI_Comm comm) {
    int size;
    MPI_Recv(&size, 1, MPI_INT, rank, tag, comm, MSI);
    v.resize(size);
    MPI_Recv(v.data(), size, GetMpiScal(), rank, tag, comm, MSI);
  }
  static void Send(const Block& b, int rank, MPI_Comm comm) {
    // bounding box
    Send(b.l, rank, comm);
    Send(b.u, rank, comm);
    // stencil
    int stsize = b.stencil.size();
    MPI_Send(&stsize, 1, MPI_INT, rank, tag, comm);
    if (stsize > 0) {
      MPI_Send(b.stencil[0].data(), stsize * dim, MPI_INT, rank, tag, comm);
    }
    Send(*b.a, rank, comm);
    Send(*b.r, rank, comm);
    Send(*b.x, rank, comm);
  }
  static void SendBlocks(
      const std::vector<Block>& bb, int rank, MPI_Comm comm) {
    int nb = bb.size();
    MPI_Send(&nb, 1, MPI_INT, rank, tag, comm);
    for (int i = 0; i < nb; ++i) {
      Send(bb[i], rank, comm);
    }
  }
  static void Recv(BlockBuffer& buf, int rank, MPI_Comm comm) {
    Block b = buf.Get();
    // bounding box
    Recv(b.l, rank, comm);
    Recv(b.u, rank, comm);
    // stencil
    int stsize;
    MPI_Recv(&stsize, 1, MPI_INT, rank, tag, comm, MSI);
    if (stsize > 0) {
      b.stencil.resize(stsize);
      MPI_Recv(
          b.stencil[0].data(), stsize * dim, MPI_INT, rank, tag, comm, MSI);
    }
    Recv(*b.a, rank, comm);
    Recv(*b.r, rank, comm);
    Recv(*b.x, rank, comm);
    buf.Update(b);
  }
  static std::vector<BlockBuffer> RecvBlocks(int rank, MPI_Comm comm) {
    int nb;
    MPI_Recv(&nb, 1, MPI_INT, rank, tag, comm, MSI);
    std::vector<BlockBuffer> vbuf(nb);
    for (auto& buf : vbuf) {
      Recv(buf, rank, comm);
    }
    return vbuf;
  }
  static void Send(Cmd cmd, int rank) {
    int a = static_cast<int>(cmd);
    MPI_Send(&a, 1, MPI_INT, rank, 1, state.commsub);
  }
  static void Send(std::string s, int rank) {
    Send(GetCmd(s), rank);
  }
  static void Recv(Cmd& cmd, int rank, MPI_Comm comm) {
    int a;
    Recv(a, rank, comm);
    cmd = static_cast<Cmd>(a);
  }
  static void RecvYield(int& a, int rank, MPI_Comm comm) {
    MPI_Request req;
    MPI_Irecv(&a, 1, MPI_INT, rank, tag, comm, &req);
    if (1) {
      while (true) {
        int flag = 0;
        MPI_Test(&req, &flag, MSI);
        if (flag) {
          break;
        }
        // sched_yield();
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
      }
    } else {
      MPI_Wait(&req, MSI);
    }
  }
  static void RecvYield(Cmd& cmd, int rank, MPI_Comm comm) {
    int a;
    RecvYield(a, rank, comm);
    cmd = static_cast<Cmd>(a);
  }

  // Returns partitioning.
  // partidx[ranksub] contains indices in bb assigned to ranksub
  static std::vector<std::vector<int>> GetPartIdx(
      const std::vector<Block>& bb, int sizesub) {
    std::vector<std::vector<int>> partidx;
    const int nb = bb.size();
    assert(nb % sizesub == 0);
    int i = 0;
    for (int ranksub = 0; ranksub < sizesub; ++ranksub) {
      std::vector<int> p;
      for (int j = 0; j < nb / sizesub; ++j) {
        p.push_back(i++);
      }
      partidx.push_back(p);
    }
    return partidx;
  }

  Imp(const std::vector<Block>& bb, MIdx gs, MIdx per)
      : bbarg_(bb), partidx_(GetPartIdx(bbarg_, state.sizesub)) {
    static int next_id = 0; // instance id
    id_ = next_id;
    ++next_id;

    Notify(state.commsub);
    for (int rank = 1; rank < state.sizesub; ++rank) {
      Send(Cmd::construct, rank);
      Send(id_, rank, state.commsub);
      DEB(std::cerr << "send construct " << EV(id_) << EV(rank) << GetPart(rank)
                    << std::endl;)
      SendBlocks(GetPart(rank), rank, state.commsub);
      Send(gs, rank, state.commsub);
      Send(per, rank, state.commsub);
    }
    DEB(std::cerr << "self construct " << EV(id_) << GetPart(0) << std::endl;)
    state.minst.emplace(
        std::piecewise_construct, std::forward_as_tuple(id_),
        std::forward_as_tuple(GetBlockBuffers(GetPart(0)), gs, per));
  }
  ~Imp() {
    Notify(state.commsub);
    for (int rank = 1; rank < state.sizesub; ++rank) {
      Send(Cmd::destruct, rank);
      Send(id_, rank, state.commsub);
    }
    state.minst.erase(id_);
    DEB(std::cerr << "self destruct" << std::endl;)
  }
  void Update() {
    Notify(state.commsub);
    for (int rank = 1; rank < state.sizesub; ++rank) {
      Send(Cmd::update, rank);
      Send(id_, rank, state.commsub);
      SendBlocks(GetPart(rank), rank, state.commsub);
    }
    auto& inst = state.minst.at(id_);
    inst.Update(GetPart(0));
    DEB(std::cerr << "self update" << std::endl;)
    inst.GetHypre().Update();
  }
  Scal GetResidual() const {
    auto& inst = state.minst.at(id_);
    return inst.GetHypre().GetResidual();
  }
  int GetIter() const {
    auto& inst = state.minst.at(id_);
    return inst.GetHypre().GetIter();
  }
  void Solve(Scal tol, int print, std::string solver, int maxiter) {
    Notify(state.commsub);
    auto comm = state.commsub;
    for (int rank = 1; rank < state.sizesub; ++rank) {
      Send(Cmd::solve, rank);
      Send(id_, rank, comm);
      Send(tol, rank, comm);
      Send(print, rank, comm);
      Send(solver, rank, comm);
      Send(maxiter, rank, comm);
    }
    auto& inst = state.minst.at(id_);
    DEB(std::cerr << "self solve" << EV(solver) << std::endl;)
    inst.GetHypre().Solve(tol, print, solver, maxiter);
    for (int rank = 1; rank < state.sizesub; ++rank) {
      auto bb = GetPart(rank);
      auto vbuf = RecvBlocks(rank, comm);
      assert(bb.size() == vbuf.size());
      for (size_t i = 0; i < bb.size(); ++i) {
        auto& b = bb[i];
        auto& s = vbuf[i].Get();
        (*b.x) = (*s.x);
      }
    }
    auto bb = GetPart(0);
    auto ss = inst.GetBlocks(); // source
    assert(bb.size() == ss.size());
    for (size_t i = 0; i < bb.size(); ++i) {
      auto& b = bb[i];
      auto& s = ss[i];
      (*b.x) = (*s.x);
    }
  }

  int id_; // instance id
  std::vector<Block> bbarg_; // blocks passed to constructor
  // partitioning
  // partidx_[ranksub] contains indices in bbarg_ assigned to ranksub
  std::vector<std::vector<int>> partidx_;
};

HypreSub::Imp::ServerState HypreSub::Imp::state;

void HypreSub::InitServer(MPI_Comm comm, MPI_Comm commsub) {
  Imp::InitServer(comm, commsub);
}

void HypreSub::RunServer(Sampler& samp) {
  Imp::RunServer(samp);
}

void HypreSub::StopServer() {
  Imp::StopServer();
}

void HypreSub::Send(std::string s, int rank) {
  Imp::Send(s, rank);
}

void HypreSub::Send(const Block& b, int rank) {
  Imp::Send(b, rank, Imp::state.commsub);
}

void HypreSub::Send(const std::vector<Block>& bb, int rank) {
  Imp::SendBlocks(bb, rank, Imp::state.commsub);
}

HypreSub::HypreSub(MPI_Comm, const std::vector<Block>& bb, MIdx gs, MIdx per)
    : imp(new Imp(bb, gs, per)) {}

void HypreSub::Update() {
  imp->Update();
}

void HypreSub::Solve(Scal tol, int print, std::string solver, int maxiter) {
  imp->Solve(tol, print, solver, maxiter);
}

auto HypreSub::GetResidual() const -> Scal {
  return imp->GetResidual();
}

int HypreSub::GetIter() const {
  return imp->GetIter();
}

HypreSub::~HypreSub() {}
