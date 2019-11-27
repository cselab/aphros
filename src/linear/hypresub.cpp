#include <iostream>
#include <stdexcept>
#include <utility>
#include <map>

#include "hypresub.h"

static std::ostream& operator<<(
    std::ostream& out, const HypreSub::MIdx& v) {
  out << "(";
  for (auto& a : v) {
    out << a << ",";
  }
  out << ")";
  return out;
}

template <class T>
static std::ostream& operator<<(
    std::ostream& out, const std::vector<T>& v) {
  out << "[";
  for (auto& a : v) {
    out << a << " ";
  }
  out << "]";
  return out;
}


static std::ostream& operator<<(std::ostream& out, const HypreSub::Block& b) {
  out
    << "b.l=" << b.l
    << " b.u=" << b.u
    << " b.st=" << b.st
    << " b.a.size()=" << b.a->size()
    ;
  return out;
}

#define EV(x) " " << (#x) << "=" << (x)

struct HypreSub::Imp {
  enum class Cmd {
      construct, destruct, update, solve, get_residual, get_iter, exit};
  struct BlockBuffer {
    BlockBuffer() = default;
    BlockBuffer(Block b0)
        : b(b0), a(*b0.a), r(*b0.r), x(*b0.x) {
      b.a = &a;
      b.r = &r;
      b.x = &x;
    }
    // disable copy to keep pointers a,r,x in b valid
    BlockBuffer(const BlockBuffer&) = delete;
    Block b;
    std::vector<Scal> a;
    std::vector<Scal> r;
    std::vector<Scal> x;
  };
  struct Instance {
    Instance() = delete;
    Instance(const Instance&) = delete;
    Instance(std::vector<BlockBuffer>&& vbuf0, MIdx gs, MIdx per)
        : vbuf(std::move(vbuf0))
        , hypre(state.comm, GetBlocks(vbuf), gs, per) {
      std::cout
          << "Instance()"
          << EV(state.rank)
          << EV(state.ranksub)
          << EV(vbuf.size())
          << EV(gs)
          << EV(per)
          << std::endl;
    }
    ~Instance() {
      std::cout
          << "~Instance()"
          << EV(state.rank)
          << EV(state.ranksub)
          << EV(vbuf.size())
          << std::endl;
    }
    std::vector<BlockBuffer> vbuf;
    Hypre hypre;
  };
  struct ServerState {
    std::vector<int> v;
    MPI_Comm comm;
    MPI_Comm commsub;
    int rank;
    int ranksub;
    std::map<int, Instance> minst;
  };

  static ServerState state;
  static constexpr MPI_Datatype MPI_SCAL =
      (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);
  static constexpr int tag = 1;
  static constexpr auto MSI = MPI_STATUS_IGNORE;

  static std::vector<Block> GetBlocks(const std::vector<BlockBuffer>& vbuf) {
    std::vector<Block> bb;
    for (auto& buf : vbuf) {
      bb.push_back(buf.b);
    }
    return bb;
  }
  static Cmd GetCmd(std::string s) {
    #define GETCMD(x) if (s == #x) return Cmd::x;
    GETCMD(construct);
    GETCMD(destruct);
    GETCMD(update);
    GETCMD(solve);
    GETCMD(get_residual);
    GETCMD(get_iter);
    GETCMD(exit);
    throw std::runtime_error("GetCmd: unknown name '" + s + "'");
  }
  static std::string GetString(Cmd cmd) {
    #define GETSTR(x) if (cmd == Cmd::x) return #x;
    GETSTR(construct);
    GETSTR(destruct);
    GETSTR(update);
    GETSTR(solve);
    GETSTR(get_residual);
    GETSTR(get_iter);
    GETSTR(exit);
    throw std::runtime_error("GetString: unknown cmd");
  }
  static void InitServer(MPI_Comm comm, MPI_Comm commsub) {
    state.comm = comm;
    state.commsub = commsub;
    MPI_Comm_rank(comm, &state.rank);
    MPI_Comm_rank(commsub, &state.ranksub);
  }
  static void RunServer() {
    while (true) {
      Cmd cmd;
      int root = 0;
      Recv(cmd, root);
      auto comm = state.comm;
      if (cmd == Cmd::construct) {
        auto vbuf = RecvBlocks(root, comm);
        MIdx gs;
        MIdx per;
        int id;
        Recv(gs, root, comm);
        Recv(per, root, comm);
        Recv(id, root, comm);
        std::cout << "recv construct" << EV(id) << std::endl;
        state.minst.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(id),
            std::forward_as_tuple(std::move(vbuf), gs, per));
      } else if (cmd == Cmd::destruct) {
        int id;
        Recv(id, root, comm);
        std::cout << "recv destruct" << EV(id) << std::endl;
        state.minst.erase(id);
      } else if (cmd == Cmd::exit) {
        break;
      }
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
  static void Send(const std::vector<Scal>& v, int rank, MPI_Comm comm) {
    int size = v.size();
    MPI_Send(&size, 1, MPI_INT, rank, tag, comm);
    MPI_Send(v.data(), size, MPI_SCAL, rank, tag, comm);
  }
  static void Recv(std::vector<Scal>& v, int rank, MPI_Comm comm) {
    int size;
    MPI_Recv(&size, 1, MPI_INT, rank, tag, comm, MSI);
    v.resize(size);
    MPI_Recv(v.data(), size, MPI_SCAL, rank, tag, comm, MSI);
  }
  static void Send(const Block& b, int rank, MPI_Comm comm) {
    // bounding box
    Send(b.l, rank, comm);
    Send(b.u, rank, comm);
    // stencil
    int stsize = b.st.size();
    MPI_Send(&stsize, 1, MPI_INT, rank, tag, comm);
    if (stsize > 0) {
      MPI_Send(b.st[0].data(), stsize * dim, MPI_INT, rank, tag, comm);
    }
    Send(*b.a, rank, comm);
    Send(*b.r, rank, comm);
    Send(*b.x, rank, comm);
  }
  static void SendBlocks(const std::vector<Block>& bb,
                         int rank, MPI_Comm comm) {
    int nb = bb.size();
    MPI_Send(&nb, 1, MPI_INT, rank, tag, comm);
    for (int i = 0; i < nb; ++i) {
      Send(bb[i], rank, comm);
    }
  }
  static void Recv(Block& b, int rank, MPI_Comm comm) {
    // bounding box
    Recv(b.l, rank, comm);
    Recv(b.u, rank, comm);
    // stencil
    int stsize;
    MPI_Recv(&stsize, 1, MPI_INT, rank, tag, comm, MSI);
    if (stsize > 0) {
      b.st.resize(stsize);
      MPI_Recv(b.st[0].data(), stsize * dim, MPI_INT, rank, tag, comm, MSI);
    }
    Recv(*b.a, rank, comm);
    Recv(*b.r, rank, comm);
    Recv(*b.x, rank, comm);
  }
  static std::vector<BlockBuffer> RecvBlocks(int rank, MPI_Comm comm) {
    int nb;
    MPI_Recv(&nb, 1, MPI_INT, rank, tag, comm, MSI);
    std::vector<BlockBuffer> vbuf(nb);
    for (int i = 0; i < nb; ++i) {
      BlockBuffer& buf = vbuf[i];
      buf.b.a = &buf.a;
      buf.b.r = &buf.r;
      buf.b.x = &buf.x;
      Recv(buf.b, rank, comm);
    }
    return vbuf;
  }
  static void Send(Cmd cmd, int rank) {
    int a = static_cast<int>(cmd);
    MPI_Send(&a, 1, MPI_INT, rank, 1, state.comm);
  }
  static void Send(std::string s, int rank) {
    Send(GetCmd(s), rank);
  }
  static void Recv(Cmd& cmd, int rank) {
    int a;
    MPI_Recv(&a, 1, MPI_INT, rank, 1, state.comm, MSI);
    cmd = static_cast<Cmd>(a);
  }

  Imp(const std::vector<Block>& bb, MIdx gs, MIdx per) {
    static int next_id = 0; // instance id
    id_ = next_id;
    ++next_id;
    for (auto rank : {1, 2}) {
      std::vector<Block> bbl = {bb[rank]};
      Send(Cmd::construct, rank);
      SendBlocks(bbl, rank, state.comm);
      Send(gs, rank, state.comm);
      Send(per, rank, state.comm);
      Send(id_, rank, state.comm);
    }
    std::vector<BlockBuffer> vbuf(1);
    vbuf[0] = BlockBuffer(bb[0]);
    state.minst.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(id_),
        std::forward_as_tuple(std::move(vbuf), gs, per));
  }
  ~Imp() {
    for (auto rank : {1, 2}) {
      Send(Cmd::destruct, rank);
      Send(id_, rank, state.comm);
      state.minst.erase(id_);
    }
  }

  int id_;
};

HypreSub::Imp::ServerState HypreSub::Imp::state;

void HypreSub::InitServer(MPI_Comm comm, MPI_Comm commsub) {
  Imp::InitServer(comm, commsub);
}

void HypreSub::RunServer() {
  Imp::RunServer();
}

void HypreSub::Send(std::string s, int rank) {
  Imp::Send(s, rank);
}

void HypreSub::Send(const Block& b, int rank) {
  Imp::Send(b, rank, Imp::state.comm);
}

void HypreSub::Send(const std::vector<Block>& bb, int rank) {
  Imp::SendBlocks(bb, rank, Imp::state.comm);
}

HypreSub::HypreSub(MPI_Comm, const std::vector<Block>& bb, MIdx gs, MIdx per)
    : imp(new Imp(bb, gs, per))
{}

HypreSub::~HypreSub() {}
