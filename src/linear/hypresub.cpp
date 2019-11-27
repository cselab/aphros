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

struct HypreSub::Imp {
  struct ServerState {
    std::vector<int> v;
    MPI_Comm comm;
    MPI_Comm commsub;
    int rank;
    int ranksub;
  };
  enum class Cmd {
      construct, destruct, update, solve, get_residual, get_iter, exit};
  struct BlockBuffer {
    BlockBuffer() = default;
    // disable copy to keep pointers a,r,x in b valid
    BlockBuffer(const BlockBuffer&) = delete;
    Block b;
    std::vector<Scal> a;
    std::vector<Scal> r;
    std::vector<Scal> x;
  };
  struct Instance {
    Instance(int nb) : vbuf(nb) {}
    std::vector<BlockBuffer> vbuf;
  };

  static ServerState state;
  static std::map<int, Instance> minst;
  static constexpr MPI_Datatype MPI_SCAL =
      (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);
  static constexpr int tag = 1;
  static constexpr auto MSI = MPI_STATUS_IGNORE;

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
      std::cout
        << "recv"
        << " rank=" << state.rank
        << " ranksub=" << state.ranksub
        << " " << GetString(cmd) << std::endl;
      if (cmd == Cmd::construct) {
        auto inst = RecvBlocks(root, comm);
        MIdx gs;
        std::array<bool, dim> per;
        Recv(gs, root, comm);
        Recv(per, root, comm);
        std::cout
          << "recv construct"
          << " rank=" << state.rank
          << " ranksub=" << state.ranksub
          << " " << gs << " " << std::endl;
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
  static void Send(const std::array<bool, dim>& d, int rank, MPI_Comm comm) {
    MIdx w;
    for (size_t i = 0; i < dim; ++i) {
      w[i] = d[i];
    }
    Send(w, rank, comm);
  }
  static void Recv(std::array<bool, dim>& d, int rank, MPI_Comm comm) {
    MIdx w;
    Recv(w, rank, comm);
    for (size_t i = 0; i < dim; ++i) {
      d[i] = w[i];
    }
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
      std::cout
          << "send i=" << i
          << " block=" << bb[i]
          << std::endl;
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
  static Instance RecvBlocks(int rank, MPI_Comm comm) {
    auto MSI = MPI_STATUS_IGNORE;
    int nb;
    MPI_Recv(&nb, 1, MPI_INT, rank, tag, comm, MSI);
    Instance inst(nb);
    for (int i = 0; i < nb; ++i) {
      BlockBuffer& buf = inst.vbuf[i];
      buf.b.a = &buf.a;
      buf.b.r = &buf.r;
      buf.b.x = &buf.x;
      Recv(buf.b, rank, comm);
      std::cout
          << "recv i=" << i
          << " block=" << buf.b
          << std::endl;
    }
    return inst;
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

  Imp(const std::vector<Block>& bb, MIdx gs, std::array<bool, dim> per) {
    for (auto rank : {1, 2}) {
      std::vector<Block> bbl = {bb[rank - 1], bb[0]};
      Send(Cmd::construct, rank);
      SendBlocks(bbl, rank, state.comm);
      Send(gs, rank, state.comm);
      Send(per, rank, state.comm);
    }
  }
};

HypreSub::Imp::ServerState HypreSub::Imp::state;
std::map<int, HypreSub::Imp::Instance> HypreSub::Imp::minst;

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

HypreSub::HypreSub(MPI_Comm, const std::vector<Block>& bb,
                   MIdx gs, std::array<bool, dim> per)
    : imp(new Imp(bb, gs, per))
{}

HypreSub::~HypreSub() {}
