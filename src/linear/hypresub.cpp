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
  static void Execute(Cmd cmd) {
    std::cout << GetString(cmd) + " received" << std::endl;
  }
  static void InitServer(MPI_Comm comm, MPI_Comm commsub) {
    state.comm = comm;
    state.commsub = commsub;
    MPI_Comm_rank(comm, &state.rank);
    MPI_Comm_rank(commsub, &state.ranksub);
  }
  static void RunServer() {
    while (true) {
      int a;
      MPI_Recv(&a, 1, MPI_INT, 0, 1, state.comm, MPI_STATUS_IGNORE);
      auto cmd = static_cast<Cmd>(a);
      std::cout 
        << "recv"
        << " rank=" << state.rank
        << " ranksub=" << state.ranksub
        << " " << GetString(cmd) << std::endl;
      if (cmd == Cmd::construct) {
        auto inst = RecvBlocks(0, state.comm);
      } else if (cmd == Cmd::exit) {
        break;
      }
    }
  }
  static void SendVector(const std::vector<Scal>& v,
                         int rank, MPI_Comm comm) {
    int size = v.size();
    MPI_Send(&size, 1, MPI_INT, rank, tag, comm);
    MPI_Send(v.data(), size, MPI_SCAL, rank, tag, comm);
  }
  static void RecvVector(std::vector<Scal>& v,
                         int rank, MPI_Comm comm) {
    int size;
    MPI_Recv(&size, 1, MPI_INT, rank, tag, comm, MSI);
    v.resize(size);
    MPI_Recv(v.data(), size, MPI_SCAL, rank, tag, comm, MSI);
  }
  static void SendBlock(const Block& b, int rank, MPI_Comm comm) {
    // bounding box
    MPI_Send(b.l.data(), dim, MPI_INT, rank, tag, comm);
    MPI_Send(b.u.data(), dim, MPI_INT, rank, tag, comm);
    // stencil
    int stsize = b.st.size();
    MPI_Send(&stsize, 1, MPI_INT, rank, tag, comm);
    if (stsize > 0) {
      MPI_Send(b.st[0].data(), stsize * dim, MPI_INT, rank, tag, comm);
    }
    SendVector(*b.a, rank, comm);
    SendVector(*b.r, rank, comm);
    SendVector(*b.x, rank, comm);
  }
  static void SendBlocks(const std::vector<Block>& bb,
                         int rank, MPI_Comm comm) {
    int nb = bb.size();
    MPI_Send(&nb, 1, MPI_INT, rank, tag, comm);
    for (int i = 0; i < nb; ++i) {
      SendBlock(bb[i], rank, comm);
      std::cout
          << "send i=" << i
          << " block=" << bb[i]
          << std::endl;
    }
  }
  static void RecvBlock(Block& b, int rank, MPI_Comm comm) {
    // bounding box
    MPI_Recv(b.l.data(), dim, MPI_INT, rank, tag, comm, MSI);
    MPI_Recv(b.u.data(), dim, MPI_INT, rank, tag, comm, MSI);
    // stencil
    int stsize;
    MPI_Recv(&stsize, 1, MPI_INT, rank, tag, comm, MSI);
    if (stsize > 0) {
      b.st.resize(stsize);
      MPI_Recv(b.st[0].data(), stsize * dim, MPI_INT, rank, tag, comm, MSI);
    }
    // matrix
    RecvVector(*b.a, rank, comm);
    RecvVector(*b.r, rank, comm);
    RecvVector(*b.x, rank, comm);
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
      RecvBlock(buf.b, rank, comm);
      std::cout
          << "recv i=" << i
          << " block=" << buf.b
          << std::endl;
    }
    return inst;
  }
  static void Send(std::string s, int rank) {
    int a = static_cast<int>(GetCmd(s));
    MPI_Send(&a, 1, MPI_INT, rank, 1, state.comm);
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
  Imp::SendBlock(b, rank, Imp::state.comm);
}

void HypreSub::Send(const std::vector<Block>& bb, int rank) {
  Imp::SendBlocks(bb, rank, Imp::state.comm);
}
