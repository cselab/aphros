#include <iostream>
#include <stdexcept>

#include "hypresub.h"

static std::ostream& operator<<(
    std::ostream& out, const HypreSub::MIdx& v) {
  for (auto& a : v) {
    out << a << " ";
  }
  return out;
}


static std::ostream& operator<<(std::ostream& out, const HypreSub::Block& b) {
  out
    << "b.l=" << b.l
    << " b.u=" << b.u
    << " b.st.size()=" << b.st.size()
    << " b.a.size()=" << b.a->size()
    << std::endl;
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

  static ServerState state;
  static constexpr MPI_Datatype MPI_SCAL =
      (sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT);

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
        std::vector<Scal> a;
        std::vector<Scal> r;
        std::vector<Scal> x;
        Block b;
        b.a = &a;
        b.r = &r;
        b.x = &x;

        RecvBlock(b, 0, state.comm);
        std::cout << b << std::endl;
      } else if (cmd == Cmd::exit) {
        break;
      }
    }
  }
  static void SendBlock(const Block& b, int rank, MPI_Comm comm) {
    const int tag = 1;
    // bounding box
    MPI_Send(b.l.data(), dim, MPI_INT, rank, tag, comm);
    MPI_Send(b.u.data(), dim, MPI_INT, rank, tag, comm);
    // stencil
    int stsize = b.st.size();
    MPI_Send(&stsize, 1, MPI_INT, rank, tag, comm);
    if (stsize > 0) {
      MPI_Send(b.st[0].data(), stsize * dim, MPI_INT, rank, tag, comm);
    }
    // matrix
    int asize = b.a->size();
    MPI_Send(&asize, 1, MPI_INT, rank, tag, comm);
    MPI_Send(b.a->data(), asize, MPI_SCAL, rank, tag, comm);
  }
  static void RecvBlock(Block& b, int rank, MPI_Comm comm) {
    auto MSI = MPI_STATUS_IGNORE;
    const int tag = 1;
    // bounding box
    MPI_Recv(b.l.data(), dim, MPI_INT, rank, tag, comm, MSI);
    MPI_Recv(b.u.data(), dim, MPI_INT, rank, tag, comm, MSI);
    // stencil
    int stsize = b.st.size();
    MPI_Recv(&stsize, 1, MPI_INT, rank, tag, comm, MSI);
    if (stsize > 0) {
      b.st.resize(stsize);
      MPI_Recv(b.st[0].data(), stsize * dim, MPI_INT, rank, tag, comm, MSI);
    }
    // matrix
    int asize = b.a->size();
    MPI_Recv(&asize, 1, MPI_INT, rank, tag, comm, MSI);
    b.a->resize(asize);
    MPI_Recv(b.a->data(), asize, MPI_SCAL, rank, tag, comm, MSI);
  }
  static void Send(std::string s, int rank) {
    int a = static_cast<int>(GetCmd(s));
    MPI_Send(&a, 1, MPI_INT, rank, 1, state.comm);
  }
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
  Imp::SendBlock(b, rank, Imp::state.comm);
}

