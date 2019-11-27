#include <iostream>
#include <stdexcept>

#include "hypresub.h"

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
      if (cmd == Cmd::exit) {
        break;
      }
    }
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
