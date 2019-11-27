#include <iostream>
#include <stdexcept>

#include "hypresub.h"

struct HypreSub::Imp {
  struct ServerState {
    std::vector<int> v;
  };
  enum class Cmd {
    construct, destruct, update, solve, get_residual, get_iter, exit};

  static ServerState server;

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
};

HypreSub::Imp::ServerState HypreSub::Imp::server;

void HypreSub::InitServer(MPI_Comm comm, MPI_Comm commsub) {
  (void) comm;
  (void) commsub;

  std::cout << "init" << std::endl;
  Imp::server.v.resize(5);
}

void HypreSub::RunServer() {
  std::cout << Imp::server.v.size() << std::endl;
}

void HypreSub::Send(std::string cmd) {
  std::cout << cmd << std::endl;
  Imp::Execute(Imp::GetCmd(cmd));
}
