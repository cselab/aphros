#include <iostream>

#include "hypresub.h"

struct HypreSub::Imp {
  struct ServerState {
    std::vector<int> v;
  };
  enum class Cmd {
    construct, destruct, update, solve, get_residual, get_iter};

  static ServerState server;
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
