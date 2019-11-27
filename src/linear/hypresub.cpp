#include <iostream>

#include "hypresub.h"

struct HypreSub::ServerState {
  std::vector<int> v;
};

HypreSub::ServerState HypreSub::server;

void HypreSub::InitServer(MPI_Comm comm, MPI_Comm commsub) {
  (void) comm;
  (void) commsub;

  std::cout << "init" << std::endl;
  server.v.resize(5);
}
void HypreSub::RunServer() {
  std::cout << server.v.size() << std::endl;
}
