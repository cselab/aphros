#include "hypresub.h"

struct HypreSub::ServerState {
};

void HypreSub::InitServer(MPI_Comm comm, MPI_Comm commsub) {
  (void) comm;
  (void) commsub;
}
void HypreSub::RunServer() {}
