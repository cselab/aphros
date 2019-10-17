#undef NDEBUG
#include <iostream>
#include <vector>
#include <mpi.h>
#include <cassert>

#include "util/metrics.h"

void Run() {
  int sc, rank; // size of communicator
  auto comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &sc);
  MPI_Comm_rank(comm, &rank);
  std::vector<char> r;
  MPI_Barrier(comm);
  for (auto i = 0; i < 100; ++i) {
    r = std::vector<char>(100, i * rank);
    int s = r.size(); // size local
    if (rank == 0) {
      std::vector<int> ss(sc); // size of r on all ranks

      MPI_Gather(&s, 1, MPI_INT, 
                 ss.data(), 1, MPI_INT, 
                 0, comm);

      int sa = 0; // size all
      std::vector<int> oo = {0}; // offsets
      for (auto& q : ss) {
        sa += q;
        oo.push_back(oo.back() + q);
      }
      oo.pop_back();
      assert(ss.size() == oo.size());

      std::vector<char> ra(sa); // result all

      SingleTimer st;
      MPI_Gatherv(r.data(), r.size(), MPI_CHAR,
                  ra.data(), ss.data(), oo.data(), MPI_CHAR,
                  0, comm);
      auto t = st.GetSeconds();
      if (t > 0.01) {
        std::cerr << "i=" << i << ", timer=" << t << std::endl;
      }
    } else {
      MPI_Gather(&s, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm);
      MPI_Gatherv(r.data(), r.size(), MPI_CHAR,
                  nullptr, nullptr, nullptr, MPI_CHAR,
                  0, comm);
    }
  }
  MPI_Barrier(comm);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  Run();
  MPI_Finalize();
}

