#undef NDEBUG
#include <iostream>
#include <vector>
#include <mpi.h>
#include <cassert>

#include "util/metrics.h"

auto comm = MPI_COMM_WORLD;
int sc, rank;

void Gather() {
  std::vector<char> r;
  for (auto i = 0; i < 100; ++i) {
    r = std::vector<char>(100, rank + i);
    int s = r.size(); // size local
    int sm = 0; // size max

    MPI_Allreduce(&s, &sm, 1, MPI_INT, MPI_MAX, comm);

    r.resize(sm);

    if (rank == 0) {
      std::vector<char> ra(sm * sc); // result all

      SingleTimer st;
      MPI_Gather(r.data(), r.size(), MPI_CHAR,
                 ra.data(), ra.size(), MPI_CHAR, 0, comm);
      auto t = st.GetSeconds();
      if (t > 0.01) {
        std::cerr << __func__ << ":i=" << i << ", timer=" << t << std::endl;
      }
    } else {
      MPI_Gather(r.data(), r.size(), MPI_CHAR,
                 nullptr, 0, MPI_CHAR, 0, comm);
    }
  }
}

void Gatherv() {
  std::vector<char> r;
  for (auto i = 0; i < 100; ++i) {
    r = std::vector<char>(100, rank + i);
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
        std::cerr << __func__ << ":i=" << i << ", timer=" << t << std::endl;
      }
    } else {
      MPI_Gather(&s, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm);
      MPI_Gatherv(r.data(), r.size(), MPI_CHAR,
                  nullptr, nullptr, nullptr, MPI_CHAR,
                  0, comm);
    }
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &sc);
  MPI_Comm_rank(comm, &rank);

  Gatherv();
  Gather();
  MPI_Finalize();
}

