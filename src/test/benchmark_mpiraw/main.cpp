// Created by Petr Karnakov on 17.10.2019
// Copyright 2019 ETH Zurich

#undef NDEBUG
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "util/mpi.h"
#include "util/timer.h"

auto comm = MPI_COMM_WORLD;
int sc, rank;

void Gather() {
  std::vector<char> r;
  std::vector<double> tt;
  for (auto i = 0; i < 100; ++i) {
    r = std::vector<char>(100, rank + i);
    int s = r.size(); // size local
    int sm = 0; // size max

    MPI_Allreduce(&s, &sm, 1, MPI_INT, MPI_MAX, comm);

    r.resize(sm);

    if (rank == 0) {
      std::vector<char> ra(sm * sc); // result all

      SingleTimer st;
      MPI_Gather(
          r.data(), r.size(), MPI_CHAR, ra.data(), r.size(), MPI_CHAR, 0, comm);
      auto t = st.GetSeconds();
      tt.push_back(t);
    } else {
      MPI_Gather(r.data(), r.size(), MPI_CHAR, nullptr, 0, MPI_CHAR, 0, comm);
    }
  }

  std::sort(tt.begin(), tt.end());
  for (size_t i = tt.size(); i + 10 > tt.size() && i > 0;) {
    --i;
    std::cerr << __func__ << ": time=" << tt[i] << std::endl;
  }
}

void Gatherv() {
  std::vector<char> r;
  std::vector<double> tt;
  for (auto i = 0; i < 20; ++i) {
    r = std::vector<char>(100, rank + i);
    int s = r.size(); // size local
    if (rank == 0) {
      std::vector<int> ss(sc); // size of r on all ranks

      MPI_Gather(&s, 1, MPI_INT, ss.data(), 1, MPI_INT, 0, comm);

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
      MPI_Gatherv(
          r.data(), r.size(), MPI_CHAR, ra.data(), ss.data(), oo.data(),
          MPI_CHAR, 0, comm);
      auto t = st.GetSeconds();
      tt.push_back(t);
    } else {
      MPI_Gather(&s, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm);
      MPI_Gatherv(
          r.data(), r.size(), MPI_CHAR, nullptr, nullptr, nullptr, MPI_CHAR, 0,
          comm);
    }
  }
  std::sort(tt.begin(), tt.end());
  for (size_t i = tt.size(); i + 10 > tt.size() && i > 0;) {
    --i;
    std::cerr << __func__ << ": time=" << tt[i] << std::endl;
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  std::cerr << std::scientific;

  MPI_Comm_size(comm, &sc);
  MPI_Comm_rank(comm, &rank);

  Gatherv();
  Gather();
  MPI_Finalize();
}
