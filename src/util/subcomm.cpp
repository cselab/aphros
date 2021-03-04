// Created by Petr Karnakov on 28.11.2019
// Copyright 2019 ETH Zurich

#include <errno.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <sched.h>
#include <unistd.h>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "logger.h"
#include "subcomm.h"
#include "sysinfo.h"

#define EV(x) (#x) << "=" << (x) << " "
#define EVV(x) std::cerr << (#x) << "=" << (x) << std::endl;
#define EVVS(s, x) s << (#x) << "=" << (x) << std::endl;

template <class T>
static std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
  out << "[";
  for (auto& a : v) {
    out << a << " ";
  }
  out << "]";
  return out;
}

void PrintFromRoot(const std::string& str, MPI_Comm comm) {
  int rank, commsize;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &commsize);
  int cnt = str.size();

  if (rank == 0) {
    std::vector<int> cnts(commsize);
    MPI_Gather(&cnt, 1, MPI_INT, cnts.data(), 1, MPI_INT, 0, comm);

    std::vector<int> displs(commsize + 1);
    displs[0] = 0;
    for (int i = 0; i < commsize; ++i) {
      displs[i + 1] = displs[i] + cnts[i];
    }

    std::vector<char> buf(displs[commsize]);
    MPI_Gatherv(
        str.data(), cnt, MPI_CHAR, buf.data(), cnts.data(), displs.data(),
        MPI_CHAR, 0, comm);

    for (auto c : buf) {
      std::cout << c;
    }
  } else {
    MPI_Gather(&cnt, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm);
    MPI_Gatherv(
        str.data(), cnt, MPI_CHAR, nullptr, nullptr, nullptr, MPI_CHAR, 0,
        comm);
  }
}

#define _SETAFFINITY(tid, cpu_map)                                             \
  do {                                                                         \
    cpu_set_t cpuset;                                                          \
    CPU_ZERO(&cpuset);                                                         \
    CPU_SET((cpu_map)[(tid)], &cpuset);                                        \
    if (0 != sched_setaffinity(0, sizeof(cpu_set_t), &cpuset)) {               \
      fprintf(stderr, "Can not set affinity for thread %d errno=%d\n", (tid)); \
    }                                                                          \
  } while (0)

#define NCHAR_HOST 512
struct Affinity {
  int node_ID;
  int core_ID;
  char hostname[NCHAR_HOST];
};

static std::ostream& operator<<(std::ostream& out, const Affinity& a) {
  out << EV(a.node_ID) << " " << EV(a.core_ID) << " " << EV(a.hostname);
  return out;
}

// Moves current thread to given cpu
void SetAffinity(int cpu) {
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(cpu, &cpuset);
  if (0 != sched_setaffinity(0, sizeof(cpu_set_t), &cpuset)) {
    fprintf(stderr, "Can not move to cpu=%d errno=%d\n", cpu, errno);
  }
}

Affinity GetAffinity() {
  Affinity ret;
  unsigned long a, d, c;
  __asm__ volatile("rdtscp" : "=a"(a), "=d"(d), "=c"(c));
  ret.node_ID = (c & 0xFFF000) >> 12;
  ret.core_ID = c & 0xFFF;
  gethostname(ret.hostname, NCHAR_HOST);
  return ret;
}

// Example for two nodes each with two threads.
// Ranks:
// comm_world:  | 0 1 | 2 3 |
// comm_omp:    | 0 1 | 0 1 |
// comm_master: | 0 - | 1 - |
void SubComm(MPI_Comm& comm_world, MPI_Comm& comm_omp, MPI_Comm& comm_master) {
  comm_world = MPI_COMM_WORLD;
  MPI_Group omp_master_group_;

  ////////////////////////////////////////////////////////////////////////////
  // 1. Find color of ranks in comm_world that belong to same socket and
  // create split communicator stored in comm_omp
  // 2. Create new communicator comm_master with only root ranks in comm_omp

  // 1.
  Affinity a = GetAffinity();
  const size_t host_ID = std::hash<std::string>{}(std::string(a.hostname));
  // const size_t c_hash = host_ID ^ (static_cast<size_t>(a.node_ID) << 1);
  const size_t c_hash = host_ID;

  // construct color and split communicator
  int hypre_size, hypre_rank;
  MPI_Comm_size(comm_world, &hypre_size);
  MPI_Comm_rank(comm_world, &hypre_rank);
  std::vector<size_t> all_hashes(hypre_size);
  MPI_Allgather(
      &c_hash, 1, MPI_UINT64_T, all_hashes.data(), 1, MPI_UINT64_T, comm_world);
  std::set<size_t> unique_nodes(all_hashes.begin(), all_hashes.end());
  std::map<size_t, int> cmap;
  int color = 0;
  for (auto ID : unique_nodes) {
    cmap[ID] = color++;
  }
  color = cmap[c_hash];
  MPI_Comm_split(comm_world, color, hypre_rank, &comm_omp);

  // check if hyper-threads are available
  int omp_size;
  MPI_Comm_size(comm_omp, &omp_size);
  std::vector<int> thread_affinity(omp_size);
  std::iota(thread_affinity.begin(), thread_affinity.end(), 0);
  std::vector<int> mpi_affinity(omp_size);
  MPI_Allgather(
      &a.core_ID, 1, MPI_INT, mpi_affinity.data(), 1, MPI_INT, comm_omp);
  /*
  if (sysinfo::HasHyperthreads()) {
    for (size_t i = 0; i < thread_affinity.size(); ++i) {
      const int mpi_core = mpi_affinity[i];
      thread_affinity[i] = (mpi_core < omp_size) ? omp_size + mpi_core
        : mpi_core - omp_size;
    }
  } else {
    thread_affinity = mpi_affinity;
  }
  */
  thread_affinity = mpi_affinity;

#ifdef _OPENMP
  fassert(
      omp_get_max_threads() <= omp_size,
      "Not enough MPI processes: " + NAMEVALUE(omp_size) + " but " +
          NAMEVALUE(omp_get_max_threads()));
#endif

  // 2.
  int omp_rank;
  MPI_Comm_rank(comm_omp, &omp_rank);

#ifdef _OPENMP
  if (omp_rank == 0) {
#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      SetAffinity(thread_affinity[tid]);
    }
  } else {
    SetAffinity(thread_affinity[omp_rank]);
  }
#endif

  // create comm_master
  if (0 != omp_rank) {
    omp_rank = -omp_rank; // if I am not root in comm_omp
  } else {
    omp_rank = hypre_rank;
  }
  std::vector<int> omp2hypre(hypre_size);
  MPI_Allgather(
      &omp_rank, 1, MPI_INT, omp2hypre.data(), 1, MPI_INT, comm_world);
  std::set<int> m2h(omp2hypre.begin(), omp2hypre.end());
  omp2hypre.clear();
  for (auto r : m2h) {
    if (r >= 0) {
      omp2hypre.push_back(r);
    }
  }
  std::sort(omp2hypre.begin(), omp2hypre.end());
  MPI_Group hypre_group;
  MPI_Comm_group(comm_world, &hypre_group);
  MPI_Group_incl(
      hypre_group, static_cast<int>(omp2hypre.size()), omp2hypre.data(),
      &omp_master_group_);
  MPI_Comm_create(comm_world, omp_master_group_, &comm_master);
  MPI_Group_free(&hypre_group);
}

void PrintStats(MPI_Comm comm_world, MPI_Comm comm_omp, MPI_Comm comm_master) {
  int size_world, rank_world;
  int size_omp, rank_omp;

  MPI_Comm_size(comm_world, &size_world);
  MPI_Comm_rank(comm_world, &rank_world);
  MPI_Comm_size(comm_omp, &size_omp);
  MPI_Comm_rank(comm_omp, &rank_omp);

  Affinity a = GetAffinity();

  std::vector<int> affall(size_omp);
  MPI_Allgather(&a.core_ID, 1, MPI_INT, affall.data(), 1, MPI_INT, comm_omp);

  {
    std::stringstream s;
    s << EV(size_world) << EV(rank_world) << std::endl;
    s << EV(size_omp) << EV(rank_omp) << std::endl;
    EVVS(s, GetAffinity());
    s << std::endl;
    PrintFromRoot(s.str(), comm_world);
  }

  {
    std::stringstream s;
    if (rank_omp == 0) {
      int size_master, rank_master;
      MPI_Comm_size(comm_master, &size_master);
      MPI_Comm_rank(comm_master, &rank_master);
      s << EV(size_master) << EV(rank_master) << std::endl;
      EVVS(s, GetAffinity());
      std::vector<int> affinity = affall;
      std::sort(affinity.begin(), affinity.end());
      EVVS(s, affinity);
      s << std::endl;
    }
    PrintFromRoot(s.str(), comm_world);
  }
}
