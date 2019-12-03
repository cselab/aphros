#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <mpi.h>
#include <numeric>
#include <omp.h>
#include <sched.h>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>
#include <iostream>
#include <errno.h>

#include "subcomm.h"
#include "sysinfo.h"

#define EV(x) (#x) << "=" << (x) << " "
#define EVV(x) std::cerr << (#x) << "=" << (x) << std::endl;

template <class T>
static std::ostream& operator<<(
    std::ostream& out, const std::vector<T>& v) {
  out << "[";
  for (auto& a : v) {
    out << a << " ";
  }
  out << "]";
  return out;
}

#define _SETAFFINITY(tid, cpu_map)                                             \
    do {                                                                       \
        cpu_set_t cpuset;                                                      \
        CPU_ZERO(&cpuset);                                                     \
        CPU_SET((cpu_map)[(tid)], &cpuset);                                    \
        if (0 != sched_setaffinity(0, sizeof(cpu_set_t), &cpuset)) {           \
            fprintf(stderr, "Can not set affinity for thread %d errno=%d\n",   \
                (tid));                                                        \
        }                                                                      \
    } while (0)


#define NCHAR_HOST 512
struct Affinity {
    int node_ID;
    int core_ID;
    char hostname[NCHAR_HOST];
};

// Moves current thread to given cpu
void SetAffinity(int cpu) {
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(cpu, &cpuset);
  if (0 != sched_setaffinity(0, sizeof(cpu_set_t), &cpuset)) {
      fprintf(stderr, "Can not move to cpu=%d errno=%d\n", cpu, errno);
  }
}

Affinity GetAffinity()
{
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
void SubComm(
    MPI_Comm& comm_world, MPI_Comm& comm_omp, MPI_Comm& comm_master) {
  comm_world = MPI_COMM_WORLD;
  MPI_Group omp_master_group_;

  ////////////////////////////////////////////////////////////////////////////
  // 1. Find color of ranks in comm_world that belong to same socket and
  // create split communicator stored in comm_omp
  // 2. Create new communicator comm_master with only root ranks in comm_omp

  // 1.
  Affinity a = GetAffinity();
  const size_t host_ID = std::hash<std::string>{}(std::string(a.hostname));
  const size_t c_hash = host_ID ^ (static_cast<size_t>(a.node_ID) << 1);

  // construct color and split communicator
  int hypre_size, hypre_rank;
  MPI_Comm_size(comm_world, &hypre_size);
  MPI_Comm_rank(comm_world, &hypre_rank);
  std::vector<size_t> all_hashes(hypre_size);
  MPI_Allgather(&c_hash,
      1,
      MPI_UINT64_T,
      all_hashes.data(),
      1,
      MPI_UINT64_T,
      comm_world);
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
  if (sysinfo::HasHyperthreads()) {
    for (size_t i = 0; i < thread_affinity.size(); ++i) {
      const int mpi_core = mpi_affinity[i];
      thread_affinity[i] = (mpi_core < omp_size) ? omp_size + mpi_core
        : mpi_core - omp_size;
    }
  }

  // 2.
  int omp_rank;
  MPI_Comm_rank(comm_omp, &omp_rank);
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

  MPI_Group_incl(hypre_group,
                 static_cast<int>(omp2hypre.size()),
                 omp2hypre.data(),
                 &omp_master_group_);
  MPI_Comm_create(comm_world, omp_master_group_, &comm_master);
  MPI_Group_free(&hypre_group);

#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    SetAffinity(thread_affinity[tid]);
  }
}

void PrintStats(MPI_Comm comm_world, MPI_Comm comm_omp, MPI_Comm comm_master) {
  int size_world, rank_world;
  int size_omp, rank_omp;

  MPI_Comm_size(comm_world, &size_world);
  MPI_Comm_rank(comm_world, &rank_world);
  MPI_Comm_size(comm_omp, &size_omp);
  MPI_Comm_rank(comm_omp, &rank_omp);

  if (rank_omp == 0) {
    int size_master, rank_master;
    MPI_Comm_size(comm_master, &size_master);
    MPI_Comm_rank(comm_master, &rank_master);
    std::cout
        << EV(size_world) << EV(rank_world) << std::endl
        << EV(size_omp) << EV(rank_omp) << std::endl
        << EV(size_master) << EV(rank_master) << std::endl
        << std::endl;
  } else {
    std::cout
        << EV(size_world) << EV(rank_world) << std::endl
        << EV(size_omp) << EV(rank_omp) << std::endl
        << std::endl;
  }
}
