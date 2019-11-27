#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <mpi.h>
#include <numeric>
#include <omp.h>
#include <pthread.h>
#include <sched.h>
#include <set>
#include <sstream>
#include <stdio.h>
#include <string>
#include <thread>
#include <unistd.h>
#include <utility>
#include <vector>

#define _USE_OMP_ 1
#if 1 == _USE_OMP_
#warning "USING OPENMP THREAD MODEL"
#else
#warning "USING C++ THREAD MODEL"
#endif /* _USE_OMP_ */

#define NCHAR_HOST 512
struct Affinity {
    int node_ID;
    int core_ID;
    char hostname[NCHAR_HOST];
};

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

// available in util/sysinfo.h
bool HasHyperthreads() {
  bool has_ht = false;
  std::ifstream f("/proc/cpuinfo");
  f >> std::skipws;
  while (f) {
    std::string s;
    f >> s;
    if (s == "flags") {
        break;
    }
  }
  std::string line;
  std::getline(f, line);
  std::istringstream is(line);
  while (!is.eof()) {
      std::string s;
      is >> s;
      if (s == "ht") {
          has_ht = true;
          break;
      }
  }
  return has_ht;
}

int main(int argc, char *argv[])
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    MPI_Comm comm_hypre_ = MPI_COMM_WORLD;
    MPI_Comm comm_, comm_omp_;
    MPI_Group omp_master_group_;
    bool omp_master_;

    ////////////////////////////////////////////////////////////////////////////
    // 1. Find color of ranks in comm_hypre_ that belong to same socket and
    // create split communicator stored in comm_omp_
    // 2. Create new communicator with only root ranks in comm_omp_, store in
    // temporarily in comm_

    // 1.
    Affinity a = GetAffinity();
    const size_t host_ID = std::hash<std::string>{}(std::string(a.hostname));
    const size_t c_hash = host_ID ^ (static_cast<size_t>(a.node_ID) << 1);

    // construct color and split communicator
    int hypre_size, hypre_rank;
    MPI_Comm_size(comm_hypre_, &hypre_size);
    MPI_Comm_rank(comm_hypre_, &hypre_rank);
    std::vector<size_t> all_hashes(hypre_size);
    MPI_Allgather(&c_hash,
                  1,
                  MPI_UINT64_T,
                  all_hashes.data(),
                  1,
                  MPI_UINT64_T,
                  comm_hypre_);
    std::set<size_t> unique_nodes(all_hashes.begin(), all_hashes.end());
    std::map<size_t, int> cmap;
    int color = 0;
    for (auto ID : unique_nodes) {
        cmap[ID] = color++;
    }
    color = cmap[c_hash];
    MPI_Comm_split(comm_hypre_, color, hypre_rank, &comm_omp_);

    // check if hyper-threads are available
    int omp_size;
    MPI_Comm_size(comm_omp_, &omp_size);
    std::vector<int> thread_affinity(omp_size);
    std::iota(thread_affinity.begin(), thread_affinity.end(), 0);
    if (HasHyperthreads()) {
        std::vector<int> mpi_affinity(omp_size);
        MPI_Allgather(
            &a.core_ID, 1, MPI_INT, mpi_affinity.data(), 1, MPI_INT, comm_omp_);
        for (size_t i = 0; i < thread_affinity.size(); ++i) {
            const int mpi_core = mpi_affinity[i];
            thread_affinity[i] = (mpi_core < omp_size) ? omp_size + mpi_core
                                                       : mpi_core - omp_size;
        }
    }

    // 2.
    int omp_rank;
    MPI_Comm_rank(comm_omp_, &omp_rank);
    if (0 != omp_rank) {
        omp_rank = -omp_rank; // if I am not root in comm_omp_
        omp_master_ = false;
    } else {
        omp_rank = hypre_rank;
        omp_master_ = true;
    }
    std::vector<int> omp2hypre(hypre_size);
    MPI_Allgather(
        &omp_rank, 1, MPI_INT, omp2hypre.data(), 1, MPI_INT, comm_hypre_);
    std::set<int> m2h(omp2hypre.begin(), omp2hypre.end());
    omp2hypre.clear();
    for (auto r : m2h) {
        if (r >= 0) {
            omp2hypre.push_back(r);
        }
    }
    std::sort(omp2hypre.begin(), omp2hypre.end());

    MPI_Group hypre_group;
    MPI_Comm_group(comm_hypre_, &hypre_group);

    MPI_Group_incl(hypre_group,
                   static_cast<int>(omp2hypre.size()),
                   omp2hypre.data(),
                   &omp_master_group_);
    MPI_Comm_create(comm_hypre_, omp_master_group_, &comm_);
    MPI_Group_free(&hypre_group);
    if (omp_master_) {
        omp_rank = 0;
    } else {
        omp_rank = -omp_rank;
    }

    {
        int comm_rank, comm_size;
        comm_rank = -1;
        comm_size = -1;
        if (omp_master_) {
            MPI_Comm_size(comm_, &comm_size);
            MPI_Comm_rank(comm_, &comm_rank);
        }
        printf(
            "HYBRID: host:%s; comm_hypre_:%d/%d; comm_omp_:%d/%d; comm_:%d/%d; "
            "node_ID=%lu; core_ID=%lu\n",
            a.hostname,
            hypre_rank,
            hypre_size,
            omp_rank,
            omp_size,
            comm_rank,
            comm_size,
            a.node_ID,
            a.core_ID);
        fflush(stdout);
    }
    MPI_Barrier(comm_hypre_);
    int is_root;
    MPI_Comm_rank(comm_hypre_, &is_root);
    is_root = (0 == is_root);
    if (is_root) {
        printf("MPI_Init_thread provided=%d\n", provided);
        fflush(stdout);
    }

    ////////////////////////////////////////////////////////////////////////////
    // do some work
    if (omp_master_) { // only subcomm ranks must do this
        int subrank, subsize, omp_size;
        MPI_Comm_rank(comm_, &subrank);
        MPI_Comm_size(comm_, &subsize);
        MPI_Comm_size(comm_omp_, &omp_size);
        printf("Root on comm_omp_: comm_hypre_:%d/%d; comm_omp_:%d/%d; "
               "comm_:%d/%d\n",
               hypre_rank,
               hypre_size,
               omp_rank,
               omp_size,
               subrank,
               subsize);
    }
    MPI_Barrier(comm_hypre_);
    if (is_root) {
        fflush(stdout);
    }

    ////////////////////////////////////////////////////////////////////////////
    // hybrid
    int command;
    if (omp_master_) { // omp root rank enters hybrid region
        int omp_size;
        MPI_Comm_size(comm_omp_, &omp_size);
#if 1 == _USE_OMP_
#pragma omp parallel num_threads(omp_size)
        {
            // each thread configures its own affinity
            const int tid = omp_get_thread_num();
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(thread_affinity[tid], &cpuset);
            if (0 != sched_setaffinity(0 /* = calling thread */,
                                       sizeof(cpu_set_t),
                                       &cpuset)) {
                fprintf(stderr, "Can not set affinity for thread %d\n", tid);
            }

            const int nthreads = omp_get_num_threads();
            Affinity a = GetAffinity();
#pragma omp critical
            printf("Thread %d/%d: host:%s; node_ID:%d; core_ID:%d\n",
                   tid,
                   nthreads,
                   a.hostname,
                   a.node_ID,
                   a.core_ID);
        }
#else
        std::vector<std::thread> threads(omp_size);
        for (int i = 0; i < static_cast<int>(threads.size()); ++i) {
            threads[i] = std::thread([i, omp_size] {
                // each thread configures its own affinity
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(thread_affinity[i], &cpuset);
                if (0 != sched_setaffinity(0 /* = calling thread */,
                                           sizeof(cpu_set_t),
                                           &cpuset)) {
                    fprintf(stderr, "Can not set affinity for thread %d\n", i);
                }

                Affinity a = GetAffinity();
                printf("Thread %d/%d: host:%s; node_ID:%d; core_ID:%d\n",
                       i,
                       omp_size,
                       a.hostname,
                       a.node_ID,
                       a.core_ID);
            });
            // XXX: [fabianw@mavt.ethz.ch; 2019-11-27] does not work
            // cpu_set_t cpuset;
            // CPU_ZERO(&cpuset);
            // CPU_SET(thread_affinity[i], &cpuset);
            // if (0 != pthread_setaffinity_np(threads[i].native_handle(),
            //                                 sizeof(cpu_set_t),
            //                                 &cpuset)) {
            //     fprintf(stderr, "Can not set affinity for thread %d\n", i);
            // }
        }
        for (auto &t : threads) {
            t.join();
        }
#endif /* _USE_OMP_ */
        command = 0;
        MPI_Bcast(&command, 1, MPI_INT, 0, comm_omp_);
    } else { // others wait for command
        MPI_Bcast(&command, 1, MPI_INT, 0, comm_omp_);
    }

    MPI_Finalize();
    return 0;
}
