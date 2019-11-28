#include "subcomm.h"
#include "sysinfo.h"

// Example for two nodes each with two threads.
// Ranks:
// comm_world:  | 0 1 | 2 3 |
// comm_omp:    | 0 1 | 0 1 |
// comm_master: | 0 - | 1 - |


void SubComm(
    MPI_Comm& comm_world, MPI_Comm& comm_omp, MPI_Comm& comm_master) {
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
      _SETAFFINITY(tid, thread_affinity);

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
          _SETAFFINITY(i, thread_affinity);

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

  ////////////////////////////////////////////////////////////////////////////
  // OpenMP strong scaling
  command = -1;
  if (omp_master_) { // omp root rank enters hybrid region
    int omp_size;
    MPI_Comm_size(comm_omp_, &omp_size);
#if 1 == _USE_OMP_
    std::vector<int> nthreads = {1, 2, 4, 8, 12};
    std::vector<double> result;
    std::vector<std::pair<double, double>> avg_time; // inner, outer
    for (auto threads : nthreads) {
      const size_t nsteps = 10000000;
      const double step = 1.0 / static_cast<double>(nsteps);
      double sum = 0.0;
      double t_inner = 0.0;
      const double t0_outer = omp_get_wtime();
#pragma omp parallel num_threads(threads) reduction(+ : sum, t_inner)
      {
        // each thread configures its own affinity
        const int tid = omp_get_thread_num();
        _SETAFFINITY(tid, thread_affinity);
        const int nthreads = omp_get_num_threads();

        const double t0_inner = omp_get_wtime();
#pragma omp for
        for (size_t i = 0; i < nsteps; ++i) {
          const double x = (i + 0.5) * step;
          sum += 4.0 / (1.0 + x * x);
        }
        t_inner = omp_get_wtime() - t0_inner;
      }
      const double t_outer = omp_get_wtime() - t0_outer;
      result.push_back(sum * step);
      avg_time.push_back(std::make_pair(t_inner / threads, t_outer));
    }
    int comm_rank, comm_size;
    MPI_Comm_rank(comm_, &comm_rank);
    MPI_Comm_size(comm_, &comm_size);
    for (size_t i = 0; i < nthreads.size(); ++i) {
      const int num_threads = nthreads[i];
      const double t0i = avg_time[0].first;
      const double ti = avg_time[i].first;
      const double to = avg_time[i].second;
      printf("omp_master %d/%d:\tnthreads:%d;\tspeedup:%4.2f;\t"
          "efficiency:%5.1f%%;\toverhead:%5.1f%%;\tresult=%.4e\n",
          comm_rank,
          comm_size,
          num_threads,
          t0i / ti,
          (t0i / ti) / num_threads * 100.0,
          (to - ti) / to * 100.0,
          result[i]);
    }
#endif /* _USE_OMP_ */
    command = 0;
    MPI_Bcast(&command, 1, MPI_INT, 0, comm_omp_);
  } else { // others wait for command
    MPI_Bcast(&command, 1, MPI_INT, 0, comm_omp_);
  }

  MPI_Finalize();
}
