// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#ifdef _OPENMP
#include <omp.h>
#endif
#include <set>

#include "distrsolver.h"
#include "linear/hypresub.h"
#include "parse/argparse.h"
#include "util/git.h"
#include "util/logger.h"
#include "util/subcomm.h"

static void RunKernelOpenMP(
    MPI_Comm comm_world, MPI_Comm comm_omp, MPI_Comm comm_master,
    std::function<void(MPI_Comm, Vars&)> kernel, Vars& var) {
  int rank_omp;
  MPI_Comm_rank(comm_omp, &rank_omp);

  // XXX moving exception handlers to RunMpi() (around RunMpi0())
  // is not possible with MPI since the destructor of MpiWrapper
  // calls MPI_Finalize() so in case of an exception the program freezes.
  try {
    Histogram hist(comm_world, "runkernelOMP", var.Int["histogram"]);
    HypreSub::InitServer(comm_world, comm_omp);
    if (rank_omp == 0) {
      kernel(comm_master, var);
      HypreSub::StopServer();
    } else {
      HypreSub::RunServer(hist);
    }
  } catch (const std::exception& e) {
    std::cerr << FILELINE + "\nabort after throwing exception\n"
              << e.what() << '\n';
    std::terminate();
  }
}

MpiWrapper::MpiWrapper(int* argc, const char*** argv) : comm_(MPI_COMM_WORLD) {
#ifdef _OPENMP
  omp_set_dynamic(0);
#endif
  int required = MPI_THREAD_FUNNELED;
  int provided;
  MPICALL(MPI_Init_thread(argc, (char***)argv, required, &provided));
  fassert_equal(required, provided, ", mismatch in thread support level");
}

MpiWrapper::~MpiWrapper() {
  MPI_Finalize();
}

MPI_Comm MpiWrapper::GetComm() const {
  return comm_;
}

int MpiWrapper::GetCommSize(MPI_Comm comm) {
  int size;
  MPICALL(MPI_Comm_size(comm, &size));
  return size;
}

int MpiWrapper::GetCommRank(MPI_Comm comm) {
  int rank;
  MPICALL(MPI_Comm_rank(comm, &rank));
  return rank;
}

int MpiWrapper::GetCommSize() const {
  return GetCommSize(comm_);
}

int MpiWrapper::GetCommRank() const {
  return GetCommRank(comm_);
}

bool MpiWrapper::IsRoot(MPI_Comm comm) {
  return GetCommRank(comm) == 0;
}

bool MpiWrapper::IsRoot() const {
  return IsRoot(comm_);
}

int RunMpi0(
    int argc, const char** argv, std::function<void(MPI_Comm, Vars&)> kernel) {
  MpiWrapper mpi(&argc, &argv);
  const int rank = mpi.GetCommRank();
  bool isroot = (!rank);

  const auto args = [&argc, &argv, &isroot]() {
    ArgumentParser parser("Distributed solver", isroot);
    parser.AddSwitch({"--verbose", "-v"}).Help("Print initial configuration");
    parser.AddSwitch("--confverbose")
        .Help("Add 'set int verbose 1' to configuration");
    parser.AddSwitch({"--version"}).Help("Print version");
    parser.AddVariable<std::string>("config", "a.conf")
        .Help("Path to configuration file");
    return parser.ParseArgs(argc, argv);
  }();
  if (const int* p = args.Int.Find("EXIT")) {
    return *p;
  }

  const bool verbose = args.Int["verbose"];
  const std::string config = args.String["config"];
  if (args.Int["version"] && isroot) {
    std::cerr << "aphros " << GetGitRev() << "\nmsg: " << GetGitMsg()
              << "\ndiff: " << GetGitDiff() << '\n';
  }

  if (verbose && isroot) {
    std::cerr << "Loading config from '" << config << "'" << std::endl;
  }

  Vars var; // parameter storage
  Parser ip(var); // parser
  ip.ParseFile(config);
  if (args.Int["confverbose"]) {
    var.Int.Set("verbose", 1);
  }

  // Print vars on root
  if (verbose && isroot) {
    std::cerr << "\n=== config begin ===\n";
    ip.PrintAll(std::cerr);
    std::cerr << "=== config end ===\n\n";
  }

  const std::string backend = var.String["backend"];

  if (backend == "local") {
    MPI_Comm comm;
    MPICALL(MPI_Comm_split(mpi.GetComm(), rank, rank, &comm));
    if (rank == 0) {
      RunKernelOpenMP(comm, comm, comm, kernel, var);
    }
  } else {
    bool openmp = var.Int["openmp"];
    if (openmp) {
      MPI_Comm comm_world;
      MPI_Comm comm_omp;
      MPI_Comm comm_master;
      SubComm(comm_world, comm_omp, comm_master);
      if (var.Int["verbose_openmp"]) {
        PrintStats(comm_world, comm_omp, comm_master);
      }
      RunKernelOpenMP(comm_world, comm_omp, comm_master, kernel, var);
    } else {
      MPI_Comm comm = mpi.GetComm();
      MPI_Comm comm_omp;
      MPICALL(MPI_Comm_split(comm, rank, rank, &comm_omp));
      RunKernelOpenMP(comm, comm_omp, comm, kernel, var);
    }
  }

  if (isroot) {
    if (var.Int("verbose_conf_reads", 0)) {
      auto print_reads = [&var](const auto& map) {
        for (auto it = map.cbegin(); it != map.cend(); ++it) {
          const auto key = it->first;
          std::cout << map.GetReads(key) << ' ' << map.GetTypeName() << ' '
                    << key << '\n';
        }
      };
      std::cout << "Number of accesses to configuration variables\n";
      var.ForEachMap(print_reads);
    }
    if (var.Int("verbose_conf_unused", 0)) {
      const std::string path = var.String["conf_unused_ignore_path"];
      std::set<std::string> ignore;
      if (path != "") {
        Vars vign;
        Parser parser(vign);
        std::ifstream f(path);
        parser.ParseStream(f);
        vign.ForEachMap([&ignore](const auto& map) {
          for (auto it = map.cbegin(); it != map.cend(); ++it) {
            ignore.insert(it->first);
          }
        });
      }
      std::cout << "Unused configuration variables:\n";
      var.ForEachMap([&var, &ignore](const auto& map) {
        for (auto it = map.cbegin(); it != map.cend(); ++it) {
          const auto key = it->first;
          if (map.GetReads(key) == 0 && !ignore.count(key)) {
            std::cout << map.GetTypeName() << ' ' << key << '\n';
          }
        }
      });
    }
  }
  return 0;
}

int RunMpi(
    int argc, const char** argv, std::function<void(MPI_Comm, Vars&)> kernel) {
  return RunMpi0(argc, argv, kernel);
}
