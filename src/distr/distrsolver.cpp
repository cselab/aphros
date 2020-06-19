// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#ifdef _OPENMP
#include <omp.h>
#endif

#include <set>

#include "distrsolver.h"
#include "linear/hypresub.h"
#include "util/git.h"
#include "util/subcomm.h"

static void RunKernelOpenMP(
    MPI_Comm comm_world, MPI_Comm comm_omp, MPI_Comm comm_master,
    std::function<void(MPI_Comm, Vars&)> kernel, Vars& var) {
  int rank_omp;
  MPI_Comm_rank(comm_omp, &rank_omp);

  Histogram hist(comm_world, "runkernelOMP", var.Int["histogram"]);
  HypreSub::InitServer(comm_world, comm_omp);
  if (rank_omp == 0) {
    kernel(comm_master, var);
    HypreSub::StopServer();
  } else {
    HypreSub::RunServer(hist);
  }
}

int RunMpi0(
    int argc, const char** argv, std::function<void(MPI_Comm, Vars&)> kernel) {
#ifdef _OPENMP
  omp_set_dynamic(0);
#endif
  char string[MPI_MAX_ERROR_STRING];
  int errorcode;
  int prov;
  int resultlen;
  if ((errorcode = MPI_Init_thread(
           &argc, (char***)&argv, MPI_THREAD_MULTIPLE, &prov)) != MPI_SUCCESS) {
    MPI_Error_string(errorcode, string, &resultlen);
    throw std::runtime_error(FILELINE + ": mpi failed: " + string);
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool isroot = (!rank);

  const std::set<std::string> novalue_args = {
      "-v", "--verbose", "-h", "--help", "--version",
  };
  const auto known_args = novalue_args;
  const auto args = ParseArgs(argc, argv, novalue_args);
  const std::map<std::string, std::string> oargs = args.first; // optional
  const std::vector<std::string> pargs = args.second; // positional

  auto exit_usage = [&argv, isroot](int status) {
    if (isroot) {
      std::cerr << "usage: " << argv[0]
                << " [--version] [-v|--verbose] [a.conf]\n";
    }
    return status;
  };

  if (oargs.count("-h") || oargs.count("--help")) {
    return exit_usage(0);
  }

  const bool verbose = oargs.count("-v") || oargs.count("--verbose");

  std::string filename = "a.conf";
  if (pargs.size() == 0) {
    // nop
  } else if (pargs.size() == 1) {
    filename = pargs[0];
  } else {
    std::cerr << "too many positional arguments\n";
    return exit_usage(1);
  }

  auto check_known_args = [&oargs, isroot]() {
    const std::set<std::string> known = {
        "-v", "--verbose", "-h", "--help", "--version",
    };
    bool pass = true;
    for (auto p : oargs) {
      if (!known.count(p.first)) {
        pass = false;
        if (isroot) {
          std::cerr << "unrecognized option: " << p.first << '\n';
        }
      }
    }
    return pass;
  };
  if (!check_known_args()) {
    return exit_usage(1);
  }

  if (oargs.count("--version") && isroot) {
    std::cerr << "aphros " << GetGitRev() << "\nmsg: " << GetGitMsg()
              << "\ndiff: " << GetGitDiff() << '\n';
  }

  if (verbose && isroot) {
    std::cerr << "Loading config from '" << filename << "'" << std::endl;
  }

  Vars var; // parameter storage
  Parser ip(var); // parser

  ip.ParseFile(filename);

  // Print vars on root
  if (verbose && isroot) {
    std::cerr << "\n=== config begin ===\n";
    ip.PrintAll(std::cerr);
    std::cerr << "=== config end ===\n\n";
  }

  std::string be = var.String["backend"];

  if (be == "local") {
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &comm);
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
      MPI_Comm comm = MPI_COMM_WORLD;
      MPI_Comm comm_omp;
      MPI_Comm_split(comm, rank, rank, &comm_omp);
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

  MPI_Finalize();
  return 0;
}

int RunMpi(
    int argc, const char** argv, std::function<void(MPI_Comm, Vars&)> kernel) {
  try {
    return RunMpi0(argc, argv, kernel);
  } catch (const std::runtime_error& e) {
    const int status = 1;
    std::cerr << //
        "\nterminate called after throwing an instance of "
        "'std::runtime_error'\n"
              << e.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, status);
    return status;
  }
}
