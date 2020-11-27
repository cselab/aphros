// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#ifdef _OPENMP
#include <omp.h>
#endif
#include <atomic>
#include <mutex>
#include <set>
#include <sstream>

#include "distrsolver.h"
#include "parse/argparse.h"
#include "util/git.h"
#include "util/logger.h"
#include "util/subcomm.h"

MpiWrapper::MpiWrapper(int* argc, const char*** argv, MPI_Comm comm)
    : comm_(comm) {
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
    int argc, const char** argv,
    const std::function<void(MPI_Comm, Vars&)>& kernel) {
  MpiWrapper mpi(&argc, &argv);
  const int rank = mpi.GetCommRank();
  bool isroot = (!rank);

  const auto args = [&argc, &argv, &isroot]() {
    ArgumentParser parser("Distributed solver", isroot);
    parser.AddSwitch({"--verbose", "-v"}).Help("Print initial configuration");
    parser.AddSwitch("--config_verbose")
        .Help("Append configuration with 'set int verbose 1'");
    parser.AddVariable<std::string>("--extra", "")
        .Help("Extra configuration (commands `set ...`)");
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

  std::map<std::string, std::atomic<int>> var_reads;
  std::mutex var_reads_mutex;
  auto hook_read = [&var_reads, &var_reads_mutex](const std::string& key) {
    auto it = var_reads.find(key);
    if (it == var_reads.end()) {
      std::lock_guard<std::mutex> guard(var_reads_mutex);
      it = var_reads.emplace(key, 0).first;
    }
    ++it->second;
  };

  Vars var(hook_read); // parameter storage
  Parser ip(var); // parser
  ip.ParseFile(config);
  if (args.Int["config_verbose"]) {
    var.Int.Set("verbose", 1);
  }
  const auto extra = args.String["extra"];
  if (extra.length()) {
    std::stringstream buf(extra);
    ip.ParseStream(buf);
  }

  const std::string backend = var.String["backend"];
  if (backend == "local") {
    fassert_equal(
        mpi.GetCommSize(), 1, "\nBackend 'local' requires a single rank.\n");
  }

  try {
    kernel(mpi.GetComm(), var);
  } catch (const std::exception& e) {
    std::cerr << FILELINE + "\nabort after throwing exception\n"
              << e.what() << '\n';
    std::terminate();
  } catch (...) {
    std::cerr << FILELINE + "\nabort after unknown exception\n";
    std::terminate();
  }

  if (isroot) {
    if (var.Int("verbose_conf_reads", 0)) {
      auto print_reads = [&var, &var_reads](const auto& map) {
        for (auto it = map.cbegin(); it != map.cend(); ++it) {
          const auto key = it->first;
          std::cout << (var_reads.count(key) ? var_reads.at(key).load() : 0)
                    << ' ' << map.GetTypeName() << ' ' << key << '\n';
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
      var.ForEachMap([&var, &ignore, &var_reads](const auto& map) {
        for (auto it = map.cbegin(); it != map.cend(); ++it) {
          const auto key = it->first;
          if (!var_reads.count(key) && !ignore.count(key)) {
            std::cout << map.GetTypeName() << ' ' << key << '\n';
          }
        }
      });
    }
  }
  return 0;
}

int RunMpi(
    int argc, const char** argv,
    const std::function<void(MPI_Comm, Vars&)>& kernel) {
  return RunMpi0(argc, argv, kernel);
}
