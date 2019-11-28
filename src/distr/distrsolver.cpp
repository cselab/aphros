#include "distrsolver.h"
#include "util/git.h"
#include "util/subcomm.h"
#include "linear/hypresub.h"

static void RunKernelOpenMP(
    MPI_Comm comm_world, MPI_Comm comm_omp, MPI_Comm comm_master,
    std::function<void(MPI_Comm, Vars&)> kernel, Vars& var) {
  int rank_omp;
  MPI_Comm_rank(comm_omp, &rank_omp);

  HypreSub::InitServer(comm_world, comm_omp);
  if (rank_omp == 0) {
    kernel(comm_master, var);
    HypreSub::StopServer();
  } else {
    HypreSub::RunServer();
  }
}

static void RunKernel(
    MPI_Comm comm, std::function<void(MPI_Comm, Vars&)> kernel, Vars& var) {
  kernel(comm, var);
}

int RunMpi(int argc, const char ** argv,
           std::function<void(MPI_Comm, Vars&)> kernel) {
  int prov;
  MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &prov);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool isroot = (!rank);

  Vars var;   // parameter storage
  Parser ip(var); // parser

  std::string fn = "a.conf";
  if (argc == 1) {
    // nop
  } else if (argc == 2) {
    fn = argv[1];
  } else {
    if (isroot) {
      std::cerr << "usage: " << argv[0] << " [a.conf]" << std::endl;
    }
    return 1;
  }

  if (isroot) {
    std::cerr << "hydro " << GetGitRev() << std::endl;
    std::cerr << "msg: " << GetGitMsg() << std::endl;
    std::cerr << "diff: " << GetGitDiff() << std::endl;
    std::cerr << "Loading config from '" << fn << "'" << std::endl;
  }

  std::ifstream f(fn);  // config file
  // Read file and run all commands
  ip.RunAll(f);

  // Print vars on root
  if (isroot) {
    std::cerr << "\n=== config begin ===" << std::endl;
    ip.PrintAll(std::cerr);
    std::cerr << "=== config end ===\n" << std::endl;
  }

  std::string be = var.String["backend"];

  if (be == "local") {
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &comm);
    if (rank == 0) {
      RunKernel(comm, kernel, var);
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
      RunKernel(MPI_COMM_WORLD, kernel, var);
    }
  }

  MPI_Finalize();
  return 0;
}
