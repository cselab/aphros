#include <mpi.h>
#include <sstream>

#include "geom/mesh.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

#include "local.ipp"

template <class KF>
std::unique_ptr<DistrMesh<KF>> CreateLocal(
    MPI_Comm comm, KF& kf, Vars& var) {
  return std::unique_ptr<DistrMesh<KF>>(new Local<KF>(comm, kf, var));
}

using M = MeshStructured<double, 3>;
using KF = KernelMeshFactory<M>;

template std::unique_ptr<DistrMesh<KF>> CreateLocal<KF>(
    MPI_Comm comm, KF& kf, Vars& var);
