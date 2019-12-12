#include <mpi.h>
#include <sstream>

#include "geom/mesh.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

#include "cubism.ipp"

// XXX: removing 'static' leads to symbol collision with cubismnc.cpp
template <size_t bx, size_t by, size_t bz, class KF>
static void Try(
    MPI_Comm comm, KF& kf, Vars& var, std::unique_ptr<DistrMesh<KF>>& r) {
  if (r) return;
  if (var.Int["bsx"] == bx && var.Int["bsy"] == by &&
      (var.Int["bsz"] == bz || (bz == 2 && var.Int["bsz"] == 1))) {
    using M = typename KF::M;
    using Scal = typename M::Scal;
    using Par = GPar<Scal, bx, by, bz, 8>;
    r.reset(new Cubism<Par, KF>(comm, kf, var));
  }
}

template <class KF>
std::unique_ptr<DistrMesh<KF>> CreateCubism(
    MPI_Comm comm, KF& kf, Vars& var) {
  std::unique_ptr<DistrMesh<KF>> r;
  // 3D
  Try<32, 32, 32>(comm, kf, var, r);
  Try<16, 16, 16>(comm, kf, var, r);
  Try<8, 8, 8>(comm, kf, var, r);
  // 2D
  Try<32, 32, 2>(comm, kf, var, r);
  Try<16, 16, 2>(comm, kf, var, r);
  Try<8, 8, 2>(comm, kf, var, r);
  if (!r) {
    std::stringstream ss;
    ss << __func__ << ": no instance with "
       << "bs=" << var.Int["bsx"] << " " << var.Int["bsy"] << " "
       << var.Int["bsz"] << " hl=" << var.Int["hl"];
    throw std::runtime_error(ss.str());
  }
  return r;
}

using M = MeshStructured<double, 3>;
using KF = KernelMeshFactory<M>;

template std::unique_ptr<DistrMesh<KF>> CreateCubism<KF>(
    MPI_Comm, KF&, Vars&);
