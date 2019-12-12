#include <mpi.h>
#include <sstream>

#include "geom/mesh.h"
#include "kernel/kernel.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

#include "cubismnc.ipp"

// XXX: removing 'static' leads to symbol collision with cubism.cpp
template <size_t bx, size_t by, size_t bz, class KF>
static void Try(
    MPI_Comm comm, KF& kf, Vars& var, std::unique_ptr<DistrMesh<KF>>& r) {
  using M = typename KF::M;
  using Scal = typename M::Scal;
  if (r) return;
  // using Par0 = GPar<Scal, bx, by, bz, 0>; // 0 halos on either side
  // using Par1 = GPar<Scal, bx, by, bz, 1>; // 1 halos on either side
  using Par2 = GPar<Scal, bx, by, bz, 2>; // 2 halos on either side
  // using Par3 = GPar<Scal, bx, by, bz, 3>; // 3 halos on either side

  if (var.Int["bsx"] == bx && var.Int["bsy"] == by && var.Int["bsz"] == bz) {
    switch (var.Int["hl"]) {
      // case 0:
      //  r.reset(new Cubismnc<Par0, KF>(comm, kf, var));
      // case 1:
      //  r.reset(new Cubismnc<Par1, KF>(comm, kf, var));
      case 2:
        r.reset(new Cubismnc<Par2, KF>(comm, kf, var));
      // case 3:
      //  r.reset(new Cubismnc<Par3, KF>(comm, kf, var));
      default:
        break;
    }
  }
}

template <class KF>
std::unique_ptr<DistrMesh<KF>> CreateCubismnc(
    MPI_Comm comm, KF& kf, Vars& var) {
  std::unique_ptr<DistrMesh<KF>> r;
  // 3D
  Try<32, 32, 32>(comm, kf, var, r);
  Try<16, 16, 16>(comm, kf, var, r);
  Try<8, 8, 8>(comm, kf, var, r);
  // 2D
  Try<32, 32, 1>(comm, kf, var, r);
  Try<16, 16, 1>(comm, kf, var, r);
  Try<8, 8, 1>(comm, kf, var, r);
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

template std::unique_ptr<DistrMesh<KF>> CreateCubismnc<KF>(
    MPI_Comm, KF&, Vars&);

// vim: expandtab:smarttab:sw=2:ts=2
