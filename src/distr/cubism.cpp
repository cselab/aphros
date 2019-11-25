// vim: expandtab:smarttab:sw=2:ts=2
#include <mpi.h>

#include "cubism.ipp"
#include "geom/mesh.h"
#include "kernel/kernel.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

template <size_t bx, size_t by, size_t bz, class KF>
Distr* Try(MPI_Comm comm, KernelFactory& kf, Vars& par) {
  using M = typename KF::M;
  using Scal = typename M::Scal;
  //using Par0 = GPar<Scal, bx, by, bz, 0>; // 0 halos on either side
  //using Par1 = GPar<Scal, bx, by, bz, 1>; // 1 halos on either side
  using Par2 = GPar<Scal, bx, by, bz, 2>; // 2 halos on either side
  //using Par3 = GPar<Scal, bx, by, bz, 3>; // 3 halos on either side

  // Check block size
  if (par.Int["bsx"] == bx &&
      par.Int["bsy"] == by &&
      par.Int["bsz"] == bz) {
    // Check kernel
    if (KF* kfd = dynamic_cast<KF*>(&kf)) {
      switch (par.Int["hl"]) {
        //case 0:
        //  return new Cubism<Par0, KF>(comm, *kfd, par);
        //case 1:
        //  return new Cubism<Par1, KF>(comm, *kfd, par);
        case 2:
          return new Cubism<Par2, KF>(comm, *kfd, par);
        //case 3:
        //  return new Cubism<Par3, KF>(comm, *kfd, par);
        default:
          break;
      }
    }
  }

  return nullptr;
}

Distr* CreateCubism(MPI_Comm comm, KernelFactory& kf, Vars& par) {
  Distr* r = nullptr;
  using KF = KernelMeshFactory<MeshStructured<double, 3>>;
  // 3D
  if (!r) r = Try<32, 32, 32, KF>(comm, kf, par);
  if (!r) r = Try<16, 16, 16, KF>(comm, kf, par);
  if (!r) r = Try<8, 8, 8, KF>(comm, kf, par);
  // 2D
  if (!r) r = Try<32, 32, 1, KF>(comm, kf, par);
  if (!r) r = Try<16, 16, 1, KF>(comm, kf, par);
  if (!r) r = Try<8, 8, 1, KF>(comm, kf, par);
  if (!r) {
    std::cerr << "CreateCubism(): no instance with "
      << "bs=("
      << par.Int["bsx"] << ","
      << par.Int["bsy"] << ","
      << par.Int["bsz"] << "); hl=" << par.Int["hl"]
      << std::endl;
    assert(false);
  }
  return r;
}

