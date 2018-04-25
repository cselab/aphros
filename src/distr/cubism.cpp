#include <mpi.h>

#include "cubism.h"
#include "geom/mesh.h"
#include "kernel/kernel.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

template <size_t bx, size_t by, size_t bz, class KF>
Distr* Try(MPI_Comm comm, KernelFactory& kf, Vars& par) {
  using M = typename KF::M;
  using Scal = typename M::Scal;
  using Par = GPar<Scal, bx, by, bz, 8>;

  // Check block size
  if (par.Int["bsx"] == bx && 
      par.Int["bsy"] == by && 
      par.Int["bsz"] == bz) {
    // Check kernel
    if (KF* kfd = dynamic_cast<KF*>(&kf)) {
      return new Cubism<Par, KF>(comm, *kfd, par);
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
  if (!r) r = Try<32, 32, 2, KF>(comm, kf, par);
  if (!r) r = Try<16, 16, 2, KF>(comm, kf, par);
  if (!r) r = Try<8, 8, 2, KF>(comm, kf, par);
  if (!r) {
    std::cerr << "CreateCubism(): no instance with "
      << "bs=(" 
      << par.Int["bsx"] << ","
      << par.Int["bsy"] << ","
      << par.Int["bsz"] << ")"
      << std::endl;
    assert(false);
  }
  return r;
}

