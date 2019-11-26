#include <mpi.h>
#include <sstream>

#include "cubism.ipp"
#include "geom/mesh.h"
#include "kernel/kernel.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

using U = std::unique_ptr<Distr>;

template <size_t bx, size_t by, size_t bz, class KF>
U Try(MPI_Comm comm, KernelFactory& kf, Vars& var) {
  using M = typename KF::M;
  using Scal = typename M::Scal;
  using Par = GPar<Scal, bx, by, bz, 8>;

  // Check block size
  if (var.Int["bsx"] == bx &&
      var.Int["bsy"] == by &&
      (var.Int["bsz"] == bz || (bz == 2 && var.Int["bsz"] == 1))) {
    // Check kernel
    if (KF* kfd = dynamic_cast<KF*>(&kf)) {
      return U(new Cubism<Par, KF>(comm, *kfd, var));
    }
  }

  return nullptr;
}

U CreateCubism(MPI_Comm comm, KernelFactory& kf, Vars& var) {
  U r;
  using KF = KernelMeshFactory<MeshStructured<double, 3>>;
  // 3D
  if (!r) r = Try<32, 32, 32, KF>(comm, kf, var);
  if (!r) r = Try<16, 16, 16, KF>(comm, kf, var);
  if (!r) r = Try<8, 8, 8, KF>(comm, kf, var);
  // 2D
  if (!r) r = Try<32, 32, 2, KF>(comm, kf, var);
  if (!r) r = Try<16, 16, 2, KF>(comm, kf, var);
  if (!r) r = Try<8, 8, 2, KF>(comm, kf, var);
  if (!r) {
    std::stringstream ss;
    ss << __func__ << ": no instance with "
        << "bs="
        << var.Int["bsx"] << " " << var.Int["bsy"] << " " << var.Int["bsz"]
        << " hl=" << var.Int["hl"];
    throw std::runtime_error(ss.str());
  }
  return r;
}

