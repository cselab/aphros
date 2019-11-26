#include <mpi.h>
#include <sstream>

#include "local.ipp"
#include "geom/mesh.h"
#include "kernel/kernel.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

using U = std::unique_ptr<Distr>;

template <class KF>
static U Try(MPI_Comm comm, KernelFactory& kf, Vars& var) {
  if (KF* kfd = dynamic_cast<KF*>(&kf)) {
    return U(new Local<KF>(comm, *kfd, var));
  }
  return nullptr;
}

U CreateLocal(MPI_Comm comm, KernelFactory& kf, Vars& var) {
  U r;
  using KF = KernelMeshFactory<MeshStructured<double, 3>>;
  if (!r) r = Try<KF>(comm, kf, var);
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
