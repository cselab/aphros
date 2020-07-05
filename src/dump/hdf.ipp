// Created by Petr Karnakov on 04.07.2020
// Copyright 2020 ETH Zurich

#include <hdf5.h>
#include <mpi.h>

#include "hdf.h"
#include "util/logger.h"

template <class M>
void Hdf<M>::Write(
    const FieldCell<typename M::Scal>& fc, std::string path, M& m) {
  auto sem = m.GetSem();
  struct {
    std::vector<MIdx> origin;
    std::vector<MIdx> size;
    std::vector<std::vector<Scal>> data;
  } * ctx(sem);

  if (sem("gather")) {
    ctx->origin.push_back(m.GetInBlockCells().GetBegin());
    ctx->size.push_back(m.GetInBlockCells().GetSize());

    ctx->data.emplace_back();
    auto& data = ctx->data.back();
    for (auto c : m.Cells()) {
      data.push_back(fc[c]);
    }

    using OpCatM = typename M::template OpCatT<MIdx>;
    m.Reduce(std::make_shared<OpCatM>(&ctx->origin));
    m.Reduce(std::make_shared<OpCatM>(&ctx->size));
    using OpCatVT = typename M::template OpCatVT<Scal>;
    m.Reduce(std::make_shared<OpCatVT>(&ctx->data));
  }
  if (sem("write") && m.IsLead()) {
    const auto hdf_type =
        (sizeof(Scal) == 4 ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE);
    MPI_Comm comm = m.GetMpiComm();

    MPI_Barrier(comm);
    H5open();

    hid_t fapl;
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL);
    const hid_t file =
        H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);

    const size_t kComps = 1;
    const size_t nblocks = ctx->data.size();
    fassert_equal(ctx->data.size(), nblocks);
    fassert_equal(ctx->size.size(), nblocks);
    fassert_equal(ctx->origin.size(), nblocks);

    const size_t kDim = 4;

    auto harray = [](MIdx w, size_t s) -> std::array<hsize_t, kDim> {
      return {
          static_cast<hsize_t>(w[2]), static_cast<hsize_t>(w[1]),
          static_cast<hsize_t>(w[0]), s};
    };

    const auto fieldname = "data";
    const auto gsize = harray(m.GetGlobalSize(), kComps);
    const hid_t fspace = H5Screate_simple(kDim, gsize.data(), NULL);
    const hid_t dataset = H5Dcreate(
        file, fieldname, hdf_type, fspace, H5P_DEFAULT, H5P_DEFAULT,
        H5P_DEFAULT);
    for (size_t i = 0; i < nblocks; ++i) {
      const auto count = harray(ctx->size[i], kComps);
      const auto offset = harray(ctx->origin[i], 0);
      const hid_t mspace = H5Screate_simple(kDim, count.data(), NULL);
      const hid_t fapl = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(fapl, H5FD_MPIO_COLLECTIVE);

      H5Sselect_hyperslab(
          fspace, H5S_SELECT_SET, offset.data(), NULL, count.data(), NULL);
      H5Dwrite(dataset, hdf_type, mspace, fspace, fapl, ctx->data[i].data());

      H5Sclose(mspace);
      H5Pclose(fapl);
    }
    H5Sclose(fspace);
    H5Dclose(dataset);
    H5Fclose(file);
  }
  if (sem()) {}
}
