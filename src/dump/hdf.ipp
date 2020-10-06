// Created by Petr Karnakov on 04.07.2020
// Copyright 2020 ETH Zurich

#include <hdf5.h>
#include <mpi.h>
#include <fstream>

#include "hdf.h"
#include "util/logger.h"

class Fapl {
 public:
  Fapl(hid_t cls_id) : fapl_(H5Pcreate(cls_id)) {}
  ~Fapl() {
    H5Pclose(fapl_);
  }
  operator hid_t() const {
    return fapl_;
  }

 private:
  const hid_t fapl_;
};

template <class M>
void Hdf<M>::Write(
    const FieldCell<typename M::Scal>& fc, std::string path, M& m,
    std::string dname) {
  auto sem = m.GetSem();
  struct {
    std::vector<MIdx> origin;
    std::vector<MIdx> size;
    std::vector<std::vector<Scal>> data;
  } * ctx(sem);

  if (sem("gather")) {
    const auto bc = m.GetInBlockCells();
    ctx->origin.push_back(bc.GetBegin());
    ctx->size.push_back(bc.GetSize());

    ctx->data.emplace_back();
    auto& data = ctx->data.back();
    for (auto c : m.Cells()) {
      data.push_back(fc[c]);
    }

    m.GatherToLead(&ctx->origin);
    m.GatherToLead(&ctx->size);
    m.GatherToLead(&ctx->data);
  }
  if (sem("write") && m.IsLead()) {
    const auto hdf_type =
        (sizeof(Scal) == 4 ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE);
    const MPI_Comm comm = m.GetMpiComm();

    MPI_Barrier(comm);
    H5open();

    const hid_t file = [&comm, &path]() {
      const Fapl fapl(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL);
      return H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    }();

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

    const auto gsize = harray(m.GetGlobalSize(), kComps);
    const hid_t fspace = H5Screate_simple(kDim, gsize.data(), NULL);
    const hid_t dataset = H5Dcreate(
        file, dname.c_str(), hdf_type, fspace, H5P_DEFAULT, H5P_DEFAULT,
        H5P_DEFAULT);
    for (size_t i = 0; i < nblocks; ++i) {
      const auto count = harray(ctx->size[i], kComps);
      const auto offset = harray(ctx->origin[i], 0);

      H5Sselect_hyperslab(
          fspace, H5S_SELECT_SET, offset.data(), NULL, count.data(), NULL);

      const hid_t mspace = H5Screate_simple(kDim, count.data(), NULL);
      const hid_t fapl = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(fapl, H5FD_MPIO_COLLECTIVE);

      H5Dwrite(dataset, hdf_type, mspace, fspace, fapl, ctx->data[i].data());

      H5Pclose(fapl);
      H5Sclose(mspace);
    }
    H5Sclose(fspace);
    H5Dclose(dataset);
    H5Fclose(file);
  }
  if (sem()) { // XXX empty stage
  }
}

template <class M>
void Hdf<M>::Read(
    FieldCell<typename M::Scal>& fc, std::string path, M& m,
    std::string dname) {
  auto sem = m.GetSem();
  struct {
    std::vector<MIdx> origin;
    std::vector<MIdx> size;
    std::vector<Scal> data;
    std::vector<std::vector<Scal>*> dataptr;
  } * ctx(sem);

  if (sem("gather")) {
    const auto bc = m.GetInBlockCells();
    ctx->origin.push_back(bc.GetBegin());
    ctx->size.push_back(bc.GetSize());
    ctx->data.resize(bc.size());
    ctx->dataptr.push_back(&ctx->data);

    m.GatherToLead(&ctx->origin);
    m.GatherToLead(&ctx->size);
    m.GatherToLead(&ctx->dataptr);
  }
  if (sem("read") && m.IsLead()) {
    const auto hdf_type =
        (sizeof(Scal) == 4 ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE);
    const MPI_Comm comm = m.GetMpiComm();

    MPI_Barrier(comm);
    H5open();

    const hid_t file = [&comm, &path]() {
      const Fapl fapl(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL);
      return H5Fopen(path.c_str(), H5F_ACC_RDONLY, fapl);
    }();
    if (file < 0) {
      throw std::runtime_error(FILELINE + ": cannot open file '" + path + "'");
    }

    const size_t kComps = 1;
    const size_t nblocks = ctx->dataptr.size();
    fassert_equal(ctx->dataptr.size(), nblocks);
    fassert_equal(ctx->size.size(), nblocks);
    fassert_equal(ctx->origin.size(), nblocks);

    const size_t kDim = 4;

    auto harray = [](MIdx w, size_t s) -> std::array<hsize_t, kDim> {
      return {
          static_cast<hsize_t>(w[2]), static_cast<hsize_t>(w[1]),
          static_cast<hsize_t>(w[0]), s};
    };

    const hid_t dataset = H5Dopen2(file, dname.c_str(), H5P_DEFAULT);
    const hid_t fspace = H5Dget_space(dataset);

    const size_t ndims = H5Sget_simple_extent_ndims(fspace);
    fassert_equal(ndims, kDim);

    hsize_t dims[kDim];
    H5Sget_simple_extent_dims(fspace, dims, NULL);
    using V = generic::Vect<size_t, 4>;
    fassert_equal(V(dims), V(harray(m.GetGlobalSize(), kComps)));

    for (size_t i = 0; i < nblocks; ++i) {
      const auto count = harray(ctx->size[i], kComps);
      const auto offset = harray(ctx->origin[i], 0);
      const hid_t mspace = H5Screate_simple(kDim, count.data(), NULL);

      H5Sselect_hyperslab(
          fspace, H5S_SELECT_SET, offset.data(), NULL, count.data(), NULL);

      const hid_t fapl = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(fapl, H5FD_MPIO_COLLECTIVE);

      fassert_equal(ctx->dataptr[i]->size(), size_t(ctx->size[i].prod()));
      H5Dread(dataset, hdf_type, mspace, fspace, fapl, ctx->dataptr[i]->data());

      H5Pclose(fapl);
      H5Sclose(mspace);
    }
    H5Sclose(fspace);
    H5Dclose(dataset);
    H5Fclose(file);
  }
  if (sem("copy")) {
    size_t i = 0;
    fc.Reinit(m);
    for (auto c : m.Cells()) {
      fc[c] = ctx->data[i++];
    }
    m.Comm(&fc);
  }
  if (sem()) { // XXX empty stage
  }
}

template <class M>
std::vector<size_t> Hdf<M>::GetShape(std::string path, std::string dname) {
  H5open();
  const hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, fapl);
  H5Pclose(fapl);
  if (file < 0) {
    throw std::runtime_error(FILELINE + ": cannot open file '" + path + "'");
  }
  const hid_t dataset = H5Dopen2(file, dname.c_str(), H5P_DEFAULT);
  const hid_t fspace = H5Dget_space(dataset);
  const size_t ndims = H5Sget_simple_extent_ndims(fspace);
  std::vector<hsize_t> dims(ndims);
  H5Sget_simple_extent_dims(fspace, dims.data(), NULL);
  H5Sclose(fspace);
  H5Dclose(dataset);
  H5Fclose(file);
  return {dims.begin(), dims.end()};
}

template <class M>
void Hdf<M>::WriteXmf(
    std::string xmfpath, std::string name, const std::array<double, 3>& origin,
    const std::array<double, 3>& spacing, const std::array<size_t, 3>& dims,
    std::string hdfpath, std::string dname) {
  std::ofstream f(xmfpath);
  f.precision(20);
  f << "<?xml version='1.0' ?>\n";
  f << "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n";
  f << "<Xdmf Version='2.0'>\n";

  f << " <Domain>\n";
  f << "   <Grid GridType='Uniform'>\n";
  f << "     <Topology TopologyType='3DCORECTMesh' Dimensions='";
  f << dims[2] + 1 << " " << dims[1] + 1 << " " << dims[0] + 1 << "'/>\n\n";

  f << "     <Geometry GeometryType='ORIGIN_DXDYDZ'>\n";
  f << "       <DataItem Name='Origin' Dimensions='3' NumberType='Float' "
       "Precision='8' Format='XML'>\n";
  f << "         " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << "\n";
  f << "       </DataItem>\n";
  f << "       <DataItem Name='Spacing' Dimensions='3' NumberType='Float' "
       "Precision='8' Format='XML'>\n";
  f << "         " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2]
    << "\n";
  f << "       </DataItem>\n";
  f << "     </Geometry>\n\n";

  f << "     <Attribute Name='" << name
    << "' AttributeType='Scalar' Center='Cell'>\n";
  f << "       <DataItem Dimensions='";
  f << dims[2] << " " << dims[1] << " " << dims[0] << " " << 1;
  f << "' NumberType='Float' Precision='" << sizeof(Scal)
    << "' Format='HDF'>\n";
  f << "        " << hdfpath << ":/" + dname + "\n";
  f << "       </DataItem>\n";
  f << "     </Attribute>\n";
  f << "   </Grid>\n";
  f << " </Domain>\n";
  f << "</Xdmf>\n";
}

template <class M>
void Hdf<M>::WriteXmf(
    std::string xmfpath, std::string name, std::string hdfpath, const M& m,
    std::string dname) {
  const std::array<double, 3> origin = {0, 0, 0};
  const std::array<double, 3> spacing = m.GetCellSize();
  const auto dims = m.GetGlobalSize();
  WriteXmf(xmfpath, name, origin, spacing, dims, hdfpath, dname);
}
