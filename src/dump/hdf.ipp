// Created by Petr Karnakov on 04.07.2020
// Copyright 2020 ETH Zurich

#include <hdf5.h>
#include <mpi.h>
#include <fstream>

#include "hdf.h"
#include "util/logger.h"

template <class Scal>
generic::Vect<Scal, 1>& GetVect(Scal& src) {
  return reinterpret_cast<generic::Vect<Scal, 1>&>(src);
}

template <class Scal, size_t dim>
generic::Vect<Scal, dim>& GetVect(generic::Vect<Scal, dim>& src) {
  return src;
}

template <class Scal>
const generic::Vect<Scal, 1>& GetConstVect(const Scal& src) {
  return reinterpret_cast<const generic::Vect<Scal, 1>&>(src);
}

template <class Scal, size_t dim>
const generic::Vect<Scal, dim>& GetConstVect(
    const generic::Vect<Scal, dim>& src) {
  return src;
}

template <class Scal>
constexpr size_t GetDim(Scal*) {
  return 1;
}

template <class Scal, size_t dim>
constexpr size_t GetDim(generic::Vect<Scal, dim>*) {
  return dim;
}

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
template <class Field>
void Hdf<M>::Write(const Field& fc, std::string path, M& m, std::string dname) {
  auto sem = m.GetSem();
  const size_t dim = GetDim((typename Field::Value*)nullptr);
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
    data.reserve(bc.size() * dim);
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        data.push_back(GetConstVect(fc[c])[d]);
      }
    }

    m.GatherToLead(&ctx->origin);
    m.GatherToLead(&ctx->size);
    m.GatherToLead(&ctx->data);
  }
  if (sem("write") && m.IsLead()) {
    const auto hdf_type =
        (sizeof(Scal) == 4 ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE);
    const MPI_Comm comm = m.GetMpiComm();

    fassert_equal(MPI_Barrier(comm), MPI_SUCCESS, "MPI_Barrier faield");
    fassert_equal(H5open(), 0, "H5open failed");

    const hid_t file = [&comm, &path]() {
      const Fapl fapl(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL);
      int rc = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
      fassert_equal(rc >= 0, 1, "H5Fcreate failed for '" + path + "'");
      return rc;
    }();

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

    const auto gsize = harray(m.GetGlobalSize(), dim);
    const hid_t fspace = H5Screate_simple(kDim, gsize.data(), NULL);
    const hid_t dataset = H5Dcreate(
        file, dname.c_str(), hdf_type, fspace, H5P_DEFAULT, H5P_DEFAULT,
        H5P_DEFAULT);
    for (size_t i = 0; i < nblocks; ++i) {
      const auto count = harray(ctx->size[i], dim);
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
template <class Field>
void Hdf<M>::Read(Field& fc, std::string path, M& m, std::string dname) {
  auto sem = m.GetSem();
  const size_t dim = GetDim((typename Field::Value*)nullptr);
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
    ctx->data.resize(bc.size() * dim);
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
    fassert(file >= 0, "cannot open file '" + path + "'");

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
    fassert_equal(V(dims), V(harray(m.GetGlobalSize(), dim)));

    for (size_t i = 0; i < nblocks; ++i) {
      const auto count = harray(ctx->size[i], dim);
      const auto offset = harray(ctx->origin[i], 0);
      const hid_t mspace = H5Screate_simple(kDim, count.data(), NULL);

      H5Sselect_hyperslab(
          fspace, H5S_SELECT_SET, offset.data(), NULL, count.data(), NULL);

      const hid_t fapl = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(fapl, H5FD_MPIO_COLLECTIVE);

      fassert_equal(ctx->dataptr[i]->size(), size_t(ctx->size[i].prod() * dim));
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
      auto& vect = GetVect(fc[c]);
      for (size_t d = 0; d < dim; ++d) {
        vect[d] = ctx->data[i++];
      }
    }
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
  fassert(file >= 0, "cannot open file '" + path + "'");
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
    std::string xmfpath, std::string name, Vect origin, Vect spacing, MIdx dims,
    std::string hdfpath, std::string dname) {
  std::ofstream f(xmfpath);
  f.precision(20);
  f << "<?xml version='1.0' ?>\n";
  f << "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n";
  f << "<Xdmf Version='2.0'>\n";

  f << " <Domain>\n";
  f << "   <Grid Name='mesh' GridType='Uniform'>\n";
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
  WriteXmf(
      xmfpath, name, Vect(0), m.GetCellSize(), m.GetGlobalSize(), hdfpath,
      dname);
}
