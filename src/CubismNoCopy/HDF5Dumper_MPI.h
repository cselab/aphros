/*
 *  Created by Fabian Wermelinger on 12.11.2019.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <hdf5.h>
#include <array>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#ifdef _FLOAT_PRECISION_
#define HDF_REAL H5T_NATIVE_FLOAT
#else
#define HDF_REAL H5T_NATIVE_DOUBLE
#endif

#include "BlockInfo.h"

// Creates XMF file with uniform mesh and scalar data attribute
// linked to HDF field with name 'data'.
// xmfpath: path to output
// name: name of attribute
// origin,spacing,dims: mesh parameters
// hdfpath: path to existing hdf
// hdfdims: dimensions of hdf, either (nz,ny,nx,1) or (ny,nx,1)
void WriteXmf(
    std::string xmfpath, std::string name, double /*time*/,
    const std::array<double, 3>& origin, const std::array<double, 3>& spacing,
    const std::array<size_t, 3>& dims, std::string hdfpath,
    const std::array<size_t, 4>& hdfdims) {
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
  f << hdfdims[0] << " " << hdfdims[1] << " " << hdfdims[2] << " "
    << hdfdims[3];
  f << "' NumberType='Float' Precision='" << sizeof(Real)
    << "' Format='HDF'>\n";
  f << "        " << hdfpath << ":/data\n";
  f << "       </DataItem>\n";
  f << "     </Attribute>\n";
  f << "   </Grid>\n";
  f << " </Domain>\n";
  f << "</Xdmf>\n";
}

template <typename Streamer, typename TGrid>
void DumpHDF5_MPI(
    std::vector<typename Streamer::B>& blocks, const int aos_idx,
    const TGrid& grid, int /*iCounter*/, Real absTime, std::string f_name,
    std::string dump_path, std::vector<Real> origin, std::vector<Real> spacing,
    bool bXMF) {
  using B = typename TGrid::BlockType;

  int rank;
  const std::string fullname = f_name + Streamer::EXT;
  std::string filename = dump_path + "/" + fullname + ".h5";

  MPI_Comm comm = grid.getCartComm();
  MPI_Comm_rank(comm, &rank);

  bool isroot = (rank == 0);

  int coords[3];
  grid.peindex(coords);

  hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

  MPI_Barrier(comm);

  // startup file
  // H5open();
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  H5Pclose(fapl_id);

  ///////////////////////////////////////////////////////////////////////////
  // write data
  const size_t NX = grid.getResidentBlocksPerDimension(0) * B::sizeX;
  const size_t NY = grid.getResidentBlocksPerDimension(1) * B::sizeY;
  const size_t NZ = grid.getResidentBlocksPerDimension(2) * B::sizeZ;
  const size_t NCHANNELS = Streamer::NCHANNELS;

  Real* array_all = new Real[NX * NY * NZ * NCHANNELS];

  std::vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

  const size_t sX = 0;
  const size_t sY = 0;
  const size_t sZ = 0;

  const size_t eX = B::sizeX;
  const size_t eY = B::sizeY;
  const size_t eZ = B::sizeZ;

  hsize_t count[4] = {
      grid.getResidentBlocksPerDimension(2) * B::sizeZ,
      grid.getResidentBlocksPerDimension(1) * B::sizeY,
      grid.getResidentBlocksPerDimension(0) * B::sizeX, NCHANNELS};

  hsize_t dims[4] = {
      grid.getBlocksPerDimension(2) * B::sizeZ,
      grid.getBlocksPerDimension(1) * B::sizeY,
      grid.getBlocksPerDimension(0) * B::sizeX, NCHANNELS};

  if (isroot) {
    std::cout << "HDF5 dump '" << fullname << "': "
              << "one rank "
              << double(NX * NY * NZ * NCHANNELS * sizeof(Real)) / (1 << 20)
              << " MB"
              << ", total "
              << double(dims[0] * dims[1] * dims[2] * dims[3] * sizeof(Real)) /
                     (1 << 20)
              << " MB" << std::endl;
  }

  hsize_t offset[4] = {
      coords[2] * grid.getResidentBlocksPerDimension(2) * B::sizeZ,
      coords[1] * grid.getResidentBlocksPerDimension(1) * B::sizeY,
      coords[0] * grid.getResidentBlocksPerDimension(0) * B::sizeX, 0};

  assert(vInfo_local.size() == blocks.size());
  for (size_t i = 0; i < vInfo_local.size(); i++) {
    BlockInfo& info = vInfo_local[i];
    const size_t idx[3] = {
        (size_t)info.index[0], (size_t)info.index[1], (size_t)info.index[2]};
    B& b = blocks[i];
    Streamer streamer(b, aos_idx);

    for (size_t iz = sZ; iz < eZ; iz++) {
      const size_t gz = idx[2] * B::sizeZ + iz;
      for (size_t iy = sY; iy < eY; iy++) {
        const size_t gy = idx[1] * B::sizeY + iy;
        for (size_t ix = sX; ix < eX; ix++) {
          const size_t gx = idx[0] * B::sizeX + ix;
          const size_t idx = NCHANNELS * (gx + NX * (gy + NY * gz));
          assert(idx < NX * NY * NZ * NCHANNELS);
          Real* const ptr = array_all + idx;
          Real output[NCHANNELS];
          for (size_t i = 0; i < NCHANNELS; ++i)
            output[i] = 0;
          streamer.operate(ix, iy, iz, (Real*)output);
          for (size_t i = 0; i < NCHANNELS; ++i)
            ptr[i] = output[i];
        }
      }
    }
  }

  fapl_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

  fspace_id = H5Screate_simple(4, dims, NULL);
  dataset_id = H5Dcreate(
      file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT);
  H5Sclose(fspace_id);

  fspace_id = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  mspace_id = H5Screate_simple(4, count, NULL);
  H5Dwrite(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);
  H5Sclose(fspace_id);

  H5Sclose(mspace_id);
  H5Dclose(dataset_id);
  H5Pclose(fapl_id);
  H5Fclose(file_id);

  delete[] array_all;

  if (bXMF && isroot) {
    size_t mesh_dims[3] = {
        grid.getBlocksPerDimension(0) * B::sizeX,
        grid.getBlocksPerDimension(1) * B::sizeY,
        grid.getBlocksPerDimension(2) * B::sizeZ};

    const std::string xmfpath = dump_path + "/" + fullname + ".xmf";
    WriteXmf(
        xmfpath, Streamer::NAME, absTime, {origin[0], origin[1], origin[2]},
        {spacing[0], spacing[1], spacing[2]},
        {mesh_dims[0], mesh_dims[1], mesh_dims[2]}, fullname + ".h5",
        {dims[0], dims[1], dims[2], dims[3]});
  }
}

template <typename Streamer, typename TGrid>
void ReadHDF5_MPI(
    std::vector<typename Streamer::B>& blocks, const int aos_idx,
    const TGrid& grid, const std::string f_name,
    const std::string dump_path = ".") {
  typedef typename TGrid::BlockType B;

  int rank;
  char filename[256];
  herr_t status;
  hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

  MPI_Comm comm = grid.getCartComm();
  MPI_Comm_rank(comm, &rank);

  int coords[3];
  grid.peindex(coords);

  const size_t NX = grid.getResidentBlocksPerDimension(0) * B::sizeX;
  const size_t NY = grid.getResidentBlocksPerDimension(1) * B::sizeY;
  const size_t NZ = grid.getResidentBlocksPerDimension(2) * B::sizeZ;
  static const size_t NCHANNELS = Streamer::NCHANNELS;

  Real* array_all = new Real[NX * NY * NZ * NCHANNELS];

  std::vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

  static const int sX = 0;
  static const int sY = 0;
  static const int sZ = 0;

  const int eX = B::sizeX;
  const int eY = B::sizeY;
  const int eZ = B::sizeZ;

  hsize_t count[4] = {
      grid.getResidentBlocksPerDimension(2) * B::sizeZ,
      grid.getResidentBlocksPerDimension(1) * B::sizeY,
      grid.getResidentBlocksPerDimension(0) * B::sizeX, NCHANNELS};

  /*
  hsize_t dims[4] = {
      grid.getBlocksPerDimension(2)*B::sizeZ,
      grid.getBlocksPerDimension(1)*B::sizeY,
      grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};
      */

  hsize_t offset[4] = {
      coords[2] * grid.getResidentBlocksPerDimension(2) * B::sizeZ,
      coords[1] * grid.getResidentBlocksPerDimension(1) * B::sizeY,
      coords[0] * grid.getResidentBlocksPerDimension(0) * B::sizeX, 0};

  sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name.c_str());

  H5open();
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
  file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
  status = H5Pclose(fapl_id);

  dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
  fapl_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

  fspace_id = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

  mspace_id = H5Screate_simple(4, count, NULL);
  status =
      H5Dread(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

  assert(vInfo_local.size() == blocks.size());
  for (int i = 0; i < int(vInfo_local.size()); i++) {
    BlockInfo& info = vInfo_local[i];
    const int idx[3] = {info.index[0], info.index[1], info.index[2]};
    B& b = blocks[i];
    Streamer streamer(b, aos_idx);

    for (int iz = sZ; iz < eZ; iz++)
      for (int iy = sY; iy < eY; iy++)
        for (int ix = sX; ix < eX; ix++) {
          const int gx = idx[0] * B::sizeX + ix;
          const int gy = idx[1] * B::sizeY + iy;
          const int gz = idx[2] * B::sizeZ + iz;

          Real* const ptr_input =
              array_all + NCHANNELS * (gx + NX * (gy + NY * gz));
          streamer.operate(ptr_input, ix, iy, iz);
        }
  }

  status = H5Pclose(fapl_id);
  status = H5Dclose(dataset_id);
  status = H5Sclose(fspace_id);
  status = H5Sclose(mspace_id);
  status = H5Fclose(file_id);

  H5close();

  delete[] array_all;
}
