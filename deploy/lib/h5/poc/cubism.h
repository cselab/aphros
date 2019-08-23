#include <cassert>
#include <cstdio>
#include <vector>
#include <string>

#include <hdf5.h>

#ifdef _FLOAT_PRECISION_
#define HDF_REAL H5T_NATIVE_FLOAT
#else
#define HDF_REAL H5T_NATIVE_DOUBLE
#endif

#include "BlockInfo.h"
#include "MeshMap.h"

template<typename TGrid, typename Streamer>
void DumpHDF5_MPI(TGrid &grid, int /*iCounter*/, Real absTime, 
                  std::string f_name, std::string dump_path=".", 
                  bool bXMF=true) {
  using B = typename TGrid::BlockType;

  int rank;
  std::string fullname = f_name + Streamer::EXT;
  char filename[256];

  sprintf(filename, "%s/%s.h5", dump_path.c_str(), fullname.c_str());

  MPI_Comm comm = grid.getCartComm();
  MPI_Comm_rank(comm, &rank);

  bool isroot = (rank == 0);

  int coords[3];
  grid.peindex(coords);

  hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;


  ///////////////////////////////////////////////////////////////////////////
  // write mesh
  std::vector<int> mesh_dims;
  std::vector<std::string> dset_name;
  dset_name.push_back("/vx");
  dset_name.push_back("/vy");
  dset_name.push_back("/vz");

  if (isroot)
  {
    //H5open();
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    H5Pclose(fapl_id);

    for (size_t i = 0; i < 3; ++i)
    {
      MeshMap<B>& m = grid.getMeshMap(i);
      std::vector<double> vertices(m.ncells()+1, m.start());
      mesh_dims.push_back(vertices.size());

      for (size_t j = 0; j < m.ncells(); ++j)
        vertices[j+1] = vertices[j] + m.cell_width(j);

      hsize_t dim[1] = {vertices.size()};

      fspace_id = H5Screate_simple(1, dim, NULL);
      dataset_id = H5Dcreate(file_id, dset_name[i].c_str(), H5T_NATIVE_DOUBLE, 
          fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(fspace_id);

      H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
          H5P_DEFAULT, vertices.data());

      H5Dclose(dataset_id);
    }

    // shutdown h5 file
    H5Fclose(file_id);
    //H5close();
  }
  MPI_Barrier(comm);

  ///////////////////////////////////////////////////////////////////////////
  // startup file
  //H5open();
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL);
  file_id = H5Fopen(filename, H5F_ACC_RDWR, fapl_id);
  H5Pclose(fapl_id);

  ///////////////////////////////////////////////////////////////////////////
  // write data
  size_t NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
  size_t NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
  size_t NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
  size_t NCHANNELS = Streamer::NCHANNELS;

  Real * array_all = new Real[NX * NY * NZ * NCHANNELS];

  std::vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

  size_t sX = 0;
  size_t sY = 0;
  size_t sZ = 0;

  size_t eX = B::sizeX;
  size_t eY = B::sizeY;
  size_t eZ = B::sizeZ;

  hsize_t count[4] = {
    grid.getResidentBlocksPerDimension(2)*B::sizeZ,
    grid.getResidentBlocksPerDimension(1)*B::sizeY,
    grid.getResidentBlocksPerDimension(0)*B::sizeX, NCHANNELS};

  hsize_t dims[4] = {
    grid.getBlocksPerDimension(2)*B::sizeZ,
    grid.getBlocksPerDimension(1)*B::sizeY,
    grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};

  if (isroot) {
    std::cout 
        << "HDF5 dump '"
        << fullname << "': " 
        << "one rank " 
        << double(NX * NY * NZ * NCHANNELS * sizeof(Real)) 
            / (1 << 20) 
        << " MB"
        << ", total "
        << double(dims[0] * dims[1] * dims[2] * dims[3] * sizeof(Real)) 
            / (1 << 20) 
        << " MB" << std::endl;
  }

  hsize_t offset[4] = {
    coords[2]*grid.getResidentBlocksPerDimension(2)*B::sizeZ,
    coords[1]*grid.getResidentBlocksPerDimension(1)*B::sizeY,
    coords[0]*grid.getResidentBlocksPerDimension(0)*B::sizeX, 0};

  for(size_t i=0; i<vInfo_local.size(); i++)
  {
    BlockInfo& info = vInfo_local[i];
    size_t idx[3] = 
        {(size_t)info.index[0], (size_t)info.index[1], (size_t)info.index[2]};
    B & b = *(B*)info.ptrBlock;
    Streamer streamer(b);

    for(size_t iz=sZ; iz<eZ; iz++)
    {
      size_t gz = idx[2]*B::sizeZ + iz;
      for(size_t iy=sY; iy<eY; iy++)
      {
        size_t gy = idx[1]*B::sizeY + iy;
        for(size_t ix=sX; ix<eX; ix++)
        {
          size_t gx = idx[0]*B::sizeX + ix;
          size_t idx = NCHANNELS * (gx + NX * (gy + NY * gz));
          assert(idx < NX * NY * NZ * NCHANNELS);
          Real * ptr = array_all + idx;
          Real output[NCHANNELS];
          for(size_t i=0; i<NCHANNELS; ++i)
            output[i] = 0;
          streamer.operate(ix, iy, iz, (Real*)output);
          for(size_t i=0; i<NCHANNELS; ++i)
            ptr[i] = output[i];
        }
      }
    }
  }

  fapl_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

  fspace_id = H5Screate_simple(4, dims, NULL);
  dataset_id = H5Dcreate(file_id, "data", HDF_REAL, fspace_id, 
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

  delete [] array_all;
}

