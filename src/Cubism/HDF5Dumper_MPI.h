/*
 *  HDF5Dumper_MPI.h
 *  Cubism
 *
 *  Created by Babak Hejazialhosseini on 5/24/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

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
void DumpHDF5_MPI(const TGrid &grid, int /*iCounter*/, Real absTime, 
                  std::string f_name, std::string dump_path=".", 
                  bool bXMF=true) {
  using B = typename TGrid::BlockType;

  int rank;
  const std::string fullname = f_name + Streamer::EXT;
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
      const MeshMap<B>& m = grid.getMeshMap(i);
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
  const size_t NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
  const size_t NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
  const size_t NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
  const size_t NCHANNELS = Streamer::NCHANNELS;

  Real * array_all = new Real[NX * NY * NZ * NCHANNELS];

  std::vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

  const size_t sX = 0;
  const size_t sY = 0;
  const size_t sZ = 0;

  const size_t eX = B::sizeX;
  const size_t eY = B::sizeY;
  const size_t eZ = B::sizeZ;

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
        << "HDF5 dump: "
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
    const size_t idx[3] = 
        {(size_t)info.index[0], (size_t)info.index[1], (size_t)info.index[2]};
    B & b = *(B*)info.ptrBlock;
    Streamer streamer(b);

    for(size_t iz=sZ; iz<eZ; iz++)
    {
      const size_t gz = idx[2]*B::sizeZ + iz;
      for(size_t iy=sY; iy<eY; iy++)
      {
        const size_t gy = idx[1]*B::sizeY + iy;
        for(size_t ix=sX; ix<eX; ix++)
        {
          const size_t gx = idx[0]*B::sizeX + ix;
          const size_t idx = NCHANNELS * (gx + NX * (gy + NY * gz));
          assert(idx < NX * NY * NZ * NCHANNELS);
          Real * const ptr = array_all + idx;
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
  //H5close();

  delete [] array_all;

  if (bXMF && isroot)
  {
    char wrapper[256];
    sprintf(wrapper, "%s/%s.xmf", dump_path.c_str(), fullname.c_str());
    FILE *xmf = 0;
    xmf = fopen(wrapper, "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Time Value=\"%e\"/>\n\n", absTime);
    fprintf(xmf, "     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n\n", mesh_dims[2], mesh_dims[1], mesh_dims[0]);
    fprintf(xmf, "     <Geometry GeometryType=\"VxVyVz\">\n");
    fprintf(xmf, "       <DataItem Name=\"mesh_vx\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", mesh_dims[0]);
    fprintf(xmf, "        %s:/vx\n",(fullname+".h5").c_str());
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Name=\"mesh_vy\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", mesh_dims[1]);
    fprintf(xmf, "        %s:/vy\n",(fullname+".h5").c_str());
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Name=\"mesh_vz\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", mesh_dims[2]);
    fprintf(xmf, "        %s:/vz\n",(fullname+".h5").c_str());
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n\n");
    fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Cell\">\n", 
        Streamer::NAME.c_str(), Streamer::getAttributeName());
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"%lu\" Format=\"HDF\">\n",(int)dims[0], (int)dims[1], (int)dims[2], (int)dims[3], sizeof(Real));
    fprintf(xmf, "        %s:/data\n",(fullname+".h5").c_str());
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
  }
}

template<typename TGrid, typename Streamer>
void ReadHDF5_MPI(TGrid &grid, const std::string f_name, const std::string dump_path=".")
{
    typedef typename TGrid::BlockType B;

    int rank;
    char filename[256];
    herr_t status;
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

    MPI_Comm comm = grid.getCartComm();
    MPI_Comm_rank(comm, &rank);

    int coords[3];
    grid.peindex(coords);

    const size_t NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
    const size_t NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
    const size_t NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
    static const size_t NCHANNELS = Streamer::NCHANNELS;

    Real * array_all = new Real[NX * NY * NZ * NCHANNELS];

    std::vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

    static const int sX = 0;
    static const int sY = 0;
    static const int sZ = 0;

    const int eX = B::sizeX;
    const int eY = B::sizeY;
    const int eZ = B::sizeZ;

    hsize_t count[4] = {
        grid.getResidentBlocksPerDimension(2)*B::sizeZ,
        grid.getResidentBlocksPerDimension(1)*B::sizeY,
        grid.getResidentBlocksPerDimension(0)*B::sizeX, NCHANNELS};

    hsize_t dims[4] = {
        grid.getBlocksPerDimension(2)*B::sizeZ,
        grid.getBlocksPerDimension(1)*B::sizeY,
        grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};

    hsize_t offset[4] = {
        coords[2]*grid.getResidentBlocksPerDimension(2)*B::sizeZ,
        coords[1]*grid.getResidentBlocksPerDimension(1)*B::sizeY,
        coords[0]*grid.getResidentBlocksPerDimension(0)*B::sizeX, 0};

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
    status = H5Dread(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

    for(int i=0; i<vInfo_local.size(); i++)
    {
        BlockInfo& info = vInfo_local[i];
        const int idx[3] = {info.index[0], info.index[1], info.index[2]};
        B & b = *(B*)info.ptrBlock;
        Streamer streamer(b);

        for(int iz=sZ; iz<eZ; iz++)
            for(int iy=sY; iy<eY; iy++)
                for(int ix=sX; ix<eX; ix++)
                {
                    const int gx = idx[0]*B::sizeX + ix;
                    const int gy = idx[1]*B::sizeY + iy;
                    const int gz = idx[2]*B::sizeZ + iz;

                    Real * const ptr_input = array_all + NCHANNELS*(gx + NX * (gy + NY * gz));
                    streamer.operate(ptr_input, ix, iy, iz);
                }
    }

    status = H5Pclose(fapl_id);
    status = H5Dclose(dataset_id);
    status = H5Sclose(fspace_id);
    status = H5Sclose(mspace_id);
    status = H5Fclose(file_id);

    H5close();

    delete [] array_all;
}
