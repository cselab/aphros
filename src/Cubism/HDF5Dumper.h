//
//  HDF5Dumper.h
//  Cubism
//
//  Created by Babak Hejazialhosseini on 1/23/12.
//  Copyright 2011 ETH Zurich. All rights reserved.
//
#pragma once

#include <iostream>
#include <vector>
#include <string>

#ifdef _USE_HDF_
#include <hdf5.h>
#endif

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#ifdef _FLOAT_PRECISION_
#define HDF_REAL H5T_NATIVE_FLOAT
#else
#define HDF_REAL H5T_NATIVE_DOUBLE
#endif

#include "BlockInfo.h"
#include "MeshMap.h"

template<typename TGrid, typename Streamer>
void DumpHDF5(const TGrid &grid, const int iCounter, const Real absTime, const std::string f_name, const std::string dump_path=".", const bool bXMF=true)
{
#ifdef _USE_HDF_
    typedef typename TGrid::BlockType B;

    const std::string fullname = f_name + Streamer::EXT;
    char filename[256];
    sprintf(filename, "%s/%s.h5", dump_path.c_str(), fullname.c_str());

    std::vector<BlockInfo> vInfo_local = grid.getBlocksInfo();

    herr_t status;
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

    ///////////////////////////////////////////////////////////////////////////
    // startup file
    H5open();
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    status = H5Pclose(fapl_id);

    ///////////////////////////////////////////////////////////////////////////
    // write mesh
    std::vector<int> mesh_dims;
    std::vector<std::string> dset_name;
    dset_name.push_back("/vx");
    dset_name.push_back("/vy");
    dset_name.push_back("/vz");
    for (size_t i = 0; i < 3; ++i)
    {
        const MeshMap<B>& m = grid.getMeshMap(i);
        std::vector<double> vertices(m.ncells()+1, m.start());
        mesh_dims.push_back(vertices.size());

        for (int j = 0; j < m.ncells(); ++j)
            vertices[j+1] = vertices[j] + m.cell_width(j);

        hsize_t dim[1] = {vertices.size()};
        fspace_id = H5Screate_simple(1, dim, NULL);
#ifndef _ON_FERMI_
        dataset_id = H5Dcreate(file_id, dset_name[i].c_str(), H5T_NATIVE_DOUBLE, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
        dataset_id = H5Dcreate2(file_id, dset_name[i].c_str(), H5T_NATIVE_DOUBLE, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vertices.data());
        status = H5Sclose(fspace_id);
        status = H5Dclose(dataset_id);
    }

    ///////////////////////////////////////////////////////////////////////////
    // write data
    static const size_t NCHANNELS = Streamer::NCHANNELS;
    const size_t NX = grid.getBlocksPerDimension(0)*B::sizeX;
    const size_t NY = grid.getBlocksPerDimension(1)*B::sizeY;
    const size_t NZ = grid.getBlocksPerDimension(2)*B::sizeZ;

    std::cout << "Allocating " << (NX * NY * NZ * NCHANNELS * sizeof(Real))/(1024.*1024.*1024.) << " GB of HDF5 data" << std::endl;
    Real * array_all = new Real[NX * NY * NZ * NCHANNELS];

    const size_t sX = 0;
    const size_t sY = 0;
    const size_t sZ = 0;

    const size_t eX = B::sizeX;
    const size_t eY = B::sizeY;
    const size_t eZ = B::sizeZ;

    hsize_t count[4] = {
        grid.getBlocksPerDimension(2)*B::sizeZ,
        grid.getBlocksPerDimension(1)*B::sizeY,
        grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};

    hsize_t dims[4] = {
        grid.getBlocksPerDimension(2)*B::sizeZ,
        grid.getBlocksPerDimension(1)*B::sizeY,
        grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};

    hsize_t offset[4] = {0, 0, 0, 0};

#pragma omp parallel for
    for(int i=0; i<(int)vInfo_local.size(); i++)
    {
        BlockInfo& info = vInfo_local[i];
        const size_t idx[3] = {info.index[0], info.index[1], info.index[2]};
        B & b = *(B*)info.ptrBlock;
        Streamer streamer(b);

        for(size_t iz=sZ; iz<eZ; iz++)
            for(size_t iy=sY; iy<eY; iy++)
                for(size_t ix=sX; ix<eX; ix++)
                {
                    Real output[NCHANNELS];
                    for(size_t i=0; i<NCHANNELS; ++i)
                        output[i] = 0;

                    streamer.operate(ix, iy, iz, (Real*)output);

                    const size_t gx = idx[0]*B::sizeX + ix;
                    const size_t gy = idx[1]*B::sizeY + iy;
                    const size_t gz = idx[2]*B::sizeZ + iz;

                    Real * const ptr = array_all + NCHANNELS*(gx + NX * (gy + NY * gz));

                    for(size_t i=0; i<NCHANNELS; ++i)
                        ptr[i] = output[i];
                }
    }

    fapl_id = H5Pcreate(H5P_DATASET_XFER);

    fspace_id = H5Screate_simple(4, dims, NULL);
#ifndef _ON_FERMI_
    dataset_id = H5Dcreate(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    dataset_id = H5Dcreate2(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

    fspace_id = H5Dget_space(dataset_id);

    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    mspace_id = H5Screate_simple(4, count, NULL);

    status = H5Dwrite(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

    status = H5Sclose(mspace_id);
    status = H5Sclose(fspace_id);
    status = H5Dclose(dataset_id);
    status = H5Pclose(fapl_id);

    // shutdown h5 file
    status = H5Fclose(file_id);
    H5close();

    delete [] array_all;

    if (bXMF)
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
        fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"%s\" Center=\"Cell\">\n", Streamer::getAttributeName());
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",(int)dims[0], (int)dims[1], (int)dims[2], (int)dims[3], sizeof(Real));
        fprintf(xmf, "        %s:/data\n",(fullname+".h5").c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "   </Grid>\n");
        fprintf(xmf, " </Domain>\n");
        fprintf(xmf, "</Xdmf>\n");
        fclose(xmf);
    }
#else
#warning USE OF HDF WAS DISABLED AT COMPILE TIME
#endif
}


template<typename TGrid, typename Streamer>
void ReadHDF5(TGrid &grid, const std::string f_name, const std::string read_path=".")
{
#ifdef _USE_HDF_
    typedef typename TGrid::BlockType B;

    char filename[256];
    herr_t status;
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

    const int NX = grid.getBlocksPerDimension(0)*B::sizeX;
    const int NY = grid.getBlocksPerDimension(1)*B::sizeY;
    const int NZ = grid.getBlocksPerDimension(2)*B::sizeZ;
    static const int NCHANNELS = Streamer::NCHANNELS;

    Real * array_all = new Real[NX * NY * NZ * NCHANNELS];

    std::vector<BlockInfo> vInfo_local = grid.getBlocksInfo();

    const int sX = 0;
    const int sY = 0;
    const int sZ = 0;

    const int eX = B::sizeX;
    const int eY = B::sizeY;
    const int eZ = B::sizeZ;

    hsize_t count[4] = {
        grid.getBlocksPerDimension(2)*B::sizeZ,
        grid.getBlocksPerDimension(1)*B::sizeY,
        grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};

    hsize_t dims[4] = {
        grid.getBlocksPerDimension(2)*B::sizeZ,
        grid.getBlocksPerDimension(1)*B::sizeY,
        grid.getBlocksPerDimension(0)*B::sizeX, NCHANNELS};

    hsize_t offset[4] = {0, 0, 0, 0};

    sprintf(filename, "%s/%s.h5", read_path.c_str(), f_name.c_str());

    H5open();
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
    status = H5Pclose(fapl_id);

    dataset_id = H5Dopen2(file_id, "data", H5P_DEFAULT);
    fapl_id = H5Pcreate(H5P_DATASET_XFER);

    fspace_id = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    mspace_id = H5Screate_simple(4, count, NULL);
    status = H5Dread(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);

#pragma omp parallel for
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
#else
#warning USE OF HDF WAS DISABLED AT COMPILE TIME
#endif
}
