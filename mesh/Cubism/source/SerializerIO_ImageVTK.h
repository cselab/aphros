/*
 *  SerializerIO_ImageVTK.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/28/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <vtkPoints.h> 
#include <vtkCell.h>
#include <vtkImageData.h>
#include <vtkImageNoiseSource.h>
#include <vtkFloatArray.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include "SerializerIO.h"

template<typename GridType, typename Streamer>
class SerializerIO_ImageVTK
{
	// Typedefs
	typedef typename GridType::BlockType TBlock;
	typedef typename TBlock::ElementType TElement;
	
public:
	
	void Write(GridType & inputGrid, string fileName, Streamer streamer = Streamer())
	{
		static const int BX = TBlock::sizeX;
		static const int BY = TBlock::sizeY;
		static const int BZ = TBlock::sizeZ;
		
		const int NBX = inputGrid.getBlocksPerDimension(0);
		const int NBY = inputGrid.getBlocksPerDimension(1);
		const int NBZ = inputGrid.getBlocksPerDimension(2);
		
		const int NX = BX*NBX;
		const int NY = BY*NBY;
		const int NZ = BZ*NBZ;
		const int NC = streamer.channels;
		
		vtkImageData * img = vtkImageData::New();
		
		img->SetExtent(0,NX-1,0,NY-1,0,NZ-1);
		img->SetDimensions(NX, NY, NZ);
		img->SetScalarTypeToFloat();
		img->SetNumberOfScalarComponents(NC);
		
		img->AllocateScalars();
		img->SetSpacing(1./NX, 1./NX, 1./NX);
		img->SetOrigin(0,0,0);
		
		for(int ibz = 0; ibz<NBZ; ibz++)
			for(int iby = 0; iby<NBY; iby++)
				for(int ibx = 0; ibx<NBX; ibx++)
				{
					const TBlock& block = inputGrid(ibx, iby, ibz);
					
					for(int iz=0; iz<BZ; iz++)
						for(int iy=0; iy<BY; iy++)
							for(int ix=0; ix<BX; ix++)
							{
								Real output[NC];
								
								streamer.operate(block.data[iz][iy][ix], output);
								
								const int gx = ix + ibx*BX;
								const int gy = iy + iby*BY;
								const int gz = iz + ibz*BZ;
								
								for(int ic=0; ic<NC; ic++)
									img->SetScalarComponentFromFloat(gx, gy, gz, ic, output[ic]);
							}
				}
		
		vtkXMLImageDataWriter * writer = vtkXMLImageDataWriter::New();
		writer->SetFileName(fileName.c_str());
		writer->SetInput(img);
		writer->Write();
		
		writer->Delete();
		img->Delete();
	}
	
	template<typename TLab>
	void WriteLabs(GridType & inputGrid, string fileName, const Real time=0, Streamer streamer = Streamer())
	{
		const vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		
		static const int BX = TBlock::sizeX;
		static const int BY = TBlock::sizeY;
		static const int BZ = TBlock::sizeZ;
		
		const int NBX = inputGrid.getBlocksPerDimension(0);
		const int NBY = inputGrid.getBlocksPerDimension(1);
		const int NBZ = inputGrid.getBlocksPerDimension(2);
		
		const int NX = BX*NBX;
		const int NY = BY*NBY;
		const int NZ = BZ*NBZ;
		const int NC = streamer.channels;
		
		vtkImageData * img = vtkImageData::New();
		
		img->SetExtent(0,NX+5,0,NY+5,0,NZ-1);
		img->SetDimensions(NX+6, NY+6, NZ);
		img->SetScalarTypeToFloat();
		img->SetNumberOfScalarComponents(NC);
		
		img->AllocateScalars();
		img->SetSpacing(1./NX, 1./NX, 1./NX);
		img->SetOrigin(0,0,0);
		
		TLab lab;
		int steStart[3] ={-3,-3,0};
		int steEnd[3] = {4,4,1};
		lab.prepare(inputGrid, steStart, steEnd);
		
		for(int i=0; i<vInfo.size(); i++)
		{
			const BlockInfo& info = vInfo[i];
			lab.load(info,time);
			
			for(int iz=0; iz<BZ; iz++)
				for(int iy=(info.index[1]==0? steStart[1]:0); iy<(info.index[1]==inputGrid.getBlocksPerDimension(1)-1? BY+steEnd[1]-1:BY); iy++)
					for(int ix=(info.index[0]==0? steStart[0]:0); ix<(info.index[0]==inputGrid.getBlocksPerDimension(0)-1? BX+steEnd[0]-1:BX); ix++)
					{
						Real output[NC];
						
						streamer.operate(lab(ix,iy,iz), output);
						
						const int gx = ix-steStart[0] + info.index[0]*BX;
						const int gy = iy-steStart[1] + info.index[1]*BY;
						const int gz = iz + info.index[2]*BZ;
						
						for(int ic=0; ic<NC; ic++)
							img->SetScalarComponentFromFloat(gx, gy, gz, ic, output[ic]);
					}
			
		}
		
		vtkXMLImageDataWriter * writer = vtkXMLImageDataWriter::New();
		writer->SetFileName(fileName.c_str());
		writer->SetInput(img);
		writer->Write();
		
		writer->Delete();
		img->Delete();
	}
};
