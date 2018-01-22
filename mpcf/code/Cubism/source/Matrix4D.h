/*
 *  Matrix4D.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 12/5/06.
 *  Copyright 2006 ETH Zurich. All rights reserved. 0
 *
 */
#pragma once
#include <assert.h>
#include <stdio.h>
#include "Matrix3D.h"

template <class DataType, bool bPrimitiveType, template <typename T> class allocator>  
class Matrix4D
{
private:
	DataType * m_pData;
	unsigned int m_vSize[4];
	unsigned int m_nElements;
	unsigned int m_nElementsPerSlice;
	unsigned int m_nElementsPerCube;
	
	void _Setup(unsigned int nSizeX, unsigned int nSizeY, unsigned int nSizeZ, unsigned int nSizeW);  
	
public:
	Matrix4D(FILE * f, bool bFileIsBigEndian=false);
	Matrix4D(FILE * f, bool bFileIsBigEndian=false, bool bFrom64bits=false, int bytes_padding=0);
	Matrix4D(unsigned int nSizeX, unsigned int nSizeY, unsigned int nSizeZ, unsigned int nSizeW);
	~Matrix4D() { delete [] m_pData; }
	
	DataType & Access(unsigned int ix, unsigned int iy, unsigned int iz, unsigned int iw) const;
	DataType & LinAccess(unsigned int i) const;
	DataType * getData() const { return m_pData; } 
	
	unsigned int * getSize();
	unsigned int getNumberOfElements() const;
	
	void operator= (DataType val);
	
	void Serialize(FILE * f);
	void Deserialize(FILE *f, bool bSwapBytes=false);
	void Deserialize(FILE *f, bool bSwapBytes=false, bool bFrom64bites=false, int bytes_padding=0);
};

template <class DataType, bool bPrimitiveType, template <typename T> class allocator>  void Matrix4D<DataType, bPrimitiveType, allocator>::_Setup(unsigned int nSizeX, unsigned int nSizeY, unsigned int nSizeZ, unsigned int nSizeW)
{
	m_vSize[0] = nSizeX;
	m_vSize[1] = nSizeY;
	m_vSize[2] = nSizeZ;
	m_vSize[3] = nSizeW;
	
	m_nElementsPerSlice = nSizeX*nSizeY;
	m_nElementsPerCube = m_nElementsPerSlice*nSizeZ;
	
	m_nElements = nSizeX*nSizeY*nSizeZ*nSizeW;
	
	if (m_pData != NULL) delete m_pData;
	m_pData = NULL;
	m_pData = new DataType[m_nElements];
	assert(m_pData !=NULL);
}

template <class DataType, bool bPrimitiveType, template <typename T> class allocator>  Matrix4D<DataType, bPrimitiveType, allocator>::Matrix4D(unsigned int nSizeX, unsigned int nSizeY, unsigned int nSizeZ, unsigned int nSizeW):
m_pData(NULL),
m_nElements(0),
m_nElementsPerSlice(0),
m_nElementsPerCube(0)
{
	_Setup(nSizeX,nSizeY,nSizeZ,nSizeW);
}


template <class DataType, bool bPrimitiveType, template <typename T> class allocator>  Matrix4D<DataType, bPrimitiveType, allocator>::Matrix4D(FILE * f, bool bSwapBytes):
m_pData(NULL),
m_nElements(0),
m_nElementsPerSlice(0),
m_nElementsPerCube(0)
{
	Deserialize(f,bSwapBytes);
}

template <class DataType, bool bPrimitiveType, template <typename T> class allocator>  Matrix4D<DataType, bPrimitiveType, allocator>::Matrix4D(FILE * f, bool bSwapBytes, bool bFrom64bits, int bytes_padding):
m_pData(NULL),
m_nElements(0),
m_nElementsPerSlice(0),
m_nElementsPerCube(0)
{
	Deserialize(f,bSwapBytes,bFrom64bits,bytes_padding);
}


template <class DataType, bool bPrimitiveType, template <typename T> class allocator> inline DataType & Matrix4D<DataType, bPrimitiveType, allocator>::Access(unsigned int ix, unsigned int iy, unsigned int iz, unsigned int iw) const
{
#ifdef _DEBUG
	if (!(ix<m_vSize[0]) || !(iy<m_vSize[1]) || !(iz<m_vSize[2]) || !(iw<m_vSize[3])) printf("excpetion %d %d %d %d\n", ix,iy,iz, iw);
	assert(ix<m_vSize[0]);
	assert(iy<m_vSize[1]);
	assert(iz<m_vSize[2]);
	assert(iw<m_vSize[3]);
#endif
	
	return m_pData[iw*m_nElementsPerCube + iz*m_nElementsPerSlice + iy*m_vSize[0] + ix];
}

template <class DataType, bool bPrimitiveType, template <typename T> class allocator> inline DataType & Matrix4D<DataType, bPrimitiveType, allocator>::LinAccess(unsigned int i) const
{
#ifdef _DEBUG
	assert(i<m_nElements);
#endif
	return m_pData[i];
}

template <class DataType, bool bPrimitiveType, template <typename T> class allocator> inline unsigned int  Matrix4D<DataType, bPrimitiveType, allocator>::getNumberOfElements() const
{
	return m_nElements;
}

template <class DataType, bool bPrimitiveType, template <typename T> class allocator> inline unsigned int *  Matrix4D<DataType, bPrimitiveType, allocator>::getSize() 
{
	return m_vSize;
}

template <class DataType, bool bPrimitiveType, template <typename T> class allocator> inline void  Matrix4D<DataType, bPrimitiveType, allocator>::operator= (DataType val)
{
	for(unsigned int i=0;i<m_nElements;i++) m_pData[i] = val;
}


template <class DataType, bool bPrimitiveType, template <typename T> class allocator> void  Matrix4D<DataType, bPrimitiveType, allocator>::Serialize(FILE * f)
{
	fwrite((void*) this, sizeof(Matrix4D<DataType, bPrimitiveType, allocator>), 1, f);
	fwrite((void*) m_pData, sizeof(DataType), m_nElements, f);
}

template <class DataType, bool bPrimitiveType, template <typename T> class allocator> void  Matrix4D<DataType, bPrimitiveType, allocator>::Deserialize(FILE *f, bool bSwapBytes)
{
	{
		const unsigned int nIntMembers = 1+4+1+1+1; 
		unsigned int buf[nIntMembers];
		fread((void *)buf, 4,nIntMembers, f);
		m_pData = NULL;
		m_vSize[0] = buf[1];
		m_vSize[1] = buf[2];
		m_vSize[2] = buf[3];
		m_vSize[3] = buf[4];
		m_nElements = buf[5];
		m_nElementsPerSlice = buf[6];
		m_nElementsPerCube = buf[7]; 
	}
	
	m_pData = NULL;
	_Setup(m_vSize[0], m_vSize[1], m_vSize[2], m_vSize[3]);
	
	if (bSwapBytes)
	{
		unsigned char * buf = new unsigned char [sizeof(DataType)* m_nElements];
		fread((void*)buf, sizeof(DataType), m_nElements, f);
		delete [] buf;
	}
	else
		fread((void*) m_pData, sizeof(DataType), m_nElements, f);
}

template <class DataType, bool bPrimitiveType, template <typename T> class allocator> void  Matrix4D<DataType, bPrimitiveType, allocator>::Deserialize(FILE *f, bool bSwapBytes, bool bFrom64bits, int bytes_padding)
{
	int mysize  =-1;
	if (bFrom64bits || sizeof(int *) != 8)
		fread(&mysize, sizeof(int), 1, f);
	
	fread((void*)this, sizeof(Matrix4D<DataType, bPrimitiveType, allocator>), 1, f);
	if (bytes_padding>0) 
	{
		char pad[bytes_padding];
		fread(&pad, sizeof(char), bytes_padding, f);
	}
	
	m_pData = NULL;
	_Setup(m_vSize[0], m_vSize[1], m_vSize[2], m_vSize[3]);
	
	if (bSwapBytes)
	{
		unsigned char * buf = new unsigned char [sizeof(DataType)* m_nElements];
		fread((void*)buf, sizeof(DataType), m_nElements, f);
		SwapBytes<DataType>(buf, m_nElements*sizeof(DataType));
		memcpy((void*) m_pData, buf, sizeof(DataType)*m_nElements);
		delete [] buf;
	}
	else
		fread((void*) m_pData, sizeof(DataType), m_nElements, f);	
}

