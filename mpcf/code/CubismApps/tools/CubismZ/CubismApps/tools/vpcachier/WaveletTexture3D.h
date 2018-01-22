/*
 *  WaveletTexture3D.h
 *  
 *
 *  Created by Diego Rossinelli on 3/28/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <string>
#include <vector>
#include <sstream>

using namespace std;

#include <mpi.h>

#include "../../MPCFnode/source/WaveletCompressor.h"

#ifndef MYASSERT
//MACRO TAKEN FROM http://stackoverflow.com/questions/3767869/adding-message-to-assert
#   define MYASSERT(condition, message) \
do { \
if (! (condition)) { \
std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
<< " line " << __LINE__ << ": " << message << std::endl; \
std::exit(EXIT_FAILURE); \
} \
} while (false)
#endif

typedef WaveletCompressorGeneric_zlib<_VOXELS_, float> TextureCompressor;

#if 1
	float swapfloat(const float inFloat)
        {
		//if (!doswapping) return inFloat;

                float retVal;
                char *floatToConvert = (char *) &inFloat;
                char *returnFloat = (char *) &retVal;

                // swap the bytes into a temporary buffer
                returnFloat[0] = floatToConvert[3];
                returnFloat[1] = floatToConvert[2];
                returnFloat[2] = floatToConvert[1];
                returnFloat[3] = floatToConvert[0];

                return retVal;
        }

	int swapint(const int inInt)
        {
                //if (!doswapping) return inInt;

                int retVal;
                char *intToConvert = (char *) &inInt;
                char *returnInt = (char *) &retVal;

                // swap the bytes into a temporary buffer
                returnInt[0] = intToConvert[3];
                returnInt[1] = intToConvert[2];
                returnInt[2] = intToConvert[1];
                returnInt[3] = intToConvert[0];

                return retVal;
        }

        long swaplong(const long inLong)
        {
                //if (!doswapping) return inLong;

                long retVal;
                char *longToConvert = (char *) &inLong;
                char *returnLong = (char *) &retVal;

                if (sizeof(long) != 8) exit(1);

                // swap the bytes into a temporary buffer
                returnLong[0] = longToConvert[7];
                returnLong[1] = longToConvert[6];
                returnLong[2] = longToConvert[5];
                returnLong[3] = longToConvert[4];
                returnLong[4] = longToConvert[3];
                returnLong[5] = longToConvert[2];
                returnLong[6] = longToConvert[1];
                returnLong[7] = longToConvert[0];

                return retVal;
        }

        size_t swapsize(const size_t inSize)
        {
                //if (!doswapping) return inSize;

                size_t retVal;
                char *sizeToConvert = (char *) &inSize;
                char *returnSize = (char *) &retVal;

                if (sizeof(size_t) != 8) exit(1);

                // swap the bytes into a temporary buffer
                returnSize[0] = sizeToConvert[7];
                returnSize[1] = sizeToConvert[6];
                returnSize[2] = sizeToConvert[5];
                returnSize[3] = sizeToConvert[4];
                returnSize[4] = sizeToConvert[3];
                returnSize[5] = sizeToConvert[2];
                returnSize[6] = sizeToConvert[1];
                returnSize[7] = sizeToConvert[0];

                return retVal;
        }

#endif
struct WaveletTexture3D
{
	struct Geometry
	{
		float pos[3], size[3]; 
		float texcoordstart[3],texcoordend[3];
		
		template<int dim>
		void setup(const int gstart, const int gend, const int ghost1side, const double gridspacing)
		{			
			assert(gend - gstart == _VOXELS_);
			
			pos[dim] = (gstart + ghost1side) * gridspacing;
			size[dim] = (gend - gstart - 2 * ghost1side) * gridspacing;
			
			const double voxelsize = 1. / _VOXELS_;

			texcoordstart[dim] = ghost1side * voxelsize;
			texcoordend[dim] = (_VOXELS_ - ghost1side) * voxelsize; 
		}		
		
		void print()
		{
			for(int i = 0; i < 3; ++i)
				printf("%d: pos: %f size: %f, tex: %f - %f\n", i, pos[i], size[i], texcoordstart[i], texcoordend[i]);
		}
		
	} geometry;
	
	TextureCompressor wavcomp;

	WaveletsOnInterval::FwtAp (& data())[_VOXELS_][_VOXELS_][_VOXELS_] { return wavcomp.uncompressed_data(); }
	
	void compress(const float threshold, const bool halffloat, const unsigned char *& compresseddata, size_t& nbytes, int wtype)
	{
		nbytes = wavcomp.compress(threshold, halffloat, wtype);
		compresseddata = (unsigned char *) wavcomp.compressed_data();
	}	

	void compress(const float threshold, const bool halffloat, const unsigned char *& compresseddata, size_t& nbytes, bool swap, int wtype)
	{
		nbytes = wavcomp.compress(threshold, halffloat, swap, wtype);
		compresseddata = (unsigned char *) wavcomp.compressed_data();
	}	
	
	void * compression_buffer() { return wavcomp.compressed_data(); }
	
	void decompress(const bool halffloat, const size_t nbytes, int wtype) { wavcomp.decompress(halffloat, nbytes, wtype); }
};

class WaveletTexture3D_Collection
{
protected:
	
#ifdef __bgq__
	struct CompressedTexData { WaveletTexture3D::Geometry geometry; size_t start; int nbytes; } __attribute__((packed));
#else
#pragma pack (push)
#pragma pack (1)
	struct CompressedTexData { WaveletTexture3D::Geometry geometry; size_t start; int nbytes; };
#pragma pack (pop)
#endif		
	string path, header, data_title;
	size_t lutfile_start, bufferfile_start, buffertitle_start, buffer_size;
	
	int xtextures, ytextures, ztextures, ntextures;
	float wavelet_threshold;
	bool halffloat;
	bool swap;
	int wtype_read, wtype_write;
	
	FILE * myfile;
	
	//this data is used when reading the file
	vector<CompressedTexData> metadata;
	
	void swapGeometry(WaveletTexture3D::Geometry &g)
	{
		for (int i = 0; i < 3; i++) g.pos[i] = swapfloat(g.pos[i]);
		for (int i = 0; i < 3; i++) g.size[i] = swapfloat(g.size[i]);
		for (int i = 0; i < 3; i++) g.texcoordstart[i] = swapfloat(g.texcoordstart[i]);
		for (int i = 0; i < 3; i++) g.texcoordend[i] = swapfloat(g.texcoordend[i]);
	}

	void swapCTD(CompressedTexData &ctd)
	{
		swapGeometry(ctd.geometry);
		ctd.start = swapsize(ctd.start);
		ctd.nbytes = swapint(ctd.nbytes);
	}

	void _setup()
	{
		{
			std::stringstream ss;
			
			ss << "\n==============START-ASCI-HEADER==============\n";
			
			{
				int one = 1;
				bool isone = *(char *)(&one);
				
				if (!swap)
					ss << "Endianess: " << (isone ? "little" : "big") << "\n";
				else 
					ss << "Endianess: " << (!isone ? "little" : "big") << "\n";
			}
			
			ss << "Voxelsperdimension: " << _VOXELS_ << "\n";
			ss << "Textures: " << xtextures << " x "  << ytextures << " x " << ztextures  << "\n";
			ss << "HalfFloat: " << (this->halffloat ? "yes" : "no") << "\n";
			ss << "Wavelets: " << WaveletsOnInterval::ChosenWavelets_GetName(wtype_write) << "\n";
			ss << "Threshold: " << wavelet_threshold << "\n";
			ss << "Encoder: " << "zlib" << "\n";
			ss << "SizeofCompressedTexData: " << sizeof(CompressedTexData) << "\n";
			ss << "==============START-BINARY-LUT==============\n";
			
			header = ss.str();
		}
		
		data_title = "==============START-BINARY-DATA==============\n";
		
		lutfile_start = header.size();
		buffertitle_start = lutfile_start + ntextures * sizeof(CompressedTexData);
		bufferfile_start = buffertitle_start + data_title.size();
	}
	
public:
	
	WaveletTexture3D_Collection(const string path, int wtype_read, bool openfile = true): path(path), wtype_read(wtype_read), buffer_size(0), myfile(NULL)
	{
		_setup();
		
		if (openfile)
		{
			myfile = fopen(path.c_str(), "rb");
			assert(myfile);
			
			char tmp[1024];
			fgets(tmp, sizeof(tmp), myfile);
			fgets(tmp, sizeof(tmp), myfile);
			assert(string(tmp) == "==============START-ASCI-HEADER==============\n");
			
			fscanf(myfile, "Endianess:  %s\n", tmp);
			printf("Endianess: <%s>\n", tmp);
			//MYASSERT(string(tmp) == "little", "ooops Endianess is not little!\n");
			
			int voxels = -1;
			fscanf(myfile, "Voxelsperdimension:  %d\n", &voxels);
			printf("Voxelsperdimension: <%d>\n", voxels);
			MYASSERT(voxels == _VOXELS_, "Ooops voxels " << voxels << " instead of " << _VOXELS_ << "\n");

			fscanf(myfile, "Textures: %d x %d x %d\n", &xtextures, &ytextures, &ztextures);
			ntextures = xtextures * ytextures * ztextures;
			printf("Textures: %d x %d x %d\n", xtextures, ytextures, ztextures);

			fscanf(myfile, "HalfFloat: %s\n", tmp);
			printf("HalfFloat: <%s>\n", tmp);
			this->halffloat = (string(tmp) == "yes");
			
			fscanf(myfile, "Wavelets: %s\n", tmp);
			printf("Wavelets: <%s>\n", tmp);
			MYASSERT(tmp == string(WaveletsOnInterval::ChosenWavelets_GetName(wtype_read)),
					 "\nATTENZIONE:\nWavelets in the file is " << tmp << 
					 " and i have " << WaveletsOnInterval::ChosenWavelets_GetName(wtype_read) << "\n");
			
			fscanf(myfile, "Threshold: %e\n", &wavelet_threshold);
			printf("Threshold: <%e>\n", wavelet_threshold);
			
			fscanf(myfile, "Encoder: %s\n", tmp);
			printf("Encoder: <%s>\n", tmp);
			MYASSERT(tmp == string("zlib"),
					 "\nATTENZIONE:\nWavelets in the file is " << tmp << 
					 " and i have zlib.\n");
			
			int sizeofstruct = -1;	
			fscanf(myfile, "SizeofCompressedTexData:  %d\n", &sizeofstruct);
			printf("SizeofCompressedTexData: <%d>\n", sizeofstruct);
			MYASSERT(sizeof(CompressedTexData) == sizeofstruct,
					 "\nATTENZIONE:\nWavelets in the file is " << sizeofstruct << 
					 " and i have" << sizeof(CompressedTexData) << "\n");

			fgets(tmp, sizeof(tmp), myfile);
			
			assert(string("==============START-BINARY-LUT==============\n") == string(tmp));
			printf("==============END ASCI-HEADER==============\n\n");	
			
			metadata.resize(ntextures);
			fread(&metadata.front(), sizeof(CompressedTexData), ntextures, myfile);
			
			printf("geometry is:\n");
			for(int i=0; i<ntextures; ++i)
			{
				printf("%d: %zd %d\n", i, metadata[i].start, metadata[i].nbytes);
				metadata[i].geometry.print();
			}
		}
	}

	//void set_wtype_read(int wtype) { wtype_read = wtype; } 
	//void set_wtype_write(int wtype) { wtype_write = wtype; } 
	
	int get_xtextures() const { return xtextures; }
	int get_ytextures() const { return ytextures; }
	int get_ztextures() const { return ztextures; }
	int get_ntextures() const { return ntextures; }
	
	WaveletTexture3D_Collection(const string path, const int xtextures, const int ytextures, const int ztextures,
								const float wavelet_threshold, const bool halffloat, bool swap, int wtype_read, int wtype_write, bool openfile = true): 
	path(path), xtextures(xtextures), ytextures(ytextures), ztextures(ztextures), ntextures(xtextures * ytextures * ztextures),
	wavelet_threshold(wavelet_threshold), halffloat(halffloat), swap(swap), wtype_read(wtype_read), wtype_write(wtype_write), buffer_size(0), myfile(NULL)
	{				
		_setup();
		
		if (openfile)
		{
			myfile = fopen(path.c_str(), "wb");
			assert(myfile);
		
			fwrite(header.c_str(), sizeof(char), header.size(), myfile);	// text
			fseek(myfile, buffertitle_start, SEEK_SET);
			fwrite(data_title.c_str(), sizeof(char), data_title.size(), myfile);	// text				
		}
	}
	
	virtual ~WaveletTexture3D_Collection() 
	{				
		if (myfile)
		{
			fclose(myfile);
			
			myfile = NULL;
			
			printf("Terminating...closing file (%.2f kB) now.\n", (bufferfile_start + buffer_size)/1024.);
		}
	}
	
	virtual void write(const int ix, const int iy, const int iz, WaveletTexture3D& texture)
	{
		//check that we are not totally nuts
		assert(ix >= 0 && ix < xtextures);
		assert(iy >= 0 && iy < ytextures);
		assert(iz >= 0 && iz < ztextures);
		
		//compress the texture
		const unsigned char * ptr = NULL;
		size_t nbytes = 0;
		texture.compress(wavelet_threshold, halffloat, ptr, nbytes, wtype_write);
		
		//spit some output for now
		{
			const size_t uncompressedbytes = sizeof(Real) * _VOXELS_ * _VOXELS_ * _VOXELS_;
			
			printf("(w) Texture: %d %d %d, CR: %.1fX\n", 
				   ix, iy, iz, uncompressedbytes * 1. / nbytes, wavelet_threshold);	
		}
		
		//allocate a region in the file
		const size_t myoffset = bufferfile_start + buffer_size;
		buffer_size += nbytes;
		
		//write lut
		CompressedTexData entry = { texture.geometry, myoffset, nbytes };
		const size_t mylutoffset = sizeof(entry) * (ix + xtextures * (iy + ytextures * iz));
		fseek(myfile, lutfile_start + mylutoffset, SEEK_SET);
		if (swap) swapCTD(entry);
		fwrite(&entry, sizeof(entry), 1, myfile);
		
		//write data
		assert(ptr != NULL);
		assert(nbytes != 0);
		fseek(myfile, myoffset, SEEK_SET);
		fwrite(ptr, sizeof(unsigned char), nbytes, myfile);	// compressed buffer
	}
	
	virtual void read(const int index, WaveletTexture3D& texture, bool onlygeometry = false) const
	{
		//check that we are not totally nuts
		assert(index >= 0 && index < metadata.size());
	
		//recover the geometry
		texture.geometry = metadata[index].geometry;
		
		if (onlygeometry) return;
		
		const size_t nbytes = metadata[index].nbytes;		
		fseek(myfile, metadata[index].start, SEEK_SET);
		fread(texture.compression_buffer(), sizeof(unsigned char), nbytes, myfile);
		
		//decompress the texture
		texture.decompress(halffloat, nbytes, wtype_read);
		
		
		//spit some output
		{
			const size_t uncompressedbytes = sizeof(float) * _VOXELS_ * _VOXELS_ * _VOXELS_;
			
			printf("(r) Texture-id %d. CR: %.1fX\n", index, uncompressedbytes * 1. / nbytes);	
		}
	}
	
	virtual void read(const int ix, const int iy, const int iz, WaveletTexture3D& texture, bool onlygeometry = false) const
	{
		//check that we are not totally nuts
		assert(ix >= 0 && ix < xtextures);
		assert(iy >= 0 && iy < ytextures);
		assert(iz >= 0 && iz < ztextures);		
		
		const int myentry = ix + xtextures * (iy + ytextures * iz);
		assert(myentry >= 0 && myentry < metadata.size());

		read(myentry, texture, onlygeometry);
	}
};

class WaveletTexture3D_CollectionMPI: public WaveletTexture3D_Collection
{	
	size_t * file_offset;
	
	MPI::Win rmawindow;
	MPI::File myfile;
	const MPI::Intracomm& mycomm;

public:
	
	WaveletTexture3D_CollectionMPI(const MPI::Intracomm& comm, 
								   const string path, const int xtextures, const int ytextures, const int ztextures,
								   const float wavelet_threshold, const bool halffloat, bool swap, int wtype_read, int wtype_write): 
	mycomm(comm), WaveletTexture3D_Collection(path, xtextures, ytextures, ztextures, wavelet_threshold, halffloat, swap, wtype_read, wtype_write, false), file_offset(NULL)
	{		
		const int mygid = comm.Get_rank();

		//file setup
		{			
			myfile = MPI::File::Open(mycomm, path.c_str(), MPI::MODE_CREATE | MPI::MODE_WRONLY, MPI::INFO_NULL);
						
			if (mygid == 0)
			{
				myfile.Write_at(0, header.c_str(), header.size(), MPI::CHAR);
				myfile.Write_at(buffertitle_start, data_title.c_str(), data_title.size(), MPI::CHAR);				
			}			
		}
		
		//one sided communication setup
		{
			file_offset = (size_t *) MPI::Alloc_mem(sizeof(size_t), MPI::INFO_NULL);
			*file_offset = bufferfile_start; 
			assert(*file_offset != 0);
			
			//if (mygid == 0)
			//	printf("at the beginning my offset was %zd\n", *file_offset);
			
			rmawindow = MPI::Win::Create(file_offset, sizeof(size_t), sizeof(size_t), MPI::INFO_NULL, mycomm);
		}
	}
	
	~WaveletTexture3D_CollectionMPI() 
	{		
		rmawindow.Free();
		
		if (mycomm.Get_rank() == 0)
			printf("Terminating...closing file (%.2f kB) now.\n", *file_offset/1024.);
		
		MPI::Free_mem(file_offset);
		
		myfile.Close(); 		
	}
		
	void write(const int ix, const int iy, const int iz, WaveletTexture3D& texture)
	{
		//check that we are not totally nuts
		assert(ix >= 0 && ix < xtextures);
		assert(iy >= 0 && iy < ytextures);
		assert(iz >= 0 && iz < ztextures);

		//compress the texture
		const unsigned char * ptr = NULL;
		size_t nbytes = 0;
		texture.compress(wavelet_threshold, halffloat, ptr, nbytes, wtype_write);

		//spit some output for now
		{
			const size_t uncompressedbytes = sizeof(float) * _VOXELS_ * _VOXELS_ * _VOXELS_;
			
			printf("Texture: %d %d %d, CR: %.1fX\n", ix, iy, iz, uncompressedbytes * 1. / nbytes);	
		}

		//obtain file offset
		size_t myoffset = 0;
		{
			rmawindow.Lock(MPI::LOCK_EXCLUSIVE, 0, 0);
			rmawindow.Get(&myoffset, 1, MPI_UINT64_T, 0, 0, 1, MPI_UINT64_T); 
			rmawindow.Accumulate(&nbytes, 1, MPI_UINT64_T, 0, 0, 1, MPI_UINT64_T, MPI::SUM);
			rmawindow.Unlock(0);
		}
		assert(myoffset != 0);

		//write lut
		CompressedTexData entry = { texture.geometry, myoffset, nbytes };
		
		const size_t mylutoffset = sizeof(entry) * (ix + xtextures * (iy + ytextures * iz));
		if (swap) swapCTD(entry);
		myfile.Write_at(lutfile_start + mylutoffset, &entry, sizeof(entry), MPI::CHAR);

		//write data
		assert(ptr != NULL);
		assert(nbytes != 0);
		myfile.Write_at(myoffset, ptr, nbytes, MPI::CHAR);	// compressed buffer
	}
};
