/*
 *  Reader_WaveletCompression.h
 *
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Extended by Panos Hadjidoukas.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */
#ifndef READER_WAVELETCOMPRESSION_H_00000000
#define READER_WAVELETCOMPRESSION_H_00000000

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <string>
#include <vector>
#include <iostream>

#include <mpi.h>

using namespace std;

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#include "../../Compressor/source/WaveletCompressor.h"

#include "../../Compressor/source/WaveletSerializationTypes.h"
#include "../../Compressor/source/CompressionEncoders0.h"
#include "../../Compressor/source/FullWaveletTransform.h"

//MACRO TAKEN FROM http://stackoverflow.com/questions/3767869/adding-message-to-assert
#   define MYASSERT(condition, message) \
do { \
if (! (condition)) { \
std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
<< " line " << __LINE__ << ": " << message << std::endl; \
std::exit(EXIT_FAILURE); \
} \
} while (false)

/* peh: support of different endianness */

/* data structures where bytes swapping is applied
struct BlockMetadata { int idcompression, subid, ix, iy, iz; }  __attribute__((packed));
struct HeaderLUT { size_t aggregate_bytes; int nchunks; }  __attribute__((packed));
struct CompressedBlock{ size_t start, extent; int subid; }  __attribute__((packed));
*/

class Reader_WaveletCompression0
{
protected:
	string path;
	bool doswapping; // peh: support of different endianness
	int wtype;	// peh:

	size_t global_header_displacement;
	int miniheader_bytes;
	int NBLOCKS;
	int totalbpd[3], bpd[3];
	bool halffloat;

	vector<CompressedBlock> idx2chunk;

	int _id(int ix, int iy, int iz) const
	{
		assert(ix >= 0 && ix < totalbpd[0]);
		assert(iy >= 0 && iy < totalbpd[1]);
		assert(iz >= 0 && iz < totalbpd[2]);

		return ix + totalbpd[0] * ( iy + totalbpd[1] * iz );
	}

	// peh: BGQ <-> x86_64
	void swapbytes(unsigned char *mem, int nbytes)
	{
		if (!doswapping) return;

		unsigned char buf[8];
		for (int i = 0; i < nbytes; i++) buf[i] = mem[i];
		for (int i = 0; i < nbytes; i++) mem[nbytes-i-1] = buf[i];
	}

	float swapfloat(const float inFloat)
	{
		if (!doswapping) return inFloat;

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
		if (!doswapping) return inInt;

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
		if (!doswapping) return inLong;

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

	void swapBM(struct BlockMetadata &bm)
	{
		if (!doswapping) return;

		bm.idcompression = swapint(bm.idcompression);
		bm.subid = swapint(bm.subid);
		bm.ix = swapint(bm.ix);
		bm.iy = swapint(bm.iy);
		bm.iz = swapint(bm.iz);
	}

	void swapHL(struct HeaderLUT &hl)
	{
		if (!doswapping) return;

		hl.aggregate_bytes = swaplong(hl.aggregate_bytes);
		hl.nchunks = swapint(hl.nchunks);
	}

	void swapCB(struct CompressedBlock &cb)
	{
		if (!doswapping) return;

		cb.start = swaplong(cb.start);
		cb.extent = swaplong(cb.extent);
		cb.subid = swapint(cb.subid);
	}

	int myendianness()
	{
		int one = 1;
		bool isone = *(char *)(&one);

		return isone ? 0: 1;
	}

public:

	Reader_WaveletCompression0(const string path, bool doswapping, int wtype): NBLOCKS(-1), global_header_displacement(-1), path(path), doswapping(doswapping), wtype(wtype) { }

	virtual void load_file()
	{
		for(int i = 0; i < 3; ++i)
			totalbpd[i] = -1;

		for(int i = 0; i < 3; ++i)
			bpd[i] = -1;

		string binaryocean_title = "\n==============START-BINARY-OCEAN==============\n";

		this->miniheader_bytes = sizeof(size_t) + binaryocean_title.size();

		vector<BlockMetadata> metablocks;
		vector<size_t> lutchunks;
		vector<int> sizechunks;

		{
			FILE * file = fopen(path.c_str(), "rb");

			MYASSERT(file, "\nAAATTENZIONE:\nOooops could not open the file. Path: " << path);

			//reading the header and mini header
			{
				size_t header_displacement = -1;
				fread(&header_displacement, sizeof(size_t), 1, file);
				header_displacement = swaplong(header_displacement);

				fseek(file, header_displacement, SEEK_SET);
				global_header_displacement = header_displacement;

				char buf[1024];
				fgets(buf, sizeof(buf), file);
				fgets(buf, sizeof(buf), file);

				printf("\n%s", buf);
				assert(string("==============START-ASCI-HEADER==============\n") == string(buf));

				fscanf(file, "Endianess:  %s\n", buf);
				printf("Endianess: <%s>\n", buf);
//				assert(string(buf) == "little");
//				int fileendianness = (string(buf) == "little");
//				printf("fileendianness = %d\n", fileendianness);
//				printf("myendianness() = %d\n", myendianness());
//				doswapping = (myendianness() != fileendianness);
//				printf("Swapping is <%d>\n", doswapping);

				int sizeofreal = -1;
				fscanf(file, "sizeofReal:  %d\n", &sizeofreal);
				printf("sizeofReal: <%d>\n", sizeofreal);

				MYASSERT(sizeof(Real) == sizeofreal,
						 "\nATTENZIONE:\nSizeof(Real) in the file is " << sizeofreal << " which is wrong\n");

				int sizeofsize_t = -1;
				fscanf(file, "sizeofsize_t:  %d\n", &sizeofsize_t);
				printf("sizeofsize_t: <%d>\n", sizeofsize_t);

				MYASSERT(sizeof(size_t) == sizeofsize_t,
						 "\nATTENZIONE:\nSizeof(size_t) in the file is " << sizeofsize_t << " which is wrong\n");

				int sizeofblockmetadata = -1;
				fscanf(file, "sizeofBlockMetadata:  %d\n", &sizeofblockmetadata);
				printf("sizeofBlockMetadata: <%d>\n", sizeofblockmetadata);

				MYASSERT(sizeof(BlockMetadata) == sizeofblockmetadata,
						 "\nATTENZIONE:\nSizeof(sizeofblockmetadata) in the file is " << sizeofblockmetadata << " which is wrong\n");

				int sizeofheaderlut = -1;
				fscanf(file, "sizeofHeaderLUT:  %d\n", &sizeofheaderlut);
				printf("sizeofHeaderLUT: <%d>\n", sizeofheaderlut);

				MYASSERT(sizeof(HeaderLUT) == sizeofheaderlut,
						 "\nATTENZIONE:\nSizeof(HeaderLUT) in the file is " << sizeofheaderlut << " which is wrong\n");

				int sizeofcompressedblock = -1;
				fscanf(file, "sizeofCompressedBlock:  %d\n", &sizeofcompressedblock);
				printf("sizeofCompressedBlock: <%d>\n", sizeofcompressedblock);

				MYASSERT(sizeof(CompressedBlock) == sizeofcompressedblock,
						 "\nATTENZIONE:\nSizeof(sizeofCompressedBlock) in the file is " << sizeofcompressedblock << " which is wrong\n");

				int bsize = -1;
				fscanf(file, "Blocksize: %d\n", &bsize);
				printf("Blocksize: <%d>\n", bsize);

				MYASSERT(bsize == _BLOCKSIZE_,
						 "\nATTENZIONE:\nBlocksize in the file is " << bsize <<
						 " and i have " << _BLOCKSIZE_ << "\n");

				fscanf(file, "Blocks: %d x %d x %d\n", totalbpd, totalbpd + 1, totalbpd + 2);
				printf("Blocks: %d x %d x %d\n", totalbpd[0], totalbpd[1], totalbpd[2]);

				float myxextent = -1, myyextent = -1, myzextent = -1;
				fscanf(file, "Extent: %f %f %f\n", &myxextent,  &myyextent, &myzextent);
				printf("Extent: <%f> x <%f> x <%f>\n", myxextent,  myyextent, myzextent);

				fscanf(file, "SubdomainBlocks: %d x %d x %d\n", bpd, bpd + 1, bpd + 2);
				printf("SubdomainBlocks: <%d x %d x %d>\n", bpd[0], bpd[1], bpd[2]);

				fscanf(file, "HalfFloat: %s\n", buf);
				printf("HalfFloat: <%s>\n", buf);
				this->halffloat = (string(buf) == "yes");

				fscanf(file, "Wavelets: %s\n", buf);
				printf("Wavelets: <%s>\n", buf);





















				MYASSERT(buf == string("none"),
						"\nATTENZIONE:\nWavelets in the file is " << buf <<
						" and i have " << "none"  << "\n");


				float mythreshold = -1;
				fscanf(file, "WaveletThreshold: %f\n", &mythreshold);
				printf("WaveletThreshold: <%f>\n", mythreshold);

				fscanf(file, "Encoder: %s\n", buf);
				printf("Encoder: <%s>\n", buf);









				MYASSERT(buf == string("none"),
						 "\nATTENZIONE:\nEncoder in the file is " << buf <<
						 " and i have none.\n");


				fgets(buf, sizeof(buf), file);

				assert(string("==============START-BINARY-METABLOCKS==============\n") == string(buf));
				printf("==============END ASCI-HEADER==============\n\n");
				NBLOCKS = totalbpd[0] * totalbpd[1] * totalbpd[2];

				//printf("Blocks: %d -> %dx%dx%d -> subdomains of %dx%dx%d\n",
				//	NBLOCKS, totalbpd[0], totalbpd[1], totalbpd[2], bpd[0], bpd[1], bpd[2]);
			}

			//reading the binary lut
			{
				metablocks.resize(NBLOCKS);

				for(int i = 0; i < NBLOCKS; ++i)
				{
					BlockMetadata entry;
					fread(&entry, sizeof(entry), 1, file);
					swapBM(entry);
					//printf("reading metablock %d -> %d %d %d  cid %d\n", i, entry.ix, entry.iy, entry.iz, entry.idcompression);
					assert(entry.idcompression >= 0 && entry.idcompression < bpd[0] * bpd[1] * bpd[2]);
					metablocks[i] = entry;
				}
			}

			//reading the lut header
			{
				char buf[1024];

				fgetc(file);
				fgets(buf, sizeof(buf), file);

				assert(string("==============START-BINARY-LUT==============\n") == string(buf));

				bool done = false;

				size_t base = miniheader_bytes;

				const int BPS = bpd[0] * bpd[1] * bpd[2];
				assert(NBLOCKS % BPS == 0);
				const int SUBDOMAINS = NBLOCKS / BPS;

				vector<HeaderLUT> headerluts(SUBDOMAINS); //oh mamma mia
				fread(&headerluts.front(), sizeof(HeaderLUT), SUBDOMAINS, file);
				{
				HeaderLUT *hl = headerluts.data();
				for (int h = 0; h < SUBDOMAINS; h++) swapHL(hl[h]);
				}

				{
					int c = fgetc(file);

					do
					{
						printf("shouldnt be here! 0x%x\n", c);
						//abort();
						c = fgetc(file);
					}
					while (! feof(file) );
				}
				//assert(feof(file));

				for(int s = 0, currblock = 0; s < SUBDOMAINS; ++s)
				{
					const int nglobalchunks = lutchunks.size();

					const int nchunks = headerluts[s].nchunks;
					const size_t myamount = headerluts[s].aggregate_bytes;
					const size_t lutstart = base + myamount - sizeof(size_t) * nchunks;

					fseek(file, lutstart, SEEK_SET);
					vector<size_t> mylut(nchunks);
					fread(&mylut.front(), sizeof(size_t), nchunks, file);
					{
					size_t *ml = mylut.data();
					for (int n = 0; n < nchunks; n++) ml[n] = swaplong(ml[n]);
					}

					for(int i=0; i< nchunks; ++i)
						assert(mylut[i] < myamount);

					for(int i=1; i< nchunks; ++i)
						assert(mylut[i-1] < mylut[i]);

					//compute the chunk sizes
					{
						for(int i = 0; i < mylut.size()-1; ++i)
						{
							const int mysize = mylut[i+1] - mylut[i];
							assert(mysize > 0);
							sizechunks.push_back(mysize);
						}

						const size_t mysize = (myamount - sizeof(size_t) * nchunks) - mylut[mylut.size() - 1];

						assert(mysize > 0);
						sizechunks.push_back(mysize);
					}

					for(int i = 0; i < mylut.size(); ++i)
					{
						assert(mylut[i] < myamount);
						mylut[i] += base;
					}

					assert(myamount > 0);
					base += myamount;
					assert(base <= global_header_displacement);

					//compute the base for this blocks
					for(int i = 0; i < BPS; ++i, ++currblock)
						metablocks[currblock].idcompression += nglobalchunks;

					lutchunks.insert(lutchunks.end(), mylut.begin(), mylut.end());
				}

				assert(base == global_header_displacement);

				lutchunks.push_back(base);
			}

			fclose(file);
		}

		idx2chunk.resize(NBLOCKS);

		for(int i = 0; i < NBLOCKS ; ++i)
		{
			BlockMetadata entry = metablocks[i];

			assert(entry.idcompression >= 0);
			assert(entry.idcompression < lutchunks.size()-1);

			size_t start_address = lutchunks[entry.idcompression];
			assert(sizechunks.size() > entry.idcompression);
			size_t end_address = start_address + sizechunks[entry.idcompression] ;
			assert(sizechunks[entry.idcompression] > 0);
			assert (end_address > start_address);
			assert (end_address <= lutchunks[entry.idcompression + 1]);

			assert(start_address < end_address);
			assert(end_address <= global_header_displacement);
			assert( start_address < global_header_displacement );

			CompressedBlock compressedblock = { start_address, end_address - start_address, entry.subid };

			idx2chunk[_id(entry.ix, entry.iy, entry.iz)] = compressedblock;
		}

		const bool verbose = true;

		if (verbose)
		{
			const size_t size_idx2chunk = idx2chunk.size() * sizeof(CompressedBlock);
			const double footprint_mb =  size_idx2chunk / 1024. / 1024.;
			printf("the header data is taking %.2f MB\n", footprint_mb);
		}
	}

	int xblocks() { return totalbpd[0]; }
	int yblocks() { return totalbpd[1]; }
	int zblocks() { return totalbpd[2]; }

	void load_block(int ix, int iy, int iz, Real MYBLOCK[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_])
	{
		//printf("trying to load <%s>\n", path.c_str());
		FILE * f = fopen(path.c_str(), "rb");

		assert(f);

		CompressedBlock compressedchunk = idx2chunk[_id(ix, iy, iz)];

		size_t start = compressedchunk.start;

		assert(start >= miniheader_bytes);
		assert(start < global_header_displacement);
		assert(start + compressedchunk.extent <= global_header_displacement);

		vector<unsigned char> compressedbuf(compressedchunk.extent);
		fseek(f, compressedchunk.start, SEEK_SET);
		fread(&compressedbuf.front(), compressedchunk.extent, 1, f);

		assert(!feof(f));

		//static vector<unsigned char> waveletbuf(2 << 21);	// 21: 4MB, 22: 8MB, 28: 512MB
		static vector<unsigned char> waveletbuf(2 << 28);	// 21: 4MB, 22: 8MB, 28: 512MB
		const size_t decompressedbytes = zdecompress0(&compressedbuf.front(), compressedbuf.size(), &waveletbuf.front(), waveletbuf.size());

		int readbytes = 0;
		for(int i = 0; i<compressedchunk.subid; ++i)
		{
			int nbytes = * (int *) & waveletbuf[readbytes];
			nbytes = swapint(nbytes);
			readbytes += sizeof(int);
			readbytes += nbytes;

			assert(readbytes <= decompressedbytes);
		}

		{
			int nbytes = *(int *)&waveletbuf[readbytes];
			nbytes = swapint(nbytes);
			readbytes += sizeof(int);
			assert(readbytes <= decompressedbytes);
			//printf("decompressing %d bytes...\n", nbytes);
			WaveletCompressor compressor;

			{ // swapping
			enum
			{
				BS3 = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_,
				BITSETSIZE = (BS3 + 7) / 8
			};

			unsigned char *buf = & waveletbuf[readbytes];
			for (int i = BITSETSIZE; i < nbytes; i+=4)
				swapbytes(buf+i, 4);
			}

			memcpy(compressor.compressed_data(), &waveletbuf[readbytes], nbytes);
			readbytes += nbytes;

			compressor.decompress(halffloat, nbytes, wtype, MYBLOCK);
		}

		fclose(f);
	}

	float load_block2(int ix, int iy, int iz, Real MYBLOCK[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_])
	{
		float zratio1, zratio2;

		//printf("trying to load <%s>\n", path.c_str());
		FILE * f = fopen(path.c_str(), "rb");

		assert(f);

		CompressedBlock compressedchunk = idx2chunk[_id(ix, iy, iz)];

		size_t start = compressedchunk.start;

		assert(start >= miniheader_bytes);
		assert(start < global_header_displacement);
		assert(start + compressedchunk.extent <= global_header_displacement);

		vector<unsigned char> compressedbuf(compressedchunk.extent);
		fseek(f, compressedchunk.start, SEEK_SET);
		fread(&compressedbuf.front(), compressedchunk.extent, 1, f);

		assert(!feof(f));

		size_t zz_bytes = compressedbuf.size();
		//static vector<unsigned char> waveletbuf(2 << 21);	// 21: 4MB, 22: 8MB, 28: 512MB
		static vector<unsigned char> waveletbuf(2 << 28);
		const size_t decompressedbytes = zdecompress0(&compressedbuf.front(), compressedbuf.size(), &waveletbuf.front(), waveletbuf.size());
		zratio1 = (1.0*decompressedbytes)/zz_bytes;
#if defined(VERBOSE)
		printf("zdecompressed %d bytes to %d bytes...(%.2lf)\n", zz_bytes, decompressedbytes, zratio1);
#endif
		int readbytes = 0;
		for(int i = 0; i<compressedchunk.subid; ++i)
		{
			int nbytes = * (int *) & waveletbuf[readbytes];
			nbytes = swapint(nbytes);
			readbytes += sizeof(int);
			readbytes += nbytes;

			assert(readbytes <= decompressedbytes);
		}

		{
			int nbytes = *(int *)&waveletbuf[readbytes];
			nbytes = swapint(nbytes);
			readbytes += sizeof(int);
			assert(readbytes <= decompressedbytes);
#if defined(VERBOSE)
			printf("wavelet decompressing %d bytes...\n", nbytes);
#endif
			WaveletCompressor compressor;

			{ // swapping
			enum
			{
				BS3 = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_,
				BITSETSIZE = (BS3 + 7) / 8
			};

			unsigned char *buf = & waveletbuf[readbytes];
			for (int i = BITSETSIZE; i < nbytes; i+=4)
				swapbytes(buf+i, 4);
			}
			memcpy(compressor.compressed_data(), &waveletbuf[readbytes], nbytes);
			readbytes += nbytes;

                        memcpy((void *) MYBLOCK, (void *)compressor.compressed_data(), nbytes);

			const int BS3 = (_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_)*sizeof(Real);
			zratio2 = (1.0*BS3)/nbytes;
#if defined(VERBOSE)
			printf("decompressed %d bytes to %d bytes...(%.2lf)\n", nbytes, BS3, zratio2);
#endif
		}

		fclose(f);

		return zratio1*zratio2;
	}
};

class Reader_WaveletCompressionMPI0: public Reader_WaveletCompression0
{
	const MPI_Comm comm;

public:

	Reader_WaveletCompressionMPI0(const MPI_Comm comm, const string path, int swapbytes, int wtype):
	Reader_WaveletCompression0(path,swapbytes,wtype), comm(comm)
	{

	}

	virtual void load_file()
	{
		int myrank;
        MPI_Comm_rank(comm, &myrank);

		if (myrank == 0)
			Reader_WaveletCompression0::load_file();

		//propagate primitive type data members
		{
			MPI_Bcast(&global_header_displacement, sizeof(global_header_displacement), MPI_CHAR, 0, comm);
			MPI_Bcast(&miniheader_bytes, sizeof(miniheader_bytes), MPI_CHAR, 0, comm);
			MPI_Bcast(&NBLOCKS, sizeof(NBLOCKS), MPI_CHAR, 0, comm);
			MPI_Bcast(totalbpd, sizeof(totalbpd), MPI_CHAR, 0, comm);
			MPI_Bcast(bpd, sizeof(bpd), MPI_CHAR, 0, comm);
			MPI_Bcast(&halffloat, sizeof(halffloat), MPI_CHAR, 0, comm);
			MPI_Bcast(&doswapping, sizeof(doswapping), MPI_CHAR, 0, comm);
		}

		size_t nentries = idx2chunk.size();

		MPI_Bcast(&nentries, sizeof(nentries), MPI_CHAR, 0, comm);

		if (myrank)
			idx2chunk.resize(nentries);

		const size_t nbytes = nentries * sizeof(CompressedBlock);
		char * const entries = (char *)&idx2chunk.front();

		MPI_Bcast(entries, nbytes, MPI_CHAR, 0, comm);
	}
};

#endif /* READER_WAVELETCOMPRESSION_H_00000000 */
