/*
 *  SerializerIO_WaveletCompression_MPI_Simple.h
 *  ... te lo do io il simple, porco porco
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <typeinfo>
#include <sstream>
#include <numeric>
#include <omp.h>
#include <Timer.h>


using namespace std;

#include "WaveletSerializationTypes.h"
#include "FullWaveletTransform.h"
#include "WaveletCompressor.h"

#include "CompressionEncoders.h"
#define		_ASYNC_IO_	1	// peh: hardcoded option, not in active code though
//#define	_WRITE_AT_ALL_	1	// peh:

//#define _USE_MORTON_
#if defined(_USE_MORTON_)
unsigned int Morton_3D_Encode_5bit(uint index1,uint index2,uint index3 )
{ // pack 3 5-bitindicesinto a 15-bit Morton code
        index1 &= 0x0000001f;
        index2 &= 0x0000001f;
        index3 &= 0x0000001f;
        index1 *= 0x01041041;
        index2 *= 0x01041041;
        index3 *= 0x01041041;
        index1 &= 0x10204081;
        index2 &= 0x10204081;
        index3 &= 0x10204081;
        index1 *= 0x00011111;
        index2 *= 0x00011111;
        index3 *= 0x00011111;
        index1 &= 0x12490000;
        index2 &= 0x12490000;
        index3 &= 0x12490000;
        return( (index1 >> 16 ) | (index2 >> 15 ) | (index3 >> 14 ) );
}
#endif

#if defined(_USE_ZEROBITS_)
void float_zero_bits(unsigned int *ul, int n)
{
        unsigned int _ul = *ul;

        if (n == 0)
                _ul = _ul;
        else if (n == 4)
                _ul = _ul & 0xfffffff0;
        else if (n == 8)
                _ul = _ul & 0xffffff00;
        else if (n == 12)
                _ul = _ul & 0xfffff000;
        else if (n == 16)
                _ul = _ul & 0xffff0000;

        *ul = _ul;
}
#endif

template<typename GridType, typename IterativeStreamer>
class SerializerIO_WaveletCompression_MPI_SimpleBlocking
{
protected:

	enum
	{
		NCHANNELS = IterativeStreamer::channels,
		NPTS = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_
	};

	enum //some considerations about the per-thread working set
	{
	  DESIREDMEM = (4 * 1024) * 1024,
	  ENTRYSIZE = sizeof(WaveletCompressor) + sizeof(int),
	  ENTRIES_CANDIDATE = DESIREDMEM / ENTRYSIZE,
	  ENTRIES = ENTRIES_CANDIDATE ? ENTRIES_CANDIDATE : 1,
	  BUFFERSIZE = ENTRIES * ENTRYSIZE, //sonnentanz
	  ALERT = (ENTRIES - 1) * ENTRYSIZE //got rain?
	};

	struct CompressionBuffer
	{
	  BlockMetadata hotblocks[ENTRIES];
	  unsigned char compressedbuffer[BUFFERSIZE];
	};

	struct TimingInfo { float total, fwt, encoding; };

	string binaryocean_title, binarylut_title, header;

	vector< BlockMetadata > myblockindices; //tells in which compressed chunk is any block, nblocks
	HeaderLUT lutheader;
	vector< size_t > lut_compression; //tells the number of compressed chunk, and where do they start, nchunks + 2
	vector< unsigned char > allmydata; //buffer with the compressed data
	size_t written_bytes, pending_writes, completed_writes;

	Real threshold;
	bool halffloat, verbosity;
	int wtype_read, wtype_write;	// peh

	vector< float > workload_total, workload_fwt, workload_encode; //per-thread cpu time for imbalance insight for fwt and encoding
	vector<CompressionBuffer> workbuffer; //per-thread compression buffer

	float _encode_and_flush(unsigned char inputbuffer[], int& bufsize, const int maxsize, BlockMetadata metablocks[], int& nblocks)
	{
		//0. setup
		//1. compress the data with zlib, obtain zptr, zbytes
		//2. obtain an offset from allmydata -> dstoffset
		//3. obtain a new entry in lut_compression -> idcompression
		//4. copy the [zptr,zptr+zbytes] into in allmydata, starting from dstoffset
		//5. for all blocks, set the myblockindices[blockids[i]] to a valid state
		//6. set nblocks to zero

		Timer timer; timer.start();
		const unsigned char * const zptr = inputbuffer;
		size_t zbytes = bufsize;
		int dstoffset = -1;
		int idcompression = -1;

		//1.
		/* ZLIB/LZ4/LZF/FPC/FPZIP/DUMMY ENCODING (LOSSLESS COMPRESSION) */
		{
			/* TODO: replace the following code with: zbytes = zcompress(inputbuffer, bufsize, maxsize); */
			z_stream myzstream = {0};
			//deflateInit(&myzstream, Z_DEFAULT_COMPRESSION); // Z_BEST_COMPRESSION
			deflateInit(&myzstream, Z_BEST_COMPRESSION);

			unsigned mah = maxsize;
			int err = deflate_inplace(&myzstream, inputbuffer, bufsize, &mah);

			assert(err == Z_OK);

			zbytes = myzstream.total_out;

			// fixed memory leak, Fabian (08 Oct 2014)
			err = deflateEnd(&myzstream);
			assert(err == Z_OK);
		}

		//2-3.
#pragma omp critical
		{
			dstoffset = written_bytes;
			written_bytes += zbytes;

			//exception: we have to resize allmydata
			if (written_bytes > allmydata.size())
			{
				//yield-based wait until writes complete
				while (pending_writes != completed_writes) sched_yield();

				//safely resize
				allmydata.resize(written_bytes);
			}

			idcompression = lut_compression.size();
			lut_compression.push_back(dstoffset);

			++pending_writes;
		}

		//4.
		assert(allmydata.size() >= written_bytes);
		memcpy(&allmydata.front() + dstoffset, zptr, zbytes);

#pragma omp atomic
		++completed_writes;

		//5.
		for(int i = 0; i < nblocks; ++i)
		{
			const int entry = metablocks[i].idcompression;
			assert(entry >= 0 && entry < myblockindices.size());

			myblockindices[entry] = metablocks[i];
			myblockindices[entry].idcompression = idcompression;
			myblockindices[entry].subid = i;
		}

		//6.
		bufsize = 0;
		nblocks = 0;

		return timer.stop();
	}


#define APPLY_MORTON()													\
	do														\
	{														\
	assert(_BLOCKSIZE_==32);											\
	Real tmp[_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_];									\
	for(int iz=0; iz<FluidBlock::sizeZ; iz++)									\
	for(int iy=0; iy<FluidBlock::sizeY; iy++)									\
	for(int ix=0; ix<FluidBlock::sizeX; ix++)									\
		tmp[Morton_3D_Encode_5bit(iz,iy,ix)] =  mysoabuffer[ix + _BLOCKSIZE_ * (iy + _BLOCKSIZE_ * iz)];	\
	memcpy(mysoabuffer, tmp, _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*sizeof(Real));					\
	} while (0);

	template<int channel>
	void _compress(const vector<BlockInfo>& vInfo, const int NBLOCKS, IterativeStreamer streamer)
	{
		int compress_threads = omp_get_max_threads();	// peh: recheck this on BGQ
		if (compress_threads < 1) compress_threads = 1;

#pragma omp parallel num_threads(compress_threads)
		{
		  const int tid = omp_get_thread_num();

		  CompressionBuffer & mybuf = workbuffer[tid];

			int mybytes = 0, myhotblocks = 0;

			float tfwt = 0, tencode = 0;
			Timer timer;
			timer.start();

#pragma omp for
			for(int i = 0; i < NBLOCKS; ++i)
			{
				Timer tw; tw.start();

				//wavelet compression
				{
					FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;

					WaveletCompressor compressor;

					Real * const mysoabuffer = &compressor.uncompressed_data()[0][0][0];
#if 0
					for(int iz=0; iz<FluidBlock::sizeZ; iz++)
						for(int iy=0; iy<FluidBlock::sizeY; iy++)
							for(int ix=0; ix<FluidBlock::sizeX; ix++)
								mysoabuffer[ix + _BLOCKSIZE_ * (iy + _BLOCKSIZE_ * iz)] = streamer.template operate<channel>(b(ix, iy, iz));
#else
					if(streamer.name() == "StreamerGridPointIterative")
					{
					for(int iz=0; iz<FluidBlock::sizeZ; iz++)
						for(int iy=0; iy<FluidBlock::sizeY; iy++)
							for(int ix=0; ix<FluidBlock::sizeX; ix++)
								mysoabuffer[ix + _BLOCKSIZE_ * (iy + _BLOCKSIZE_ * iz)] = streamer.template operate<channel>(b(ix, iy, iz));
					}
					else
					{
					IterativeStreamer mystreamer(b);
					for(int iz=0; iz<FluidBlock::sizeZ; iz++)
						for(int iy=0; iy<FluidBlock::sizeY; iy++)
							for(int ix=0; ix<FluidBlock::sizeX; ix++)
								mysoabuffer[ix + _BLOCKSIZE_ * (iy + _BLOCKSIZE_ * iz)] = mystreamer.operate(ix, iy, iz);
					}
#endif

					//wavelet digestion
#if defined(_USE_WAVZ_) /* WAVELET COMPRESSION */
					const int nbytes = (int)compressor.compress(this->threshold, this->halffloat, this->wtype_write);
					memcpy(mybuf.compressedbuffer + mybytes, &nbytes, sizeof(nbytes));
					mybytes += sizeof(nbytes);

					memcpy(mybuf.compressedbuffer + mybytes, compressor.compressed_data(), sizeof(unsigned char) * nbytes);

#elif defined(_USE_FPZIP_)
					const int inbytes = FluidBlock::sizeX * FluidBlock::sizeY * FluidBlock::sizeZ * sizeof(Real);
					int nbytes;

#if defined(_USE_MORTON_)
					APPLY_MORTON();
#endif

					int layout[4] = {_BLOCKSIZE_, _BLOCKSIZE_, _BLOCKSIZE_, 1};
					int fpzip_prec = (int)this->threshold;
					fpz_compress3D((void *)mysoabuffer, inbytes, layout, (void *) compressor.compressed_data(), (unsigned int *)&nbytes, (sizeof(Real)==4)?1:0, fpzip_prec);

					//printf("fpz_compress3D: from %d to %d bytes\n", inbytes, nbytes);

					memcpy(mybuf.compressedbuffer + mybytes, &nbytes, sizeof(nbytes));
					mybytes += sizeof(nbytes);
					memcpy(mybuf.compressedbuffer + mybytes, compressor.compressed_data(), sizeof(unsigned char) * nbytes);

#elif defined(_USE_DRAIN_)
					const int inbytes = FluidBlock::sizeX * FluidBlock::sizeY * FluidBlock::sizeZ * sizeof(Real);
					int nbytes;
#if defined(_USE_MORTON_)
					APPLY_MORTON();
#endif
					int drain_level = (int)this->threshold;
					char drain_level_value[16];
					sprintf(drain_level_value, "%d", drain_level);
					if(getenv("DRAIN_LEVEL") == NULL) setenv("DRAIN_LEVEL",drain_level_value, 1));
					int layout[4] = {_BLOCKSIZE_, _BLOCKSIZE_, _BLOCKSIZE_, 1};
					nbytes = drain_3df_buffer((float * const)mysoabuffer, layout[0], layout[1], layout[2], (unsigned char *)compressor.compressed_data());
#if VERBOSE
					printf("draind_3db_buffer: from %d to %d bytes\n", inbytes, nbytes);
#endif

					memcpy(mybuf.compressedbuffer + mybytes, &nbytes, sizeof(nbytes));
					mybytes += sizeof(nbytes);
					memcpy(mybuf.compressedbuffer + mybytes, compressor.compressed_data(), sizeof(unsigned char) * nbytes);

#elif defined(_USE_ZFP_)
					const int inbytes = FluidBlock::sizeX * FluidBlock::sizeY * FluidBlock::sizeZ * sizeof(Real);
					int nbytes;
#if defined(_USE_MORTON_)
					APPLY_MORTON();
#endif
//					double zfp_acc = 0.0;
//					if(getenv("ZFP_ACC")) zfp_acc = atof(getenv("ZFP_ACC"));
					double zfp_acc = this->threshold;
					int is_float = 1;
					int layout[4] = {_BLOCKSIZE_, _BLOCKSIZE_, _BLOCKSIZE_, 1};
					size_t nbytes_zfp;
					int status = zfp_compress_buffer(mysoabuffer, layout[0], layout[1], layout[2], zfp_acc, is_float, (unsigned char *)compressor.compressed_data(), &nbytes_zfp);
					nbytes = nbytes_zfp;
#if VERBOSE
					printf("zfp_compress status = %d, from %d to %d bytes = %d\n", status, inbytes, nbytes);
#endif

					memcpy(mybuf.compressedbuffer + mybytes, &nbytes, sizeof(nbytes));
					mybytes += sizeof(nbytes);
					memcpy(mybuf.compressedbuffer + mybytes, compressor.compressed_data(), sizeof(unsigned char) * nbytes);

#elif defined(_USE_SZ_)
					const int inbytes = FluidBlock::sizeX * FluidBlock::sizeY * FluidBlock::sizeZ * sizeof(Real);
					int nbytes;
#if defined(_USE_MORTON_)
					APPLY_MORTON();
#endif
					double sz_abs_acc = 0.0;
					double sz_rel_acc = 0.0;
					if(getenv("SZ_ABS_ACC")) sz_abs_acc = atof(getenv("SZ_ABS_ACC"));
					//SZ_ABS_ACC=$PARAM
					sz_abs_acc = (double) this->threshold;

					int layout[4] = {_BLOCKSIZE_, _BLOCKSIZE_, _BLOCKSIZE_, 1};

					int *bytes_sz = (int *)malloc(sizeof(int));
					unsigned char *compressed_sz = SZ_compress_args(SZ_FLOAT, (unsigned char *)mysoabuffer, bytes_sz, ABS, sz_abs_acc, sz_rel_acc, 0, 0, layout[2], layout[1], layout[0]);
					nbytes = *bytes_sz;
					memcpy(compressor.compressed_data(), compressed_sz, nbytes);
					free(bytes_sz);
					free(compressed_sz);

#if VERBOSE
					printf("SZ_compress_args: from %d to %d bytes\n", inbytes, nbytes);
#endif

					memcpy(mybuf.compressedbuffer + mybytes, &nbytes, sizeof(nbytes));
					mybytes += sizeof(nbytes);
					memcpy(mybuf.compressedbuffer + mybytes, compressor.compressed_data(), sizeof(unsigned char) * nbytes);

#elif defined(_USE_ISA_)
					const int inbytes = FluidBlock::sizeX * FluidBlock::sizeY * FluidBlock::sizeZ * sizeof(Real);
					int nbytes;
#if defined(_USE_MORTON_)
					APPLY_MORTON();
#endif
					double isa_rate = (double) this->threshold;
					int is_float = 1;
					int layout[4] = {_BLOCKSIZE_, _BLOCKSIZE_, _BLOCKSIZE_, 1};
					int status = isa_compress_buffer((char *)mysoabuffer, layout[0], layout[1], layout[2], isa_rate, is_float, (char *)compressor.compressed_data(), &nbytes);
#if VERBOSE
					printf("isa_compressed_buffer: from %d to %d bytes\n", inbytes, nbytes);
#endif

					memcpy(mybuf.compressedbuffer + mybytes, &nbytes, sizeof(nbytes));
					mybytes += sizeof(nbytes);
					memcpy(mybuf.compressedbuffer + mybytes, compressor.compressed_data(), sizeof(unsigned char) * nbytes);

#else /* NO COMPRESSION */
//#error	"Compression method for blocks has not been defined!"
					const int inbytes = FluidBlock::sizeX * FluidBlock::sizeY * FluidBlock::sizeZ * sizeof(Real);
					int nbytes = inbytes;
					memcpy(mybuf.compressedbuffer + mybytes, &nbytes, sizeof(nbytes));
					mybytes += sizeof(nbytes);

#if defined(_USE_ZEROBITS_)
					// set some bits to zero
					for(int iz=0; iz<FluidBlock::sizeZ; iz++)
					for(int iy=0; iy<FluidBlock::sizeY; iy++)
					for(int ix=0; ix<FluidBlock::sizeX; ix++)
						float_zero_bits((unsigned int *)&mysoabuffer[ix + _BLOCKSIZE_ * (iy + _BLOCKSIZE_ * iz)], _ZEROBITS_);
#endif

#if defined(_USE_SHUFFLE_)
#if !defined(_USE_MORTON_)
					shuffle((char *)mysoabuffer, nbytes, sizeof(Real));
#else
					assert(_BLOCKSIZE_==32);
					Real tmp[_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_];
					for(int iz=0; iz<FluidBlock::sizeZ; iz++)
					for(int iy=0; iy<FluidBlock::sizeY; iy++)
					for(int ix=0; ix<FluidBlock::sizeX; ix++)
						tmp[Morton_3D_Encode_5bit(ix,iy,iz)] =  mysoabuffer[ix + _BLOCKSIZE_ * (iy + _BLOCKSIZE_ * iz)];
					memcpy(mysoabuffer, tmp, _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*sizeof(Real));
                                        shuffle((char *)mysoabuffer, nbytes, sizeof(Real));
#endif
#endif
					memcpy(mybuf.compressedbuffer + mybytes, mysoabuffer, sizeof(unsigned char) * nbytes);
#endif
					mybytes += nbytes;
				}

				tfwt += tw.stop();

				//building the meta data
				{
					BlockMetadata curr = { i, myhotblocks, vInfo[i].index[0], vInfo[i].index[1], vInfo[i].index[2]};
					mybuf.hotblocks[myhotblocks] = curr;
					myhotblocks++;
				}

				if (mybytes >= ALERT || myhotblocks >= ENTRIES)
					tencode = _encode_and_flush(mybuf.compressedbuffer, mybytes, BUFFERSIZE, mybuf.hotblocks, myhotblocks);
			}

			if (mybytes > 0)
				tencode = _encode_and_flush(mybuf.compressedbuffer, mybytes, BUFFERSIZE, mybuf. hotblocks, myhotblocks);

			workload_total[tid] = timer.stop();
			workload_fwt[tid] = tfwt;
			workload_encode[tid] = tencode;
		}
	}

	virtual void _to_file(const MPI_Comm mycomm, const string fileName)
	{
		int mygid;
		int nranks;
		MPI_Comm_rank(mycomm, &mygid);
		MPI_Comm_size(mycomm, &nranks);

		MPI_Info myfileinfo;
		MPI_Info_create(&myfileinfo);

		string key("access_style");
		string val("write_once");
		MPI_Info_set(myfileinfo, (char*)key.c_str(), (char*)val.c_str());

		MPI_File myfile;
#if defined(_WRITE_AT_ALL_)
		MPI_File_open(mycomm, (char*)fileName.c_str(),  MPI_MODE_WRONLY | MPI_MODE_CREATE, myfileinfo, &myfile);
#else
		MPI_File_open(MPI_COMM_SELF, (char*)fileName.c_str(),  MPI_MODE_WRONLY | MPI_MODE_CREATE, myfileinfo, &myfile);
#endif
		size_t current_displacement = 0;

		//write the mini-header
		{
			size_t blank_address = -1;

			const int miniheader_bytes = sizeof(blank_address) + binaryocean_title.size();

			current_displacement += miniheader_bytes;
		}

		//write the buffer - alias the binary ocean
		{
			// apparently this suck:
			// myfile.Seek_shared(current_displacement, MPI_SEEK_SET);
			// myfile.Write_ordered(&allmydata.front(), written_bytes, MPI_CHAR);
			// current_displacement = myfile.Get_position_shared();
			// so here we do it manually. so nice!

			size_t myfileoffset = 0;
			MPI_Exscan(&written_bytes, &myfileoffset, 1, MPI_UINT64_T, MPI_SUM, mycomm);

			if (mygid == 0)
				myfileoffset = 0;

			MPI_Status status;
#if defined(_WRITE_AT_ALL_)
			MPI_File_write_at_all(myfile, current_displacement + myfileoffset, &allmydata.front(), written_bytes, MPI_CHAR, &status);
#else
			MPI_File_write_at(myfile, current_displacement + myfileoffset, &allmydata.front(), written_bytes, MPI_CHAR, &status);
#endif

			//here we update current_displacement by broadcasting the total written bytes from rankid = nranks -1
			size_t total_written_bytes = myfileoffset + written_bytes;
			MPI_Bcast(&total_written_bytes, 1, MPI_UINT64_T, nranks - 1, mycomm);

			current_displacement += total_written_bytes;
		}

		//go back at the blank address and fill it with the displacement
		if (mygid == 0)
		{
			MPI_Status status;
			MPI_File_write(myfile, &current_displacement, sizeof(current_displacement), MPI_CHAR, &status);
			MPI_File_write(myfile, (void*)binaryocean_title.c_str(), binaryocean_title.size(), MPI_CHAR, &status);
		}

		//write the header
		{
			const size_t header_bytes = header.size();

			MPI_Status status;
			if (mygid == 0)
				MPI_File_write_at(myfile, current_displacement, (void*)header.c_str(), header_bytes, MPI_CHAR, &status);

			current_displacement += header_bytes;
		}

		//write block metadata
		{
			const int metadata_bytes = myblockindices.size() * sizeof(BlockMetadata);

			MPI_Status status;
#if defined(_WRITE_AT_ALL_)
			MPI_File_write_at_all(myfile, current_displacement + mygid * metadata_bytes, &myblockindices.front(), metadata_bytes, MPI_CHAR, &status);
#else
			MPI_File_write_at(myfile, current_displacement + mygid * metadata_bytes, &myblockindices.front(), metadata_bytes, MPI_CHAR, &status);
#endif
			current_displacement += metadata_bytes * nranks;
		}

		//write the lut title
		{
			const int title_bytes = binarylut_title.size();

			MPI_Status status;
			if (mygid == 0)
				MPI_File_write_at(myfile, current_displacement, (void*)binarylut_title.c_str(), title_bytes, MPI_CHAR, &status);

			current_displacement += title_bytes;
		}

		//write the local buffer entries
		{
			assert(lut_compression.size() == 0);

			const int lutheader_bytes = sizeof(lutheader);

			MPI_Status status;
#if defined(_WRITE_AT_ALL_)
			MPI_File_write_at_all(myfile, current_displacement + mygid * lutheader_bytes, &lutheader, lutheader_bytes, MPI_CHAR, &status);
#else
			MPI_File_write_at(myfile, current_displacement + mygid * lutheader_bytes, &lutheader, lutheader_bytes, MPI_CHAR, &status);
#endif
		}

		MPI_Info_free(&myfileinfo);
		MPI_File_close(&myfile); //bon voila tu vois ou quoi
	}

	float _profile_report(const char * const workload_name, vector<float>& workload, const MPI_Comm mycomm, bool isroot)
	{
		float tmin = *std::min_element(workload.begin(), workload.end());
		float tmax = *std::max_element(workload.begin(), workload.end());
		float tsum = std::accumulate(workload.begin(), workload.end(), 0.f);

		MPI_Reduce(isroot ? MPI_IN_PLACE : &tmin, &tmin, 1, MPI_FLOAT, MPI_MIN, 0, mycomm);
		MPI_Reduce(isroot ? MPI_IN_PLACE : &tsum, &tsum, 1, MPI_FLOAT, MPI_SUM, 0, mycomm);
		MPI_Reduce(isroot ? MPI_IN_PLACE : &tmax, &tmax, 1, MPI_FLOAT, MPI_MAX, 0, mycomm);

		int comm_size;
		MPI_Comm_size(mycomm, &comm_size);
		tsum /= comm_size * workload.size();

		if (isroot)
			printf("TLP %-10s min:%.3e s avg:%.3e s max:%.3e s imb:%.0f%%\n",
				   workload_name, tmin, tsum, tmax, (tmax - tmin) / tsum * 100);

		return tsum;
	}

	template<int channel>
	void _write(GridType & inputGrid, string fileName, IterativeStreamer streamer)
	{
		const vector<BlockInfo> infos = inputGrid.getBlocksInfo();
		const int NBLOCKS = infos.size();

		//prepare the headers
		{
			this->binaryocean_title = "\n==============START-BINARY-OCEAN==============\n";
			this->binarylut_title = "\n==============START-BINARY-LUT==============\n";

			{
				const int xtotalbpd = inputGrid.getBlocksPerDimension(0);
				const int ytotalbpd = inputGrid.getBlocksPerDimension(1);
				const int ztotalbpd = inputGrid.getBlocksPerDimension(2);

				const int xbpd = inputGrid.getResidentBlocksPerDimension(0);
				const int ybpd = inputGrid.getResidentBlocksPerDimension(1);
				const int zbpd = inputGrid.getResidentBlocksPerDimension(2);

				const double xExtent = inputGrid.getH()*xtotalbpd*_BLOCKSIZE_;
				const double yExtent = inputGrid.getH()*ytotalbpd*_BLOCKSIZE_;
				const double zExtent = inputGrid.getH()*ztotalbpd*_BLOCKSIZE_;

				std::stringstream ss;

				ss << "\n==============START-ASCI-HEADER==============\n";

				{
					int one = 1;
					bool isone = *(char *)(&one);

					ss << "Endianess: " << (isone ? "little" : "big") << "\n";
				}

				ss << "sizeofReal: " << sizeof(Real) << "\n";
				ss << "sizeofsize_t: " << sizeof(size_t) << "\n";
				ss << "sizeofBlockMetadata: " << sizeof(BlockMetadata) << "\n";
				ss << "sizeofHeaderLUT: " << sizeof(HeaderLUT) << "\n";
				ss << "sizeofCompressedBlock: " << sizeof(CompressedBlock) << "\n";
				ss << "Blocksize: " << _BLOCKSIZE_ << "\n";
				ss << "Blocks: " << xtotalbpd << " x "  << ytotalbpd << " x " << ztotalbpd  << "\n";
				ss << "Extent: " << xExtent << " " << yExtent << " " << zExtent << "\n";
				ss << "SubdomainBlocks: " << xbpd << " x "  << ybpd << " x " << zbpd  << "\n";
				ss << "HalfFloat: " << (this->halffloat ? "yes" : "no") << "\n";
#if defined(_USE_WAVZ_)
				ss << "Wavelets: " << WaveletsOnInterval::ChosenWavelets_GetName(this->wtype_write) << "\n";
#elif defined(_USE_FPZIP_)
				ss << "Wavelets: " << "fpzip" << "\n";
#elif defined(_USE_DRAIN_)
				ss << "Wavelets: " << "drain" << "\n";
#elif defined(_USE_ZFP_)
				ss << "Wavelets: " << "zfp" << "\n";
#elif defined(_USE_SZ_)
				ss << "Wavelets: " << "sz" << "\n";
#elif defined(_USE_ISA_)
				ss << "Wavelets: " << "isa" << "\n";
#elif defined(_USE_SHUFFLE_)
				ss << "Wavelets: " << "shuffle" << "\n";
#else
				ss << "Wavelets: " << "none" << "\n";
#endif
				ss << "WaveletThreshold: " << threshold << "\n";
#if defined(_USE_ZLIB_)
				ss << "Encoder: " << "zlib" << "\n";
#elif defined(_USE_LZ4_)
                                ss << "Encoder: " << "lz4" << "\n";
#elif defined(_USE_LZMA_)
                                ss << "Encoder: " << "lzma" << "\n";
#elif defined(_USE_ZSTD_)
                                ss << "Encoder: " << "zstd" << "\n";
#elif defined(_USE_ZSTD0_)
                                ss << "Encoder: " << "zstd0" << "\n";
#elif defined(_USE_BLOSC_)
                                ss << "Encoder: " << "blosc" << "\n";
#else
                                ss << "Encoder: " << "none" << "\n";
#endif
				ss << "==============START-BINARY-METABLOCKS==============\n";

				this->header = ss.str();
			}
		}

		//compress my data, prepare for serialization
		{
			written_bytes = pending_writes = completed_writes = 0;

			if (allmydata.size() == 0)
			{
				const size_t speculated_compression_rate = 10;
				allmydata.resize(NBLOCKS * sizeof(Real) * NPTS);
			}

			myblockindices.clear();
			myblockindices.resize(NBLOCKS);

			lut_compression.clear();

			_compress<channel>(infos, infos.size(), streamer);

			//manipulate the file data (allmydata, lut_compression, myblockindices)
			//so that they are file-friendly
			{
				const int nchunks = lut_compression.size();
				const size_t extrabytes = lut_compression.size() * sizeof(size_t);
				const char * const lut_ptr = (char *)&lut_compression.front();

				allmydata.insert(allmydata.begin() + written_bytes, lut_ptr, lut_ptr + extrabytes);
				lut_compression.clear();

				HeaderLUT newvalue = { written_bytes + extrabytes, nchunks };
				lutheader = newvalue;

				written_bytes += extrabytes;
			}
		}

		const MPI_Comm mycomm = inputGrid.getCartComm();
		int mygid;
		int comm_size;
		MPI_Comm_rank(mycomm, &mygid);
		MPI_Comm_size(mycomm, &comm_size);

		//write into the file
		Timer timer; timer.start();
		_to_file(mycomm, fileName);
		vector<float> workload_file(1, timer.stop());

		//just a report now
		if (verbosity)
		{
			size_t aggregate_written_bytes = -1;

			MPI_Reduce(&written_bytes, &aggregate_written_bytes, 1, MPI_UINT64_T, MPI_SUM, 0, mycomm);
			const bool isroot = mygid == 0;
			if (mygid == 0)
				printf("Channel %d: %.2f kB, wavelet-threshold: %.1e, compr. rate: %.2f\n",
				   channel, aggregate_written_bytes/1024.,
				   threshold, NPTS * sizeof(Real) * NBLOCKS * comm_size / (float) aggregate_written_bytes);

			const float tavgcompr = _profile_report("Compr", workload_total, mycomm, isroot);
			const float tavgfwt =_profile_report("FWT+decim", workload_fwt, mycomm, isroot);
			const float tavgenc =_profile_report("Encoding", workload_encode, mycomm, isroot);
			const float tavgio =_profile_report("FileIO", workload_file, mycomm, isroot);
			const float toverall = tavgio + tavgcompr;

			if (isroot)
			{
				printf("Time distribution: %+5s:%.0f%% %+5s:%.0f%% %+5s:%.0f%% %+5s:%.0f%%\n",
					   "FWT", tavgfwt / toverall * 100,
					   "ENC", tavgenc / toverall * 100,
					   "IO", tavgio / toverall * 100,
					   "Other",  (tavgcompr - tavgfwt - tavgenc)/ toverall * 100);
			}
		}
	}

	void _read(string path)
	{
		//THE FIRST PART IS SEQUENTIAL
		//THE SECOND ONE IS RANDOM ACCESS

		size_t global_header_displacement = -1;
		int NBLOCKS = -1;
		int totalbpd[3] = {-1, -1, -1};
		int bpd[3] = { -1, -1, -1};
		string binaryocean_title = "\n==============START-BINARY-OCEAN==============\n";
		const int miniheader_bytes = sizeof(size_t) + binaryocean_title.size();

		vector<BlockMetadata> metablocks;

		//random access data structures: meta2subchunk, lutchunks;
		map<int, map<int, map<int, CompressedBlock > > > meta2subchunk;
		vector<size_t> lutchunks;

		{
			FILE * file = fopen(path.c_str(), "rb");

			assert(file);

			//reading the header and mini header
			{
				size_t header_displacement = -1;
				fread(&header_displacement, sizeof(size_t), 1, file);

				fseek(file, header_displacement, SEEK_SET);
				global_header_displacement = header_displacement;

				char buf[1024];
				fgets(buf, sizeof(buf), file);
				fgets(buf, sizeof(buf), file);
				assert(string("==============START-ASCI-HEADER==============\n") == string(buf));

				fscanf(file, "Endianess:  %s\n", buf);
				{
					int one = 1;
					bool isone = *(char *)(&one);

					if (isone)
						assert(string(buf) == "little");
					else
						assert(string(buf) == "big");
				}


				int sizeofreal = -1;
				fscanf(file, "sizeofReal:  %d\n", &sizeofreal);
				assert(sizeof(Real) == sizeofreal);

				int bsize = -1;
				fscanf(file, "Blocksize: %d\n", &bsize);
				assert(bsize == _BLOCKSIZE_);
				fscanf(file, "Blocks: %d x %d x %d\n", totalbpd, totalbpd + 1, totalbpd + 2);
				fscanf(file, "SubdomainBlocks: %d x %d x %d\n", bpd, bpd + 1, bpd + 2);

				fscanf(file, "HalfFloat: %s\n", buf);
				this->halffloat = (string(buf) == "yes");

				fscanf(file, "Wavelets: %s\n", buf);
				assert(buf == string(WaveletsOnInterval::ChosenWavelets_GetName(this->wtype_read)));

				fscanf(file, "Encoder: %s\n", buf);
#if defined(_USE_ZLIB_)
				assert(buf == string("zlib"));
#else	/* _USE_LZ4_ */
				assert(buf == string("lz4"));
#endif


				fgets(buf, sizeof(buf), file);
				assert(string("==============START-BINARY-METABLOCKS==============\n") == string(buf));

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
					//printf("reading metablock %d -> %d %d %d  cid %d\n", i, entry.ix, entry.iy, entry.iz, entry.idcompression);
					assert(entry.idcompression >= 0 && entry.idcompression < bpd[0] * bpd[1] * bpd[2]);
					metablocks[i] = entry;
				}
			}

			//reading the compression lut
			{
				char buf[1024];
				fgetc(file);//buf, sizeof(buf), file);
				fgets(buf, sizeof(buf), file);
				assert(string("==============START-BINARY-LUT==============\n") == string(buf));


				//printf("reading compression lut..at location %d\n", ftell(file));

				bool done = false;

				size_t base = miniheader_bytes;

				const int BPS = bpd[0] * bpd[1] * bpd[2];
				assert(NBLOCKS % BPS == 0);
				const int SUBDOMAINS = NBLOCKS / BPS;

				vector<HeaderLUT> headerluts(SUBDOMAINS); //oh mamma mia
				fread(&headerluts.front(), sizeof(HeaderLUT), SUBDOMAINS, file);

				for(int s = 0, currblock = 0; s < SUBDOMAINS; ++s)
				{
					const int nglobalchunks = lutchunks.size();

					assert(!feof(file));

					/* this is history
					 size_t nchunks = -1;
					 fread(&nchunks, sizeof(nchunks), 1, file);
					 //printf("this buffer has %d entries\n", nchunks);

					 vector<size_t> mylut(nchunks + 1);
					 fread(&mylut.front(), sizeof(size_t), mylut.size(), file);
					 const size_t myamount = mylut.back();
					 printf("my amount is %d\n", myamount);
					 mylut.pop_back();
					 */

					const int nchunks = headerluts[s].nchunks;
					const size_t myamount = headerluts[s].aggregate_bytes;
					const size_t lutstart = base + myamount - sizeof(size_t) * nchunks;
					//printf("my header lut is %d %d (0x%x 0x%x)\n", headerluts[s].aggregate_bytes, headerluts[s].nchunks,
					//	   headerluts[s].aggregate_bytes, headerluts[s].nchunks);
					//read the lut
					fseek(file, lutstart, SEEK_SET);
					vector<size_t> mylut(nchunks);
					fread(&mylut.front(), sizeof(size_t), nchunks, file);

					//for(int i=0; i< nchunks; ++i)
					//	printf("reading %d) 0x%x\n", i, mylut[i]);

					for(int i=0; i< nchunks; ++i)
						assert(mylut[i] < myamount);

					for(int i=1; i< nchunks; ++i)
						assert(mylut[i-1] < mylut[i]);

					//printf("goodu p to here\n");

					//bump the lut with the base
					//for(int i=0; i< nchunks; ++i)
					//	mylut[i] += base;

					//exit(0);
					//compute the global positioning of the compressed chunks within the file
					for(int i = 0; i < mylut.size(); ++i)
					{
						assert(mylut[i] < myamount);
						mylut[i] += base;
					}

					assert(myamount > 0);
					base += myamount;
					//printf("new base is 0x%x\n", base);
					//printf("global_header_displacement is 0x%x\n", global_header_displacement);
					assert(base <= global_header_displacement);

					//compute the base for this blocks
					for(int i = 0; i < BPS; ++i, ++currblock)
						metablocks[currblock].idcompression += nglobalchunks;

					lutchunks.insert(lutchunks.end(), mylut.begin(), mylut.end());
				}

				//printf("minheader takes %d bytes\n", miniheader_bytes);
				//printf("my base is now 0x%x whereas my header is at 0x%x -> discrepancy is %d bytes\n", base, global_header_displacement, global_header_displacement - base);
				assert(base == global_header_displacement);
				lutchunks.push_back(base);

				{
					int c = fgetc(file);

					do
					{
						//printf("shouldnt be here! 0x%x\n", c);
						c = fgetc(file);
					}
					while (! feof(file) );
				}
			}

			fclose(file);
		}

		for(int i = 0; i < NBLOCKS ; ++i)
		{
			BlockMetadata entry = metablocks[i];

			assert(entry.idcompression >= 0);
			assert(entry.idcompression < lutchunks.size()-1);

			size_t start_address = lutchunks[entry.idcompression];
			size_t end_address = lutchunks[entry.idcompression + 1];

			assert(start_address < end_address);
			assert(end_address <= global_header_displacement);
			assert( start_address < global_header_displacement );

			CompressedBlock compressedblock = { start_address, end_address - start_address, entry.subid };

			meta2subchunk[entry.iz][entry.iy][entry.ix] = compressedblock;
		}

		//at this point we can access every block
		{
			//lets try with block 7, 3, 5
			FILE * f = fopen(path.c_str(), "rb");
			assert(f);

			const int ix = 0;
			const int iy = 2;
			const int iz = 0;

			assert(ix >= 0 && ix < totalbpd[0]);
			assert(iy >= 0 && iy < totalbpd[1]);
			assert(iz >= 0 && iz < totalbpd[2]);

			assert(meta2subchunk.find(iz) != meta2subchunk.end());
			assert(meta2subchunk[iz].find(iy) != meta2subchunk[iz].end());
			assert(meta2subchunk[iz][iy].find(ix) != meta2subchunk[iz][iy].end());

			CompressedBlock compressedchunk = meta2subchunk[iz][iy][ix];

			size_t start = compressedchunk.start;
			assert(start >= miniheader_bytes);
			assert(start < global_header_displacement);
			assert(start + compressedchunk.extent < global_header_displacement);

			vector<unsigned char> compressedbuf(compressedchunk.extent);
			fseek(f, compressedchunk.start, SEEK_SET);
			fread(&compressedbuf.front(), compressedchunk.extent, 1, f);
			assert(!feof(f));
			//printf("my steart: 0x%x extent is %d my subid %d\n" , compressedchunk.start, compressedchunk.extent, compressedchunk.subid);


			vector<unsigned char> waveletbuf(4 << 20);
			const size_t decompressedbytes = zdecompress(&compressedbuf.front(), compressedbuf.size(), &waveletbuf.front(), waveletbuf.size());
			//printf("decompressed bytes is %d\n", decompressedbytes);
			int readbytes = 0;
			for(int i = 0; i<compressedchunk.subid; ++i)
			{
				int nbytes = *(int *)&waveletbuf[readbytes];
				readbytes += sizeof(int);
				assert(readbytes <= decompressedbytes);
				//printf("scanning nbytes...%d\n", nbytes);
				readbytes += nbytes;
				assert(readbytes <= decompressedbytes);

			}

			Real MYBLOCK[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];

			{
				int nbytes = *(int *)&waveletbuf[readbytes];
				readbytes += sizeof(int);
				assert(readbytes <= decompressedbytes);

				//printf("OK MY BYTES ARE: %d\n", nbytes);

				WaveletCompressor compressor;
				memcpy(compressor.compressed_data(), &waveletbuf[readbytes], nbytes);
				readbytes += nbytes;

				compressor.decompress(halffloat, nbytes, wtype_read, MYBLOCK);
			}

			printf("OK FINAL TEST: THE DATA\n");
			for(int iz = 0; iz< _BLOCKSIZE_; ++iz)
				for(int iy = 0; iy< _BLOCKSIZE_; ++iy)
					for(int ix = 0; ix< _BLOCKSIZE_; ++ix)
						printf("%d %d %d: %e\n", ix, iy, iz, MYBLOCK[iz][iy][ix]);

			fclose(f);
		}
	}

public:

	void set_threshold(const Real threshold) { this->threshold = threshold; }

	void set_wtype_write(const int wtype) { this->wtype_write = wtype; }
	void set_wtype_read(const int wtype) { this->wtype_read = wtype; }

	void float16() { halffloat = true; }

	void verbose() { verbosity = true; }

	SerializerIO_WaveletCompression_MPI_SimpleBlocking():
	threshold(0), halffloat(false), verbosity(false),
	workload_total(omp_get_max_threads()), workload_fwt(omp_get_max_threads()), workload_encode(omp_get_max_threads()),
	workbuffer(omp_get_max_threads()),
	written_bytes(0), pending_writes(0)
	{
		wtype_write = 1;	// peh
		wtype_read = 1;		// peh
	}

	template< int channel >
	void Write(GridType & inputGrid, string fileName, IterativeStreamer streamer = IterativeStreamer())
	{
		std::stringstream ss;
		ss << "." << streamer.name() << ".channel"  << channel;

		_write<channel>(inputGrid, fileName + ss.str(), streamer);
	}

	void Read(string fileName, IterativeStreamer streamer = IterativeStreamer())
	{
		for(int channel = 0; channel < NCHANNELS; ++channel)
		{
			std::stringstream ss;
			ss << "." << streamer.name() << ".channel"  << channel;

			_read(fileName + ss.str());
		}
	}
};

template<typename GridType, typename IterativeStreamer>
class SerializerIO_WaveletCompression_MPI_Simple : public SerializerIO_WaveletCompression_MPI_SimpleBlocking<GridType, IterativeStreamer>
{
	size_t nofcalls;
	vector<MPI_Request> pending_requests;
	MPI_File myopenfile;

	void _wait_all_quiet()
	{
		//wait for pending requests
		assert(pending_requests.size() > 0); //bah ouaie ou quois, tu vois bon voila
		MPI_Waitall(pending_requests.size(), &pending_requests.front(), MPI_STATUSES_IGNORE);
		pending_requests.clear();

		//close the split collective io
		MPI_Status status;
		MPI_File_write_ordered_end(myopenfile, &this->allmydata.front(), &status);

		//e buonanotte
		MPI_File_close(&myopenfile);
	}

	//et bon, la on ce lance le fleurs
	void _to_file(const MPI_Comm mycomm, const string fileName)
	{
		if (nofcalls)
			_wait_all_quiet();

		int mygid;
		int nranks;
		MPI_Comm_rank(mycomm, &mygid);
		MPI_Comm_size(mycomm, &nranks);

		MPI_Info myfileinfo;
		MPI_Info_create(&myfileinfo);

		string key("access_style");
		string val("write_once");
		MPI_Info_set(myfileinfo, (char*)key.c_str(), (char*)val.c_str());

		MPI_File_open(mycomm, (char*)fileName.c_str(),  MPI_MODE_WRONLY | MPI_MODE_CREATE, myfileinfo, &myopenfile);

		size_t current_displacement = 0;

		//write the mini-header
		{
			size_t blank_address = -1;

			const int miniheader_bytes = sizeof(blank_address) + this->binaryocean_title.size();

			current_displacement += miniheader_bytes;
		}

		//write the buffer - alias the binary ocean
		{
			MPI_File_seek_shared(myopenfile, current_displacement, MPI_SEEK_SET);

			MPI_File_write_ordered_begin(myopenfile, &this->allmydata.front(), this->written_bytes, MPI_CHAR);

			MPI_File_get_position_shared(myopenfile, &current_displacement);
		}

		//go back at the blank address and fill it with the displacement
		if (mygid == 0)
		{
#if defined(_ASYNC_IO_)
			MPI_Request reqA, reqB;
			MPI_File_iwrite(myopenfile, &current_displacement, sizeof(current_displacement), MPI_CHAR, &reqA);
			MPI_File_iwrite(myopenfile, (void*)(this->binaryocean_title.c_str()), this->binaryocean_title.size(), MPI_CHAR, &reqB);
			pending_requests.push_back(reqA);
			pending_requests.push_back(reqB);
#else
			MPI_Status status;
			MPI_File_write(myopenfile, &current_displacement, sizeof(current_displacement), MPI_CHAR, &status) ;
			MPI_File_write(myopenfile, (void*)(this->binaryocean_title.c_str()), this->binaryocean_title.size(), MPI_CHAR, &status) ;
#endif
		}

		//write the header
		{
			const size_t header_bytes = this->header.size();

			if (mygid == 0)
            		{
#if defined(_ASYNC_IO_)
				MPI_Request req;
                		MPI_File_iwrite_at(myopenfile, current_displacement, (void*)(this->header.c_str()), header_bytes, MPI_CHAR, &req);
				pending_requests.push_back(req);
#else
                		MPI_Status status;
                		MPI_File_write_at(myopenfile, current_displacement, (void*)(this->header.c_str()), header_bytes, MPI_CHAR, &status);
#endif
			}

			current_displacement += header_bytes;
		}

		//write block metadata
		{
			const int metadata_bytes = this->myblockindices.size() * sizeof(BlockMetadata);
#if defined(_ASYNC_IO_)
			MPI_Request req;
			MPI_File_iwrite_at(myopenfile, current_displacement + mygid * metadata_bytes, &this->myblockindices.front(), metadata_bytes, MPI_CHAR, &req);
			pending_requests.push_back(req);
#else
			MPI_Status status;
			MPI_File_write_at(myopenfile, current_displacement + mygid * metadata_bytes, &this->myblockindices.front(), metadata_bytes, MPI_CHAR, &status);
#endif
			current_displacement += metadata_bytes * nranks;
		}

		//write the lut title
		{
			const int title_bytes = this->binarylut_title.size();

			if (mygid == 0)
            		{
#if defined(_ASYNC_IO_)
				MPI_Request req;
				MPI_File_iwrite_at(myopenfile, current_displacement, (void*)(this->binarylut_title.c_str()), title_bytes, MPI_CHAR, &req);
				pending_requests.push_back(req);
#else
				MPI_Status status;
				MPI_File_write_at(myopenfile, current_displacement, (void*)(this->binarylut_title.c_str()), title_bytes, MPI_CHAR, &status);
#endif
            		}

			current_displacement += title_bytes;
		}

		//write the local buffer entries
		{
			assert(this->lut_compression.size() == 0);

			const int lutheader_bytes = sizeof(this->lutheader);

#if defined(_ASYNC_IO_)
			MPI_Request req;
			MPI_File_iwrite_at(myopenfile, current_displacement + mygid * lutheader_bytes, &this->lutheader, lutheader_bytes, MPI_CHAR, &req);
			pending_requests.push_back(req);
#else
			MPI_Status status;
			MPI_File_write_at(myopenfile, current_displacement + mygid * lutheader_bytes, &this->lutheader, lutheader_bytes, MPI_CHAR, &status) ;
#endif
		}

#if defined(_ASYNC_IO_)
		++nofcalls;
#endif
	}

public:
	SerializerIO_WaveletCompression_MPI_Simple():
	SerializerIO_WaveletCompression_MPI_SimpleBlocking<GridType, IterativeStreamer>(),
	nofcalls(0)
	{
	}

	void force_close() { _wait_all_quiet(); }

	~SerializerIO_WaveletCompression_MPI_Simple()
	{
		if (nofcalls)
			_wait_all_quiet();
	}
};
