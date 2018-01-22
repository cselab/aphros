/*
 *  WaveletSerializationTypes.h
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */
#pragma once

struct BlockMetadata { int idcompression, subid, ix, iy, iz; }  __attribute__((packed));
struct HeaderLUT { size_t aggregate_bytes; int nchunks; }  __attribute__((packed));
struct CompressedBlock{ size_t start, extent; int subid; }  __attribute__((packed));

