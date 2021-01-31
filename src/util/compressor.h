// Created by Fabian Wermelinger on 03.12.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <limits.h>
#include <cassert>

#include "macros.h"
#if USEFLAG(FPZIP)
#include <fpzip.h>
#endif

namespace compression {
// Envelope with compression meta data
struct Envelope {
  void* buf;
  size_t cbytes; // compressed bytes
  size_t ubytes; // uncompressed bytes
};

// Base interface for floating point compressors
class FPCompressor {
 public:
  FPCompressor(const bool c) : compress_(c) {}
  virtual ~FPCompressor() {}

  FPCompressor(const FPCompressor& c) = delete;
  FPCompressor& operator=(const FPCompressor& c) = delete;

  virtual Envelope GetEnvelope() const = 0; // return current envelope
  virtual Envelope Compress() = 0; // encode
  virtual size_t Decompress() = 0; // decode
  bool IsActive() const {
    return compress_;
  }

 protected:
  bool compress_;
};

template <typename Scal>
class PassThrough : public FPCompressor {
 public:
  PassThrough(Scal* data, const size_t N, const bool = false, const size_t = 0)
      : FPCompressor(false) {
    const size_t dbytes = N * sizeof(Scal); // data bytes
    env_.buf = static_cast<void*>(data);
    env_.cbytes = dbytes;
    env_.ubytes = dbytes;
  }

  Envelope GetEnvelope() const override {
    return env_;
  }

  Envelope Compress() override {
    return env_;
  }

  size_t Decompress() override {
    return 0;
  }

 private:
  Envelope env_; // compression meta data
};

#if USEFLAG(FPZIP)
// Floating point compressor based on the FPZIP library
template <typename Scal>
class FPZIP : public FPCompressor {
 public:
  FPZIP(
      Scal* data, const size_t N, const bool active = true,
      const size_t thresh = 1024 /* byte */)
      : FPCompressor(N * sizeof(Scal) > thresh), data_(data), N_(N) {
    type_ = (4 == sizeof(Scal)) ? FPZIP_TYPE_FLOAT : FPZIP_TYPE_DOUBLE;
    const size_t dsize =
        (type_ == FPZIP_TYPE_FLOAT ? sizeof(float) : sizeof(double));
    prec_ = CHAR_BIT * dsize; // lossless compression
    const size_t dbytes = N_ * dsize; // data bytes
    if (!active) {
      this->compress_ = false;
    }
    if (this->compress_) {
      env_.buf = new unsigned char[dbytes];
      env_.cbytes = 0;
      env_.ubytes = dbytes;
    } else {
      env_.buf = static_cast<void*>(data_);
      env_.cbytes = dbytes;
      env_.ubytes = dbytes;
    }
  }
  ~FPZIP() {
    if (this->compress_) {
      delete[] static_cast<unsigned char*>(env_.buf);
    }
  }

  Envelope GetEnvelope() const override {
    return env_;
  }

  Envelope Compress() override {
    if (this->compress_) {
      fpz_ = fpzip_write_to_buffer(env_.buf, env_.ubytes);
      fpz_->type = type_;
      fpz_->prec = prec_;
      fpz_->nx = N_;
      fpz_->ny = 1;
      fpz_->nz = 1;
      fpz_->nf = 1;
      env_.cbytes = fpzip_write(fpz_, data_);
      fpzip_write_close(fpz_);
    }
    return env_;
  }

  size_t Decompress() override {
    if (!this->compress_) {
      return 0;
    }
    fpz_ = fpzip_read_from_buffer(env_.buf);
    fpz_->type = type_;
    fpz_->prec = prec_;
    fpz_->nx = N_;
    fpz_->ny = 1;
    fpz_->nz = 1;
    fpz_->nf = 1;
    const size_t cbytes = fpzip_read(fpz_, data_);
    fpzip_read_close(fpz_);
    return cbytes;
  }

 private:
  Envelope env_; // compression meta data
  Scal* data_; // pointer to data
  const size_t N_;
  int type_, prec_;
  FPZ* fpz_; // meta data for compression stream
};
#endif
} // namespace compression
