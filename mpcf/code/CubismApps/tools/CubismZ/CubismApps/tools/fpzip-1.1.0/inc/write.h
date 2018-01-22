#ifndef FPZIP_WRITE_H
#define FPZIP_WRITE_H

#define subsize(T, n) (CHAR_BIT * sizeof(T) * (n) / 32)

// file writer for compressed data
#if FPZIP_BLOCK_SIZE > 1
class RCfileencoder : public RCencoder {
public:
  RCfileencoder(FILE* file) : RCencoder(), file(file), count(0), size(0) {}
  ~RCfileencoder() { flush(); }
  void putbyte(unsigned byte)
  {
    if (size == FPZIP_BLOCK_SIZE)
      flush();
    buffer[size++] = (unsigned char)byte;
  }
  void flush()
  {
    if (fwrite(buffer, 1, size, file) != size)
      error = true;
    else
      count += size;
    size = 0;
  }
  size_t bytes() const { return count; }
private:
  FILE* file;
  size_t count;
  size_t size;
  unsigned char buffer[FPZIP_BLOCK_SIZE];
};
#else
class RCfileencoder : public RCencoder {
public:
  RCfileencoder(FILE* file) : RCencoder(), file(file), count(0) {}
  void putbyte(unsigned byte)
  {
    if (fputc(byte, file) == EOF)
      error = true;
    else
      count++;
  }
  size_t bytes() const { return count; }
private:
  FILE* file;
  size_t count;
};
#endif

// memory writer for compressed data
class RCmemencoder : public RCencoder {
public:
  RCmemencoder(void* buffer, size_t size) : RCencoder(), ptr((unsigned char*)buffer), begin(ptr), end(ptr + size) {}
  void putbyte(unsigned byte)
  {
    if (ptr == end) {
      error = true;
      fpzip_errno = fpzipErrorBufferOverflow;
    }
    else
      *ptr++ = (unsigned char)byte;
  }
  size_t bytes() const { return ptr - begin; }
private:
  unsigned char* ptr;
  const unsigned char* const begin;
  const unsigned char* const end;
};

#endif
