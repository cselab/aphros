#ifndef FPZIP_READ_H
#define FPZIP_READ_H

#define subsize(T, n) (CHAR_BIT * sizeof(T) * (n) / 32)

// file reader for compressed data
#if FPZIP_BLOCK_SIZE > 1
class RCfiledecoder : public RCdecoder {
public:
  RCfiledecoder(FILE* file) : RCdecoder(), file(file), count(0), index(0), size(0) {}
  unsigned getbyte()
  {
    if (index == size) {
      size = fread(buffer, 1, FPZIP_BLOCK_SIZE, file);
      if (!size) {
        size = 1;
        error = true;
      }
      else
        count += size;
      index = 0;
    }
    return buffer[index++];
  }
  size_t bytes() const { return count; }
private:
  FILE* file;
  size_t count;
  size_t index;
  size_t size;
  unsigned char buffer[FPZIP_BLOCK_SIZE];
};
#else
class RCfiledecoder : public RCdecoder {
public:
  RCfiledecoder(FILE* file) : RCdecoder(), file(file), count(0) {}
  unsigned getbyte()
  {
    int byte = fgetc(file);
    if (byte == EOF)
      error = true;
    else
      count++;
    return byte;
  }
  size_t bytes() const { return count; }
private:
  FILE* file;
  size_t count;
};
#endif

class RCmemdecoder : public RCdecoder {
public:
  RCmemdecoder(const void* buffer) : RCdecoder(), ptr((const unsigned char*)buffer), begin(ptr) {}
  unsigned getbyte() { return *ptr++; }
  size_t bytes() const { return ptr - begin; }
private:
  const unsigned char* ptr;
  const unsigned char* const begin;
};

#endif
