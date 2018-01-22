#include "fpzip.h"

fpzipError fpzip_errno;

const char* fpzip_errstr[] = {
  "success",
  "cannot read stream",
  "cannot write stream",
  "not an fpz stream",
  "fpz format version not supported",
  "precision not supported",
  "memory buffer overflow",
};
