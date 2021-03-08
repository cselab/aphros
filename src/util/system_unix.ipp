#define POSIX_C_SOURCE 200112L
#define _XOPEN_SOURCE 500
#include <stdlib.h>
#include <libgen.h>
#include <unistd.h>

char *SystemBaseName(char *path) {
     return basename(path);
}

char *SystemDirName(char *path) {
     return dirname(path);
}

char *SystemMakeDir(char *path, int parent) {
  // TODO
  return NULL;
}

char *SystemRealPath(const char *path, char *resolved) {
     return realpath(path, resolved);
}

int SystemGetHostName(char *name, size_t len) {
  return gethostname(name, len);
}

int SystemHasHyperthreads(void) {
  // TODO
  return 0;
}

int SystemIsDir(char *path) {
  // TODO
  return 0;
}

int SystemIsFile(char *path) {
  // TODO
  return 0;
}

size_t SystemGetMem(void) {
  // TODO
  return 0;
}
