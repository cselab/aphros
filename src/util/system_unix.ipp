#define POSIX_C_SOURCE 200112L
#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <sys/stat.h>
#include <unistd.h>

char* SystemBaseName(char* path) {
  return basename(path);
}

char* SystemDirName(char* path) {
  return dirname(path);
}

int SystemMakeDir(char* path, int parent) {
  char cmd[2048];
  int rc;

  if (rc)
    ;
  else
    ; /* USED */
  if (parent)
    snprintf(cmd, sizeof cmd - 1, "mkdir -p -- '%s'", path);
  else
    snprintf(cmd, sizeof cmd - 1, "mkdir -- '%s'", path);
  rc = system(cmd);
  return SystemIsDir(path);
}

char* SystemRealPath(const char* path, char* resolved) {
  return realpath(path, resolved);
}

int SystemGetHostName(char* name, size_t len) {
  return gethostname(name, len);
}

int SystemHasHyperthreads(void) {
  // TODO
  return 0;
}

int SystemIsDir(char* path) {
  struct stat info;
  return (stat(path, &info) == 0) && (info.st_mode & S_IFDIR);
}

int SystemIsFile(char* path) {
  struct stat info;
  return (stat(path, &info) == 0) && (info.st_mode & (S_IFREG | S_IFLNK));
}

size_t SystemGetMem(void) {
  // TODO
  return 0;
}
