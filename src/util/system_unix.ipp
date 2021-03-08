#define POSIX_C_SOURCE 200112L
#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

int SystemBaseName(char* path, char *name) {
  char *ans;
  if ((ans = basename(path)) == NULL)
    return 1;
  strcpy(name, ans);
  return 0;
}

int SystemDirName(char *path, char *dir) {
  strcpy(dir, dirname(path));
  return 0;
}

int SystemMakeDir(char* path, int parent) {
  char cmd[2048];
  int rc;

  (void)rc;
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

int SystemJoin(const char *a, const char *b, char *c) {
  return 0;
}

size_t SystemGetMem(void) {
  // TODO
  return 0;
}

int SystemSplitExt(const char *a, char *base, char *ext) {
  return 0;
}
