// Created by Sergey Litvinov on 08.03.2021
// Copyright 2021 ETH Zurich

#include <stdlib.h>
int SystemBaseName(const char* path, char *fname) {
  char drive[_MAX_DRIVE];
  char dir[_MAX_DIR];
  char ext[_MAX_EXT];
  _splitpath(path, drive, dir, fname, ext);
  return 0;
}

int SystemDirName(const char* path, char* dir) {
  // TODO
  char drive[_MAX_DRIVE];
  char fname[_MAX_FNAME];
  char ext[_MAX_EXT];
  _splitpath(path, drive, dir, fname, ext);
  return 0;
}

int SystemMakeDir(const char* path, int parent) {
  // TODO
  char cmd[2048];
  int rc;

  (void)rc;
  strcpy(cmd, "mkdir \"");
  strcat(cmd, path);
  strcat(cmd, "\"");
  rc = system(cmd);
  return SystemIsDir(path);
}

char* SystemRealPath(const char* path, char* resolved) {
  // TODO
  strcpy(resolved, path);
  return resolved;
}

int SystemGetHostName(char* name, size_t len) {
  //return gethostname(name, len);
  name[0] = 'w'; name[1] = 'i'; name[2] = 'n'; name[3] = '\0';
  return 0;
}

int SystemHasHyperthreads(void) {
  // TODO
  return 0;
}

int SystemIsDir(const char* path) {
  // TODO
  return 0;
}

int SystemIsFile(const char* path) {
  // TODO
  return 0;
}

int SystemJoin(const char *a, const char *b, char *c) {
  // TODO
  return 0;
}

size_t SystemGetMem(void) {
  // TODO
  return 0;
}

int SystemSplitExt(const char *a, char *base, char *ext) {
  // TODO
  return 0;
}
