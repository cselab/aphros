char* SystemBaseName(char* path) {
  сhar drive[_MAX_DRIVE];
  char dir[_MAX_DIR];
  char fname[_MAX_FNAME];
  char ext[_MAX_EXT];
  _splitpath(path, drive, dir, fname, ext);
  // TODO
  return NULL;
}

char* SystemDirName(char* path) {
  // TODO
  return NULL;
}

int SystemMakeDir(char* path, int parent) {
  // TODO
  return 0;
}

char* SystemRealPath(const char* path, char* resolved) {
  // TODO
  return NULL;
}

int SystemGetHostName(char* name, size_t len) {
  // TODO
  return 0;
}

int SystemHasHyperthreads(void) {
  // TODO
  return 0;
}

int SystemIsDir(char* path) {
  // TODO
  return 0;
}

int SystemIsFile(char* path) {
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
