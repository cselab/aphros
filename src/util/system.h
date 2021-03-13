// Created by Sergey Litvinov on 08.03.2021
// Copyright 2021 ETH Zurich

#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
int SystemBaseName(const char*, char*);
int SystemDirName(const char*, char*);
char* SystemRealPath(const char*, char* resolved);
int SystemGetHostName(char*, size_t size);
int SystemHasHyperthreads(void);
int SystemIsDir(const char*);
int SystemIsFile(const char*);
int SystemJoin(const char*, const char*, char*);
int SystemMakeDir(const char*, int parent);
int SystemSplitExt(const char*, char*, char*);
size_t SystemGetMem(void);

#ifdef __cplusplus
}
#endif
