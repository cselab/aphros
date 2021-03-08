#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
int SystemBaseName(char*, char*);
int SystemDirName(char*, char*);
char* SystemRealPath(const char*, char* resolved);
int SystemGetHostName(char*, size_t size);
int SystemHasHyperthreads(void);
int SystemIsDir(char*);
int SystemIsFile(char*);
int SystemJoin(const char *, const char *, char *);
int SystemMakeDir(char*, int parent);
int SystemSplitExt(const char *, char *, char *);
size_t SystemGetMem(void);

#ifdef __cplusplus
}
#endif
