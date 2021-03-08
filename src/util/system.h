#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
char *SystemBaseName(char*);
char *SystemDirName(char*);
char *SystemMakeDir(char*, int parent);
char *SystemRealPath(const char *, char *resolved);
int SystemGetHostName(char *, size_t size);
int SystemHasHyperthreads(void);
int SystemIsDir(char*);
int SystemIsFile(char*);
size_t SystemGetMem(void);

#ifdef __cplusplus
}
#endif
