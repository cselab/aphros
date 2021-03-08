#ifdef __cplusplus
extern "C" {
#endif

char *SystemRealPath(const char *, char *resolved);
char *SystemDirName(char*);
char *SystemBaseName(char*);
char *SystemMakeDir(char*, int parent);
int SystemIsFile(char*);
int SystemIsDir(char*);
int SystemGetHostName(char *, size_t size);
size_t SystemGetMem(void);
int SystemHasHyperthreads(void);

#ifdef __cplusplus
}
#endif
