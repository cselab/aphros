int err_print(const char* fmt, ...);
void err_exit(int);

#define ERR(x)                                              \
  do {                                                      \
    fprintf(stderr, "%s: %s:%d: ", me, __FILE__, __LINE__); \
    err_print x;                                            \
    fputs("\n", stderr);                                    \
    err_exit(2);                                            \
  } while (0)

#define MSG(x)                                      \
  do {                                              \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    err_print x;                                    \
    fputs("\n", stderr);                            \
  } while (0)
