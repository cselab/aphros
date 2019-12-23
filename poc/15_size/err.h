#define ERR(x)                                              \
  do {                                                      \
    fprintf(stderr, "%s: %s:%d: ", me, __FILE__, __LINE__); \
    err_print x;                                            \
    fputs("\n", stderr);                                    \
    exit(2);                                                \
  } while (0)

#define MSG(x)                                      \
  do {                                              \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    err_print x;                                    \
    fputs("\n", stderr);                            \
  } while (0)

static int err_print(const char* fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  return 0;
}
