#define WARN(x)                            \
  do {                                     \
    dprint("%s:%d: ", __FILE__, __LINE__); \
    dprint x;                              \
    dprint("\n");                          \
  } while (0)
static int dprint(const char* fmt, ...) {
  int r, rank;
  va_list ap;
  static FILE* f = NULL;
  static char n[4092];
  if (f == NULL) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    sprintf(n, "log-%d", rank);
    f = fopen(n, "w");
  } else {
    f = fopen(n, "a");
  }
  va_start(ap, fmt);
  r = vfprintf(f, fmt, ap);
  va_end(ap);
  fclose(f);
  return r;
}
