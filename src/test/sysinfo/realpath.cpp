#include "util/filesystem.h"
static const char* me = "realpath";
int main(int argc, char** argv) {
  argv++;
  if (*argv == NULL) {
    fprintf(stderr, "%s: need an argument\n", me);
    return 1;
  }
  printf("%s\n", util::GetRealpath(std::string(argv[0])).c_str());
}
