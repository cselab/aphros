#include "util/filesystem.h"
static const char* me = "splitext";
int main(int argc, char** argv) {
  std::string a;
  std::string b;
  std::array<std::string, 2> ans;

  argv++;
  if (*argv == NULL) {
    fprintf(stderr, "%s: need an argument\n", me);
    return 1;
  }
  ans = util::SplitExt(std::string(argv[0]));
  printf("'%s' '%s'\n", ans[0].c_str(), ans[1].c_str());
}
