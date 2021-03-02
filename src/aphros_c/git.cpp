#include "util/git.h"
#include "aphros_c.h"

const char* aphros_GetGitRev(void) {
  return GetGitRev();
}

const char* aphros_GetLogo(void) {
  return GetLogo();
}
