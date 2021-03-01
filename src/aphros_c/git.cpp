#include "util/git.h"

const char* aphros_get_git_rev(void) {
  return GetGitRev();
}

const char* aphros_get_logo(void) {
  return GetLogo();
}
