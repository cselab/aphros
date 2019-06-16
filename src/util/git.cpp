#include "git.h"

const char* GetGitRev() {
  static const char s[] = GITREV;
  return s;
}

const char* GetGitMsg() {
  static const char s[] = GITMSG;
  return s;
}

const char* GetGitDiff() {
  static const char s[] = GITDIFF;
  return s;
}
