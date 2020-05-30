// Created by Petr Karnakov on 16.06.2019
// Copyright 2019 ETH Zurich

#include "git.h"

extern const char* kGitRev;
extern const char* kGitMsg;
extern const char* kGitDiff;

const char* GetGitRev() {
  return kGitRev;
}

const char* GetGitMsg() {
  return kGitMsg;
}

const char* GetGitDiff() {
  return kGitDiff;
}
