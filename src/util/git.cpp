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

const char* GetLogo() {
  return
      R"EOF(/*************************************************\
|    O   ____ O  O  O  _  O  O   O   O    O    O  |
| O   O /    \  o o o | |   o  O   O   O   o o  o |
|  O   /  /\ |   o o  | |  o  o   o  o   o    o   |
|   o /  / | |  ____  | |__   _____ _____  ______ |
| o  /  /  | | / __ \ |  _ \ |  __//  _  |/  ___/ |
| o /  /___| || |  | || | | || |   | | | ||_|___  |
|  /  /    | || |__| || | | || |   | |_| | ___| | |
| /__/  O  |_||  ___/ |_| |_||_|   |____/ /_____/ |
|  o   o  o   | |    o           o               o|
|   o    O    |_| o    O  o  O    o   o   O   O   |
|                                                 |
\*************************************************/
)EOF";
}
