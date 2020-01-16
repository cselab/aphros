// Created by Sergey Litvinov on 24.10.2019
// Copyright 2019 ETH Zurich

enum { INFO_VERT, INFO_TRI };
struct Info {
  long nv, nt;
  char* name[999];
  int size[999];
} info_read(FILE*, *info);
