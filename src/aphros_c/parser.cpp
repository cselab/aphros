// Created by Sergey Litvinov on 01.03.2021
// Copyright 2021 ETH Zurich

#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "aphros_c.h"

#include "parse/parser.h"
#include "util/logger.h"

void ErrorHandler(int, const char* str) {
  fputs(str, stderr);
  fputs("\n", stderr);
}

struct aphros_Parser {
  Vars* par;
  Parser* parser;
  char* string;
  double* vect;
  int status;
};

struct aphros_Parser* aphros_ParserFileIni(const char* path) {
  struct aphros_Parser* q;

  aphros_SetErrorHandler(ErrorHandler);
  if ((q = (struct aphros_Parser*)malloc(sizeof *q)) == NULL) return NULL;
  q->string = NULL;
  q->vect = NULL;
  q->status = 0;
  try {
    q->par = new Vars();
  } catch (const std::bad_alloc& e) {
    return NULL;
  }
  try {
    q->parser = new Parser(*q->par);
  } catch (const std::bad_alloc& e) {
    delete q->par;
    return NULL;
  }
  try {
    q->parser->ParseFile(path);
  } catch (const std::runtime_error& e) {
    delete q->parser;
    delete q->par;
    return NULL;
  }
  return q;
}

int aphros_ParserPrintVars(struct aphros_Parser* q) {
  q->parser->PrintVars(std::cout);
  return 0;
}

int aphros_ParserGetInt(struct aphros_Parser* q, const char* name) {
  try {
    return q->par->Int[name];
  } catch (const std::runtime_error& e) {
    q->status = 1;
    return INT_MAX;
  }
}

double aphros_ParserGetDouble(struct aphros_Parser* q, const char* name) {
  try {
    return q->par->Double[name];
  } catch (const std::runtime_error& e) {
    q->status = 1;
    return DBL_MAX;
  }
}

char* aphros_ParserGetString(struct aphros_Parser* q, const char* name) {
  std::string s;

  try {
    s = q->par->String[name];
  } catch (const std::runtime_error& e) {
    q->status = 1;
    return NULL;
  }
  free(q->string);
  if ((q->string = strdup(s.c_str())) == NULL) {
    fprintf(stderr, "%s:%d: strdup failed\n", __FILE__, __LINE__);
    q->status = 2;
    return NULL;
  }
  return q->string;
}

double* aphros_ParserGetVect(
    struct aphros_Parser* q, const char* name, /**/ int* size) {
  int i;
  std::vector<double> s;

  try {
    s = q->par->Vect[name];
  } catch (const std::runtime_error& e) {
    q->status = 1;
    return NULL;
  }
  free(q->vect);
  if ((q->vect = (double*)malloc(s.size() * sizeof *q->vect)) == NULL) {
    fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
    q->status = 2;
    return NULL;
  }
  for (i = 0; i < (int)s.size(); i++)
    q->vect[i] = s.at(i);
  *size = s.size();
  return q->vect;
}

int aphros_ParserStatus(struct aphros_Parser* q) {
  int status;

  status = q->status;
  q->status = 0;
  return status;
}

int aphros_ParserFin(struct aphros_Parser* q) {
  delete q->par;
  delete q->parser;
  free(q->string);
  free(q->vect);
  free(q);
  return 0;
}
