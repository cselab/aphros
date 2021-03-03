#include "parse/parser.h"
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "aphros_c.h"

struct aphros_Parser {
  Vars* par;
  Parser* parser;
  char* string;
  int status;
};

struct aphros_Parser* aphros_ParserFileIni(const char* path) {
  struct aphros_Parser* q;
  if ((q = (struct aphros_Parser*)malloc(sizeof *q)) == NULL) return NULL;
  q->string = NULL;
  q->status = 0;
  q->par = new Vars();
  q->parser = new Parser(*q->par);
  try {
    q->parser->ParseFile(path);
  } catch (const std::runtime_error& e) {
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
    if (q->string == NULL)
      free(q->string);
    if ((q->string = strdup(s.c_str())) == NULL) {
      fprintf(stderr, "%s:%d: strdup failed\n", __FILE__, __LINE__);
      q->status = 2;
      return NULL;
    }
    return q->string;
  } catch (const std::runtime_error& e) {
    q->status = 1;
    return NULL;
  }
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
  free(q);
  return 0;
}
