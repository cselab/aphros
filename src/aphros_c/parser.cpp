#include "aphros_c.h"
#include <iostream>
#include "parse/parser.h"

struct aphros_Parser {
  Vars *par;
  Parser *parser;
  int status;
};

struct aphros_Parser *aphros_parser_file_ini(const char *path) {
  struct aphros_Parser *q;
  if ((q = (struct aphros_Parser*)malloc(sizeof *q)) == NULL)
    return NULL;
  q->status = 0;
  q->par = new Vars();
  q->parser = new Parser(*q->par);
  try {
    q->parser->ParseFile(path);
  } catch (const std::runtime_error &e) {
    fprintf(stderr, "%s\n", e.what());
    return NULL;
  }
  return q;
}

int aphros_parser_print_vars(struct aphros_Parser *q) {
  q->parser->PrintVars(std::cout);
  return 0;
}

int aphros_parser_int(struct aphros_Parser *q, const char *name) {
  try {
    return q->par->Int[name];
  } catch (const std::runtime_error &e) {
    fprintf(stderr, "%s\n", e.what());
    q->status = 1;
    return 0;
  }
}

int aphros_parser_status(struct aphros_Parser *q) {
  int status;
  status = q->status;
  q->status = 0;
  return status;
}

int aphros_parser_fin(struct aphros_Parser *q)
{
  delete q->par;
  delete q->parser;
  free(q);
  return 0;
}
