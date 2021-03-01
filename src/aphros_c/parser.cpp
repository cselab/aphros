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
  q->par = new Vars();
  q->parser = new Parser(*q->par);
  q->parser->ParseFile(path);
  q->status = 0;
  return q;
}

int aphros_parser_print_vars(struct aphros_Parser *q) {
  q->parser->PrintVars(std::cout);
  return 0;
}

int aphros_parser_int(struct aphros_Parser *q, const char *name) {
  return 0;
}

int aphros_parser_status(struct aphros_Parser *q) {
  return q->status;
}

int aphros_parser_fin(struct aphros_Parser *q)
{
  delete q->par;
  delete q->parser;
  free(q);
  return 0;
}
