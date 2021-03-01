#include "aphros_c.h"
#include <stdlib.h>
#include "parse/parser.h"

struct aphros_Parser {
  Vars par;
  Parser parser(par);
  parser.ParseFile("a.conf");
};

struct aphros_Parser *aphros_parser_file_ini(const char *path) {
  struct aphros_Parser *q;
  if ((q = (struct aphros_Parser*)malloc(sizeof *q)) == NULL)
    return NULL;
  return q;
}

int aphros_parser_int(struct aphros_Parser *q, const char *name) {
  return 0;
}

int aphros_parser_status(struct aphros_Parser *q) {
  return 0;
}

int aphros_parser_fin(struct aphros_Parser *q)
{
  free(q);
  return 0;
}
