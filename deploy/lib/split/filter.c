#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <table.h>
#include <vtk.h>

enum { N = 999 };
static char me[] = "filter";

#include "util.h"

#define USED(x) \
  if (x)        \
    ;           \
  else {        \
  }

static void usg() {
  fprintf(stderr, "%s -f field -p prefix [key ..] -- [vtk ..]\n", me);
  exit(1);
}

int main(int argc, char** argv) {
  char* Field;
  char output[N];
  char* path;
  char* Prefix;
  char* Volume;
  FILE* file;
  float* field;
  int cl;
  int* flag;
  int i;
  int nt;
  int tmp;
  struct Table* table;
  struct VTK* vtk;

  USED(argc);
  Prefix = Field = Volume = NULL;
  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
      case 'h':
        usg();
        break;
      case 'f':
        argv++;
        if ((Field = *argv) == NULL) {
          fprintf(stderr, "%s: -f needs an argument\n", me);
          exit(2);
        }
        break;
      case 'p':
        argv++;
        if ((Prefix = *argv) == NULL) {
          fprintf(stderr, "%s: -p needs an argument\n", me);
          exit(2);
        }
        break;
      case '-':
        argv++;
        goto end;
      default:
        fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
        exit(1);
    }
end:
  if (Field == NULL) {
    fprintf(stderr, "%s: field (-f) is not given\n", me);
    exit(1);
  }
  if (Prefix == NULL) {
    fprintf(stderr, "%s: prefix (-p) is not given\n", me);
    exit(1);
  }
  table = table_ini(999);
  for (;;) {
    if (*argv == NULL) {
      fprintf(stderr, "%s: missing '--' in arguments\n", me);
      exit(2);
    }
    if (util_eq(*argv, "--")) break;
    if (sscanf(*argv, "%d", &cl) != 1) {
      fprintf(stderr, "%s: not a number '%s'\n", me, *argv);
      exit(2);
    }
    table_put(table, cl, 0);
    argv++;
  }

  while ((path = *++argv) != NULL) {
    file = fopen(path, "r");
    if (file == NULL) {
      fprintf(stderr, "%s: cannot read '%s'\n", me, path);
      exit(2);
    }
    vtk = vtk_read(file);
    if (vtk == NULL) {
      fprintf(stderr, "%s: fail to parse vtk file '%s'\n", me, path);
      exit(2);
    }
    fclose(file);
    nt = vtk_nt(vtk);
    field = vtk_data(vtk, Field);
    if (field == NULL) {
      fprintf(stderr, "%s: no field '%s' in vtk file\n", me, Field);
      exit(2);
    }
    i = vtk_index(vtk, Field);
    if (vtk->type[i] != VTK_FLOAT) {
      fprintf(stderr, "%s: field '%s' is not float\n", me, Field);
      exit(2);
    }

    flag = malloc(nt * sizeof(*flag));
    if (flag == NULL) {
      fprintf(stderr, "%s: malloc failed (nt = %d)\n", me, nt);
      exit(2);
    }
    for (i = 0; i < nt; i++) {
      flag[i] = (table_get(table, (int)field[i], &tmp) == TABLE_EMPY);
    }
    vtk_remove_tri(vtk, flag);
    vtk_remove_orphan(vtk);

    util_name(Prefix, path, output);
    if ((file = fopen(output, "w")) == NULL) {
      fprintf(stderr, "%s: fail to write to '%s'\n", me, output);
      exit(2);
    }
    vtk_write(vtk, file);
    vtk_fin(vtk);
    fclose(file);
    free(flag);
  }
  table_fin(table);
}
