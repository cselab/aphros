#include <stdio.h>
#include <stdlib.h>
#include <quadtree.h>

static const char *me = "fig";
enum {SIZE = 999, MAX_NODE = 999};

static void
usg(void)
{
  fprintf(stderr, "%s points query\n", me);
  exit(1);
}

static void found(struct Node *node);

int
main(int argc, char **argv)
{
  FILE *f;
  int i;
  int n;
  struct Node *node[MAX_NODE];
  double x;
  double y;
  double L;
  double R;
  double B;
  double T;
  char line[SIZE];

  while (*++argv != NULL && argv[0][0] == '-')
    switch (argv[0][1]) {
    case 'h':
      usg();
      break;
    default:
      fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
      exit(2);
    }
  if (argv[0] == NULL || argv[1] == NULL) {
    fprintf(stderr, "%s: need two arguments\n", me);
    exit(1);
  }

  if ((f = fopen(argv[0], "r")) == NULL) {
    fprintf(stderr, "%s: fail to open '%s'\n", me, argv[0]);
    exit(1);
  }
  i = 0;
  while (fgets(line, SIZE, f) != NULL) {
    if (sscanf(line, "%lf %lf\n", &x, &y) != 2) {
      fprintf(stderr, "%s: fail to parse '%s'\n", me, line);
      exit(1);
    }
    node[i] = node_ini(x, y, i);
    if (i > 0)
      insert(node[i], node[0]);
    i++;
  }
  n = i;
  fclose(f);

  if ((f = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "%s: fail to open '%s' in '%s'\n", me, argv[1], argv[0]);
    exit(1);
  }
  while (fgets(line, SIZE, f) != NULL) {
    if (sscanf(line, "%lf %lf %lf %lf\n", &L, &R, &B, &T) != 4) {
      fprintf(stderr, "%s: fail to parse '%s' in '%s'\n", me, line, argv[1]);
      exit(1);
    }
    regionsearch(node[0], L, R, B, T, found);
  }
  fclose(f);

  for (i = 0; i < n; i++)
    free(node[i]);
}

static void
found(struct Node *node)
{
  printf("%g %g %d\n", node->x, node->y, node->data);
}
