#include <stdio.h>
#include <stdlib.h>
#include <octree.h>

static const char *me = "main";
enum {SIZE = 999};
enum {X, Y, Z};

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
  int cap;
  struct Node **node;
  double x;
  double y;
  double z;
  double lo[3];
  double hi[3];
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

  cap = 1;
  node = malloc(cap * sizeof(*node));
  i = 0;
  while (fgets(line, SIZE, f) != NULL) {
    if (sscanf(line, "%lf %lf %lf\n", &x, &y, &z) != 3) {
      fprintf(stderr, "%s: fail to parse '%s' in '%s'\n", me, line, argv[0]);
      exit(1);
    }
    if (i >= cap) {
      cap *= 2;
      if ((node = realloc(node, cap * sizeof(*node))) == NULL) {
	fprintf(stderr, "%s: realloc failed\n", me);
	exit(1);
      }
    }
    node[i] = node_ini(x, y, z, i);
    if (i > 0)
      insert(node[i], node[0]);
    i++;
  }
  n = i;
  fclose(f);

  if ((f = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "%s: fail to open '%s'\n", me, argv[1]);
    exit(1);
  }
  while (fgets(line, SIZE, f) != NULL) {
    if (sscanf(line, "%lf %lf %lf %lf %lf %lf\n", &lo[X], &hi[X], &lo[Y], &hi[Y], &lo[Z], &hi[Z]) != 6) {
      fprintf(stderr, "%s: fail to parse '%s' in '%s'\n", me, line, argv[1]);
      exit(1);
    }
    regionsearch(node[0], lo, hi, found);
  }
  fclose(f);

  for (i = 0; i < n; i++)
    free(node[i]);
  free(node);
}

static void
found(struct Node *node)
{
  printf("%g %g %g %d\n", node->x, node->y, node->z, node->data);
}
