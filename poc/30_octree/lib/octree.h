struct Node {
  double x;
  double y;
  double z;
  int data;
  struct Node *elm[8];
};
int insert(struct Node *, struct Node *);
int regionsearch(struct Node *, const double lo[3], const double hi[3], void found(struct Node*));
struct Node *node_ini(double x, double y, double z, int data);
