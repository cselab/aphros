struct Node {
  double x;
  double y;
  int data;
  struct Node *elm[4];
};
int insert(struct Node *, struct Node *);
int regionsearch(struct Node *, double L, double R, double B, double T, void found(struct Node*));
struct Node *node_ini(double x, double y, int data);
