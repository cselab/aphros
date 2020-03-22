struct Bbox;
int bbox_ini(struct Bbox**);
int bbox_fin(struct Bbox*);
int bbox_update(struct Bbox*, int, const double *);
int bbox_inside(struct Bbox*, const double[3]);
double bbox_zhi(struct Bbox*);
