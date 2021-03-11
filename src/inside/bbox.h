// Created by Sergey Litvinov on 31.01.2021
// Copyright 2021 ETH Zurich

struct Bbox;
int bbox_ini(struct Bbox**);
int bbox_fin(struct Bbox*);
int bbox_update(struct Bbox*, int, const double*);
int bbox_inside(struct Bbox*, const double[3]);
double bbox_zhi(struct Bbox*);
int bbox_lo(struct Bbox*, const double**);
int bbox_hi(struct Bbox*, const double**);
