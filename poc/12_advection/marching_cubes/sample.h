// Created by Sergey Litvinov on 29.08.2019
// Copyright 2019 ETH Zurich

struct Sample;
struct Sample* sample_ini(void);
void sample_fin(struct Sample*);
double sample_f(struct Sample*, double[3]);
int sample_time(struct Sample*, double);
int sample_next(struct Sample*);
int sample_inc(struct Sample*);
int sample_dec(struct Sample*);
struct Sample;
struct Sample* sample_ini(void);
void sample_fin(struct Sample*);
double sample_f(struct Sample*, double[3]);
int sample_time(struct Sample*, double);
int sample_next(struct Sample*);
int sample_inc(struct Sample*);
int sample_dec(struct Sample*);
