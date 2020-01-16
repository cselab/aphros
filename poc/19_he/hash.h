// Created by Sergey Litvinov on 18.11.2019
// Copyright 2019 ETH Zurich

#define T HeHash

typedef struct T T;

int he_hash_ini(int n, T**);

/* hash[i, j] = v */
int he_hash_set(T*, int i, int j, int v);

/* ans = hash[i, j], ans = -1 if unset  */
int he_hash_get(T*, int i, int j);
int he_hash_fin(T*);

#undef T
