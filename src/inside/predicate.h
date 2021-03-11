// Created by Sergey Litvinov on 31.01.2021
// Copyright 2021 ETH Zurich

int predicate_ini(void);
/* does a ray "de" cross a triangl "abc"? */
int predicate_ray(
    const double d[3], const double e[3], const double a[3], const double b[3],
    const double c[3]);
