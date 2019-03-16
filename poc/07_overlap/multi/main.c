#include <stdio.h>
#include <stdlib.h>

typedef double Scal;

static int M = 200;
static int N = 200;

struct Circle { Scal x, y, r; };
static int f(Scal x, Scal y, void *p) {
    Scal r;
    struct Circle *c = p;
    x -= c->x;
    y -= c->y;
    r = c->r;
    return x*x + y*y < r*r;
}

static Scal g(Scal x, Scal y, Scal u, Scal v, void *P) {
    int i, j;
    Scal p, q, dx, dy;
    long s, a, b;

    s = 0;
    dx = (u - x)/M;
    dy = (v - y)/N;
    for (i = 0; i < M + 1; i++) {
        a = (i == 0 || i == M) ? 1 : 2;
        p = x + i*dx;
        for (j = 0; j < N + 1; j++) {
            b = (j == 0 || j == N) ? 1 : 2;
            q = y + j*dy;
            s += a * b * f(p, q, P);
        }
    }
    return s*dx*dy/4;
}

int main(void) {
    Scal x, y;
    struct Circle c;
    c.x = 0;
    c.y = 0;
    c.r = 0.25;

    printf("%e\n", g(-0.5, -0.5, 0.5, 0.5, &c));
}
