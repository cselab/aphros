#include <stdio.h>
#include <stdlib.h>

typedef double Scal;

struct Circle { Scal x, y, r; };
static int f(Scal x, Scal y, void *p) {
    Scal r;
    struct Circle *c = p;
    x -= c->x;
    y -= c->y;
    r = c->r;
    fprintf(stderr, ": %g %g %g\n", x, y, r);
    return x*x + y*y < r*r;
}

enum {IN, OUT, MIXED};
static int g(Scal x, Scal y, Scal u, Scal v, void *p) {
    int a, b, c, d;
    a = f(x, y, p);
    b = f(x, v, p);
    c = f(u, y, p);
    d = f(u, v, p);
    if (a && b && c && d) {
        fprintf(stderr, "in\n");
        return IN;
    } else if (!a && !b && !c && !d) {
        fprintf(stderr, "out\n");
        return OUT;
    } else {
        fprintf(stderr, "mixed\n");
        return MIXED;
    }
}

static Scal h(int l, Scal x, Scal y, Scal u, Scal v, void *P) {
#   define hh(x, y, u, v) (h(l, (x), (y), (u), (v), P))
    int s;
    Scal p, q;

    if (x >= u || y >= v) {
        fprintf(stderr, "error: %d %g %g %g %g\n", l, x, y, u, v);
        exit(2);
    }
    if (l <= 0)
        return 0;

    l--;
    s = g(x, y, u, v, P);
    switch (s) {
    case IN: return (u - x) * (v - y);
    case OUT: return 0;
    case MIXED:
        p = (x + y)/2;
        q = (u + v)/2;
        return \
            hh(x, y, p, q) +
            hh(p, y, u, q) +
            hh(x, q, p, v) +
            hh(p, q, u, v);
    }
#   undef hh
}

int main(void) {
    Scal x, y;
    struct Circle c;
    c.x = 0.5;
    c.y = 0;
    c.r = 0.4;

    printf("%g\n", h(10, -0.5, -0.5, 0.5, 0.5, &c));
}
