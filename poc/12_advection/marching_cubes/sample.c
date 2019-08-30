#include <stdlib.h>
#include <tgmath.h>
#include "sample.h"
struct Vec {
    double x;
    double y;
    double z;
};
struct Sample {
    int type;
    struct Vec point[3];
    double t;
    double level;
};
static double Sample1(struct Sample *, double, double, double);
static double Sample2(struct Sample *, double, double, double);
static double Sample3(struct Sample *, double, double, double);
enum { SAMPLE1, SAMPLE2, SAMPLE3, EMPTY };

static double (*sample[3]) (struct Sample *, double, double, double) = {
Sample1, Sample2, Sample3};

static double
Sample1(struct Sample *q, double x, double y, double z)
{
    double res = 0.0;
    double dx, dy, dz;
    struct Vec *point;

    point = q->point;
    dx = x - point[0].x;
    dy = y - point[0].y;
    dz = z - point[0].z;
    res += 0.5 / (dx * dx + dy * dy + dz * dz);
    dx = x - point[1].x;
    dy = y - point[1].y;
    dz = z - point[1].z;
    res += 1.0 / (dx * dx + dy * dy + dz * dz);
    dx = x - point[2].x;
    dy = y - point[2].y;
    dz = z - point[2].z;
    res += 1.5 / (dx * dx + dy * dy + dz * dz);
    return res;
}

static double
Sample2(struct Sample *q, double x, double y, double z)
{
    double res = 0.0;
    double dx, dy, dz;
    struct Vec *point;

    point = q->point;
    dx = x - point[0].x;
    dy = y - point[0].y;
    res += 0.5 / (dx * dx + dy * dy);
    dx = x - point[1].x;
    dz = z - point[1].z;
    res += 0.75 / (dx * dx + dz * dz);
    dy = y - point[2].y;
    dz = z - point[2].z;
    res += 1.0 / (dy * dy + dz * dz);
    return res;
}

static double
Sample3(struct Sample *q, double x, double y, double z)
{
    double res;
    double Height;
    double t;

    t = q->t;
    Height =
	20.0 * (t + sqrt((0.5 - x) * (0.5 - x) + (0.5 - y) * (0.5 - y)));
    Height = 1.5 + 0.1 * (sinf(Height) + cosf(Height));
    res = (Height - z) * 50.0;
    return res;
}

struct Sample *
sample_ini(void)
{
    struct Sample *q;

    q = malloc(sizeof(*q));
    q->type = SAMPLE1;
    q->level = 48.0;
    return q;
}

void
sample_fin(struct Sample *q)
{
    free(q);
}

double
sample_f(struct Sample *q, double x, double y, double z)
{
    double f, level;

    f = sample[q->type] (q, x, y, z);
    level = q->level;
    return f - level;
}

int
sample_time(struct Sample *q, double t)
{
    double Offset;
    struct Vec *point;
    int i;

    point = q->point;
    for (i = 0; i < 3; i++) {
	point[i].x = 0.5;
	point[i].y = 0.5;
	point[i].z = 0.5;
    }
    Offset = 1.0 + sin(t);
    point[0].x *= Offset;
    point[1].y *= Offset;
    point[2].z *= Offset;
    q->t = t;
    return 0;
}

int
sample_next(struct Sample *q)
{
    q->type++;
    if (q->type == EMPTY)
	q->type = SAMPLE1;
    return 0;
}

int
sample_inc(struct Sample *q)
{
    if (q->level < 1000)
	q->level *= 1.1;
    return 0;
}

int
sample_dec(struct Sample *q)
{
    if (q->level > 1)
	q->level /= 1.1;
    return 0;
}
