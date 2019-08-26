#include <stdio.h>
#include <assert.h>

enum {EMPTY = -1};
enum {FREE = -1};
#define SIZE(a) (int)(sizeof(a)/sizeof(*(a)))
static int root[10*10*10];
static int a[10*10*10];
static int lbl[10*10*10];

void
make_(int v)
{
    root[v] = v;
}

int
find_(int v)
{
    if (v == root[v])
	return v;
    return root[v] = find_(root[v]);
}

void
union_(int a, int b)
{
    a = find_(a);
    b = find_(b);
    if (a != b)
	root[b] = a;
}

int
color(int n, /**/ int *pcnt, int *a)
{
    enum {I, J, K};
    int i, j, k, m, N;
    int u, v, w;
    int idx, jdx;
    int x, cnt;
    int d[][3] = {
	{0, 0, 1},
	{0, 1, 0},
	{0, 1, 1},
	{1, 0, 0},
	{1, 0, 1},
	{1, 1, 0},
	{1, 1, 1},
    };
    N = n*n*n;
    for (i = 0; i < N; i++) {
	make_(i);
	lbl[i] = FREE;
    }
    for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	    for (k = 0; k < n; k++)
		for (m = 0; m < SIZE(d); m++) {
		    if ((u = i + d[m][I]) >= n) continue;
		    if ((v = j + d[m][J]) >= n) continue;
		    if ((w = k + d[m][K]) >= n) continue;
		    idx = i*n*n + j*n + k;
		    jdx = u*n*n + v*n + w;
		    assert(idx < N);
		    assert(jdx < N);
		    if (a[idx] == EMPTY) continue;
		    if (a[jdx] == EMPTY) continue;
		    union_(idx, jdx);
		}
    cnt = 0;
    for (i = 0; i < N; i++) {
	x = find_(i);
	assert(x < N);
	if (lbl[x] == FREE)
	    lbl[x] = cnt++;
	a[i] = lbl[x];
    }
    *pcnt = cnt;
    return 0;
}

int
main()
{
    int n, cnt;
    n = 3;
    color(n, &cnt, a);
    fprintf(stderr, "cnt: %d\n", cnt);
}
