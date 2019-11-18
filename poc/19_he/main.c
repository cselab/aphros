#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "err.h"
#include "memory.h"
#include "hash.h"
#include "he.h"

static char me[] = "he";

#define T HeRead
enum {SIZE = 999};
enum {MAGIC = 42};
enum {HE_OK = 0};

struct T {
    int nv, nt, ne, nh;
    int *nxt, *flp;
    int *ver, *tri, *edg;
    int *hdg_ver, *hdg_edg, *hdg_tri;
    int magic;
};

static char *util_fgets(char *, FILE *);
static int util_eq(const char *, const char *);

enum { OK, FAIL };
static int
valid(int n, int *a, int *b)
{
    int i, i0;

    for (i = 0; i < n; i++) {
        i0 = b[a[i]];
        if (i != i0) {
            MSG(("i=%d a[i]=%d i0=%d ", i, a[i], i0));
            return FAIL;
        }
    }
    return OK;
}

static int
alloc_hdg(int nh, T * q)
{
    int i;

    MALLOC(nh, &q->nxt);
    MALLOC(nh, &q->flp);
    MALLOC(nh, &q->ver);
    MALLOC(nh, &q->tri);
    MALLOC(nh, &q->edg);
    for (i = 0; i < nh; i++)
        q->edg[i] = -1;
    return HE_OK;
}

static int
alloc_ver(int nv, T * q)
{
    MALLOC(nv, &q->hdg_ver);
    return HE_OK;
}

static int
alloc_edg(int ne, T * q)
{
    MALLOC(ne, &q->hdg_edg);
    return HE_OK;
}

static int
alloc_tri(int nt, T * q)
{
    MALLOC(nt, &q->hdg_tri);
    return HE_OK;
}

static int
alloc(int nv, int ne, int nt, int nh, T * q)
{
    alloc_hdg(nh, q);
    alloc_ver(nv, q);
    alloc_edg(ne, q);
    alloc_tri(nt, q);
    return HE_OK;
}

static int
afree(T * q)
{
    FREE(q->nxt);
    FREE(q->flp);
    FREE(q->ver);
    FREE(q->tri);
    FREE(q->edg);
    FREE(q->hdg_ver);
    FREE(q->hdg_edg);
    FREE(q->hdg_tri);
    return HE_OK;
}

static void
setup_edg(T * q, HeHash * hdg, int i, int j, int *pe)
{
    int h, f, e;

    h = he_hash_get(hdg, i, j);
    if (h == -1)
        return;
    e = *pe;

    q->edg[h] = e;
    if ((f = he_hash_get(hdg, j, i)) != -1)
        q->edg[f] = e;

    he_hash_set(hdg, i, j, -1); /* counted it */
    he_hash_set(hdg, j, i, -1);

    e++;
    *pe = e;
}

static void
setup_hdg(T * q, HeHash * hdg, int h, int t, int i, int j)
{
    he_hash_set(hdg, i, j, h);
    q->hdg_tri[t] = q->hdg_ver[i] = h;
    q->ver[h] = i;
    q->tri[h] = t;
}

static void
setup_flip(T * q, HeHash * hdg, int i, int j)
{
    int h, f;

    h = he_hash_get(hdg, i, j);
    f = he_hash_get(hdg, j, i);
    q->flp[h] = f;
}

int
he_read_tri_ini(int nv, int nt, int *tri0, T ** pq)
{
    int ne, nh;
    int e, t, i, j, k;
    int h, hi, hj, hk;
    int *tri;
    T *q;
    HeHash *hdg;                /* will maps [i, j] to half-edg */

    MALLOC(1, &q);

    alloc_ver(nv, q);
    alloc_tri(nt, q);
    alloc_hdg(nh = 3 * nt, q);
    he_hash_ini(nv, &hdg);
    for (tri = tri0, t = h = 0; t < nt; t++) {
        i = *tri++;
        j = *tri++;
        k = *tri++;
        assert(i < nv);
        assert(j < nv);
        assert(k < nv);
        setup_hdg(q, hdg, hi = h++, t, i, j);
        setup_hdg(q, hdg, hj = h++, t, j, k);
        setup_hdg(q, hdg, hk = h++, t, k, i);
        q->nxt[hi] = hj;
        q->nxt[hj] = hk;
        q->nxt[hk] = hi;
    }

    for (tri = tri0, t = h = 0; t < nt; t++) {
        i = *tri++;
        j = *tri++;
        k = *tri++;
        setup_flip(q, hdg, i, j);
        setup_flip(q, hdg, j, k);
        setup_flip(q, hdg, k, i);
    }

    for (tri = tri0, ne = t = h = 0; t < nt; t++) {
        i = *tri++;
        j = *tri++;
        k = *tri++;
        setup_edg(q, hdg, i, j, &ne);
        setup_edg(q, hdg, j, k, &ne);
        setup_edg(q, hdg, k, i, &ne);
    }
    he_hash_fin(hdg);

    alloc_edg(ne, q);
    for (h = 0; h < nh; h++) {
        e = q->edg[h];
        if (e != -1)
            q->hdg_edg[e] = h;
    }

    q->nv = nv;
    q->ne = ne;
    q->nt = nt;
    q->nh = nh;
    q->magic = MAGIC;
    *pq = q;
    return HE_OK;
}

int
he_read_ini(const char *path, T ** pq)
{
    T *q;
    FILE *f;
    char line[SIZE];
    int i, cnt;
    int nv, nt, ne, nh;

#define NXT() if (util_fgets(line, f) == NULL)  \
	ERR(("unexpected EOF in '%s'", path))
    MALLOC(1, &q);
    f = fopen(path, "r");
    if (f == NULL)
        ERR(("fail to open '%s'", path));
    NXT();
    if (!util_eq(line, "HE"))
        ERR(("'%s' is not a he file", path));
    NXT();
    cnt = sscanf(line, "%d %d %d %d", &nv, &ne, &nt, &nh);
    if (cnt != 4)
        ERR(("'%s' != [nv nt ne nh] in '%s'", line, path));
    if (nv <= 0 || nt <= 0 || ne <= 0 || nh <= 0)
        ERR(("wrong sizes '%s' in '%s'", line, path));

    alloc(nv, ne, nt, nh, /**/ q);

    for (i = 0; i < nh; i++) {
        NXT();
        cnt = sscanf(line, "%d %d  %d %d %d",
                     &q->nxt[i], &q->flp[i], &q->ver[i], &q->edg[i],
                     &q->tri[i]);
        if (cnt != 5)
            ERR(("wrong half-edg line '%s' in '%s'", line, path));
    }
    for (i = 0; i < nv; i++) {
        NXT();
        cnt = sscanf(line, "%d", &q->hdg_ver[i]);
        if (cnt != 1)
            ERR(("wrong ver line '%s' in '%s'", line, path));
    }
    for (i = 0; i < ne; i++) {
        NXT();
        cnt = sscanf(line, "%d", &q->hdg_edg[i]);
        if (cnt != 1)
            ERR(("wrong edg line '%s' in '%s'", line, path));
    }
    for (i = 0; i < nt; i++) {
        NXT();
        cnt = sscanf(line, "%d", &q->hdg_tri[i]);
        if (cnt != 1)
            ERR(("wrong tri line '%s' in '%s'", line, path));
    }
    if (fclose(f) != 0)
        ERR(("fail to close file '%s'", path));

    q->nv = nv;
    q->nt = nt;
    q->ne = ne;
    q->nh = nh;
    q->magic = MAGIC;

    if (valid(nv, q->hdg_ver, q->ver) != OK)
        ERR(("invalid ver references"));
    if (valid(ne, q->hdg_edg, q->edg) != OK)
        ERR(("invalid edg references"));
    if (valid(nt, q->hdg_tri, q->tri) != OK)
        ERR(("invalid tri references"));

    *pq = q;
    return HE_OK;
}

int
he_read_fin(T * q)
{
    if (q->magic != MAGIC)
        ERR(("wrong fin() call"));
    afree(q);
    FREE(q);
    return HE_OK;
}

int
he_read_info(T * q, FILE * f)
{
    int r;
    int nv, nt, ne, nh;
    int *nxt, *flp, *ver, *tri, *edg;
    int *hdg_ver, *hdg_edg, *hdg_tri;

    nv = q->nv;
    nt = q->nt;
    ne = q->ne;
    nh = q->nh;
    nxt = q->nxt;
    flp = q->flp;
    ver = q->ver;
    tri = q->tri;
    edg = q->edg;
    hdg_ver = q->hdg_ver;
    hdg_edg = q->hdg_edg;
    hdg_tri = q->hdg_tri;

    r = fprintf(f, "%d %d %d %d\n", nv, ne, nt, nh);
    if (r <= 0)
        ERR(("fprintf() failed"));

    fprintf(f, "[nh=%d lines]\n", nh);
    fprintf(f, "%d %d %d %d %d\n", nxt[0], flp[0], ver[0], tri[0], edg[0]);
    fputs("...\n", f);
    fprintf(f, "%d %d %d %d %d\n",
            nxt[nh - 1], flp[nh - 1], ver[nh - 1], tri[nh - 1],
            edg[nh - 1]);
    fprintf(f, "[nv=%d lines]\n", nv);
    fprintf(f, "%d\n", hdg_ver[0]);
    fputs("...\n", f);
    fprintf(f, "%d\n", hdg_ver[nv - 1]);

    fprintf(f, "[ne=%d lines]\n", ne);
    fprintf(f, "%d\n", hdg_edg[0]);
    fputs("...\n", f);
    fprintf(f, "%d\n", hdg_edg[ne - 1]);

    fprintf(f, "[nt=%d lines]\n", nt);
    fprintf(f, "%d\n", hdg_tri[0]);
    fputs("...\n", f);
    fprintf(f, "%d\n", hdg_tri[nt - 1]);

    return HE_OK;
}

int
he_read_nv(T * q)
{
    return q->nv;
}

int
he_read_nt(T * q)
{
    return q->nt;
}

int
he_read_ne(T * q)
{
    return q->ne;
}

int
he_read_nh(T * q)
{
    return q->nh;
}

int
he_read_nxt(T * q, int **p)
{
    *p = q->nxt;
    return HE_OK;
}

int
he_read_flp(T * q, int **p)
{
    *p = q->flp;
    return HE_OK;
}

int
he_read_ver(T * q, int **p)
{
    *p = q->ver;
    return HE_OK;
}

int
he_read_tri(T * q, int **p)
{
    *p = q->tri;
    return HE_OK;
}

int
he_read_edg(T * q, int **p)
{
    *p = q->edg;
    return HE_OK;
}

int
he_read_hdg_ver(T * q, int **p)
{
    *p = q->hdg_ver;
    return HE_OK;
}

int
he_read_hdg_edg(T * q, int **p)
{
    *p = q->hdg_edg;
    return HE_OK;
}

int
he_read_hdg_tri(T * q, int **p)
{
    *p = q->hdg_tri;
    return HE_OK;
}

static char *
util_fgets(char *s, FILE * stream)
{
    char *c;

    if (fgets(s, SIZE, stream) == NULL)
        return NULL;
    if ((c = strchr(s, '\n')) != NULL)
        *c = '\0';
    return s;
}

int
util_eq(const char *a, const char *b)
{
    return strncmp(a, b, SIZE) == 0;
}
