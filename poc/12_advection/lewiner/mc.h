struct MarchLewiner;
struct MarchLewiner *march_lewiner_ini(int, int, int);
int march_lewiner_apply(struct MarchLewiner *, double *, int *nv,
			double **vert, int *nt, int **tri);
int march_lewiner_fin(struct MarchLewiner *);
