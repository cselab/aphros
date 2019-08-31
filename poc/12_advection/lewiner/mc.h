void MarchingCubes(int, int, int, double *);
void run();
void writeObj();

struct MarchLewiner;
struct MarchLewiner *march_ini(int[3]);
int march_apply(struct MarchLewiner *, double *, /**/ int *nv,
		double **vert, int *nt, int **tri);
int march_fin(struct MarchLewiner *);
