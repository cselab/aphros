typedef unsigned char uchar;
typedef signed char schar;
typedef struct {
    double x, y, z;
} Vertex;
typedef struct {
    int v1, v2, v3;
} Triangle;
struct MarchingCubes {
    MarchingCubes(int size_x, int size_y, int size_z);
    ~MarchingCubes();
    double get_data(int i, int j, int k);
    void set_data(double, int, int, int);
    void clean_temps();
    void clean_all();
    void run();
    void process_cube();
    bool test_face(schar face);
    bool test_interior(schar s);

    void compute_intersection_points();
    void add_triangle(const char *trig, char n, int v12 = -1);
    void test_vertex_addition();
    int add_x_vertex();
    int add_y_vertex();
    int add_z_vertex();
    int add_c_vertex();

    int get_x_vert(int i, int j, int k);
    int get_y_vert(int i, int j, int k);
    int get_z_vert(int i, int j, int k);

    void set_x_vert(int val, int i, int j, int k);
    void set_y_vert(int val, int i, int j, int k);
    void set_z_vert(int val, int i, int j, int k);

    int size_x;
    int size_y;
    int size_z;
    double *data;
    int *x_verts;
    int *y_verts;
    int *z_verts;
    int nverts;
    int ntrigs;
    int Nverts;
    int Ntrigs;
    Vertex *vertices;
    Triangle *triangles;
    int i;
    int j;
    int k;
    double cube[8];
    uchar lut_entry;
    uchar Case;
    uchar config;
    uchar subconfig;
};

void writeObj(MarchingCubes *);
