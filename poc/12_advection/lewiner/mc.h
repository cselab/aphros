typedef unsigned char uchar;
typedef signed char schar;
typedef float real;
typedef struct {
    real x, y, z;
} Vertex;
typedef struct {
    int v1, v2, v3;
} Triangle;
struct MarchingCubes {
MarchingCubes(int size_x, int size_y, int size_z);
    ~MarchingCubes();
    inline real get_data(const int i, const int j, const int k) const {
	return _data[i + j * _size_x + k * _size_x * _size_y];
    }
    inline void set_data(const real val, const int i, const int j,
			 const int k) {
	_data[i + j * _size_x + k * _size_x * _size_y] = val;
    }
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
    inline int get_x_vert( int i,  int j,  int k)  {
	return _x_verts[i + j * _size_x + k * _size_x * _size_y];
    }
    inline int get_y_vert( int i,  int j,  int k)  {
	return _y_verts[i + j * _size_x + k * _size_x * _size_y];
    }
    inline int get_z_vert( int i,  int j,  int k)  {
	return _z_verts[i + j * _size_x + k * _size_x * _size_y];
    }
    inline void set_x_vert( int val,  int i,  int j,
			    int k) {
	_x_verts[i + j * _size_x + k * _size_x * _size_y] = val;
    }
    inline void set_y_vert( int val,  int i,  int j,
			    int k) {
	_y_verts[i + j * _size_x + k * _size_x * _size_y] = val;
    }
    inline void set_z_vert( int val,  int i,  int j,
			    int k) {
	_z_verts[i + j * _size_x + k * _size_x * _size_y] = val;
    }
    int _size_x;
    int _size_y;
    int _size_z;
    real *_data;
    int *_x_verts;
    int *_y_verts;
    int *_z_verts;
    int _nverts;
    int _ntrigs;
    int _Nverts;
    int _Ntrigs;
    Vertex *_vertices;
    Triangle *_triangles;
    int _i;
    int _j;
    int _k;
    real _cube[8];
    uchar _lut_entry;
    uchar _case;
    uchar _config;
    uchar _subconfig;
};

void writeObj(MarchingCubes*);
