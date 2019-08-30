#include <stdio.h>
#include <tgmath.h>
#include <GL/glut.h>
#include "sample.h"
#include "lib.h"
#include "march.h"
#define	USED(x)		if(x);else{}

struct Sample *sample;
static double
f(double x, double y, double z, void *p)
{
    struct Sample *sample;

    sample = p;
    return sample_f(sample, x, y, z);
}

static int
normal(double x, double y, double z, void *p)
{
    double u, v, w;
    USED(p);
    GetColor(x, y, z, &u, &v, &w);
    glColor3f(u, v, w);
    glNormal3f(x, y, z);
    return 0;
}

static int
vertex(double x, double y, double z, void *p)
{
    USED(p);
    glVertex3f(x, y, z);
    return 0;
}

static float AmbientGreen[] = { 0.00, 0.25, 0.00, 1.00 };
static float AmbientBlue[] = { 0.00, 0.00, 0.25, 1.00 };
static float DiffuseGreen[] = { 0.00, 0.75, 0.00, 1.00 };
static float DiffuseBlue[] = { 0.00, 0.00, 0.75, 1.00 };
static float SpecularWhite[] = { 1.00, 1.00, 1.00, 1.00 };

int n = 16;
double h;
int PolygonMode = GL_FILL;
int Spin = 0;
int Move = 0;
int Light = 0;

static void
Resize(int Width, int Height)
{
    double Aspect, HalfWorldSize = (1.4142135623730950488016887242097 / 2);

    glViewport(0, 0, Width, Height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (Width <= Height) {
	Aspect = (double) Height / (double) Width;
	glOrtho(-HalfWorldSize, HalfWorldSize, -HalfWorldSize * Aspect,
		HalfWorldSize * Aspect, -10 * HalfWorldSize,
		10 * HalfWorldSize);
    } else {
	Aspect = (double) Width / (double) Height;
	glOrtho(-HalfWorldSize * Aspect, HalfWorldSize * Aspect,
		-HalfWorldSize, HalfWorldSize, -10 * HalfWorldSize,
		10 * HalfWorldSize);
    }

    glMatrixMode(GL_MODELVIEW);
}

static int
(*march0) (struct March*) = march_cube;

static void
Keyboard(unsigned char Key, int i, int j)
{
    USED(i);
    USED(j);
    switch (Key) {
    case 'w':
	if (PolygonMode == GL_LINE) {
	    PolygonMode = GL_FILL;
	} else {
	    PolygonMode = GL_LINE;
	}
	glPolygonMode(GL_FRONT_AND_BACK, PolygonMode);
	break;
    case '+':
    case '=':
	++n;
	h = 1.0 / n;
	break;
    case '-':
	if (n > 1) {
	    --n;
	    h = 1.0 / n;
	}
	break;
    case 'c':
	march0 = (march0 == march_cube) ? march_tetrahedron : march_cube;
	break;
    case 's':
	sample_next(sample);
	break;
    case 'l':
	if (Light)
	    glDisable(GL_LIGHTING);
	else
	    glEnable(GL_LIGHTING);
	Light = !Light;
    }
}

static void
Special(int Key, int i, int j)
{
    USED(i);
    USED(j);
    switch (Key) {
    case GLUT_KEY_PAGE_UP:
	sample_inc(sample);
	break;
    case GLUT_KEY_PAGE_DOWN:
	sample_dec(sample);
	break;
    case GLUT_KEY_HOME:
	Spin = !Spin;
	break;
    case GLUT_KEY_END:
	Move = !Move;
	break;
    }
}

static void
Idle(void)
{
    glutPostRedisplay();
}

static void
DrawScene(void)
{
    struct March march;
    static double Pitch = 0.0;
    static double Yaw = 0.0;
    static double Time = 0.0;
    march.f = f;
    march.fdata = sample;
    march.size[0] = n;
    march.size[1] = n;
    march.size[2] = n;
    march.spacing = h;
    march.normal = normal;
    march.vertex = vertex;
    march.cdata = NULL;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();

    if (Spin) {
	Pitch += 4.0;
	Yaw += 2.5;
    }
    if (Move) {
	Time += 0.025;
    }

    sample_time(sample, Time);

    glTranslatef(0.0, 0.0, -1.0);
    glRotatef(-Pitch, 1.0, 0.0, 0.0);
    glRotatef(0.0, 0.0, 1.0, 0.0);
    glRotatef(Yaw, 0.0, 0.0, 1.0);

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glColor3f(1.0, 1.0, 1.0);
    glutWireCube(1.0);
    glPopAttrib();

    glPushMatrix();
    glTranslatef(-0.5, -0.5, -0.5);
    glBegin(GL_TRIANGLES);
    march0(&march);
    glEnd();
    glPopMatrix();
    glPopMatrix();

    glutSwapBuffers();
}

int
main(int argc, char **argv)
{
    float PropertiesAmbient[] = { 0.50, 0.50, 0.50, 1.00 };
    float PropertiesDiffuse[] = { 0.75, 0.75, 0.75, 1.00 };
    float PropertiesSpecular[] = { 1.00, 1.00, 1.00, 1.00 };

    int Width = 640;
    int Height = 480;

    h = 1.0 / n;
    sample = sample_ini();
    glutInit(&argc, argv);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(Width, Height);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    glutCreateWindow("Marching Cubes");
    glutDisplayFunc(DrawScene);
    glutIdleFunc(Idle);
    glutReshapeFunc(Resize);
    glutKeyboardFunc(Keyboard);
    glutSpecialFunc(Special);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClearDepth(1.0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, PolygonMode);

    glLightfv(GL_LIGHT0, GL_AMBIENT, PropertiesAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, PropertiesDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, PropertiesSpecular);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);

    glEnable(GL_LIGHT0);

    glMaterialfv(GL_BACK, GL_AMBIENT, AmbientGreen);
    glMaterialfv(GL_BACK, GL_DIFFUSE, DiffuseGreen);
    glMaterialfv(GL_FRONT, GL_AMBIENT, AmbientBlue);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, DiffuseBlue);
    glMaterialfv(GL_FRONT, GL_SPECULAR, SpecularWhite);
    glMaterialf(GL_FRONT, GL_SHININESS, 25.0);

    Resize(Width, Height);

    PrintHelp();
    glutMainLoop();
}
