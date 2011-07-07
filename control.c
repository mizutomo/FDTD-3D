#include <stdlib.h>
#include "fdtd.h"
#include <GL/glut.h>

double rx=15, rz=30, sz=2;
int fm=0;
double mx, my, ms;

double dt;
int step, last_step;
int last_step;
double stop_time = 100e-9;

void finalize();

void display()
{
  int i, j, k;
  double x, y, z;
  double val;

  glClearColor(0, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  gluPerspective(60, 1, 1, 100);
  gluLookAt(0, -70, 0, 0, 0, 0, 0, 0, 1);
  glRotated(rx, 1, 0, 0);
  glRotated(-rz, 0, 0, 1);
  glPointSize(sz);
  glBegin(GL_POINTS);

  for (i = 0; i < 100; i++) {
    for (j = 0; j < 100; j++) {
      for (k = 0; k < 100; k++) {
        x = (double)(i - 100/2);
        y = (double)(j - 100/2);
        z = (double)(k - 100/2);
        val = eabs(i, j, k);
        if (val < 0.001) {
          continue;
        } else if (val > 0.3) {
          if (val > 1) val = 1;
          glColor4d(1, 1, 1, val);
        } else if (val > 0.1) {
          glColor4d(1, 0, 0, 2*val);
        } else if (val > 0.03) {
          glColor4d(1, 1, 0, 4*val);
        } else if (val > 0.01) {
          glColor4d(0, 1, 0, 8*val);
        } else if (val > 0.003) {
          glColor4d(0, 1, 1, 16*val);
        } else if (val >= 0.001) {
          glColor4d(0, 0, 1, 32*val);
        }
        glVertex3d(x, y, z);
      }
    }
  }

  glEnd();
  glutSwapBuffers();
}

void mouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    if (button == GLUT_LEFT_BUTTON) {
      if (fm) return;
      fm = 1;
      mx = rz + (double)x / 2;
      my = rx - (double)y / 2;
    } else if (button == GLUT_RIGHT_BUTTON) {
      if (fm) return;
      fm = 2;
      ms = sz + (double)y / 10;
    }
  } else {
    fm = 0;
  }
}

void motion(int x, int y)
{
  if (fm == 1) {
    rz = mx - (double)x / 2;
    rx = my + (double)y / 2;
    glutPostRedisplay();
  } else if (fm == 2) {
    sz = ms - (double)y / 10;
  }
}

void idle()
{
  int i;

  for (i = 0; i < 10; i++) {
    if (step >= last_step) {
      finalize();
    } else {
      calc_fdtd(step, last_step);
      step++;
    }
  }

  glutPostRedisplay();
}

void resize(int w, int h)
{
  glViewport(0, 0, w, h);
}

void finalize()
{
  close_output_files();
  exit(0);
}

int main_loop()
{
  double dt;
  int step, last_step;
  double stop_time = 100e-9;

  init();
  dt = get_delta_t();
  open_output_files();

  last_step = (int)(stop_time / dt);
  for (step = 0; step <= last_step; step++) {
    calc_fdtd(step, last_step);
  }

  return 0;
}

int main(int argc, char* argv[])
{
  init();
  dt = get_delta_t();
  last_step = (int)(stop_time / dt);
  open_output_files();

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutCreateWindow(argv[1]);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  glutReshapeFunc(resize);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutMainLoop();

  close_output_files();
  return 0;
}
