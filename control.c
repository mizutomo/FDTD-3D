#include "fdtd.h"
#include <GL/glut.h>

double rx=15, rz=30, sz=2;
int fm=0;
double mx, my, ms;

void display()
{
  glClearColor(0, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  gluPerspective(60, 1, 1, 100);
  gluLookAt(0, -70, 0, 0, 0, 0, 0, 0, 1);
  glRotated(rx, 1, 0, 0);
  glRotated(-rz, 0, 0, 1);
  glPointSize(sz);
  glBegin(GL_POINTS);

  //koko

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
  //koko
  glutPostRedisplay();
}

void resize(int w, int h)
{
  glViewport(0, 0, w, h);
}

int main()
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

  close_output_files();

  return 0;
}
