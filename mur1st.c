#include <stdio.h>
#include "common.h"

void update_xy(ex_t*** ex, ey_t*** ey, double* dz, int nx, int ny, int nz, double dt)
{
  int x, y;
  double c = LIGHT * dt;

  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny+1; y++) {
      ex[x][y][0].val = ex[x][y][1].hist
        + (c - dz[0]) / (c + dz[0]) * (ex[x][y][1].val - ex[x][y][0].hist);
      ex[x][y][nz].val = ex[x][y][nz-1].hist
        + (c - dz[nz-1]) / (c + dz[nz-1]) * (ex[x][y][nz-1].val - ex[x][y][nz].hist);
    }
  }

  for (y = 0; y < ny; y++) {
    for (x = 0; x < nx+1; x++) {
      ey[x][y][0].val = ey[x][y][1].hist
        + (c - dz[0]) / (c + dz[0]) * (ey[x][y][1].val - ey[x][y][0].hist);
      ey[x][y][nz].val = ey[x][y][nz-1].hist
        + (c - dz[nz-1]) / (c + dz[nz-1]) * (ey[x][y][nz-1].val - ey[x][y][nz].hist);
    }
  }
}

void update_yz(ey_t*** ey, ez_t*** ez, double* dx, int nx, int ny, int nz, double dt)
{
  int y, z;
  double c = LIGHT * dt;

  for (y = 0; y < ny; y++) {
    for (z = 0; z < nz+1; z++) {
      ey[0][y][z].val = ey[1][y][z].hist
        + (c - dx[0]) / (c + dx[0]) * (ey[1][y][z].val - ey[0][y][z].hist);
      ey[nx][y][z].val = ey[nx-1][y][z].hist
        + (c - dx[nx-1]) / (c + dx[nx-1]) * (ey[nx-1][y][z].val - ey[nx][y][z].hist);
    }
  }

  for (z = 0; z < nz; z++) {
    for (y = 0; y < ny+1; y++) {
      ez[0][y][z].val = ez[1][y][z].hist
        + (c - dx[0]) / (c + dx[0]) * (ez[1][y][z].val - ez[0][y][z].hist);
      ez[nx][y][z].val = ez[nx-1][y][z].hist
        + (c - dx[nx-1]) / (c + dx[nx-1]) * (ez[nx-1][y][z].val - ez[nx][y][z].hist);
    }
  }
}

void update_zx(ez_t*** ez, ex_t*** ex, double* dy, int nx, int ny, int nz, double dt)
{
  int z, x;
  double c = LIGHT * dt;

  for (z = 0; z < nz; z++) {
    for (x = 0; x < nx+1; x++) {
      ez[x][0][z].val = ez[x][1][z].hist
        + (c - dy[0]) / (c + dy[0]) * (ez[x][1][z].val - ez[x][0][z].hist);
      ez[x][ny][z].val = ez[x][ny-1][z].hist
        + (c - dy[ny-1]) / (c + dy[ny-1]) * (ez[x][ny-1][z].val - ez[x][ny][z].hist);
    }
  }

  for (x = 0; x < nx; x++) {
    for (z = 0; z < nz+1; z++) {
      ex[x][0][z].val = ex[x][1][z].hist
        + (c - dy[0]) / (c + dy[0]) * (ex[x][1][z].val - ex[x][0][z].hist);
      ex[x][ny][z].val = ex[x][ny-1][z].hist
        + (c - dy[ny-1]) / (c + dy[ny-1]) * (ex[x][ny-1][z].val - ex[x][ny][z].hist);
    }
  }
}

void ex_hist_update(ex_t*** ex, int nx, int ny, int nz)
{
	int x, y, z;

	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny+1; y++) {
      for (z = 0; z < nz+1; z++) {
        ex[x][y][z].hist = ex[x][y][z].val;
      } 
		}
	}
}

void ey_hist_update(ey_t*** ey, int nx, int ny, int nz)
{
	int x, y, z;

	for (x = 0; x < nx+1; x++) {
		for (y = 0; y < ny; y++) {
      for (z = 0; z < nz+1; z++) {
        ey[x][y][z].hist = ey[x][y][z].val;
      }
		}
	}
}

void ez_hist_update(ez_t*** ez, int nx, int ny, int nz)
{
	int x, y, z;

	for (x = 0; x < nx+1; x++) {
		for (y = 0; y < ny+1; y++) {
      for (z = 0; z < nz; z++) {
        ez[x][y][z].hist = ez[x][y][z].val;
      }
		}
	}
}
