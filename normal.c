#include "common.h"

void calc_mg_hx(hx_t*** hx, ey_t*** ey, ez_t*** ez, int nx, int ny, int nz)
{
  int x, y, z;

#pragma omp parallel for private(x, y, z)
  for (x = 0; x < nx+1; x++) {
    for (y = 0; y < ny; y++) {
      for (z = 0; z < nz; z++) {
        hx[x][y][z].val = hx[x][y][z].val
          - hx[x][y][z].c_hxly * (ez[x][y+1][z].val - ez[x][y][z].val)
          + hx[x][y][z].c_hxlz * (ey[x][y][z+1].val - ey[x][y][z].val);
      }
    }
  }
}

void calc_mg_hy(hy_t*** hy, ez_t*** ez, ex_t*** ex, int nx, int ny, int nz)
{
  int x, y, z;

#pragma omp parallel for private(x, y, z)
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny+1; y++) {
      for (z = 0; z < nz; z++) {
        hy[x][y][z].val = hy[x][y][z].val
          - hy[x][y][z].c_hylz * (ex[x][y][z+1].val - ex[x][y][z].val)
          + hy[x][y][z].c_hylx * (ez[x+1][y][z].val - ez[x][y][z].val);
      }
    }
  }
}

void calc_mg_hz(hz_t*** hz, ex_t*** ex, ey_t*** ey, int nx, int ny, int nz)
{
	int x, y, z;

#pragma omp parallel for private(x, y, z)
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny; y++) {
      for (z = 0; z < nz+1; z++) {
        hz[x][y][z].val = hz[x][y][z].val
          - hz[x][y][z].c_hzlx * (ey[x+1][y][z].val - ey[x][y][z].val)
          + hz[x][y][z].c_hzly * (ex[x][y+1][z].val - ex[x][y][z].val);
      }
		}
	}
}

void calc_mg_ex(ex_t*** ex, hy_t*** hy, hz_t*** hz, int nx, int ny, int nz)
{
	int x, y, z;
	
#pragma omp parallel for private(x, y, z)
	for (x = 0; x < nx; x++) {
		for (y = 1; y < ny; y++) {
      for (z = 1; z < nz; z++) {
        ex[x][y][z].val = ex[x][y][z].c_ex * ex[x][y][z].val
          + ex[x][y][z].c_exly * (hz[x][y][z].val - hz[x][y-1][z].val)
          - ex[x][y][z].c_exlz * (hy[x][y][z].val - hy[x][y][z-1].val);
      }
		}
	}
}

void calc_mg_ey(ey_t*** ey, hz_t*** hz, hx_t*** hx, int nx, int ny, int nz)
{
	int x, y, z;
	
#pragma omp parallel for private(x, y, z)
	for (x = 1; x < nx; x++) {
		for (y = 0; y < ny; y++) {
      for (z = 1; z < nz; z++) {
        ey[x][y][z].val = ey[x][y][z].c_ey * ey[x][y][z].val
          + ey[x][y][z].c_eylz * (hx[x][y][z].val - hx[x][y][z-1].val)
          - ey[x][y][z].c_eylx * (hz[x][y][z].val - hz[x-1][y][z].val);
      }
		}
	}
}

void calc_mg_ez(ez_t*** ez, hx_t*** hx, hy_t*** hy, int nx, int ny, int nz)
{
	int x, y, z;
	
#pragma omp parallel for private(x, y, z)
	for (x = 1; x < nx; x++) {
		for (y = 1; y < ny; y++) {
      for (z = 0; z < nz; z++) {
        ez[x][y][z].val = ez[x][y][z].c_ez * ez[x][y][z].val
          + ez[x][y][z].c_ezlx * (hy[x][y][z].val - hy[x-1][y][z].val)
          - ez[x][y][z].c_ezly * (hx[x][y][z].val - hx[x][y-1][z].val);
      }
		}
	}
}

void std_initialize_ex(ex_t*** ex, mat_t*** mat, double* dy, double* dz, int nx, int ny, int nz, double dt)
{
	int x, y, z;
	double common;

#pragma omp parallel for private(x, y, z, common)
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny+1; y++) {
      for (z = 0; z < nz+1; z++) {
        common = (mat[x][y][z].sig * dt) / (2 * mat[x][y][z].eps);
        ex[x][y][z].c_ex   = (1 - common) / (1 + common);
        ex[x][y][z].c_exly = (dt / mat[x][y][z].eps) / (1 + common) * (1 / dy[y]);
        ex[x][y][z].c_exlz = (dt / mat[x][y][z].eps) / (1 + common) * (1 / dz[z]);
      }
		}
	}
}

void std_initialize_ey(ey_t*** ey, mat_t*** mat, double* dz, double* dx, int nx, int ny, int nz, double dt)
{
	int x, y, z;
	double common;

#pragma omp parallel for private(x, y, z, common)
	for (x = 0; x < nx+1; x++) {
		for (y = 0; y < ny; y++) {
      for (z = 0; z < nz+1; z++) {
        common = (mat[x][y][z].sig * dt) / (2 * mat[x][y][z].eps);
        ey[x][y][z].c_ey   = (1 - common) / (1 + common);
        ey[x][y][z].c_eylz = (dt / mat[x][y][z].eps) / (1 + common) * (1 / dz[z]);
        ey[x][y][z].c_eylx = (dt / mat[x][y][z].eps) / (1 + common) * (1 / dx[x]);
      }
		}
	}
}

void std_initialize_ez(ez_t*** ez, mat_t*** mat, double* dx, double* dy, int nx, int ny, int nz, double dt)
{
	int x, y, z;
	double common;

#pragma omp parallel for private(x, y, z, common)
	for (x = 0; x < nx+1; x++) {
		for (y = 0; y < ny+1; y++) {
      for (z = 0; z < nz; z++) {
        common = (mat[x][y][z].sig * dt) / (2 * mat[x][y][z].eps);
        ez[x][y][z].c_ez   = (1 - common) / (1 + common);
        ez[x][y][z].c_ezlx = (dt / mat[x][y][z].eps) / (1 + common) * (1 / dx[x]);
        ez[x][y][z].c_ezly = (dt / mat[x][y][z].eps) / (1 + common) * (1 / dy[y]);
      }
		}
	}
}

void std_initialize_hx(hx_t*** hx, mat_t*** mat, double* dy, double* dz, int nx, int ny, int nz, double dt)
{
	int x, y, z;

#pragma omp parallel for private(x, y, z)
	for (x = 0; x < nx+1; x++) {
		for (y = 0; y < ny; y++) {
      for (z = 0; z < nz; z++) {
        hx[x][y][z].c_hxly = dt / mat[x][y][z].mu * (1 / dy[y]);
        hx[x][y][z].c_hxlz = dt / mat[x][y][z].mu * (1 / dz[z]);
      }
		}
	}
}

void std_initialize_hy(hy_t*** hy, mat_t*** mat, double* dz, double* dx, int nx, int ny, int nz, double dt)
{
	int x, y, z;

#pragma omp parallel for private(x, y, z)
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny+1; y++) {
      for (z = 0; z < nz; z++) {
        hy[x][y][z].c_hylz = dt / mat[x][y][z].mu * (1 / dz[z]);
        hy[x][y][z].c_hylx = dt / mat[x][y][z].mu * (1 / dx[x]);
      }
		}
	}
}

void std_initialize_hz(hz_t*** hz, mat_t*** mat, double* dx, double* dy, int nx, int ny, int nz, double dt)
{
	int x, y, z;

#pragma omp parallel for private(x, y, z)
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny; y++) {
      for (z = 0; z < nz+1; z++) {
        hz[x][y][z].c_hzlx = dt / mat[x][y][z].mu * (1 / dx[x]);
        hz[x][y][z].c_hzly = dt / mat[x][y][z].mu * (1 / dy[y]);
      }
		}
	}
}

