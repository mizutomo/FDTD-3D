#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <signal.h>
#include "common.h"
#include "normal.h"
#include "mur1st.h"

// グローバル変数
#define FILENUM (8)
FILE* fps[FILENUM];

double square(double x)
{
  return x * x;
}

void set_delta_map(double* d, int leng, float size)
{
  int i;

  for (i = 0; i < leng; i++) {
    d[i] = size;
  }
}

double calc_cfl_constant(float cfl, double dx, double dy, double dz)
{
  double dt;

  dt = cfl * 1.0 / (LIGHT * sqrt(square(1.0/dx) + square(1.0/dy) + square(1.0/dz)));

  return dt;
}

mat_t*** alloc_3d_array_mat(int nx, int ny, int nz)
{
  int i, j;
  mat_t*** ary;

  ary = (mat_t***)malloc(sizeof(mat_t**) * nx);
  for (i = 0; i < nx; i++) {
    ary[i] = (mat_t**)malloc(sizeof(mat_t*) * ny);
    for (j = 0; j < ny; j++) {
      ary[i][j] = (mat_t*)malloc(sizeof(mat_t) * nz);
    }
  }

  return ary;
}

ex_t*** alloc_3d_array_ex(int nx, int ny, int nz)
{
  int i, j;
  ex_t*** ary;

  ary = (ex_t***)malloc(sizeof(ex_t**) * nx);
  for (i = 0; i < nx; i++) {
    ary[i] = (ex_t**)malloc(sizeof(ex_t*) * ny);
    for (j = 0; j < ny; j++) {
      ary[i][j] = (ex_t*)malloc(sizeof(ex_t) * nz);
    }
  }

  return ary;
}

ey_t*** alloc_3d_array_ey(int nx, int ny, int nz)
{
  int i, j;
  ey_t*** ary;

  ary = (ey_t***)malloc(sizeof(ey_t**) * nx);
  for(i = 0; i < nx; i++) {
    ary[i] = (ey_t**)malloc(sizeof(ey_t*) * ny);
    for (j = 0; j < ny; j++) {
      ary[i][j] = (ey_t*)malloc(sizeof(ey_t) * nz);
    }
  }

  return ary;
}

ez_t*** alloc_3d_array_ez(int nx, int ny, int nz)
{
  int i, j;
  ez_t*** ary;

  ary = (ez_t***)malloc(sizeof(ez_t**) * nx);
  for(i = 0; i < nx; i++) {
    ary[i] = (ez_t**)malloc(sizeof(ez_t*) * ny);
    for (j = 0; j < ny; j++) {
      ary[i][j] = (ez_t*)malloc(sizeof(ez_t) * nz);
    }
  }

  return ary;
}

hx_t*** alloc_3d_array_hx(int nx, int ny, int nz)
{
  int i, j;
  hx_t*** ary;

  ary = (hx_t***)malloc(sizeof(hx_t**) * nx);
  for (i = 0; i < nx; i++) {
    ary[i] = (hx_t**)malloc(sizeof(hx_t*) * ny);
    for (j = 0; j < ny; j++) {
      ary[i][j] = (hx_t*)malloc(sizeof(hx_t) * nz);
    }
  }

  return ary;
}

hy_t*** alloc_3d_array_hy(int nx, int ny, int nz)
{
  int i, j;
  hy_t*** ary;

  ary = (hy_t***)malloc(sizeof(hy_t**) * nx);
  for (i = 0; i < nx; i++) {
    ary[i] = (hy_t**)malloc(sizeof(hy_t*) * ny);
    for (j = 0; j < ny; j++) {
      ary[i][j] = (hy_t*)malloc(sizeof(hy_t) * nz);
    }
  }

  return ary;
}

hz_t*** alloc_3d_array_hz(int nx, int ny, int nz)
{
  int i, j;
  hz_t*** ary;

  ary = (hz_t***)malloc(sizeof(hz_t**) * nx);
  for (i = 0; i < nx; i++) {
    ary[i] = (hz_t**)malloc(sizeof(hz_t*) * ny);
    for (j = 0; j < ny; j++) {
      ary[i][j] = (hz_t*)malloc(sizeof(hz_t) * nz);
    }
  }

  return ary;
}

void initialize_ex(ex_t*** ary, int nx, int ny, int nz)
{
	int i, j, k;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        ary[i][j][k].val = 0.0;
        ary[i][j][k].hist = 0.0;
      }
		}
	}
}

void initialize_ey(ey_t*** ary, int nx, int ny, int nz)
{
	int i, j, k;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        ary[i][j][k].val = 0.0;
        ary[i][j][k].hist = 0.0;
      }
		}
	}
}

void initialize_ez(ez_t*** ary, int nx, int ny, int nz)
{
	int i, j, k;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        ary[i][j][k].val = 0.0;
        ary[i][j][k].hist = 0.0;
      }
		}
	}
}

void initialize_hx(hx_t*** ary, int nx, int ny, int nz)
{
	int i, j, k;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        ary[i][j][k].val = 0.0;
      }
		}
	}
}

void initialize_hy(hy_t*** ary, int nx, int ny, int nz)
{
	int i, j, k;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        ary[i][j][k].val = 0.0;
      }
		}
	}
}

void initialize_hz(hz_t*** ary, int nx, int ny, int nz)
{
	int i, j, k;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        ary[i][j][k].val = 0.0;
      }
		}
	}
}

void initialize_mat_array(mat_t*** ary, int nx, int ny, int nz)
{
	int i, j, k;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        ary[i][j][k].eps = 1.0 * PERMITTIVITY;
        ary[i][j][k].sig = 0.0;
        ary[i][j][k].mu  = 1.0 * PERMEABILITY;
      }
		}
	}
}

void free_3d_ary(void*** ary, int nx, int ny)
{
	int i, j;

	for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      free(ary[i][j]);
    }
    free(ary[i]);
	}
	free(ary);
}

void check_opt(int argc, char** argv, float* cfl)
{
	if (argc != 2) {
		printf("Error: Invalid Arguments.\n");
		printf("Usage: fdtd cfl\n");
		exit(10);
	} else {
		*cfl = atof(argv[1]);
	}
}

void setup_grid(grid_t* grid, float cfl, int nx, int ny, int nz,
                double delta_x, double delta_y, double delta_z, double dt)
{
	mat_t*** mat;

	grid->nx = nx;
	grid->ny = ny;
  grid->nz = nz;
	grid->dx = (double*)malloc(sizeof(double) * nx+1);
	grid->dy = (double*)malloc(sizeof(double) * ny+1);
  grid->dz = (double*)malloc(sizeof(double) * nz+1);
	grid->dt = dt;

	// グリッド初期化
	set_delta_map(grid->dx, nx+1, delta_x);
	set_delta_map(grid->dy, ny+1, delta_y);
  set_delta_map(grid->dz, nz+1, delta_z);

	// CFL条件の設定
  printf("Main Grid Info\n");
	printf("  Num X  : %d\n", nx);
	printf("  Num Y  : %d\n", ny);
  printf("  Num Z  : %d\n", nz);
	printf("  Min X  : %f\n", delta_x);
	printf("  Min Y  : %f\n", delta_y);
  printf("  Min Z  : %f\n", delta_z);
	printf("  CFL    : %f\n", cfl);
	printf("  dt     : %g\n", grid->dt);

	// 材料配列の確保
	mat = alloc_3d_array_mat(nx+1, ny+1, nz+1);
	initialize_mat_array(mat, nx+1, ny+1, nz+1);

	// 電磁界配列の確保
	grid->ex = alloc_3d_array_ex(nx,   ny+1, nz+1);
	grid->ey = alloc_3d_array_ey(nx+1, ny,   nz+1);
  grid->ez = alloc_3d_array_ez(nx+1, ny+1, nz);
  grid->hx = alloc_3d_array_hx(nx+1, ny,   nz);
  grid->hy = alloc_3d_array_hy(nx,   ny+1, nz);
	grid->hz = alloc_3d_array_hz(nx,   ny,   nz+1);
	initialize_ex(grid->ex, nx,   ny+1, nz+1);
	initialize_ey(grid->ey, nx+1, ny,   nz+1);
  initialize_ez(grid->ez, nx+1, ny+1, nz);
  initialize_hx(grid->hx, nx+1, ny,   nz);
  initialize_hy(grid->hy, nx,   ny+1, nz);
  initialize_hz(grid->hz, nx,   ny,   nz+1);

	// 係数行列の計算
	std_initialize_ex(grid->ex, mat, grid->dy, grid->dz, nx, ny, nz, grid->dt);
	std_initialize_ey(grid->ey, mat, grid->dz, grid->dx, nx, ny, nz, grid->dt);
  std_initialize_ez(grid->ez, mat, grid->dx, grid->dy, nx, ny, nz, grid->dt);
	std_initialize_hx(grid->hx, mat, grid->dy, grid->dz, nx, ny, nz, grid->dt);
	std_initialize_hy(grid->hy, mat, grid->dz, grid->dx, nx, ny, nz, grid->dt);
	std_initialize_hz(grid->hz, mat, grid->dx, grid->dy, nx, ny, nz, grid->dt);

	free_3d_ary((void***)mat, nx, ny);
}

void open_output_files(FILE** fps)
{
  fps[0] = fopen("wave/fdtd_point.csv", "w");
  fps[1] = fopen("wave/fdtd_point_subgrid.csv", "w");
  fps[2] = fopen("wave/fdtd_x.csv", "w");
  fps[3] = fopen("wave/fdtd_y.csv", "w");
  fps[4] = fopen("wave/fdtd_xy.csv", "w");
  fps[5] = fopen("wave/fdtd_map.csv", "w");
  fps[6] = fopen("wave/fdtd_map_subgrid.csv", "w");
  fps[7] = fopen("wave/fdtd_energy.csv", "w");
}

void close_output_files(FILE** fps)
{
  int i;
  for (i = 0; i < FILENUM; i++) {
    fclose(fps[i]);
  }
}

void print_point_value(FILE* fp, double time, int num, ...)
{
  va_list list;
  int i;

  fprintf(fp, "%g", time);

  va_start(list, num);
  for (i = 0; i < num; i++) {
    fprintf(fp, ", %g", va_arg(list, double));
  }
  fprintf(fp, "\n");

  va_end(list);
}

void print_map_value(FILE* fp, grid_t* mg, int step)
{
	int x, y;

	fprintf(fp, "# STEP = %d\n", step);
	for (x = 0; x < mg->nx; x++) {
		for (y = 0; y < mg->ny; y++) {
      fprintf(fp, "%d, %d, %g\n", x, y, sqrt(square(mg->hz[x][y][mg->nz/2].val)));
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}

void write_waveforms(FILE** fps, int step, grid_t* mg)
{
  double current_time = step * mg->dt;

  print_point_value(fps[0], current_time, 3,
                    mg->hz[mg->nx-10][mg->ny/2][mg->nz/2].val, mg->hz[mg->nx-30][mg->ny/2][mg->nz/2].val,
                    mg->hz[10][mg->ny/2][mg->nz/2].val);

  if (step % 10 == 0) {
    print_map_value(fps[5], mg, step);
  }
}

void update_mg_h_field(grid_t* mg)
{
  calc_mg_hx(mg->hx, mg->ey, mg->ez, mg->nx, mg->ny, mg->nz);
  calc_mg_hy(mg->hy, mg->ez, mg->ex, mg->nx, mg->ny, mg->nz);
  calc_mg_hz(mg->hz, mg->ex, mg->ey, mg->nx, mg->ny, mg->nz);
}

void update_mg_e_field(grid_t* mg)
{
  // Ex&Eyの計算
  calc_mg_ex(mg->ex, mg->hy, mg->hz, mg->nx, mg->ny, mg->nz);
  calc_mg_ey(mg->ey, mg->hz, mg->hx, mg->nx, mg->ny, mg->nz);
  calc_mg_ez(mg->ez, mg->hx, mg->hy, mg->nx, mg->ny, mg->nz);
}

void update_external_boundary(grid_t* mg)
{
  update_xy(mg->ex, mg->ey, mg->dz, mg->nx, mg->ny, mg->nz, mg->dt);
  update_yz(mg->ey, mg->ez, mg->dx, mg->nx, mg->ny, mg->nz, mg->dt);
  update_zx(mg->ez, mg->ex, mg->dy, mg->nx, mg->ny, mg->nz, mg->dt);
  ex_hist_update(mg->ex, mg->nx, mg->ny, mg->nz);
  ey_hist_update(mg->ey, mg->nx, mg->ny, mg->nz);
  ez_hist_update(mg->ez, mg->nx, mg->ny, mg->nz);  
}

double calc_gauss_wave(int step, double dt)
{
	double tc = 10e-9; 
	double pw = 1e-9;
	double t_current = ((float)step - 0.5) * dt;
	double t_eff;

	t_eff = (t_current - tc) / pw;

	return 1.0 * exp(-1 * (t_eff * t_eff));	
}

void inject_stimulus(grid_t* mg, int step)
{
  double stimulus = calc_gauss_wave(step, mg->dt);
  mg->hz[(mg->nx)/2][(mg->ny)/2][(mg->nz)/2].val += stimulus * mg->dt;   // Soft Source 
}

void calc_fdtd(grid_t* mg, double stop_time)
{
  int step;
  int last_step = (int)(stop_time / mg->dt);

  open_output_files(fps);

  for (step = 0; step <= last_step; step++) {
    if (step % 100 == 0) {
      printf("%d / %d (%.2f %%)\n", step, last_step, (double)step/last_step*100);
    }

    update_mg_h_field(mg);
    inject_stimulus(mg, step);

    update_mg_e_field(mg);

    update_external_boundary(mg);

    write_waveforms(fps, step, mg);
  }

  close_output_files(fps);
}

void sig_handler_sigint(int sig) {
  close_output_files(fps);
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
  double stop_time = 100e-9;
  int nx = 100;
  int ny = 100;
  int nz = 100;
  double dx = 0.09;
  double dy = 0.09;
  double dz = 0.09;

  float cfl;

  grid_t mg;

  check_opt(argc, argv, &cfl);

  double dt = calc_cfl_constant(cfl, dx, dy, dz);

  setup_grid(&mg, cfl, nx, ny, nz, dx, dy, dz, dt);

  if (SIG_ERR == signal(SIGINT, sig_handler_sigint)) {
    printf("[Error] failed to set signal handler.\n");
    exit(EXIT_FAILURE);
  }  

  printf("FDTD Calculating...\n");
	calc_fdtd(&mg, stop_time);
	printf("Finished FDTD Calculating.\n");

  return 0;
}
