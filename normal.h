#ifndef __NORMAL_H__
#define __NORMAL_H__

void calc_mg_hx(hx_t*** hx, ey_t*** ey, ez_t*** ez, int nx, int ny, int nz);
void calc_mg_hy(hy_t*** hy, ez_t*** ez, ex_t*** ex, int nx, int ny, int nz);
void calc_mg_hz(hz_t*** hz, ex_t*** ex, ey_t*** ey, int nx, int ny, int nz);
void calc_mg_ex(ex_t*** ex, hy_t*** hy, hz_t*** hz, int nx, int ny, int nz);
void calc_mg_ey(ey_t*** ey, hz_t*** hz, hx_t*** hx, int nx, int ny, int nz);
void calc_mg_ez(ez_t*** ez, hx_t*** hx, hy_t*** hy, int nx, int ny, int nz);
void std_initialize_ex(ex_t*** ex, mat_t*** mat, double* dy, double* dz, int nx, int ny, int nz, double dt);
void std_initialize_ey(ey_t*** ey, mat_t*** mat, double* dz, double* dx, int nx, int ny, int nz, double dt);
void std_initialize_ez(ez_t*** ez, mat_t*** mat, double* dx, double* dy, int nx, int ny, int nz, double dt);
void std_initialize_hx(hx_t*** hx, mat_t*** mat, double* dy, double* dz, int nx, int ny, int nz, double dt);
void std_initialize_hy(hy_t*** hy, mat_t*** mat, double* dz, double* dx, int nx, int ny, int nz, double dt);
void std_initialize_hz(hz_t*** hz, mat_t*** mat, double* dx, double* dy, int nx, int ny, int nz, double dt);

#endif // __NORMAL_H__
