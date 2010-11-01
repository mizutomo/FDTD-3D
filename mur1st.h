#ifndef __MUR1ST_H__
#define __MUR1ST_H__

void update_xy(ex_t*** ex, ey_t*** ey, double* dz, int nx, int ny, int nz, int dt);
void update_yz(ey_t*** ey, ez_t*** ez, double* dx, int nx, int ny, int nz, int dt);
void update_zx(ez_t*** ez, ex_t*** ex, double* dy, int nx, int ny, int nz, int dt);
void ex_hist_update(ex_t*** ex, int nx, int ny, int nz);
void ey_hist_update(ey_t*** ey, int nx, int ny, int nz);
void ez_hist_update(ez_t*** ez, int nx, int ny, int nz);

#endif // __MUR1ST_H__
