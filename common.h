#ifndef __COMMON_H__
#define __COMMON_H__

#define PI            3.14159
#define LIGHT         2.99792458e8
#define PERMITTIVITY  8.85418782e-12
#define PERMEABILITY  (4.0 * PI * 1.0e-7)
#define CONDUCTIVITY  0.0

// 材料を現す構造体宣言
typedef struct {
	double eps;     // 誘電率(Relative値ではない)
	double sig;     // 導電率
	double mu;      // 透磁率(Relative値ではない)
} mat_t;

// 電磁界計算用の係数を現す構造体
typedef struct {
	double val;
	double c_hxly;
	double c_hxlz;
} hx_t;

typedef struct {
	double val;
	double c_hylz;
	double c_hylx;
} hy_t;

typedef struct {
	double val;
	double c_hzlx;
	double c_hzly;
} hz_t;

typedef struct {
	double val;
  double c_ex;
	double c_exly;
  double c_exlz;
  double hist;
} ex_t;

typedef struct {
	double val;
  double c_ey;
  double c_eylz;
	double c_eylx;
  double hist;
} ey_t;

typedef struct {
	double val;
  double c_ez;
  double c_ezlx;
	double c_ezly;
  double hist;
} ez_t;

typedef struct {
	int nx;
	int ny;
  int nz;

	double* dx;
	double* dy;
  double* dz;

	ex_t*** ex;
	ey_t*** ey;
  ez_t*** ez;
  hx_t*** hx;
  hy_t*** hy;
	hz_t*** hz;

	double dt;
} grid_t;

#endif // __COMMON_H__
