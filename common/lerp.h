#ifndef _LERP_H
#define _LERP_H

#include <stdbool.h>
#include <cmath>

bool inside_st_sz(int num_dims, const int *st, const int *sz, const float *p);
bool inside_lb_ub(int num_dims, const int *lb, const int *ub, const int *p);

float texel2D(const float *ptr, const int *sz, int x, int y);
float texel3D(const float *ptr, const int *sz, int x, int y, int z);
float texel4D(const float *ptr, const int *sz, int x, int y, int z, int t);

bool lerp2D(const float *pt, const int *gst,
    const int *gsz, const int *st, const int *sz, int num_vars, const float **ptrs, float *vars);
bool lerp3D(const float *pt, const int *gst,
    const int *gsz, const int *st, const int *sz, int num_vars, const float **ptrs, float *vars);
bool lerp4D(const float *pt, const int *gst,
    const int *gsz, const int *st, const int *sz, int num_vars, const float **ptrs, float *vars);

#endif
