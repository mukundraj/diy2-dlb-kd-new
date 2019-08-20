#ifndef _ADVECT_H
#define _ADVECT_H

#include <stdbool.h>

enum {
  TRACE_SUCC=0,
  TRACE_OUT_OF_BOUND=1,
  TRACE_NO_VALUE=2,
  TRACE_CRITICAL_POINT=3
};

double geos5_haversine(float x0, float y0, float x1, float y1);
void geos5_wrap_x(float *pt);
double geos5_x2lon(float x);
double geos5_y2lat(float y);

int trace_2D_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  );

int trace_3D_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  );

int trace_4D_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  );

int trace_3D_isabel_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  );

int trace_4D_isabel_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  );

int trace_4D_isabel_rk4(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  );

int trace_3D_geos5_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    bool wrap,
    float h
  );

int trace_4D_geos5_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    bool wrap,
    float h
  );

int velocity_4D_geos5(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    const float *X,
    float *v
  );

int trace_4D_geos5_rk4(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    bool wrap,
    float h
  );

#endif
