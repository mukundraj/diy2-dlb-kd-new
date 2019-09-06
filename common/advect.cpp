#include "advect.h"
#include "lerp.h"
#include <cstring>
#include <sstream>
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>

using namespace std;

double geos5_haversine(float x0, float y0, float x1, float y1) 
{
  double lon0 = geos5_x2lon(x0), 
         lat0 = geos5_y2lat(y0), 
         lon1 = geos5_x2lon(x1), 
         lat1 = geos5_y2lat(y1); 

  const double R = 6371000.0; // earth radius in meters. 
  double dlat = (lat1 - lat0) * M_PI / 180.0, 
         dlon = (lon1 - lon0) * M_PI / 180.0; 
  double a = pow(sin(dlat/2), 2.0) + cos(lat0*M_PI/180.0)*cos(lat1*M_PI/180.0)*pow(sin(dlon/2), 2.0); 
  double c = 2.0 * atan2(sqrt(a), sqrt(1-a)); 
#if 0
  if (std::isnan(c) || !isfinite(c)) {
    fprintf(stderr, "x0=%f, y0=%f, x1=%f, y1=%f, c=%f\n", x0, y0, x1, y1, c); 
    assert(false);
  }
#endif
  return R * c; 
}

void geos5_wrap_x(float *pt)
{
  const int geos5_domain[3] = {288, 181, 72}; 
  if (pt[0] > geos5_domain[0] - 1) pt[0] -= geos5_domain[0] - 1; 
  else if (pt[0] < 0) pt[0] += geos5_domain[0] - 1; 
}

double geos5_x2lon(float x)
{
  const int geos5_domain[3] = {288, 181, 72}; 
  double lon = x * (360.0 / (geos5_domain[0] - 1)) - 180.0; 
  lon = (lon > 180)  ? 180.0  : lon; 
  lon = (lon < -180) ? -180.0 : lon; 
  return lon; 
}

double geos5_y2lat(float y) 
{
  double lat = y - 90.0; 
  lat = (lat > 90)  ? 90.0  : lat; 
  lat = (lat < -90) ? -90.0 : lat; 
  return lat; 
}

int trace_2D_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  )
{
  if (!inside_st_sz(2, gst, gsz, X)) return TRACE_OUT_OF_BOUND;

  float v[2];
  //if (!lerp3D(X, gst, gsz, 3, vec, v))
  if (!lerp2D(X, gst, gsz, lst, lsz, 2, vec, v))
    return TRACE_OUT_OF_BOUND;

  if (v[0]*v[0] + v[1]*v[1] == 0)
    return TRACE_CRITICAL_POINT;

  for (int i = 0; i < 2; i++)
    X[i] += 0.5 * h * v[i];

  return TRACE_SUCC;
}

int trace_3D_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  )
{
  if (!inside_st_sz(3, gst, gsz, X)) return TRACE_OUT_OF_BOUND;

  float v[3];
  //if (!lerp3D(X, gst, gsz, 3, vec, v))
  if (!lerp3D(X, gst, gsz, lst, lsz, 3, vec, v))
    return TRACE_OUT_OF_BOUND;

  if (v[0]*v[0] + v[1]*v[1] + v[2]*v[2] == 0)
    return TRACE_CRITICAL_POINT;

  for (int i = 0; i < 3; i++)
    X[i] += 0.5 * h * v[i];

  return TRACE_SUCC;
}

int trace_3D_rk1_core(
    const float *clb,
    const float *cub,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  )
{
  // if (!inside_st_sz(3, gst, gsz, X)) return TRACE_OUT_OF_BOUND;
  if (!inside_clb_cub(3, clb, cub, X)) return TRACE_OUT_OF_BOUND;

  float v[3];
  //if (!lerp3D(X, gst, gsz, 3, vec, v))
  if (!lerp3D_core(X, clb, cub, lst, lsz, 3, vec, v))
    return TRACE_OUT_OF_BOUND;

  if (v[0]*v[0] + v[1]*v[1] + v[2]*v[2] == 0)
    return TRACE_CRITICAL_POINT;

  for (int i = 0; i < 3; i++)
    X[i] += 0.5 * h * v[i];

  return TRACE_SUCC;
}

int trace_4D_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  )
{
  if (!inside_st_sz(4, gst, gsz, X)) return TRACE_OUT_OF_BOUND;

  float v[3];
  //if (!lerp4D(X, gst, gsz, 3, vec, v))
  if (!lerp4D(X, gst, gsz, lst, lsz, 3, vec, v))
    return TRACE_OUT_OF_BOUND;

  X[0] = X[0] + h*v[0];
  X[1] = X[1] + h*v[1];
  X[2] = X[2] + h*v[2];
  X[3] = X[3] + h; // TODO

  return TRACE_SUCC;
}

int trace_3D_isabel_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  )
{
  if (!inside_st_sz(3, gst, gsz, X)) return TRACE_OUT_OF_BOUND;

  float grid2meter[3] = {1940.f / 500.f * 1000.f, 2004.f / 500.f * 1000.f, 19.8 / 100.f * 1000.f};

  float v[3];
  //if (!lerp3D(X, gst, gsz, 3, vec, v))
  if (!lerp3D(X, gst, gsz, lst, lsz, 3, vec, v))
    return TRACE_OUT_OF_BOUND;

  if (v[0]*v[0] + v[1]*v[1] + v[2]*v[2] == 0)
    return TRACE_CRITICAL_POINT;

  float k[3] = {h * v[0] / grid2meter[0], h * v[1] / grid2meter[1], h * v[2] / grid2meter[2]};
  for (int i = 0; i < 3; i++)
    X[i] += 0.5 * k[i];

  return TRACE_SUCC;
}

int trace_4D_isabel_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  )
{
  if (!inside_st_sz(4, gst, gsz, X)) return TRACE_OUT_OF_BOUND;

  double TIME_SCALE = 3600.f; 

  float grid2meter[3] = {1940.f / 500.f * 1000.f, 2004.f / 500.f * 1000.f, 19.8 / 100.f * 1000.f};

  if (h < 0) TIME_SCALE = -TIME_SCALE;

  double distance;
  float p0[4] = {X[0], X[1], X[2], X[3]};
  float v[3];

  if (!lerp4D(X, gst, gsz, lst, lsz, 3, vec, v))
    return TRACE_OUT_OF_BOUND;

  float mag = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (std::isnan(mag) || !std::isfinite(mag)) return TRACE_NO_VALUE;

  float k1[3] = {h * v[0] / grid2meter[0], h * v[1] / grid2meter[1], h * v[2] / grid2meter[2]};
  for (int i = 0; i < 3; i++)
    X[i] = p0[i] + 0.5 * k1[i];
  distance = sqrt(pow((p0[0] - X[0]) * grid2meter[0], 2) + pow((p0[1] - X[1]) * grid2meter[1], 2) + pow((p0[2] - X[2]) * grid2meter[2], 2));
  X[3] = p0[3] + distance / mag / TIME_SCALE;

  return TRACE_SUCC;
}

int trace_4D_isabel_rk4(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    float h
  )
{
  if (!inside_st_sz(4, gst, gsz, X)) return TRACE_OUT_OF_BOUND;

  double TIME_SCALE = 3600.f; 

  float grid2meter[3] = {1940.f / 500.f * 1000.f, 2004.f / 500.f * 1000.f, 19.8 / 100.f * 1000.f};

  if (h < 0) TIME_SCALE = -TIME_SCALE;

  double distance;
  float p0[4] = {X[0], X[1], X[2], X[3]};
  float v[3];

  // rk1
  float mag;
  while (1) {
    if (!lerp4D(X, gst, gsz, lst, lsz, 3, vec, v)) return TRACE_OUT_OF_BOUND;
    mag = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (mag < 1e-5) { // critical point
      // printf("critical point.\n");
      X[3] += 0.01;
      p0[3] += 0.01;
    } else break;
    // fprintf(stderr, "after all\n");
  }

  float k1[3] = {h * v[0] / grid2meter[0], h * v[1] / grid2meter[1], h * v[2] / grid2meter[2]};
  for (int i = 0; i < 3; ++ i) {
    X[i] = p0[i] + 0.5 * k1[i];
    if (std::isnan(X[i]) || !std::isfinite(X[i])) return TRACE_NO_VALUE;
  }
  distance = sqrt(pow((p0[0] - X[0]) * grid2meter[0], 2) + pow((p0[1] - X[1]) * grid2meter[1], 2) + pow((p0[2] - X[2]) * grid2meter[2], 2));
  if (std::isnan(distance) || !std::isfinite(distance)) return TRACE_NO_VALUE; //delta_t = 0.1;
  X[3] = p0[3] + distance / mag / TIME_SCALE;

  // rk2
  if (!lerp4D(X, gst, gsz, lst, lsz, 3, vec, v)) return TRACE_OUT_OF_BOUND;
  float k2[3] = {h * v[0] / grid2meter[0], h * v[1] / grid2meter[1], h * v[2] / grid2meter[2]};
  for (int i = 0; i < 3; ++ i) {
    X[i] = p0[i] + 0.5 * k2[i];
    if (std::isnan(X[i]) || !std::isfinite(X[i])) return TRACE_NO_VALUE;
  }
  distance = sqrt(pow((p0[0] - X[0]) * grid2meter[0], 2) + pow((p0[1] - X[1]) * grid2meter[1], 2) + pow((p0[2] - X[2]) * grid2meter[2], 2));
  if (std::isnan(distance) || !std::isfinite(distance)) return TRACE_NO_VALUE; //delta_t = 0.1;
  X[3] = p0[3] + distance / mag / TIME_SCALE;

  // rk3
  if (!lerp4D(X, gst, gsz, lst, lsz, 3, vec, v)) return TRACE_OUT_OF_BOUND;
  float k3[3] = {h * v[0] / grid2meter[0], h * v[1] / grid2meter[1], h * v[2] / grid2meter[2]};
  for (int i = 0; i < 3; ++ i) {
    X[i] = p0[i] + 0.5 * k3[i];
    if (std::isnan(X[i]) || !std::isfinite(X[i])) return TRACE_NO_VALUE;
  }
  distance = sqrt(pow((p0[0] - X[0]) * grid2meter[0], 2) + pow((p0[1] - X[1]) * grid2meter[1], 2) + pow((p0[2] - X[2]) * grid2meter[2], 2));
  if (std::isnan(distance) || !std::isfinite(distance)) return TRACE_NO_VALUE; //delta_t = 0.1;
  X[3] = p0[3] + distance / mag / TIME_SCALE;

  // rk4
  if (!lerp4D(X, gst, gsz, lst, lsz, 3, vec, v)) return TRACE_OUT_OF_BOUND;
  for (int i = 0; i < 3; ++ i) {
    X[i] = p0[i] + (k1[i] + 2.0 * (k2[i] + k3[i]) + h * v[i] / grid2meter[i]) / 6.0;
    if (std::isnan(X[i]) || !std::isfinite(X[i])) return TRACE_NO_VALUE;
  }
  distance = sqrt(pow((p0[0] - X[0]) * grid2meter[0], 2) + pow((p0[1] - X[1]) * grid2meter[1], 2) + pow((p0[2] - X[2]) * grid2meter[2], 2));
  if (std::isnan(distance) || !std::isfinite(distance)) return TRACE_NO_VALUE; //delta_t = 0.1;
  X[3] = p0[3] + distance / mag / TIME_SCALE;

  return TRACE_SUCC;
}

int trace_3D_geos5_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    bool wrap,
    float h
  )
{
  if (!inside_st_sz(3, gst, gsz, X)) return TRACE_OUT_OF_BOUND;

  const double DEG2RAD = M_PI/180.0; 

  double cos_lat = cos(geos5_y2lat(X[1]) * DEG2RAD) + 0.01; 
  //double distance; 
  //float p0[3] = {X[0], X[1], X[2]}; 
  float v[3]; 

  float vals[3];
  if (!lerp3D(X, gst, gsz, lst, lsz, 4, vec, vals))
    return TRACE_OUT_OF_BOUND;

  if (vals[0]*vals[0] + vals[1]*vals[1] + vals[2]*vals[2] == 0)
    return TRACE_CRITICAL_POINT;

  float above_p[3], below_p[3]; 
  float above_zl, below_zl; 
  memcpy(above_p, X, sizeof(float)*3); 
  memcpy(below_p, X, sizeof(float)*3); 
  above_p[2] -= 1.f; 
  below_p[2] += 1.f; 

  float zl = vals[3]; 

  bool above_succ = lerp3D(above_p, gst, gsz, lst, lsz, 1, &vec[3], &above_zl), 
       below_succ = lerp3D(below_p, gst, gsz, lst, lsz, 1, &vec[3], &below_zl); 
  if (!above_succ) 
    above_zl = zl + (zl - below_zl); 
  if (!below_succ)
    below_zl = zl - (above_zl - zl); 

  float thickness = (above_zl - below_zl)*0.5f + (zl - below_zl)*0.5f; 
  float omega = vals[2]; 

  v[0] = vals[0]; 
  v[1] = vals[1]; 
  v[2] = omega*1000.f*44*thickness/111000.f; 

  //double mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); 

  float k1[3] = {h*v[0], h*v[1]*(float)cos_lat, h*v[2]}; 
  for (int i=0; i<3; i++) { 
    X[i] += 0.5 * k1[i]; 
    //if (isnan(X[i]) || !isfinite(X[i])) return TRACE_OUT_OF_BOUNDS;
  }
  //if (!inside_domain(X)) return TRACE_OUT_OF_BOUNDS;
  if (wrap) geos5_wrap_x(X); 
  //distance = geos5_haversine(p0[0], p0[1], X[0], X[1]); 

  return TRACE_SUCC; 
}

int trace_4D_geos5_rk1(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    bool wrap,
    float h
  )
{
  if (!inside_st_sz(4, gst, gsz, X)) return TRACE_OUT_OF_BOUND;

  double TIME_SCALE = 2629744.f; // seconds in a month
  const double DEG2RAD = M_PI/180.0; 

  if (h < 0) TIME_SCALE = -TIME_SCALE;

  double cos_lat = cos(geos5_y2lat(X[1]) * DEG2RAD) + 0.01; 
  double distance; 
  float p0[4] = {X[0], X[1], X[2], X[3]}; 
  float v[4]; 

  float vals[4];
  if (!lerp4D(X, gst, gsz, lst, lsz, 4, vec, vals))
    return TRACE_OUT_OF_BOUND;

  float above_p[4], below_p[4]; 
  float above_zl, below_zl; 
  memcpy(above_p, X, sizeof(float)*4); 
  memcpy(below_p, X, sizeof(float)*4); 
  above_p[2] -= 1.f; 
  below_p[2] += 1.f; 

  float zl = vals[3]; 

  bool above_succ = lerp4D(above_p, gst, gsz, lst, lsz, 1, &vec[3], &above_zl), 
       below_succ = lerp4D(below_p, gst, gsz, lst, lsz, 1, &vec[3], &below_zl); 
  if (!above_succ) 
    above_zl = zl + (zl - below_zl); 
  if (!below_succ)
    below_zl = zl - (above_zl - zl); 

  float thickness = (above_zl - below_zl)*0.5f + (zl - below_zl)*0.5f; 
  float omega = vals[2]; 

  v[0] = vals[0]; 
  v[1] = vals[1]; 
  v[2] = omega*1000.f*44*thickness/111000.f; 
  v[3] = 0;//1.f / 2629744.f; // monthly average

  double mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); 

  float k1[3] = {h*v[0], h*v[1]*(float)cos_lat, h*v[2]}; 
  for (int i=0; i<3; i++) { 
    X[i] += 0.5 * k1[i]; 
    //if (isnan(pt[i]) || !isfinite(pt[i])) return TRACE_OUT_OF_BOUNDS;
  }
  //if (!inside_domain(pt)) return TRACE_OUT_OF_BOUNDS;
  if (wrap) geos5_wrap_x(X); 
  distance = geos5_haversine(p0[0], p0[1], X[0], X[1]); 

  if (std::isnan(distance) || !isfinite(distance))
    X[3] += 0.1; //0.00226888f; 
  else
    X[3] += distance / mag / TIME_SCALE; 

  return TRACE_SUCC; 
}

int velocity_4D_geos5(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    const float *X,
    float *v
  ) 
{
  if (!inside_st_sz(4, gst, gsz, X)) return TRACE_OUT_OF_BOUND;

  float vals[4];
  if (!lerp4D(X, gst, gsz, lst, lsz, 4, vec, vals))
    return TRACE_OUT_OF_BOUND;

  float above_p[4], below_p[4]; 
  float above_zl, below_zl; 
  memcpy(above_p, X, sizeof(float)*4); 
  memcpy(below_p, X, sizeof(float)*4); 
  above_p[2] -= 1.f; 
  below_p[2] += 1.f; 

  float zl = vals[3]; 

  bool above_succ = lerp4D(above_p, gst, gsz, lst, lsz, 1, &vec[3], &above_zl), 
       below_succ = lerp4D(below_p, gst, gsz, lst, lsz, 1, &vec[3], &below_zl); 
  if (!above_succ) 
    above_zl = zl + (zl - below_zl); 
  if (!below_succ)
    below_zl = zl - (above_zl - zl); 

  float thickness = (above_zl - below_zl)*0.5f + (zl - below_zl)*0.5f; 
  float omega = vals[2]; 

  v[0] = vals[0]; 
  v[1] = vals[1]; 
  v[2] = omega*1000.f*44*thickness/111000.f; 
  v[3] = 0;//1.f / 2629744.f; // monthly average

  return TRACE_SUCC; 
}

int trace_4D_geos5_rk4(
    const int *gst,
    const int *gsz,
    const int *lst,
    const int *lsz,
    const float **vec,
    float *X,
    bool wrap,
    float h
  )
{
  if (!inside_st_sz(4, gst, gsz, X)) return TRACE_OUT_OF_BOUND;

  double TIME_SCALE = 2629744.f; // seconds in a month
  const double DEG2RAD = M_PI/180.0; 

  if (h < 0) TIME_SCALE = -TIME_SCALE;

  int stat;
  double cos_lat; 
  double distance, delta_t=0.1;
  float p0[4] = {X[0], X[1], X[2], X[3]}; 
  float v[4]; 

  // rk1
  float mag;
  while (true) {
    stat = velocity_4D_geos5(gst, gsz, lst, lsz, vec, X, v);
    if (stat != TRACE_SUCC) return stat;
    mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (mag < 1e-5) { // critical point
      X[3] += 1e-3;
      p0[3] += 1e-3;
    } else break;
  }
  cos_lat = cos(geos5_y2lat(X[1]) * DEG2RAD) + 0.01;
  float k1[3] = {h*v[0], h*v[1]*(float)cos_lat, h*v[2]};
  for (int i=0; i<3; i++) {
    X[i] = p0[i] + 0.5 * k1[i];
    if (std::isnan(X[i]) || !std::isfinite(X[i])) return TRACE_NO_VALUE;
  }
  if (wrap) geos5_wrap_x(X);
  distance = geos5_haversine(p0[0], p0[1], X[0], X[1]);
  if (std::isnan(distance) || !std::isfinite(distance)) {
    //fprintf(stderr, "\n--- 1.revise finite time value, distance=%f, delta_t=%f ---\n", distance, delta_t);
    delta_t = 0.1; // return TRACE_OUT_OF_BOUNDS;
  } else {
    delta_t = distance / mag / TIME_SCALE;
  }
  X[3] = p0[3] + delta_t;

  // rk2
  stat = velocity_4D_geos5(gst, gsz, lst, lsz, vec, X, v);
  if (stat != TRACE_SUCC) return stat;
  cos_lat = cos(geos5_y2lat(X[1]) * DEG2RAD) + 0.01;
  float k2[3] = {h*v[0], h*v[1]*(float)cos_lat, h*v[2]};
  for (int i=0; i<3; i++) {
    X[i] = p0[i] + 0.5 * k2[i];
    if (std::isnan(X[i]) || !std::isfinite(X[i])) return TRACE_NO_VALUE;
  }
  if (wrap) geos5_wrap_x(X);
  distance = geos5_haversine(p0[0], p0[1], X[0], X[1]);
  if (std::isnan(distance) || !std::isfinite(distance)) {
    //fprintf(stderr, "\n--- 1.revise finite time value, distance=%f, delta_t=%f ---\n", distance, delta_t);
    delta_t = 0.1; // return TRACE_OUT_OF_BOUNDS;
  } else {
    delta_t = distance / mag / TIME_SCALE;
  }
  X[3] = p0[3] + delta_t;

  // rk3
  stat = velocity_4D_geos5(gst, gsz, lst, lsz, vec, X, v);
  if (stat != TRACE_SUCC) return stat;
  cos_lat = cos(geos5_y2lat(X[1]) * DEG2RAD) + 0.01;
  float k3[3] = {h*v[0], h*v[1]*(float)cos_lat, h*v[2]};
  for (int i=0; i<3; i++) {
    X[i] = p0[i] + 0.5 * k3[i];
    if (std::isnan(X[i]) || !std::isfinite(X[i])) return TRACE_NO_VALUE;
  }
  if (wrap) geos5_wrap_x(X);
  distance = geos5_haversine(p0[0], p0[1], X[0], X[1]);
  if (std::isnan(distance) || !std::isfinite(distance)) {
      //fprintf(stderr, "\n--- 1.revise finite time value, distance=%f, delta_t=%f ---\n", distance, delta_t);
    delta_t = 0.1; // return TRACE_OUT_OF_BOUNDS;
  } else {
    delta_t = distance / mag / TIME_SCALE;
  }
  X[3] = p0[3] + delta_t;

  // rk4
  stat = velocity_4D_geos5(gst, gsz, lst, lsz, vec, X, v);
  if (stat != TRACE_SUCC) return stat;
  cos_lat = cos(geos5_y2lat(X[1]) * DEG2RAD) + 0.01;
  for (int i=0; i<3; i++) {
    X[i] = p0[i] + (k1[i] + 2.0*(k2[i]+k3[i]) + h*v[i])/6.0;
    if (std::isnan(X[i]) || !std::isfinite(X[i])) return TRACE_NO_VALUE;
  }
  if (wrap) geos5_wrap_x(X);
  distance = geos5_haversine(p0[0], p0[1], X[0], X[1]);
  if (std::isnan(distance) || !std::isfinite(distance)) {
    //fprintf(stderr, "\n--- 1.revise finite time value, distance=%f, delta_t=%f ---\n", distance, delta_t);
    delta_t = 0.1; // return TRACE_OUT_OF_BOUNDS;
  } else {
    delta_t = distance / mag / TIME_SCALE;
  }
  X[3] = p0[3] + delta_t;

  return TRACE_SUCC;
}