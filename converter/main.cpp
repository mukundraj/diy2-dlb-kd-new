#include <algorithm>
#include <stdio.h>
#include <stdlib.h> 
#include <cassert>
#include <string>
#include <limits.h>
//#include <mpi.h>
#include <netcdf.h>

using namespace std;

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

void ReadStaticDataRaw(const char* fname)
{
  int retval;

  FILE * fin = fopen(fname, "rb");
  assert(fin != NULL);

  int n_dims = 3;
  size_t dims[] = { 1024, 1024, 1024, };
  int dimids[3];
  size_t n_elems = dims[0] * dims[1] * dims[2];

  const int N_VARS = 3;
  const char *VARNAMES[] = { "U", "V", "W", };
  int varid[N_VARS];

  for (int k = 0; k < N_VARS; ++k)
  {
    fprintf(stderr, "start k = %d\n", k);
    int ncid;
    string nc_filename = string("/projects/SDAV/jiangz/nek5000_/") + VARNAMES[k] + string(".nc");
    if ((retval = nc_create(nc_filename.c_str(), NC_64BIT_OFFSET, &ncid))) ERR(retval);
    
    if ((retval = nc_def_dim(ncid, "z", dims[0], &dimids[0]))) ERR(retval);
    if ((retval = nc_def_dim(ncid, "y", dims[1], &dimids[1]))) ERR(retval);
    if ((retval = nc_def_dim(ncid, "x", dims[2], &dimids[2]))) ERR(retval);
    if ((retval = nc_def_var(ncid, VARNAMES[k], NC_FLOAT, n_dims, dimids, &varid[k]))) ERR(retval);
    if ((retval = nc_enddef(ncid))) ERR(retval);

    char *buf = new char[n_elems * sizeof(float)];
#if 0
    fread(buf, sizeof(char), n_elems * sizeof(float), fin);
    //float *p = (float*)buf;
#else
    int size = sizeof(float) * N_VARS;
    //for (size_t i = 0; i < n_elems * size; i += size) {
    //  size_t offset = i + k * sizeof(float);
    //  fseek(fin, offset, SEEK_SET);
    //  fread(buf, sizeof(char), sizeof(float), fin);
    //  buf += sizeof(float);
    //}
    int count = 4;
    size_t num = n_elems * sizeof(float) * N_VARS / count;
    size_t n = 0;
    fseek(fin, 0, SEEK_SET);

    for (int i = 0; i < count; i ++) {
      char *temp_buf = new char[num];
      fread(temp_buf, sizeof(char), num, fin);

      for (size_t j = 0; j < num; j += size) {
        buf[n++] = temp_buf[j+ k * sizeof(float)];
        buf[n++] = temp_buf[j+1 + k * sizeof(float)];
        buf[n++] = temp_buf[j+2 + k * sizeof(float)];
        buf[n++] = temp_buf[j+3 + k * sizeof(float)];
      }

      delete []temp_buf;
    }
#endif 

#if 0 // for test
    for (int x = 0; x < dims[0]; ++x)
    {
      for (int y = 0; y < dims[1]; ++y)
      {
        for (int z = 0; z < dims[2]; ++z)
          *p++ = x + y + z;
      }
    }
#endif

    char tmp;
    for (size_t i = 0; i < n_elems * sizeof(float); i += 4)
    {
      char *p = &buf[i];
      tmp = p[0]; p[0] = p[3]; p[3] = tmp;
      tmp = p[1]; p[1] = p[2]; p[2] = tmp;
    }


/*
    float *f_buf = (float*)buf;
    float *r_buf = new float[n_elems];

    int dd = 1024*1024;
    for (size_t n = 0; n < n_elems; n ++) {
      int r_x = n / dd;
      int r_y = (n - r_x*dd) / 1024;
      int r_z = n - r_y*1024 - r_x*dd;

      size_t index = r_x + 1024 *(r_y + 1024*r_z);
      r_buf[index] = f_buf[n];
    }
*/


    if ((retval = nc_put_var_float(ncid, varid[k], (float*)buf))) ERR(retval);
    //if ((retval = nc_put_var_float(ncid, varid[k], r_buf))) ERR(retval);

    delete []buf;
    //delete []f_buf;
    //delete []r_buf;
    if ((retval = nc_close(ncid))) ERR(retval);

    fprintf(stderr, "end k = %d\n", k);
  }
  fclose(fin);
}

int main(int argc, char** argv)
{
  //int dimension[3] = {1024, 1024, 1024};
  string filename = "/projects/SDAV/jiangz/0.xyz";

  ReadStaticDataRaw(filename.c_str());
  return 0;
}