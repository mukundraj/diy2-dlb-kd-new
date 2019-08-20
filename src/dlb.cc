#include "dlb.h"
#include <cassert>

std::string pt2buf(const int num_dims, const float *pt)
{
  std::string buf((char*)pt, num_dims*sizeof(float));
  return buf;
}

void buf2pt(const int num_dims, const std::string &buf, float *pt)
{
  memcpy(pt, buf.data(), sizeof(float)*num_dims);
}

std::string id2str(const int num_dims, const int *id)
{
  std::string buf((char*)id, num_dims*sizeof(int));
  return buf;
}

void str2id(const int num_dims, const std::string &buf, int *id)
{
  memcpy(id, buf.data(), sizeof(int)*num_dims);
}

void fill_bounds(
    const int& ghost_block_size, 
    const int* domain_minmax, 
    const int* core, 
    int* bounds
  )
{
  float bounds_max, bounds_min;
  float add = (float)ghost_block_size/2. - (float)(core[1]-core[0]+1)/2.;
  assert (add >= 0);
  bounds_max = core[1] + add;
  bounds_min = core[0] - add;

  if (ghost_block_size > (domain_minmax[1] - domain_minmax[0] + 1)) {
    bounds_max = domain_minmax[1];
    bounds_min = domain_minmax[0];
  } else {
    if (bounds_max > domain_minmax[1]) {
      bounds_max = domain_minmax[1];
      bounds_min = bounds_max - ghost_block_size + 1;
    }

    if (bounds_min < domain_minmax[0]) {
      bounds_min = domain_minmax[0];
      bounds_max = bounds_min + ghost_block_size - 1;
    }
  }

  bounds[1] = std::floor(bounds_max);
  bounds[0] = std::floor(bounds_min);
  assert (bounds[1] >= core[1] && bounds[0] <= core[0]);
}
