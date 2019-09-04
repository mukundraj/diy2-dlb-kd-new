#include "def.h"
#include <iostream>
#include <netcdf.h>
#include <sstream>
#include <iomanip>
#include <vector>


int main(){

	const std::string filename = "/nfs/proj-tpeterka/hguo/data/nek5000.nc";

	const std::string op_filename = "/nfs/proj-tpeterka/mraj/datasets/nek5000_subset.nc";


	int ncid;
    int dimid_x, dimid_y, dimid_z, dimid_t;
    int varid_U, varid_V, varid_W;
    size_t nX, nY, nZ, cnZ, cnY, cnX;

	// Open the file for read access
    NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid));

    NC_SAFE_CALL( nc_inq_dimid(ncid, "x", &dimid_x));
    NC_SAFE_CALL( nc_inq_dimid(ncid, "y", &dimid_y));
    NC_SAFE_CALL( nc_inq_dimid(ncid, "z", &dimid_z));

    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_x, &nX));
    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_y, &nY));
    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimid_z, &nZ));

    // fprintf(stderr, "%d %d %d\n", dimid_x, dimid_y, dimid_z );

    cnX = nX - 128; cnY = nY - 128; cnZ = nZ - 128;

    NC_SAFE_CALL( nc_inq_varid(ncid, "U", &varid_U));
    NC_SAFE_CALL( nc_inq_varid(ncid, "V", &varid_V));
    NC_SAFE_CALL( nc_inq_varid(ncid, "W", &varid_W));


    std::vector<float> U, V, W;
    const size_t start_z_y_x[3] = {64, 64, 64}, size_z_y_x[3] = {cnX, cnX, cnX};
    U.resize(cnZ*cnY*cnX);
    V.resize(cnZ*cnY*cnX);
    W.resize(cnZ*cnY*cnX);

    //   NC_SAFE_CALL(nc_get_vara_float(ncid, varid_U, &U[0]));
    NC_SAFE_CALL( nc_get_vara_float(ncid, varid_U, start_z_y_x, size_z_y_x, &U[0]));
    NC_SAFE_CALL( nc_get_vara_float(ncid, varid_V, start_z_y_x, size_z_y_x, &V[0]));
    NC_SAFE_CALL( nc_get_vara_float(ncid, varid_W, start_z_y_x, size_z_y_x, &W[0]));


    NC_SAFE_CALL( nc_close(ncid));


    NC_SAFE_CALL(nc_create(op_filename.c_str(), NC_CLOBBER, &ncid));

     /* Define dimensions */
    NC_SAFE_CALL(nc_def_dim(ncid, "z", cnZ, &dimid_z));
    NC_SAFE_CALL(nc_def_dim(ncid, "y", cnY, &dimid_y));
    NC_SAFE_CALL(nc_def_dim(ncid, "x", cnX, &dimid_x));

    fprintf(stderr, "%d %d %d\n", dimid_x, dimid_y, dimid_z );
    fprintf(stderr, "%d %d %d\n", varid_U, varid_V, varid_W );
    /* Define variables */

    int ndims = 3;
    int dimids[] = {dimid_z, dimid_y, dimid_x};
    NC_SAFE_CALL(nc_def_var(ncid, "U", NC_FLOAT, ndims, dimids, &varid_U));
    NC_SAFE_CALL(nc_def_var(ncid, "V", NC_FLOAT, ndims, dimids, &varid_V));
    NC_SAFE_CALL(nc_def_var(ncid, "W", NC_FLOAT, ndims, dimids, &varid_W));

    NC_SAFE_CALL(nc_enddef(ncid));

    const size_t start[3] = {0, 0, 0}, size[3] = {cnX, cnX, cnX};
    NC_SAFE_CALL(nc_put_vara_float(ncid, varid_U, start, size, &U[0]));
    NC_SAFE_CALL(nc_put_vara_float(ncid, varid_V, start, size, &V[0]));
    NC_SAFE_CALL(nc_put_vara_float(ncid, varid_W, start, size, &W[0]));

    NC_SAFE_CALL(nc_close(ncid));

	return 0;
}