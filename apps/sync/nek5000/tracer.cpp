#include <cmath>
#include <algorithm>
#include <sstream>
#include <boost/foreach.hpp>
#include "tracer.h"
#include "common/advect.h"
#include "common/lerp.h"
#include "src/utils.hpp"
#include "src/misc.h"

const int max_trace_size = 2048; //2048;
const float stepsize = 1.0;
const int EPOCH_STEPS = 128;
const int NUM_STEPS = 32; // 50 for small data, 200 for 12GB data
// const float pred_step = 10.0;

typedef diy::RegularContinuousLink RCLink;
typedef diy::ContinuousBounds Bounds;

CSyncNekApp::CSyncNekApp()
{
}

CSyncNekApp::~CSyncNekApp()
{
}

void CSyncNekApp::init(int argc, char **argv)
{
	CPTApp_Sync::init(argc, argv);
}

void CSyncNekApp::initialize_particles(Block &b,
									   std::vector<Particle> &particles)
{
#if 1 // for k-d tree paper
	  // {64, 128, 128} -> {256, 256, 256}
	// const int stride[3] = {512, 512, 512};
	// const int stride[3] = {128, 128, 128};
	// const int stride[3] = {4, 4, 4};
	const int stride[3] = {8, 8, 8};
	const float gap[3] = {1.f / (stride[0] - 1) * (float)(domain_size()[0] - 1),
						  1.f / (stride[1] - 1) * (float)(domain_size()[1] - 1),
						  1.f / (stride[2] - 1) * (float)(domain_size()[2] - 1)};

	float i = 0, j = 0, k = 0;
	int ctr = 0;
	while (i <= domain_size()[0] - 1)
	{
		j = 0;
		k = 0;
		while (j <= domain_size()[1] - 1)
		{
			k = 0;
			while (k <= domain_size()[2] - 1)
			{
				const float idx[3] = {i, j, k};
				ctr++;
				if (is_ptinblock(b, idx))
				{
					Particle p;
					p.id = ctr;
					p.home_gid = b.gid;
					p[0] = i;
					p[1] = j;
					p[2] = k;
					particles.push_back(p);
				}
				k += gap[2];
			}
			j += gap[1];
		}
		i += gap[0];
	}
#endif

#if 0 // original
  const int stride[3] = {4, 4, 4};
  for (int i = b.lb[0]; i <= b.ub[0]; i += stride[0])
    for (int j = b.lb[1]; j <= b.ub[1]; j += stride[1])
      for (int k = b.lb[2]; k <= b.ub[2]; k += stride[2]) {
        const int idx[3] = {i, j, k}; 
        const int id = idx2id(idx);
        Particle p;
        p.home_gid = b.gid;
        p[0] = i;  p[1] = j;  p[2] = k;
        particles.push_back(p);
      }
#endif
}

void CSyncNekApp::trace_particles(Block &b,
								  std::vector<Particle> &particles,
								  std::map<int, std::vector<Particle>> &unfinished_particles,
								  std::map<int, std::vector<Particle>> &finished_particles)
{
	const float **vars = (const float **)(b.vars.data());
	int gst[4], gsz[4], lst[4], lsz[4];
	b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

	BOOST_FOREACH (Particle &p, particles)
	{
		int steps = NUM_STEPS;
		while (p.num_steps < max_trace_size)
		{
			int rtn = trace_3D_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
			if (rtn == TRACE_OUT_OF_BOUND)
				break; // out of ghost size
			add_workload();
			if (rtn == TRACE_CRITICAL_POINT || rtn == TRACE_NO_VALUE)
			{
				p.finished = true;
				break;
			}
			p.num_steps++;
			steps--;

			if (steps <= 0)
				break;
		}

		if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
			p.finished = true;

		if (p.finished)
		{
			finished_particles[p.home_gid].push_back(p);
			_local_done++;
		}
		else
		{
			const int dst_gid = bound_gid(pt2gid(p.coords)); // TODO
			unfinished_particles[dst_gid].push_back(p);
		}
	}
}

static int pt2gid_core_old(int id, const float *pt, std::vector<float> &nbr_bounds, std::vector<int> &nbr_gids, int cgid)
{

	int gid = -1;
	for (int i = 0; i < nbr_gids.size(); i++)
	{
		if (pt[0] >= nbr_bounds[i * 6 + 0] && pt[0] < nbr_bounds[i * 6 + 1] &&
			pt[1] >= nbr_bounds[i * 6 + 2] && pt[1] < nbr_bounds[i * 6 + 3] &&
			pt[2] >= nbr_bounds[i * 6 + 4] && pt[2] < nbr_bounds[i * 6 + 5])
		{

			gid = nbr_gids[i];

			// if (id==27)
			//   fprintf(stderr, "calort cgid %d (%f %f %f), %d %d, %f %f, %f %f, %f %f, gid %d\n", cgid, pt[0], pt[1], pt[2], nbr_bounds.size(), nbr_gids.size(), nbr_bounds[i*6+0], nbr_bounds[i*6+1], nbr_bounds[i*6+2], nbr_bounds[i*6+3], nbr_bounds[i*6+4], nbr_bounds[i*6+5], gid);
		}
	}

	int idx = nbr_gids.size();
	if (gid <= 0)
		if (pt[0] >= nbr_bounds[idx * 6 + 0] && pt[0] < nbr_bounds[idx * 6 + 1] &&
			pt[1] >= nbr_bounds[idx * 6 + 2] && pt[1] < nbr_bounds[idx * 6 + 3] &&
			pt[2] >= nbr_bounds[idx * 6 + 4] && pt[2] < nbr_bounds[idx * 6 + 5])
		{

			gid = cgid;
		}

	if (gid <= 0)
	{

		for (int i = 0; i < nbr_gids.size(); i++)
			fprintf(stderr, "calort cgid %d (%f %f %f), %d %d, %f %f, %f %f, %f %f\n", cgid, pt[0], pt[1], pt[2], nbr_bounds.size(), nbr_gids.size(), nbr_bounds[i * 6 + 0], nbr_bounds[i * 6 + 1], nbr_bounds[i * 6 + 2], nbr_bounds[i * 6 + 3], nbr_bounds[i * 6 + 4], nbr_bounds[i * 6 + 5]);
		// fprintf(stderr, "assert fail. gidsize %ld, (%f %f %f), %d , [%f %f] [%f %f] [%f %f], %d\n", nbr_gids.size(), pt[0], pt[1], pt[2], nbr_gids[i], nbr_gids[i*6+0], nbr_bounds[i*6+1], nbr_bounds[i*6+2], nbr_bounds[i*6+3], nbr_bounds[i*6+4], nbr_bounds[i*6+5], 77);
	}
	assert(gid > -1);

	// if (id==14)
	//          fprintf(stderr, "AAAAA pid %d cgid %d (%f %f %f), %d %d, %f %f, %f %f, %f %f gid %d\n", id, cgid, pt[0], pt[1], pt[2], nbr_bounds.size(), nbr_gids.size(), nbr_bounds[idx*6+0], nbr_bounds[idx*6+1], nbr_bounds[idx*6+2], nbr_bounds[idx*6+3], nbr_bounds[idx*6+4], nbr_bounds[idx*6+5], gid);
	return gid;
}

double distance_new(const float *pt, int i, std::vector<float> &nbr_bounds, bool flag, double &resval){
	double res = 0;

	for (int dim=0; dim<3; dim++){
		double diff = 0, d;

		d = (double) nbr_bounds[i*6+dim*2] - (double) pt[dim];
		if (d>diff) diff = d;
		d = (double) pt[dim] - (double) nbr_bounds[i*6+dim*2+1];
		if (d>diff) diff = d;

		res += diff*diff;
	}

	if (flag)
		resval = (double)i;
	else
		resval = 0;

	return sqrt(res);
}


static int pt2gid_core(Block &b, int id, const float *pt, std::vector<float> &nbr_bounds, std::vector<int> &nbr_gids, int cgid)
{

	int gid = -1;

	//   diy::RegularContinuousLink* kdtree_link = static_cast<diy::RegularContinuousLink*>(&b.cons_kdtree_link);

	//   fprintf(stderr, "BBlklinksize %d\n", b.cons_kdtree_link.size());

	// for each neighbor get distance and check if equals zero 

	int in_count = 0;
	int selectedi = -1;
	double resval=1;
	for (int i=0; i<nbr_gids.size(); i++){

		if (distance_new(pt, i, nbr_bounds, false, resval)==0){
			in_count++;
			gid = nbr_gids[i];
			selectedi = i;
		}
		
	}

	int cur_idx = (int) nbr_gids.size();
	double dist;
	dist = distance_new(pt, cur_idx, nbr_bounds, true, resval);
	if (dist == 0 ){
			in_count++;
			in_count *= -1;
			gid = cgid;
			resval = -1*resval;
	}

	// dprint("ingid %d, nbrs %ld, nbrfullsize %ld", gid, nbr_gids.size(), nbr_bounds.size());

	// dprint("ingid %d, nbrs %ld, p (%f %f %f) nbrb [%f %f, %f %f, %f %f] [%f %f, %f %f, %f %f] [%f %f, %f %f, %f %f] [%f %f, %f %f, %f %f], cgid %d, dist %f, rval %f, in_ct %d", gid, nbr_bounds.size(),pt[0], pt[1], pt[2], \
	// nbr_bounds[0], nbr_bounds[1], nbr_bounds[2], nbr_bounds[3], nbr_bounds[4], nbr_bounds[5], \ 
	// nbr_bounds[6], nbr_bounds[7], nbr_bounds[8], nbr_bounds[9], nbr_bounds[10], nbr_bounds[11], \
	// nbr_bounds[12], nbr_bounds[13], nbr_bounds[14], nbr_bounds[15], nbr_bounds[16], nbr_bounds[17], \
	// nbr_bounds[18], nbr_bounds[19], nbr_bounds[20], nbr_bounds[21], nbr_bounds[22], nbr_bounds[23], cgid, dist, resval, in_count);
	// assert(in_count==1); //(only one nbr dist should be zero) 
	return gid;
	
}

void CSyncNekApp::trace_particles_core(Block &b,
									   std::vector<Particle> &particles,
									   std::map<int, std::vector<Particle>> &unfinished_particles,
									   std::map<int, std::vector<Particle>> &finished_particles,
									   const diy::Master::ProxyWithLink &cp)
{
	const float **vars = (const float **)(b.vars.data());
	int gst[4], gsz[4], lst[4], lsz[4];
	float clb[4], cub[4]; // core_start and core_size
	b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);
	b.get_core_st_sz(num_dims(), clb, cub);

	// fprintf(stderr, "(%f %f) (%f %f) (%f %f) \n", clb[0], cub[0],  clb[1], cub[1],  clb[2], cub[2]);

	BOOST_FOREACH (Particle &p, particles)
	{

		int steps = NUM_STEPS;
		int round_steps = 0;

		while (p.num_steps < max_trace_size && !p.epoch_finished)
		{
			// int rtn = trace_3D_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
			int rtn = trace_3D_rk1_core(clb, cub, lst, lsz, vars, p.coords, stepsize);
			if (rtn == TRACE_OUT_OF_BOUND)
			{
				break; // out of core size
			}
			add_workload();
			if (rtn == TRACE_CRITICAL_POINT || rtn == TRACE_NO_VALUE)
			{
				p.finished = true;
				// _local_done_epoch ++;
				break;
			}
			p.num_steps++;
			p.num_esteps++;
			round_steps++;
			steps--;

			// if (p.id==204)
			// 	dprint("stp%d p(%f %f %f)", p.num_steps, p.coords[0], p.coords[1], p.coords[2]);

			if (p.num_esteps == EPOCH_STEPS)
			{   // && !p.epoch_finished){
				// if epoch is done, then kd-tree rebalance, else continue with same partition
				_local_done_epoch++;
				p.epoch_finished = true;
				p.num_esteps = 0;
				break;
			}

			// if (steps <= 0) break;
		}
		// dprint("particle %d ((%f %f %f)), %d %d %d", p.id, p[0], p[1], p[2], p.num_steps, p.finished, cp.gid());

		if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
		{
			p.finished = true;
			// _local_done_epoch++;
		}

		if (p.finished)
		{	
			// dprint("finished pid %d", p.id);
			finished_particles[p.home_gid].push_back(p);
			_local_done++;
			finished_status[p.id]++;
			_local_finished_in_this_round ++;
			if (!p.epoch_finished)
			{
				p.epoch_finished = true;
				_local_done_epoch++;
			}
			// fprintf(stderr, " F (pid (%d) %d %d, %d)\n", p.id,  bp.num_steps, p.finished,.gid);
		}
		else
		{

			// if (p.id==204)
			// 	dprint("piddd %d, %f %f %f, %d %d %d", p.id, p.coords[0], p.coords[1], p.coords[2], p.num_steps, p.finished, b.gid);
			// int dst_gid = bound_gid(pt2gid(p.coords)); // TODO
			// int dst_gid = bound_gid(pt2gid_core(p.id, p.coords, b.nbr_bounds, b.nbr_gids, cp.gid())); // TODO
			
			// int dst_gid = pt2gid_core_old(p.id, p.coords, b.nbr_bounds, b.nbr_gids, cp.gid()); // TODO
			int dst_gid = pt2gid_core(b, p.id, p.coords, b.nbr_bounds, b.nbr_gids, cp.gid()); // TODO

				//  RCLink*  link      = static_cast<RCLink*>(cp.link());
				//  fprintf(stderr, "size %d\n", link->size());
				//  for (int i = 0; i < link->size(); ++i)
				// { fprintf(stderr, "link->bounds().min[0] (%ld %ld)\n", link->bounds().min[0], link->bounds().max[0]);
				//   const Bounds& bounds = link->bounds(i);

				//   if (p.coords[0]>link->bounds().min[0] && p.coords[0] < link->bounds().max[0]
				//     && p.coords[1]>link->bounds().min[1] && p.coords[1] < link->bounds().max[1]
				//     && p.coords[2]>link->bounds().min[2] && p.coords[2] < link->bounds().max[2])
				//       {
				//           dst_gid = link->target(i).gid;
				//           assert(0);

				//       }

				// }

				// RCLink*  link      = static_cast<RCLink*>(&b.cons_kdtree_link);
				// Bounds neigh_bounds {0}; // neighbor block bounds
				// for (int i=0; i<link->size(); i++){
				//   neigh_bounds = link->bounds(i);
				// }

				// std::vector<int> dests;
				// std::vector<int>::iterator it = dests.begin();
				// std::insert_iterator<std::vector<int>> insert_it(dests, it);
				// utl::in(b.core_bounds, p.coords, insert_it, 1);

				// int dst_gidd = link->target(dests[0]).gid;

				// fprintf(stderr, "destssize %ld\n", dests.size());

				// if (p.num_esteps>0)
				//   fprintf(stderr, " UNF (pid (%d) [%d %d %d], %d, %d %d) (%f %f %f)\n", p.id, round_steps, p.num_esteps, p.num_steps, p.finished, b.gid,  dst_gidd, p.coords[0], p.coords[1], p.coords[2]);

			if (p.num_esteps == EPOCH_STEPS)
			{
				// if epoch is done, then kd-tree rebalance, else continue with same partition
				// _local_done_epoch++; // already incremented earlier
				p.num_esteps = 0;
			}

			unfinished_particles[dst_gid].push_back(p);

			// if (b.gid ==7){
			//   p.id=999;
			//     unfinished_particles[6].push_back(p);
			// }
		}
	}
}

#if 0

#else
void CSyncNekApp::trace_particles_kdtree(Block &b,
										 std::vector<Particle> &particles,
										 std::map<int, std::vector<Particle>> &unfinished_particles,
										 std::map<int, std::vector<Particle>> &finished_particles)
{
	const float **vars = (const float **)(b.vars.data());
	int gst[4], gsz[4], lst[4], lsz[4];
	b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

	BOOST_FOREACH (Particle &p, particles)
	{
		int steps = NUM_STEPS;
		while (p.num_steps < max_trace_size)
		{
			int rtn = trace_3D_rk1(gst, gsz, lst, lsz, vars, p.coords, stepsize);
			if (rtn == TRACE_OUT_OF_BOUND)
				break; // out of ghost size
			add_workload();
			if (rtn == TRACE_CRITICAL_POINT || rtn == TRACE_NO_VALUE)
			{
				p.finished = true;
				break;
			}
			p.num_steps++;
			steps--;

			if (steps <= 0)
				break;
		}
		// fprintf(stderr, "%d ", p.wgt);
		_pred_mismatch += abs((NUM_STEPS - steps) - p.wgt) - 1;

		if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
			p.finished = true;

		if (p.finished)
		{
			_local_done_epoch++;
			_local_done++;
		}
		else
		{

			if (p.num_steps % EPOCH_STEPS == 0)
			{
				// if epoch is done, then kd-tree rebalance, else continue with same partition
				_local_done_epoch++;
			}
			b.particles.push_back(p);
		}
	}
}

void CSyncNekApp::trace_particles_kdtree_predict(Block &b, int factor)
{
	float pred_step = factor * stepsize;
	const float **vars = (const float **)(b.vars.data());
	int gst[4], gsz[4], lst[4], lsz[4];
	b.get_ghost_load_st_sz(num_dims(), gst, gsz, lst, lsz);

	BOOST_FOREACH (Particle &p, b.particles)
	{
		int steps = NUM_STEPS;
		int orig_num_steps = p.num_steps;
		p.wgt = 1;
		float coords[4] = {p.coords[0], p.coords[1], p.coords[2], p.coords[3]};

		while (p.num_steps < max_trace_size)
		{
			int rtn = trace_3D_rk1(gst, gsz, lst, lsz, vars, p.coords, pred_step);
			if (rtn == TRACE_OUT_OF_BOUND)
				break; // out of ghost size
			// add_workload();
			if (rtn == TRACE_CRITICAL_POINT || rtn == TRACE_NO_VALUE)
			{
				// p.finished = true;
				break;
			}
			p.num_steps += pred_step;
			p.wgt += pred_step;
			steps -= pred_step;

			if (steps <= 0)
				break;
		}

		p.num_steps = orig_num_steps;
		// p.finished = false;
		p.coords[0] = coords[0];
		p.coords[1] = coords[1];
		p.coords[2] = coords[2];
		p.coords[3] = coords[3];

		//   if (!inside_domain(p.coords) || p.num_steps >= max_trace_size)
		//     p.finished = true;

		//   if (p.finished) {
		//     _local_done ++;
		//   } else {
		//     b.particles.push_back(p);
		//   }
		// fprintf(stderr, "%d ", p.wgt );
	}
}

#endif

void CSyncNekApp::trace_particles_video(Block &b,
										std::vector<Particle> &particles,
										std::vector<Particle> &particles_first,
										std::vector<Particle> &particles_second)
{
	// TODO
}