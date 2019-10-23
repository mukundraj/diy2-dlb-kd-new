#ifndef UTILS_HPP
#define UTILS_HPP

#include <math.h>
#include "diy/link.hpp"
#include <cstdio>

// This utility is the same as diy's pick.hpp, but ensures that distance computation is
// done in double precision even though the bounds are integer
//
// This is a safety measure to reduce the chance of numerical errors and should be used
namespace utl{

    // Find the distance between point `p` and box `bounds`
    template<class Point, class Bounds>
        double
        distance(const Bounds& bounds, const Point& p)
        {
            double res = 0;
            // for (int i = 0; i < p.size(); ++i)
            for (int i = 0; i < 3; ++i)
            {
                // avoids all the annoying case logic by finding
                // diff = max(bounds.min[i] - p[i], 0, p[i] - bounds.max[i])
                double diff = 0, d;

                d = (double)(bounds.min[i]) - (double)(p[i]);
                if (d > diff) diff = d;
                d = (double)(p[i]) - (double)(bounds.max[i]);
                if (d > diff) diff = d;

                // TP, 10/10/19: I reverted this back to match diy::distance()
//                 if (bounds.max[i] == p[i])
//                     res++;

                fprintf(stderr, " distance %f, (%f) [%f %f] \n", res, p[i], bounds.min[i], bounds.max[i]);

                res += diff*diff;
            }

            return sqrt(res);
        }


    // // Finds the neighbor(s) containing the target point
    // template<class Bounds, class Point, class OutIter>
    //     void
    //     in(
    //         const diy::RegularLink<Bounds>& link,   // neighbors
    //         const Point&                    p,      // target point
    //         OutIter                         out,    // insert iterator for output set of neighbors
    //         const Bounds&                   domain, // global domain bounds
    //         bool                            core)   // check against core (or bounds, if false)
    //     {
    //         Bounds neigh_bounds {0}; // neighbor block bounds

    //         // for all neighbors of this block
    //         for (int n = 0; n < link.size(); n++)
    //         {
    //             if (core)
    //                 neigh_bounds = link.core(n);
    //             else
    //                 neigh_bounds = link.bounds(n);

    //             // wrap neighbor bounds, if necessary, otherwise bounds will be unchanged
    //             wrap_bounds(neigh_bounds, link.wrap(n), domain);

    //             if (utl::distance(neigh_bounds, p) == 0)
    //                 *out++ = n;
    //         } // for all neighbors
    //     }


           // Finds the neighbor(s) containing the target point
    template<class Bounds, class Point, class OutIter>
        void
        in(
            const Bounds& neigh_bounds,             // neighbor bounds
            const Point&                    p,      // target point
            OutIter                         out,    // insert iterator for output set of neighbors
            // const Bounds&                   domain, // global domain bounds
            bool                            core)   // check against core (or bounds, if false)
        {
            // Bounds neigh_bounds {0}; // neighbor block bounds
            // fprintf(stderr, "linksizeee %d\n", link.size());

            // for all neighbors of this block
            for (int n = 0; n < 3; n++)
            {
                fprintf(stderr, "herenoww\n");
                // if (core)
                //     neigh_bounds = link.bounds(n);
                // else
                //     neigh_bounds = link.bounds(n);

                // wrap neighbor bounds, if necessary, otherwise bounds will be unchanged
                // wrap_bounds(neigh_bounds, link.wrap(n), domain);

                if (utl::distance(neigh_bounds, p) == 0)
                    *out++ = n;
            } // for all neighbors
        }
}

#endif
