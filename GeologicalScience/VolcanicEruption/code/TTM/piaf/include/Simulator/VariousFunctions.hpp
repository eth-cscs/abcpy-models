/*
Particles in Advection Field (PIAF)
Copyright (C) 2015  University of Geneva, Switzerland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef VARIOUSFUNCTIONS_H_
#define VARIOUSFUNCTIONS_H_

#include <cmath>
#include "gsl/gsl_interp.h"
#include "mpfr.h"
// #include <vector>

namespace piaf{

/* \brief compute the relative difference between two numbers */
inline double relativeDif(double x, double y){
    if( x==0.0 && y==0.0 ) return 0.0;
    return fabs(x-y)/std::max(fabs(x),fabs(y));
}

/* \brief compute the distance between two points on a plane (take into account only x and y coordinates)*/
inline double dist2D( Double3 p1, Double3 p2 ){
    return sqrt( pow( p2.x_ - p1.x_, 2.0 ) + pow( p2.y_ - p1.y_, 2.0) );
}

/* \brief compute the distance between two points*/
inline double dist3D( Double3 p1, Double3 p2 ){
    return sqrt( pow( p2.x_ - p1.x_, 2.0 ) + pow( p2.y_ - p1.y_, 2.0) + pow( p2.z_ - p1.z_, 2.0) );
}

/* \brief test if p is located inside the box defined by p1, p2 */
inline bool isInside( Double3 p, Double3 p1, Double3 p2 ){
    double startX = std::min(p1.x_, p2.x_); double startY = std::min(p1.y_, p2.y_); double startZ = std::min(p1.z_, p2.z_);
    double endX = std::max(p1.x_, p2.x_); double endY = std::max(p1.y_, p2.y_); double endZ = std::max(p1.z_, p2.z_);
    if( p.x_ >= startX && p.x_ <= endX && p.y_ >= startY && p.y_ <= endY && p.z_ >= startZ && p.z_ <= endZ ) return true;
    return false;
}

/* \brief compute the distribution of particles in number of particles per m3 */
inline std::vector<double> computeDistribution( std::vector<Particle> particles, double volume, std::vector<double> scalingFactors){
    std::vector<double> res(scalingFactors.size(), 0);

    for( auto p : particles ){
        // if( p.familyId_+1 > res.size() ) res.resize( p.familyId_+1 , 0 );
        res[p.familyId_] += scalingFactors[p.familyId_];
        // res[p.familyId_] += 1;
    }

    for( auto& v : res ) v /= volume;

    return res;
}

inline bool coordInSphere(Double3 coord, Double3 sphereCenter, double sphereRadius){
    if( dist3D(coord, sphereCenter) > sphereRadius ) return false;
    return true;
}

// precise sum of vector
inline double preciseSum(std::vector<double> v)
{
    mpfr_t res_mpfr;
    mpfr_init2 (res_mpfr, 200);
    mpfr_set_d (res_mpfr, 0.0, MPFR_RNDD);

    mpfr_t e_mpfr;
    mpfr_init2 (e_mpfr, 200);

    for (auto e : v){
        mpfr_set_d (e_mpfr, e, MPFR_RNDD);
        mpfr_add (res_mpfr, res_mpfr, e_mpfr, MPFR_RNDD);
    }

    return mpfr_get_d (res_mpfr, MPFR_RNDD);
}

// normalize v to make its sum = targetSum
inline std::vector<double> preciseNormalize(std::vector<double> v, double targetsum){
    mpfr_t e_mpfr, sum_mpfr, targetsum_mpfr, sumdivtargetsum_mpfr;
    mpfr_init2 (e_mpfr, 200);
    mpfr_init2 (sum_mpfr, 200);
    mpfr_init2 (targetsum_mpfr, 200);
    mpfr_init2 (sumdivtargetsum_mpfr, 200);
    mpfr_set_d (sum_mpfr, preciseSum(v), MPFR_RNDD);
    mpfr_set_d (targetsum_mpfr, targetsum, MPFR_RNDD);
    mpfr_div (sumdivtargetsum_mpfr, sum_mpfr, targetsum_mpfr, MPFR_RNDD);

    for (auto &e : v){
        mpfr_set_d (e_mpfr, e, MPFR_RNDD);
        mpfr_div (e_mpfr, e_mpfr, sumdivtargetsum_mpfr, MPFR_RNDD);
        e = mpfr_get_d (e_mpfr, MPFR_RNDD);
    }

    return v;
}

/* \brief create a vector of linearly spaced values, ranging from start to stop, with n values */
inline std::vector<double> createVec(double start, double stop, uint n){
    std::vector<double> res(n);
    int i = 0;
    for (double v = start; v <= stop; v += (stop - start) / (n - 1))
    {
        res[i] = v;
        i++;
    }
    if (res.back() != stop)
        res[res.size() - 1] = stop;
    return res;
}

/* \brief interpolate data points xs, ys with new number of points n */
inline std::vector<std::vector<double>> interpolate(std::vector<double> xs, std::vector<double> ys, uint n){
    auto xs_interp = createVec(xs.front(), xs.back(), n);
    std::vector<double> ys_interp(xs_interp.size());

    auto interp = gsl_interp_alloc(gsl_interp_linear, xs.size());
    gsl_interp_init(interp, xs.data(), ys.data(), xs.size());
    auto acc = gsl_interp_accel_alloc();

    for (uint i = 0; i < ys_interp.size(); i++)
    {
        ys_interp[i] = gsl_interp_eval(interp, xs.data(), ys.data(), xs_interp[i], acc);
    }

    gsl_interp_free(interp);
    gsl_interp_accel_free(acc);

    return std::vector<std::vector<double>>{xs_interp, ys_interp};
}

}

#endif
