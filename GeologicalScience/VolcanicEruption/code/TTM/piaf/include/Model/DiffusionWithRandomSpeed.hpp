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


#ifndef DIFFUSIONWITHRANDOMSPEED_HPP_
#define DIFFUSIONWITHRANDOMSPEED_HPP_

#include "../Simulator/SpeedFunctor.hpp"
#include <boost/math/constants/constants.hpp>

namespace piaf{

class DiffusionWithRandomSpeed: public SpeedFunctor{

public:

    DiffusionWithRandomSpeed( double d ): SpeedFunctor("piaf_diffusion"), d_ ( d )
	{
		speedType_ = LAGRANGIAN_DYNAMIC;
	};

	Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

	void setD( double d );

private:
	double d_;
};

inline void DiffusionWithRandomSpeed::setD( double d ){
	d_ = d;
}

inline Double3 DiffusionWithRandomSpeed::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
	//Double3 speed;

	// 2D
	//double ur = sqrt( ( 4.0 * d_ ) / dt );
	// 3D
    double ur = sqrt( ( 6.0 * d_ ) / dt ) * double(rand())/double(RAND_MAX);
	double theta = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * boost::math::constants::pi<double>();
	double phi = acos( ( double( rand() ) / double( RAND_MAX ) ) * 2.0 - 1.0 );

	// 2d
	//return { ur * cos( theta ), ur * sin( theta ), 0.0 };

	// 3d
	return { ur * cos( theta ) * sin( phi ), ur * sin( theta ) * sin( phi ), ur * cos( phi ) };

}

} // namespace piaf

#endif
