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


#ifndef DIRECTIONALDIFFUSION_HPP_
#define DIRECTIONALDIFFUSION_HPP_

#include "../Simulator/SpeedFunctor.hpp"
#include <boost/math/constants/constants.hpp>

namespace piaf{

class DirectionalDiffusion: public SpeedFunctor{

public:

	DirectionalDiffusion( double dV, double dH ): SpeedFunctor("piaf_anisotropic_diffusion"), dV_ ( dV ), dH_ ( dH )
	{
		speedType_ = LAGRANGIAN_DYNAMIC;
	};

	Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

	void setDV( double dV );
	void setDH( double dH );

private:
	double dV_;
	double dH_;
};

inline void DirectionalDiffusion::setDV( double dV ){
	dV_ = dV;
}

inline void DirectionalDiffusion::setDH( double dH ){
	dH_ = dH;
}

inline Double3 DirectionalDiffusion::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
	//Double3 speed;

	//double urV = 1.8*sqrt( ( 2.0 * dV_ ) / dt );
	double urV = sqrt( ( 2.0 * dV_ ) / dt );

	if( ( double( rand() ) / double( RAND_MAX ) ) > 0.5 ) urV = -urV;

	//urV = double( rand() ) / double( RAND_MAX ) * urV;


    //double alpha = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * boost::math::constants::pi<double>();
	//urV = urV * sin(alpha) + urV * cos(alpha);
	//urV = urV * (sin(alpha)+cos(alpha));
	//urV = urV * (sin(alpha));


	double urH = sqrt( ( 4.0 * dH_ ) / dt );
	double theta = ( double( rand() ) / double( RAND_MAX ) ) * 2.0 * boost::math::constants::pi<double>();

	return { urH * cos( theta ), urH * sin( theta ), urV };

}

} // namespace piaf

#endif
