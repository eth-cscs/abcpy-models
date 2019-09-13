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


#ifndef CONSTANTSPEED_HPP_
#define CONSTANTSPEED_HPP_

#include "../Simulator/SpeedFunctor.hpp"

namespace piaf{

class ConstantSpeed: public SpeedFunctor{

public:

	ConstantSpeed( Double3 speed ): SpeedFunctor("piaf_constant"), speed_ ( speed )
	{
		speedType_ = EULERIAN_STATIC;
	};

	Double3 operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family);

	void setSpeed( Double3 speed );

private:

	Double3 speed_;

};

inline void ConstantSpeed::setSpeed( Double3 speed ){
	speed_ = speed;
}

inline Double3 ConstantSpeed::operator()(Double3 pos, double t, double dt, std::shared_ptr<GenericParticleFamily> family){
	return speed_;
}

} // namespace piaf

#endif
