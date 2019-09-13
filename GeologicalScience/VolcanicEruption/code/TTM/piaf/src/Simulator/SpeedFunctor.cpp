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


#include "../../include/Simulator/SpeedFunctor.hpp"

namespace piaf{

/*
 * (virtual) implementation of constructor and destructor
 */

SpeedFunctor::SpeedFunctor( std::string name ):
speedType_	( LAGRANGIAN_DYNAMIC ),
isDiffusion_ (false),
name_ ( name )
{}

SpeedFunctor::~SpeedFunctor() {}

void SpeedFunctor::setSpeedType( ESpeedType sr ){
	speedType_ = sr;
}

ESpeedType SpeedFunctor::getSpeedType(){
	return speedType_;
}

bool SpeedFunctor::isDiffusion(){ return isDiffusion_; }

void SpeedFunctor::setIsDiffusion( bool id ){ isDiffusion_ = id; }

std::string SpeedFunctor::getName(){ return name_; }


} // namespace piaf
