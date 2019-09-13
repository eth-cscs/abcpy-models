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


#include "../../include/Simulator/AbstractEulerianDomain.hpp"

namespace piaf{

AbstractEulerianDomain::AbstractEulerianDomain( Double3 pos, Double3 size, double dx ):
    dx_					( dx ),
    size_				( { int( size.x_ / dx_ ) , int( size.y_ / dx_ ) , int( size.z_ / dx_ ) } ),
    position_			( pos ),
    containsParticles_	( 0 ),
    staticSpeedsUpToDate_( true )
{}

AbstractEulerianDomain::AbstractEulerianDomain():
    dx_					( 0.0 ),
    size_				( { 0 , 0 , 0 } ),
    position_			( { 0.0 , 0.0 , 0.0 } ),
    containsParticles_	( 0 ),
    staticSpeedsUpToDate_( true )
{}

AbstractEulerianDomain::~AbstractEulerianDomain() {}

double AbstractEulerianDomain::getDx(){ return dx_; }

int AbstractEulerianDomain::getXSize(){ return size_.x_; }

int AbstractEulerianDomain::getYSize(){ return size_.y_; }

int AbstractEulerianDomain::getZSize(){ return size_.z_; }

Int3 AbstractEulerianDomain::getSize(){ return size_; }

Double3 AbstractEulerianDomain::getPosition(){ return position_; }

void AbstractEulerianDomain::addSpeed( std::shared_ptr< SpeedFunctor > newSpeed ){
    staticSpeedsUpToDate_ = false;
}

void AbstractEulerianDomain::removeSpeed( std::string name ){
    staticSpeedsUpToDate_ = false;
}

void AbstractEulerianDomain::computeSpeeds( double t, double dt ){
    if( !staticSpeedsUpToDate_ ){
        computeStaticSpeeds();
        staticSpeedsUpToDate_ = true;
    }
}

void AbstractEulerianDomain::addParticleFamily( std::shared_ptr< GenericParticleFamily > newFamily ){
    staticSpeedsUpToDate_ = false;
}

} // namespace piaf
