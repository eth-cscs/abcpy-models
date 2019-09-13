/*
TEphra TRAsport Simulator (TETRAS)
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


#include "Atmosphere.hpp"
#include <math.h>

Atmosphere::Atmosphere( double tropopauseHeight, double topSteadyTemperatureHeight, double T0, double P0 ):
tropopauseHeight_			( tropopauseHeight ),
topSteadyTemperatureHeight_	( topSteadyTemperatureHeight ),
T0_							( T0 ),
P0_							( P0 ),
Ra_							( 285.0 ),
// values are given k/1000m but we need k/m
omegaT_						( ( 6.5 / 1000.0 ) ),
omegaS_						( ( 2.0 / 1000.0 ) )
{
	assert( tropopauseHeight_ <= topSteadyTemperatureHeight_ );
}

Atmosphere::~Atmosphere() {}

double Atmosphere::getTropopauseHeight(){return tropopauseHeight_;}

double Atmosphere::getTopSteadyHeight(){return topSteadyTemperatureHeight_;}

double Atmosphere::getT0(){return T0_;}

double Atmosphere::getP0(){return P0_;}

double Atmosphere::getRa(){return Ra_;}

double Atmosphere::getT(double z){
	double T;

	if( z <= tropopauseHeight_ ){
		T = T0_ - omegaT_ * z;
	}
	else if( z > tropopauseHeight_ && z <= topSteadyTemperatureHeight_ ){
		T = T0_ - omegaT_ * tropopauseHeight_;
	}
	else{
		T = T0_ - omegaT_ * tropopauseHeight_ + omegaS_ * ( z - topSteadyTemperatureHeight_ );
	}

	return T;
}

double Atmosphere::getP( double z ){
	double P;
	double exponent;
	double t1;
	double t2;
	double t3;

	exponent = 9.81 / ( Ra_ * omegaT_ );
	t1 = ( -9.81 * ( z - tropopauseHeight_ ) ) / ( Ra_ * ( T0_ - omegaT_ * tropopauseHeight_ ) );
	t2 = ( -9.81 * ( topSteadyTemperatureHeight_ - tropopauseHeight_ ) ) / ( Ra_ * ( T0_ - omegaT_ * tropopauseHeight_ ) );
	t3 = ( pow( ( ( T0_ - omegaT_ * tropopauseHeight_ + omegaS_ * ( z - topSteadyTemperatureHeight_ ) ) / ( T0_ - omegaT_ * tropopauseHeight_ ) ) , ( -9.81 / ( Ra_ * omegaS_ ) ) ) );

	if( z <= tropopauseHeight_ ){
		P = P0_ * ( pow( ( ( T0_ - omegaT_ * z ) / T0_ ) , exponent ) );
	}
	else if( z > tropopauseHeight_ && z <= topSteadyTemperatureHeight_ ){
		P = P0_ * ( pow( ( ( T0_ - omegaT_ * tropopauseHeight_ ) ) / T0_ , exponent ) ) * exp( t1 );
	}
	else{
		P = P0_ * ( pow( ( ( T0_ - omegaT_ * tropopauseHeight_ ) ) / T0_ , exponent ) ) * exp( t2 ) * t3;
	}

	return P;
}

double Atmosphere::getAlpha( double T, double P ){
	return P / ( Ra_ * T );
}
