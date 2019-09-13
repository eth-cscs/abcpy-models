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


#ifndef VOLCANICTEPHRAPARTICLEFAMILIES_H_
#define VOLCANICTEPHRAPARTICLEFAMILIES_H_

#include <math.h>
#include <boost/math/constants/constants.hpp>

//#include "../../Tools/FormattedLog/Log.hpp"
//#include "../../Simulator/Sim.hpp"

#include "piaf.hpp"

using namespace boost::math::constants;
using namespace std;

using namespace piaf;

/*!
 * \class TephraParticleFamily
 * \brief this class is used to represent families of tephra particles during an eruption
 * the needed computations are made inside the constructor and the class present the values through
 * constants public members
 * grain size (phi) and density are given as input and diameter and mass are computed from phi and density
 */

class TephraParticleFamily:public GenericParticleFamily{
public:

	/*!
	 * \fn TephraParticleFamily( int id, double grainSize, double density, double weightPercent )
	 * \brief constructor of the class TephraParticleFamily
	 * \param id : the id of the family
	 * \param grainSize : the size of the particle family phi ( phi = log2( diameter ) )
	 * \param density : the density of the particles of the family
	 * \param weightPercent : the weight percentage of this family regarding the total erupted mass
	 */
	TephraParticleFamily( int id, double grainSize, double density, double originalWeightPercent ):
		GenericParticleFamily	( id ),
		density_				( density ),
		grainSize_				( grainSize ),
		originalWeightPercent_	( originalWeightPercent ),
		diameter_				( ( pow( 2.0 , -grainSize ) ) / 1000.0 ),
		mass_					( ( 4.0 / 3.0 ) * pi<double>() * ( pow ( diameter_ / 2.0 , 3.0 ) ) * density ){

        //log<LOG_DEBUG>("particle family created - id = %1% - phi = %2% - density = %3% - wtpct = %4%\n"
        //		"                                                        diameter = %5% - mass = %6%") % id % grainSize % density % weightPercent % diameter_ % mass_;

	}

	/*! \brief density of particles */
	const double density_;
	/*! \brief grain size of particles (phi) */
	const double grainSize_;
	/*! \brief the weight percentage of this family regarding the total erupted mass */
	const double originalWeightPercent_;
	/*! \brief particle diameter in meters */
	const double diameter_;
	/*! \brief mass of particles */
	const double mass_;
};

#endif
