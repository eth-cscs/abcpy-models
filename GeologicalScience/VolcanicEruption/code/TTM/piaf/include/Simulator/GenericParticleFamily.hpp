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


#ifndef GENERICPARTICLEFAMILY_H_
#define GENERICPARTICLEFAMILY_H_

#include "SimulatorTypes.hpp"

namespace piaf{

/*!
 * \class GenericParticleFamily
 * \brief Abstract class used to represent family of particles
 */
class GenericParticleFamily {
public:

	/*!
	 * \fn GenericParticleFamily( int id )
	 * \brief Constructor of GenericParticleFamily
	 * \param id : ID of the new family, must be unique and all IDs must be contiguous to avoid problems with the rest of the code
	 */
	GenericParticleFamily( int id );

	/*!
	 * \fn ~GenericParticleFamily()
	 * \brief Destructor of GenericParticleFamily (pure virtual method)
	 */
	virtual ~GenericParticleFamily();

	/*! \brief Identify the family, this ID is used as index */
	const int familyId_;

};

} //namespace piaf

#endif /* GENERICPARTICLEFAMILY_H_ */
