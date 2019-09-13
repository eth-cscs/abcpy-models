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


#ifndef PARTICLEREPOSITORY_H_
#define PARTICLEREPOSITORY_H_

#include <memory>
#include "SimulatorTypes.hpp"
#include "vector"
#include "GenericParticleFamily.hpp"
//#include <boost/shared_ptr.hpp>
#include <iostream>

//using namespace piaf;

namespace piaf{

class ParticleRepository {

protected :

	/*! \brief The list of stored particles */
	std::vector< Particle > storedParticles_;

	std::vector< BoundaryParticle > storedBoundaryParticles_;

	/* \brief Tells if the repository must store particles */
	bool isActive_;

	/*! \brief List of Particle families */
	std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies_;

public:

	/*!
	 * \fn ParticleRepository( std::vector< std::shared_ptr< GenericParticleFamily> > particleFamilies )
	 * \brief constructor of the class ParticleRepository
	 * \param particleFamilies the list of particle families that will be stored
	 */
	ParticleRepository( std::vector< std::shared_ptr< GenericParticleFamily> > particleFamilies );

	/*!
	 * \fn virtual ~ParticleRepository()
	 * \brief destructor of the class ParticleRepository
	 */
	virtual ~ParticleRepository();

	/*!
	 * \fn depositParticle( Particle part )
	 * \brief Store a new particle in the repository
	 * \param part : the particle to store
	 */
	void depositParticle( Particle part );

	void depositParticle( BoundaryParticle part );

	/*!
	 * \fn getNumberOfParticles()
	 * \brief Get the number of particles stored in the repository
	 * \return the number of particles stored in the repository
	 */
	int getNumberOfParticles();

	/*!
	 * \fn setActive( bool state )
	 * \brief set if the repository must store particles or not
	 * \param state : the state to set (true = store, false = don't store)
	 */
	void setActive( bool state );

	/*!
	 * \fn std::vector< std::shared_ptr< GenericParticleFamily > > getParticleFamilies()
	 * \brief return the particle families stored in this repository
	 * \return list of particle families stored in this repository
	 */
	std::vector< std::shared_ptr< GenericParticleFamily > > getParticleFamilies();

	/*!
	* return the number of particles stored for each class
	*/
	virtual std::vector< int > getStoredParticlesCounts();

};

} // namespace piaf

#endif /* PARTICLEREPOSITORY_H_ */
