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


#ifndef CASIMULATOR_H_
#define CASIMULATOR_H_

#include "Domain.hpp"
#include "AbstractSimulator.hpp"
#include "GridTerrain.hpp"
#include "SimulatorTypes.hpp"
#include "boost/foreach.hpp"
#include <cstdlib>
#include <assert.h>
#include <ctime>

namespace piaf{

/*!
 * \class CASimulator
 * \brief This class define a simulator of particles in a domain under an advecting field
 * this simulator use a cellular automata representation of space, each site has a certain amount of particles
 * and each particle has a probability to jump to a neighbor site depending of the speed of the particle
 */
class CASimulator : public AbstractSimulator {
public:

	/*!
	 * \fn CASimulator( Domain *d, double dt, GridTerrain *terrain, ParticleRepository *repository )
	 * \brief Constructor of the class CASimulator
	 * \param d : The domain
	 * \param dt : The time integration step
	 * \param terrain : The terrrain
	 * \param repository : a particle repository (store particles which are outside the domain but not deposited on the ground )
	 */
	CASimulator( AbstractEulerianDomain *d, double initialTime, double dt, GridTerrain *terrain, ParticleRepository *repository, bool ddt );


	/*!
	 * \fn ~CASimulator()
	 * \brief The destructor of the class CASimulator
	 */
	virtual ~CASimulator();

	/*!
	 * \fn step()
	 * \brief Perform a step of simulation, ie move the particles according to advecting field, dx and dt
	 */
	void step( );

private :

	/*!
	 * \fn double computeProba( double speed, double dx )
	 * \brief compute the probability to jump from a site to another according to the speed ans dx (compute speed * dt / dx)
	 * and ensure that probability <= 1.0
	 * \return the probability to jump from a site to another
	 */
	double computeProba( double speed, double dx );

};

} //namespace piaf

#endif /* CASIMULATOR_H_ */
