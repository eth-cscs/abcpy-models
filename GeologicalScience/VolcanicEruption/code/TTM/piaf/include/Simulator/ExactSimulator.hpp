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


#ifndef EXACTSIMULATOR_H_
#define EXACTSIMULATOR_H_

#include "Domain.hpp"
#include "AbstractSimulator.hpp"
#include "GridTerrain.hpp"
#include "../Tools/FormattedLog/Trace.hpp"
#include <cstdlib>
#include <ctime>

namespace piaf{

/*!
 * \class ExactSimulator
 * \brief This class define a simulator of particles in a domain under an advecting field and keep exact position of particles
 * each particle remain linked to its domain site, but each particle keep its exact position and speed, so its a lagrangian on grid model
 */
class ExactSimulator : public AbstractSimulator {

public:



	/*!
	 * \fn ExactSimulator( Domain *d, double dt, GridTerrain *terrain, ParticleRepository *repository )
	 * \brief Constructor of the class ExactSimulator
	 * \param d : The domain
	 * \param dt : The time integration step
	 * \param terrain : The terrain
	 * \param repository : a particle repository (store particles which are outside the domain but not deposited on the ground )
	 */
	ExactSimulator( AbstractEulerianDomain *d, double initialTime, double dt, GridTerrain *terrain, ParticleRepository *repository, bool ddt );


	/*!
	 * \fn ~ExactSimulator()
	 * \brief The destructor of the class ExactSimulator
	 */
	virtual ~ExactSimulator();

	/*!
	 * \fn step()
	 * \brief Perform a step of simulation, ie move the particles according to advecting field, dx and dt
	 */
	void step( );

	//private:

	//double diffusionSpeed_;

};

}//namespacce piaf

#endif /* EXACTSIMULATOR_H_ */
