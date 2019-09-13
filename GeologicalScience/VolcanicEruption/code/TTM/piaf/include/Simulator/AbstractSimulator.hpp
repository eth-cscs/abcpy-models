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


#ifndef ABSTRACTSIMULATOR_H_
#define ABSTRACTSIMULATOR_H_

#include "GridTerrain.hpp"
#include "ParticleRepository.hpp"
#include "ParticleTracker.hpp"
#include "Domain.hpp"
#include <iostream>

#include <boost/foreach.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

namespace piaf{

class AbstractSimulator {
protected:

	/*! \brief Domain to simulate */
	AbstractEulerianDomain *domain_;
	/*! \brief Current time of the simulation */
	double currentTime_;
	/*! \brief Time integration step */
	double dt_;

	bool ddt_;

	double ddtValue_;
	/*! \brief The terrain to deposit particles */
	GridTerrain *terrain_;
	/* \brief repository for particles, used to verify mass conservation */
	ParticleRepository *repository_;

    ParticleTracker *partTracker_;

	int seed_;


public:

	/*!
	 * \fn Simulator( Domain *d, double dt, GridTerrain *terrain, ParticleRepository *repository, bool ddt )
	 * \brief Constructor of the class Simulator
	 * \param d : The domain
	 * \param dt : The time integration step
	 * \param terrain : The terrain
	 * \param repository : a repository to store particles falling outside the terrain
	 * \param ddt : use or not the dynamic dt functionnality
	 */
	AbstractSimulator( AbstractEulerianDomain *d, double initialTime, double dt, GridTerrain *terrain, ParticleRepository *repository, bool ddt );

	/*!
	 * \fn virtual ~Simulator()
	 * \brief destructor of the class Simulator
	 */
	virtual ~AbstractSimulator();

	void setParticleTracker( ParticleTracker *partTracker );

	double getDt();

	/*!
	 * \fn step()
	 * \brief Perform a step of simulation, ie move the particles according to advecting field, dx and dt
	 */
	virtual void step( ) = 0;

	void enableDdt();

	void disableDdt();

	void setDdtValue( double value );
	double getDdtValue( );

    virtual void initRand( int seed );

	int getSeed();

	double getTime();


};

} // namespace piaf

#endif /* SIMULATOR_H_ */
