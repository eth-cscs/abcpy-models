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


#ifndef MPIABSTRACTSIMULATOR_HPP_
#define MPIABSTRACTSIMULATOR_HPP_

#include "../AbstractSimulator.hpp"
#include "MPIParticleRepository.hpp"

#include "mpi.h"

namespace piaf {

class MPIAbstractSimulator : public AbstractSimulator {
public:
    MPIAbstractSimulator( AbstractEulerianDomain *d, double initialTime, double dt, GridTerrain *terrain, MPIParticleRepository *repository, bool ddt, MpiManager & mpiManager );
	virtual ~MPIAbstractSimulator();

    virtual void initRand( int seed );

    //virtual void stepOverlappedCommunication() = 0;
protected:
    MpiManager & mpiManager;
};

} /* namespace piaf */
#endif /* MPIABSTRACTSIMULATOR_HPP_ */
