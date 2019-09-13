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


#ifndef MPIEXACTSIMULATOR_H_
#define MPIEXACTSIMULATOR_H_

#include "MPIAbstractSimulator.hpp"
#include "MPIGridTerrain.hpp"
#include "MPIParticleRepository.hpp"
#include "MPI3DBlockCyclicDomain.hpp"
#include "MPIParticleTracker.hpp"

#include "../../Tools/FormattedLog/Trace.hpp"

#include <ctime>

namespace piaf{

class MPIExactSimulator : public MPIAbstractSimulator  {
public:
    MPIExactSimulator( AbstractEulerianDomain *d, double initialTime, double dt, GridTerrain *terrain, MPIParticleRepository *repository, bool ddt, MpiManager & mpiManager);
	virtual ~MPIExactSimulator();


	/* inherited from AbstractSimulator */
	void step( );

	/* inherited from MPIAbstractSimulator */
	/*! \brief : bug, don't use */
    //void stepOverlappedCommunication();

private:

	std::vector< Particle > trackBuffer_;

};

} // namespace piaf

#endif /* MPIEXACTSIMULATOR_H_ */
