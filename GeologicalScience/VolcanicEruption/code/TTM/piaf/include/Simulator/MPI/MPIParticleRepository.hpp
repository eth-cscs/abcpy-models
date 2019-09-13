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

#ifndef MPIPARTICLEREPOSITORY_H_
#define MPIPARTICLEREPOSITORY_H_

#include "../ParticleRepository.hpp"
#include "MPI3DGridTopology.hpp"
#include <boost/serialization/vector.hpp>
//#include <boost/mpi.hpp>


namespace piaf{

class MPIParticleRepository : public ParticleRepository {
public:

	MPIParticleRepository( std::vector< std::shared_ptr< GenericParticleFamily> > particleFamilies, MPI3DGridTopology* topology );
  virtual ~MPIParticleRepository();

  std::vector< BoundaryParticle > removeBoundaryParticles();

  /*!
  * return the number of particles stored for each class
  */
  virtual std::vector< int > getStoredParticlesCounts();

private:

//	void gather();


	MPI3DGridTopology* topology_;

};

}// namespace piaf

#endif /* MPIPARTICLEREPOSITORY_H_ */
