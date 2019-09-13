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

#ifndef MPIPARTICLETRACKER_H_
#define MPIPARTICLETRACKER_H_

#include "../ParticleTracker.hpp"
#include "MPI3DGridTopology.hpp"
//#include <boost/mpi.hpp>

namespace piaf{

class MPIParticleTracker : public ParticleTracker {
public:

	MPIParticleTracker( std::vector< std::shared_ptr< GenericParticleFamily> > particleFamilies, MPITopology* topology );

	MPIParticleTracker( std::vector< std::shared_ptr< GenericParticleFamily> > particleFamilies, MPITopology* topology, double trackinterval );

  virtual ~MPIParticleTracker();

  void addTime( double time );

  virtual void addParticle( Particle particle, double time );

  virtual std::tuple< std::vector< double >, std::vector< std::vector< Particle > > > getParticles();

  virtual std::tuple< std::vector< double >, std::vector< std::vector< Particle > > > removeParticles();


private:

	MPITopology* topology_;

};

}// namespace piaf

#endif /* MPIPARTICLETRACKER_H_ */
