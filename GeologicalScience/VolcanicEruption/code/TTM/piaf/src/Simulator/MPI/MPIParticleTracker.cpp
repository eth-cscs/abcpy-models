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

#include "../../include/Simulator/MPI/MPIParticleTracker.hpp"

namespace piaf{

  MPIParticleTracker::MPIParticleTracker( std::vector< std::shared_ptr< GenericParticleFamily> > particleFamilies, MPITopology* topology ):
  ParticleTracker( particleFamilies ),
  topology_ (topology){}

  MPIParticleTracker::MPIParticleTracker( std::vector< std::shared_ptr< GenericParticleFamily> > particleFamilies, MPITopology* topology, double trackinterval ):
  ParticleTracker( particleFamilies, trackinterval ),
  topology_ (topology){}

  MPIParticleTracker::~MPIParticleTracker(){}

  void MPIParticleTracker::addTime( double time ){
    if( times_.size()==0 || time > times_.back() ){
      times_.push_back( time );
      lastTrackTime = time;
      particles_.push_back( std::vector< Particle >() );
    }
    else if( time < times_.back() ){
      throw std::runtime_error("Error in MPIParticleTracker, particles must be inserted in increasing time order");
    }
  }

  void MPIParticleTracker::addParticle( Particle particle, double time ){
    if( times_.size() > 0 && times_.back() == time ){
      particles_.back().push_back( particle );
    }
    else{
      throw std::runtime_error("Error in MPIParticleTracker, time must be inserted with addTime()");
    }
  }

  std::tuple< std::vector< double >, std::vector< std::vector< Particle > > > MPIParticleTracker::getParticles(){
    std::vector< std::vector< Particle > > particlesRes;

    MpiManager & mpiManager = this->topology_->getMpiManager();
    for( std::vector< Particle >& particles : particles_ ){
      int localSize = particles.size()*sizeof(Particle);
      std::vector< int > localSizes( topology_->getSize() );
      mpiManager.gather<int>(&localSize, 1, &localSizes[0], 1, mpiManager.bossId());
      std::vector< int > displacements;
      int last = 0;
      for( int e : localSizes ) {
        displacements.push_back( last );
        last += e;
      }
      std::vector< Particle > incomingBuffer;
      if( topology_->getRank() == 0 ){
        int globalSize = 0;
        for( int e : localSizes ) globalSize += e;
        incomingBuffer.resize( globalSize/sizeof(Particle) );
      }
       mpiManager.gatherV((char*)particles.data() , localSize, (char*)incomingBuffer.data(), localSizes.data(), displacements.data(), mpiManager.bossId());
      particlesRes.push_back(incomingBuffer);
    }

    return make_tuple( times_, particlesRes );

  }

  std::tuple< std::vector< double >, std::vector< std::vector< Particle > > > MPIParticleTracker::removeParticles(){
    auto res = getParticles();
    clearParticles();
    return res;
  }

} // namespace piaf
