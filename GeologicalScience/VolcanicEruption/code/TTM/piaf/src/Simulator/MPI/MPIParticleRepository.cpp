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

#include "../../../include/Simulator/MPI/MPIParticleRepository.hpp"

namespace piaf {

MPIParticleRepository::MPIParticleRepository( std::vector< std::shared_ptr< GenericParticleFamily> > particleFamilies, MPI3DGridTopology* topology ) :
        ParticleRepository    ( particleFamilies ),
        topology_ ( topology ){
}
MPIParticleRepository::~MPIParticleRepository() {
}

/**
 * @Brief removes the local particles from local processes.
 */

std::vector< BoundaryParticle >   MPIParticleRepository::removeBoundaryParticles(){

        //mohamed
        MpiManager &mpiManager=this->topology_->getMpiManager();

        int localSize = storedBoundaryParticles_.size()*sizeof(BoundaryParticle);
        std::vector< int > displacements;
        std::vector< int > localSizes( topology_->getSize() );

        mpiManager.gather<int>(&localSize, 1, &localSizes[0], 1, mpiManager.bossId());
        int last = 0;
        for( int e : localSizes ) {
                displacements.push_back( last );
                last += e;
        }
        std::vector< BoundaryParticle > incomingBuffer;
        if( topology_->getRank() == 0 ) {
                int globalSize = 0;
                for( int e : localSizes ) globalSize += e;
                incomingBuffer.resize( globalSize/sizeof(BoundaryParticle) );
        }

        mpiManager.gatherV<char>( (char*) &storedBoundaryParticles_[0], localSize,  (char*) &incomingBuffer[0], &localSizes[0], &displacements[0], mpiManager.bossId());
        storedBoundaryParticles_.clear();
        return incomingBuffer;
}

std::vector< int > MPIParticleRepository::getStoredParticlesCounts(){
        MpiManager &mpiManager=this->topology_->getMpiManager();

        std::vector< int > resLoc( particleFamilies_.size(), 0 );
        std::vector< int > res( particleFamilies_.size(), 0 );
        for( auto part : storedParticles_ ) resLoc[ part.familyId_ ]++;

        mpiManager.reduceVect<int>( resLoc, res, MPI_SUM );
        return res;
}

}
