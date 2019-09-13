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


#ifndef MPISIMPLEDOMAIN_H_
#define MPISIMPLEDOMAIN_H_

#include "MPIAbstractEulerianDomain.hpp"

namespace piaf{

/*! \class MPISimpleDomain */

class MPISimpleDomain : public MPIAbstractEulerianDomain {
public:
    MPISimpleDomain( Double3 pos, Double3 size, double dx , MPITopology* topology );
    MPISimpleDomain( Double3 pos, Double3 size, double dx , MPITopology* topology, int parallelDim );
    MPISimpleDomain(  );
    virtual ~MPISimpleDomain();

    Int3 global2Local( Int3 coord );

    // inherited from AbstractDomain
    //void insertParticles( Double3 coord, int incr, int familyId );
    void insertParticles( Double3 coord, int incr, int familyId, double scaling );
    Double3 getSpeed( Int3 coord, int familyId );
    bool isInGroundContactAt( Int3 coord );
    void putParticle( Particle part, Int3 coord );
    void putParticle( Particle part );
    void init( Double3 pos, Double3 size, double dx );

    // inherited from MPIAbstractEulerianDomain
    Domain* getDomainBlockAtLocalIndex3d( Int3 addr );
    void setTopology( MPITopology* topology );

    // std::vector< piaf::Particle > getBox( Double3 p1, Double3 p2 );

protected:

    // used for asynchronous exchange
    std::vector< std::vector<char> > incomingBuffers_;
    std::vector< std::string > sendBuffers_;
    std::vector< int > recvSizes_, sendSizes_;
    std::vector < MPI_Request > requestsSend_;
    std::vector < MPI_Request > requestsRecv_;
    std::vector < MPI_Status > statusesSend_;
    std::vector < MPI_Status > statusesRecv_;



    int parallelDim_;

    void exchangeParticles();

    void requestExchangeParticles();
    void completeExchangeParticles();

    void putParticleInPacket( Particle part, Domain* lastDomainBlock );


    // inherited from AbstractDomain
    void insertParticlesInBuffer( Int3 coord, int incr, int familyId );
    void eraseAt( Int3 coord );
    void eraseBufferAt( Int3 coord );
    void putParticleInBuffer( Particle part, Int3 coord );
    void putParticleInBuffer( Particle part );

    void commonConstruct( Double3 pos, Double3 size, double dx , MPITopology* topology, int parallelDim );
};

} // namespace piaf

#endif /* MPISIMPLEDOMAIN_H_ */
