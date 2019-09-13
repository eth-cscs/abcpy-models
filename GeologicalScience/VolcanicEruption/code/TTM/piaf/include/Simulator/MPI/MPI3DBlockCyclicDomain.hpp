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


//#include "../Domain.hpp"
#include "MPIAbstractEulerianDomain.hpp"
#include "MPI3DGridTopology.hpp"
#include <boost/shared_array.hpp>
#include <boost/foreach.hpp>
#include <cmath>



#ifndef MPI3DBLOCKCYCLICDOMAIN_H_
#define MPI3DBLOCKCYCLICDOMAIN_H_

//using namespace piaf;

namespace piaf{

/*! \class MPI3DBlockCyclicDomain */

class MPI3DBlockCyclicDomain : public MPIAbstractEulerianDomain{
public:
    MPI3DBlockCyclicDomain( Double3 pos, Double3 size, double dx , MPI3DGridTopology* topology );
    MPI3DBlockCyclicDomain( Double3 pos, Double3 size, double dx , MPI3DGridTopology* topology, int nBlock );
    MPI3DBlockCyclicDomain();

    virtual ~MPI3DBlockCyclicDomain();

    // specific to MPI3DBlockCyclicDomain

    void test();
    int blockLocalIndex3d( Int3 coord );
    int blockGlobalIndex3d( Int3 coord );
    Double3 blockLocalAddr2AbsolutePosition( Int3 addr );
    Int3 absolutePosition2BlockGlobalAddr( Double3 coord );
    Int3 blockGlobalAddr2BlockLocalAddr( Int3 addr );
    Int3 domainGlobalAddr2BlockGlobalAddr( Int3 addr );
    Int3 domainGlobalAddr2DomainLocalAddr( Int3 addr );
    Int3 absolutePosition2BlockLocalAddr( Double3 coord );
    //Domain* getDomainBlockAtLocalIndex3d( Int3 addr );

    // inherited from MPIAbstractEulerianDomain
    //Int3 getGlobalSize();
    //Int3 getBlockSize();
    Domain* getDomainBlockAtLocalIndex3d( Int3 addr );
    void setTopology( MPITopology* topology );
    //MPITopology* getTopology();
    /* \brief in the MPI version, this function returns the same value for each processor, given the real state of the global domain
     * all MPI process must call it at the same time
     */
    //bool containsParticles( );

    // inherited from AbstractDomain
    //void insertParticles( Double3 coord, int incr, int familyId );
    void insertParticles( Double3 coord, int incr, int familyId, double scaling );
    //void addSpeed( std::shared_ptr< SpeedFunctor > newSpeed );
    //void addParticleFamily( std::shared_ptr< GenericParticleFamily > newFamily );
    //int getNextFamilyId();
    Double3 getSpeed( Int3 coord, int familyId );
    bool isInGroundContactAt( Int3 coord );
    void putParticle( Particle part, Int3 coord );
    void putParticle( Particle part );
    //std::vector< std::shared_ptr< GenericParticleFamily> > getParticleFamilies();
    //void setTerrain( GridTerrain *terrain );
    void init( Double3 pos, Double3 size, double dx );


    // std::vector< piaf::Particle > getBox( Double3 p1, Double3 p2 );

protected:

    // used for asynchronous exchange
    std::vector< std::vector<char> > incomingBuffers_;
    std::vector< std::string > sendBuffers_;
    std::vector< int > recvSizes_, sendSizes_;
    std::vector < MPI_Request > requestsSend_;
    std::vector < MPI_Status > statusesSend_;
    std::vector < MPI_Request > requestsRecv_;
    std::vector < MPI_Status > statusesRecv_;

    // specific to MPI3DBlockCyclicDomain
    /*! \brief tells how many times the "topology of processors" must be copied in each dimension ( used to compute block size ), ie the number of block in each dimension per process */
    const int adjustmentFactor_;
    MPI3DGridTopology* gridTopology_;

    // inherited from MPIAbstractEulerianDomain
    /*!
     * \brief put a particle in a packet, packets of particles are then sent to other processors
     * the address of the receiver is computed from the particle position, so the particle position must be a complete absolute position
     * \param part : the particule to package
     * \param lastDomainBlock : the domain block from which the particle come, allow to compute the neighbor process
     */
    void putParticleInPacket( Particle part, Domain* lastDomainBlock );
    /*!
     * \brief send packets of particles to corresponding neighbors and receive packets of particles
     * then put newly received particles in local domain blocks ( not buffer ! so be careful when combining calls of this function and swapBuffer ), all MPI process must call it at the same time, pending packets are cleaned during the process
     */
    // inherited
    void exchangeParticles();

    void requestExchangeParticles();
    void completeExchangeParticles();

    // inherited from AbstractDomain
    /* \brief pure virtual member functions from AbstractDomain */
    void insertParticlesInBuffer( Int3 coord, int incr, int familyId );
    void eraseAt( Int3 coord );
    void eraseBufferAt( Int3 coord );
    //void swapBuffer();
    //void computeSpeeds( double t, double dt );
    //void computeStaticSpeeds();
    void putParticleInBuffer( Particle part, Int3 coord );
    void putParticleInBuffer( Particle part );

    //friend class MPIExactSimulator;

    void commonConstruct( Double3 pos, Double3 size, double dx , MPI3DGridTopology* topology );

};

} // namespace piaf

#endif /* MPI3DBLOCKCYCLICDOMAIN_H_ */
