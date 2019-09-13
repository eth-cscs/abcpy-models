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


#ifndef MPIAbstractEulerianDomain_HPP_
#define MPIAbstractEulerianDomain_HPP_

#include "../Domain.hpp"
#include "../AbstractEulerianDomain.hpp"
#include "MPITopology.hpp"

namespace piaf{

	class MPIAbstractEulerianDomain : public AbstractEulerianDomain {
	public:
		MPIAbstractEulerianDomain(  Double3 pos, Double3 size, double dx , MPITopology* topology );
		MPIAbstractEulerianDomain(  );
		virtual ~MPIAbstractEulerianDomain();

        void gather();


        Domain getGatheredDomain();

		Int3 getGlobalSize();
		Int3 getBlockSize();



		MPITopology* getTopology();
		/* \brief in the MPI version, this function returns the same value for each processor, given the real state of the global domain
		* all MPI process must call it at the same time
		*/
		bool containsParticles();

		virtual Domain* getDomainBlockAtLocalIndex3d( Int3 addr ) = 0;
		virtual void setTopology( MPITopology* topology ) = 0;


		// inherited from AbstractEulerianDomain
		bool isInBounds( Int3 coord );
		void addSpeed( std::shared_ptr< SpeedFunctor > newSpeed );
        void removeSpeed( std::string name );
		void addParticleFamily( std::shared_ptr< GenericParticleFamily > newFamily );
		int getNextFamilyId();
		std::vector< std::shared_ptr< GenericParticleFamily> > getParticleFamilies();
		void setTerrain( GridTerrain *terrain );
		Int3 getNumBlockLocal();
        std::vector< double > getNumberOfparticles();
        std::shared_ptr<std::vector< Particle >> getParticles();
        std::vector< std::vector<double> > getDomainParticlesNumber();
        void insertParticlesRandomPosition ( int incr, int familyId, double scaling );
        void insertParticlesRandomPositionSphere ( int incr, int familyId, double scaling, Double3 sphereCenter, double sphereRadius );
        void insertParticlesRandomPositionSphereGaussian ( int incr, int familyId, double scaling, Double3 sphereCenter, double sphereRadius, double stdev );

        int countDomainParticles();

		std::vector< piaf::Particle > getParticlesInBox( piaf::Double3 p1, piaf::Double3 p2 );

        std::vector< std::string > getSpeedsName();

	protected:

		MPITopology* topology_;

		Int3 blockSize_;

		/*! \brief number of local blocks in each dimension */
		Int3 numBlockLocal_;



		//int *containsParticlesArray_;

		//bool *containsParticlesArray_;
		//std::vector< bool > containsParticlesArray_;

		std::vector< Domain > domainBlocks_;

		std::vector< std::vector < Particle > > particlePackets_;

        Domain gatheredDomain_;





		/*!
		* \brief send packets of particles to corresponding neighbors and receive packets of particles
		* then put newly received particles in local domain blocks ( not buffer ! so be careful when combining calls of this function and swapBuffer ), all MPI process must call it at the same time, pending packets are cleaned during the process
		*/
		virtual void exchangeParticles() = 0;

		virtual void requestExchangeParticles() = 0;

		virtual void completeExchangeParticles() = 0;

		/*!
		* \brief put a particle in a packet, packets of particles are then sent to other processors
		* the address of the receiver is computed from the particle position, so the particle position must be a complete absolute position
		* \param part : the particule to package
		* \param lastDomainBlock : the domain block from which the particle come, allow to compute the neighbor process
		*/
		virtual void putParticleInPacket( Particle part, Domain* lastDomainBlock ) = 0;

		//virtual Domain* getDomainBlockAtLocalIndex3d( Int3 addr ) = 0;

		friend class MPIExactSimulator;

		void swapBuffer();
		void computeSpeeds( double t, double dt );
		void computeSpeedsOnBorder( double t, double dt );
		void computeSpeedsExceptBorder( double t, double dt );
		void computeStaticSpeeds();

	};

} //namespace piaf

#endif /* MPIAbstractEulerianDomain_HPP_ */
