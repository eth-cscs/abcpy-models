#ifndef ABSTRACTDOMAIN_H
#define ABSTRACTDOMAIN_H

#include <vector>
#include <memory>

#include "SpeedFunctor.hpp"
#include "SimulatorTypes.hpp"
#include "GridTerrain.hpp"

namespace piaf{
    
    class AbstractDomain{

        public:

            AbstractDomain();

            virtual ~AbstractDomain();

            virtual std::vector< double > getNumberOfparticles() = 0;

            virtual void insertParticles( Double3 coord, int incr, int familyId, double scaling ) = 0;

            virtual void addSpeed( std::shared_ptr< SpeedFunctor > newSpeed );

            virtual void removeSpeed( std::string name );

            virtual void addParticleFamily( std::shared_ptr< GenericParticleFamily > newFamily );

            virtual int getNextFamilyId() = 0;

            virtual void putParticle( Particle part ) = 0;

            virtual std::vector< std::shared_ptr< GenericParticleFamily> > getParticleFamilies() = 0;

            virtual void setTerrain( GridTerrain *terrain ) = 0;

            virtual bool containsParticles( ) = 0;

		    virtual int countDomainParticles() = 0;

            virtual std::vector< piaf::Particle > getParticlesInBox( Double3 p1, Double3 p2 ) = 0;

        protected:

            virtual void computeSpeeds( double t, double dt );

            virtual std::shared_ptr<std::vector< Particle >> getParticles() = 0;

    };

}

#endif