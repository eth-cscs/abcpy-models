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


#ifndef DOMAIN_H_
#define DOMAIN_H_

#include <functional>
#include <iostream>
#include <map>
#include "SpeedFunctor.hpp"
#include "GenericParticleFamily.hpp"
#include "SimulatorTypes.hpp"
#include "AbstractEulerianDomain.hpp"
#include <algorithm> 

namespace piaf{

//class Site {
//public:
//    std::vector< Particle > particles_;
//    std::vector<double> rests_;
//};


/*!
 * \class Domain
 * \brief Class representing a domain containing particles
 */
class Domain : public AbstractEulerianDomain {
public:

    //boost::posix_time::time_duration timeDebug_;

    /*!
     * \fn Domain( Double3 pos, Double3 size, double dx )
     * \brief Constructor of the classe Domain, the domain has a size in space and a position
     * the position represents an offset regarding the coordinate 0,0,0
     * there are two types of coordinate, absolute which are real coordinate in space
     * and relative which represent a coordinate into the domain
     * in general, discrete (integer) coordinates are relative and real (float) coordinates are absolutes
     * the type of each coordinate is specified for each method
     * \param pos : position of the point 0, 0, 0 of the domain ( ( x, y ) UTM , z altitude )
     * \param size : the size of the domain (x,y,z)
     * \param dx : the spatial integration step, must be > 0 (not checked)
     */
    Domain( Double3 pos, Double3 size, double dx );

    Domain();

    /*!
     * \fn ~Domain()
     * \brief Destructor of the class Domain
     */
    ~Domain();

    /*!
     * \fn test()
     * \brief for debug
     */
    void test();

    /*!
     * \fn void computeMaxNumPart()
     * \brief compute the highest number of particles in a single site of the domain
     */
    void computeMaxNumPart();

    /*!
     * \fn int getMaxNParticles()
     * \brief Get the highest number of particles in a single site of the domain
     * \return The highest number of particles in a single site of the domain
     */
    int getMaxNParticles();

    /*!
     * \fn int getNParticleAt( Int3 coord )
     * \brief Get the number of particles at the place (x,y,z) in the domain
     * \param coord : coordinate (x,y,z)
     * \return The number of particles at the given place
     */
    int getNParticleAt( Int3 coord );

    /*!
     * \fn int getNParticleAt( Int3 coord, int familyId )
     * \brief Get the number of particles of the given family at the place (x,y) in the domain
     * \param coord : coordinate (x,y,z)
     * \param familyId : The id of the Particle family
     * \return The number of particles at the given place
     */
    int getNParticleAt( Int3 coord, int familyId );

    /*!
     * \fn int index3d( Int3 coord )
     * \brief Return the index of the domain corresponding to this relative coord (x,y,z)
     * \param coord : relative coordinate (x,y,z)
     * \return The index of the domain corresponding to this coord (x,y,z)
     */
    int index3d( Int3 coord );

    /*!
     * \fn getNParticles()
     * \brief Get the total number of particles in the domain
     * \return total number of particles in the domain
     */
    int getNParticles();

    std::vector< std::vector< Particle > > 					getDomainParticles();

    std::vector< std::string >                              getSpeedsName();

    /* \brief pure virtual member functions from AbstractEulerianDomain */
    void 													init				( Double3 pos, Double3 size, double dx );
    bool 													isInBounds			( Int3 coord );
    std::vector< std::shared_ptr< GenericParticleFamily> > 	getParticleFamilies	();
    //void 													insertParticles		( Double3 coord, int incr, int familyId );
    void 													insertParticles		( Double3 coord, int incr, int familyId, double scaling );
    void 													insertParticlesRandomPosition ( int incr, int familyId, double scaling );
    void 													insertParticlesRandomPositionSphere ( int incr, int familyId, double scaling, Double3 sphereCenter, double sphereRadius );
    void 													insertParticlesRandomPositionSphereGaussian ( int incr, int familyId, double scaling, Double3 sphereCenter, double sphereRadius, double stdev );
    void 													addSpeed			( std::shared_ptr< SpeedFunctor > newSpeed );
    void 													removeSpeed			( std::string name );
    void 													addParticleFamily	( std::shared_ptr< GenericParticleFamily > newFamily );
    int 													getNextFamilyId		();
    Double3 												getSpeed			( Int3 coord, int familyId );
    bool 													isInGroundContactAt	( Int3 coord );
    void 													putParticle			( Particle part, Int3 coord );
    void 													putParticle			( Particle part );
    void 													setTerrain			( GridTerrain *terrain );
    bool 													containsParticles	();
    std::vector< double >                                   getNumberOfparticles();
    std::shared_ptr<std::vector< Particle >>                getParticles();
    std::vector< std::vector<double> >                      getDomainParticlesNumber();

    int countDomainParticles();

    std::vector< piaf::Particle > getParticlesInBox( Double3 p1, Double3 p2 );

protected:


    friend class boost::serialization::access;

    template< class Archive >
    void serialize( Archive & ar, const unsigned int version ){

        AbstractEulerianDomain::serialize( ar, version );
        ar & domain_;
        ar & buffer_;
    }

    /*! \brief List of Particle families */
    std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies_;
    /*! \brief highest number of particles in a single site of the domain */
    int maxNParticle_;

    /*! \brief Domain itself, a vector of particles for each site */
    //boost::shared_array< std::vector< Particle > > domain_;
    std::vector< std::vector< Particle > > domain_;
    //    std::vector< Site > domain_;

    //    std::vector< std::vector< double > > rests_;

    /*! \brief Buffer, used for computations */
    //boost::shared_array< std::vector< Particle > > buffer_;
    std::vector< std::vector< Particle > > buffer_;
    /*! \brief Speed vectors of the domain, there is one speed for each family in each site (static) */
    std::vector< std::vector< Double3 > > eulerianStaticSpeeds_;
    /*! \brief Speed vectors of the domain, there is one speed for each family in each site (dynamic) */
    std::vector< std::vector< Double3 > > eulerianDynamicSpeeds_;
    /*! \brief List of speed functions applying to the domain (eulerian dynamics) */
    std::vector< std::shared_ptr< SpeedFunctor > > eulerianDynamicSpeedFunctors_;
    /*! \brief List of speed functions applying to the domain (eulerian statics) */
    std::vector< std::shared_ptr< SpeedFunctor > > eulerianStaticSpeedFunctors_;
    /*! \brief List of speed functions applying to the domain (lagrangian dynamics) */
    std::vector< std::shared_ptr< SpeedFunctor > > lagrangianDynamicSpeedFunctors_;
    /*! \brief Tells if a given site of the domain is in ground contact */
    //boost::shared_array< bool > isInGroundContact_;
    std::vector< bool > isInGroundContact_;



    /*!
     * \fn void setAt( Int3 coord, int val, int familyId )
     * \brief Set the number of particles of this Particle family in the domain at the given coordinate
     * \param coord : coordinate (x,y,z)
     * \param incr : The increment
     * \param familyId : id of the Particle family
     */
    void setAt( Int3 coord, int val, int familyId );

    /*\brief pure virtual member functions from AbstractEulerianDomain */
    void insertParticlesInBuffer( Int3 coord, int incr, int familyId );
    void eraseAt( Int3 coord );
    void eraseBufferAt( Int3 coord );
    void swapBuffer();
    void computeStaticSpeeds();
    void computeSpeeds( double t, double dt );
    void computeSpeedsOnBorder( double t, double dt );
    void computeSpeedsExceptBorder( double t, double dt );
    void putParticleInBuffer( Particle part, Int3 coord );
    void putParticleInBuffer( Particle part );


    /*! \todo : avoid friend class*/
    friend class CASimulator;
    friend class ExactSimulator;
    friend class AbstractSimulator;
    friend class MPIAbstractEulerianDomain;
    friend class MPI3DBlockCyclicDomain;
    friend class MPISimpleDomain;
    friend class AbstractEulerianDomain;
    friend class MPIExactSimulator;
};

} //namespace piaf

#endif /* DOMAIN_H_ */
