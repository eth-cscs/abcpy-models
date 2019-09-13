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


#ifndef AbstractEulerianDomain_H_
#define AbstractEulerianDomain_H_

#include "SimulatorTypes.hpp"
#include "GenericParticleFamily.hpp"
#include "GridTerrain.hpp"
#include "SpeedFunctor.hpp"
#include "VariousFunctions.hpp"
#include <boost/serialization/access.hpp>
#include <vector>
#include <memory>
#include <random>
#include "../Tools/FormattedLog/Log.hpp"

//using namespace piaf;

namespace piaf{

	class AbstractEulerianDomain {
	public:
		/*!
		* \fn AbstractEulerianDomain()
		* \brief Constructor of the classe AbstractEulerianDomain, the domain has a size in space and a position
		* the position represents an offset regarding the coordinate 0,0,0
		* there are two types of coordinate, absolute which are real coordinate in space
		* and relative which represent a coordinate into the domain
		* in general, discrete (integer) coordinates are relative and real (float) coordinates are absolutes
		* the type of each coordinate is specified for each method
		* \param pos : the position of the point 0, 0, 0 of the domain
		* \param size : the size of the domain (x,y,z)
		* \param dx : the spatial integration step, must be > 0 (not checked)
		*/
		AbstractEulerianDomain( Double3 pos, Double3 size, double dx );

		AbstractEulerianDomain();

		/*!
		* \fn ~Domain()
		* \brief Destructor of the class Domain
		*/
		virtual ~AbstractEulerianDomain();
		
		/*\fn test()
		* \brief debug
		*/
		//virtual void test() = 0;

		/*!
		* \fn double getDx()
		* \brief Get the spatial integration step
		* \return The spatial integration step
		*/
		double getDx();

		/*!
		* \fn int getXSize()
		* \brief Get the discrete horizontal size of the domain
		* \return The discrete horizontal size of the domain
		*/
		int getXSize();

		/*!
		* \fn int getYSize()
		* \brief Get the discrete vertical size of the domain
		* \return The discrete vertical size of the domain
		*/
		int getYSize();

		/*!
		* \fn int getZSize()
		* \brief Get the discrete depth size of the domain
		* \return The discrete depth size of the domain
		*/
		int getZSize();

		/*!
		* \fn Int3 getSize()
		* \brief Get the discrete size of the domain
		* \return The discrete size of the domain
		*/
		Int3 getSize();

		/*!
		* \fn Double3 getPosition()
		* \brief get the position of the domain in space ( coordinate of the point (0,0,0) of the domain)
		* \return the position of the domain
		*/
		Double3 getPosition();

		/*!
		* \fn bool equals( AbstractEulerianDomain* dom )
		* \brief test if the given domain equals the current object ( dx, size and particles ) ! causes complete copy of dom ( and gather for mpi version ), use only for tests !
		* \param dom : the other domain to compare
		* \return true if dom == curent domain, false otherwise
		*/
		//bool equals( AbstractEulerianDomain* dom );

		/*!
		* \fn std::vector< std::vector< Particle > > getDomainParticles()
        * \brief get the complete domain in the form of a vector of vector of particles, return a copy of the full domain and not a pointer to the internal data !
		* the vector must be seen as a 3d array with row major ordering ( index x,y,z located at x + y * sizex + z * sizex * sizey  )
		* \return the full domain with particles
		*/
        //virtual std::vector< std::vector< Particle > > getDomainParticles() = 0;

        /*!
         * \brief getNumbersOfparticles
         * \return the total nomber of particles for each family
         */
        virtual std::vector< double > getNumberOfparticles() = 0;

		virtual void init( Double3 pos, Double3 size, double dx ) = 0;

		/*!
		* \fn bool isInBounds( Int3 coord )
		* \brief Check if the given relative coordinate (x,y,z) lies into the domain
		* \param coord : relative coordinate (x,y,z)
		* \return True if the coordinate lie into the domain, false otherwise
		*/
		virtual bool isInBounds( Int3 coord ) = 0;

		/*!
		* \fn void insertParticles( Double3 coord, int incr, int familyId )
		* \brief Increment the number of particles of the given family at the given coordinate, if incr < 0 particles are deleted
		* \param coord : absolute coordinate (x,y,z)
		* \param incr : The increment
		* \param familyId : The id of the Particle family
		*/
		//virtual void insertParticles( Double3 coord, int incr, int familyId ) = 0;

        virtual void insertParticles( Double3 coord, int incr, int familyId, double scaling ) = 0;

        virtual void insertParticlesRandomPosition( int incr, int familyId, double scaling ) = 0;

        virtual void insertParticlesRandomPositionSphere( int incr, int familyId, double scaling, Double3 sphereCenter, double sphereRadius ) = 0;

        virtual void insertParticlesRandomPositionSphereGaussian( int incr, int familyId, double scaling, Double3 sphereCenter, double sphereRadius, double stdev ) = 0;

		/*!
		* \fn addSpeed( std::shared_ptr< SpeedFunctor > newSpeed )
		* \brief Add a new speed function to the domain
		* \param newSpeed : The new speed to take into account
		*/
		virtual void addSpeed( std::shared_ptr< SpeedFunctor > newSpeed );

        virtual void removeSpeed( std::string name );

		/*!
		* \fn addParticleFamily( std::shared_ptr< GenericParticleFamily > newFamily )
		* \brief Add a new family of particles to the domain
		* \param newFamily : The new family of particles
		*/
		virtual void addParticleFamily( std::shared_ptr< GenericParticleFamily > newFamily );

		/*!
		* \fn int getNextFamilyId()
		* \brief Return the next id to use for a new Particle family
		* \return a new id
		*/
		virtual int getNextFamilyId() = 0;

		/*!
		* \fn Double3 getSpeed( Int3 coord , int familyId )
		* \brief get the speed of a Particle lying on (x,y,z) with given family id
		* \param coord : relative coordinate (x,y,z)
		* \param familyId : id of the Particle family
		* \return The speed of the Particle
		*/
		virtual Double3 getSpeed( Int3 coord, int familyId ) = 0;


		/*!
		* \fn isInGroundContactAt( Int3 coord )
		* \brief Tells if the domain is in contact with the ground a this site
		* \param coord : relative coordinate
		* \return true if the site is in contact with the ground, false otherwise
		*/
		virtual bool isInGroundContactAt( Int3 coord ) = 0;

		/*!
		* \fn void putParticle( Particle part, Int3 coord )
		* \brief place the given particle in the given site of the domain
		* \param part : the particle
		* \param coord : the relative coordinate in the domain
		*/
		virtual void putParticle( Particle part, Int3 coord ) = 0;

		/*!
		* \fn void putParticle( Particle part )
		* \brief place the given particle in the domain ( the particle must hold its absolute position )
		* \param part : the particle
		*/
		virtual void putParticle( Particle part ) = 0;

		/*!
		* \fn std::vector< std::shared_ptr< GenericParticleFamily > > getParticleFamilies()
		* \brief return the current particle families registered in the domain
		* \return particle families
		*/
		virtual std::vector< std::shared_ptr< GenericParticleFamily> > getParticleFamilies() = 0;

		/*!
		* \fn void setTerrain( GridTerrain *terrain )
		* \brief specifies the terrain linked to the domain
		* \param terrain : the terrain used
		*/
		virtual void setTerrain( GridTerrain *terrain ) = 0;

		/*!
		* \fn bool containsParticles( )
		* \brief tells if the domain contain particles or not
		* \return true if the domain contain particles, false otherwise
		*/
		virtual bool containsParticles( ) = 0;

		virtual int countDomainParticles() = 0;

		/*!
		* \fn std::vector< piaf:Particle > getParticlesInBox( Double3 p1, Double3 p2 )
		* \brief get all the particles contained within the box defined by p1, p2
		*/
		virtual std::vector< piaf::Particle > getParticlesInBox( Double3 p1, Double3 p2 ) = 0;

        virtual std::vector< std::string > getSpeedsName() = 0;

        /*!
         * \brief getDomainParticlesNumber
         * \return a vector where each element is for a class of particle and contains a 3d array of number of particles
         */
        virtual std::vector< std::vector<double> > getDomainParticlesNumber() = 0;


	protected:

		std::function<int (int)>  testFunction;

		friend class boost::serialization::access;

		template<class Archive>
		void serialize( Archive & ar, const unsigned int version ){
			ar & dx_;
			ar & size_;
			ar & position_;
			ar & containsParticles_;
		}

		/*! \brief Spatial integration step */
		double dx_;

		/*! \brief Discrete size of the domain */
		Int3 size_;

		/*! \brief position of the point 0, 0, 0 of the domain */
		Double3 position_;

		/*! \brief tells if the domain contains particles !!! true means maybe and false means no !!! */
		int containsParticles_;

		bool staticSpeedsUpToDate_;

		/*!
		* \fn void insertParticlesInBuffer( Int3 coord, int incr, int familyId )
		* \brief Increment the number of particles of the given family in the buffer at the given coordinate, if incr < 0 particles are deleted
		* \param coord : relative coordinate (x,y,z)
		* \param incr : The increment
		* \param familyId : id of the Particle family
		*/
		virtual void insertParticlesInBuffer( Int3 coord, int incr, int familyId ) = 0;

		/* !
		* \fn void eraseAt( Int3 coord )
		* \brief Erase all particles in this site
		* \param coord : relative coordinate(x,y,z)
		*/
		virtual void eraseAt( Int3 coord ) = 0;

		virtual void eraseBufferAt( Int3 coord ) = 0;

		/*!
		* \fn void swapBuffer()
		* \brief Exchange the values of the domain pointer and the buffer pointer
		*/
		virtual void swapBuffer() = 0;

		virtual void computeStaticSpeeds() = 0;

		/*!
		* \fn void computeSpeeds( double t, double dt )
		* \brief compute the speeds at each point of the domain for each Particle family
		* \param t : The time of the simulation
		* \param dt : time integration step
		*/
		virtual void computeSpeeds( double t, double dt );

		virtual void computeSpeedsOnBorder( double t, double dt ) = 0;

		virtual void computeSpeedsExceptBorder( double t, double dt ) = 0;

		/*!
		* \fn void putParticleInBuffer( Particle part, Int3 coord )
		* \brief put a given Particle in the buffer
		* \param part : the Particle
		* \param coord : relative coordinate
		*/
		virtual void putParticleInBuffer( Particle part, Int3 coord ) = 0;

		/*!
		* \fn void putParticleInBuffer( Particle part )
		* \brief place the given particle in the buffer ( the particle must hold its absolute position )
		* \param part : the particle
		*/
		virtual void putParticleInBuffer( Particle part ) = 0;

        /*!
        * \fn std::vector< Particle > getParticles()
        * \brief get raw particles
        */
        virtual std::shared_ptr<std::vector< Particle >> getParticles() = 0;


		/*! \todo : heritage (Simulator) + avoid friend class*/
		friend class CASimulator;
		friend class ExactSimulator;
		friend class AbstractSimulator;
	};

} // namespace piaf

#endif /* AbstractEulerianDomain_H_ */
