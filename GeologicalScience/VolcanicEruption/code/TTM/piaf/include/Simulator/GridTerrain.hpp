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


#ifndef GRIDTERRAIN_H_
#define GRIDTERRAIN_H_

//#include <boost/shared_array.hpp>
//#include <boost/shared_ptr.hpp>
//#include <boost/foreach.hpp>
#include "ParticleRepository.hpp"
#include "SimulatorTypes.hpp"

//using namespace piaf;

namespace piaf{

//#include <assert.h>

/*!
 * \class GridTerrain
 * \brief Class representing the terrain with a 2D grid of heights, inherits from ParticleRepository
 */
class GridTerrain : public ParticleRepository {
protected :

	/*! \brief Spatial integration step */
	double dx_;
	/*! \brief Discrete size of the terrain (size of the array) */
	Int2 discreteSize_;
	/*! \brief UTM position */
	Double2 position_;
	/*! \brief Terrain itself, a height for each point */
	double* terrain_;
	/*! \brief The list of particles which fall outside the terrain */
	std::vector< Particle > garbage_;

public:

	/*!
	 * \fn GridTerrain( Double2 size, double dx, std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies )
	 * \brief Constructor of the classe GridTerrain
	 * \param pos : the position of the terrain ( (x, y) UTM)
	 * \param size : the size of the domain (x,y)
	 * \param dx : the spatial integration step, must be > 0 (not checked)
	 * \param particleFamilies : the list of particle families that can be deposited
	 */
	GridTerrain( Double2 pos, Double2 size, double dx, std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies );

	/*!
	 * \fn ~Domain()
	 * \brief Destructor of the class Domain
	 */
	virtual ~GridTerrain();

	/*!
	 * \fn depositParticle( Particle part )
	 * \brief Take the particle as input and deposit on the ground, the particle must hold his speed and position !
	 * \param part : the particle to deposit
	 */
	void depositParticle( Particle part );

	/*
	 * \fn isInGround( Double3 coord )
	 * \brief Test if a given coordinate is in the ground or not
	 * \param coord : the coordinate to test
	 * \return : true if the coordinate lies in the ground, false otherwise
	 */
	bool isInGround( Double3 coord );

	/*!
	 * \fn bool isInBounds( Int2 coord )
	 * \brief Check if the given coordinate (x,y) lie into the terrain
	 * \param coord : coordinate (x,y)
	 * \return True if the coordinate lies into the terrain, false otherwise
	 */
	bool isInBounds( Int2 coord );

	/*!
	 * \fn bool isInBounds( Double2 coord )
	 * \brief Check if the given coordinate (x,y) lie into the terrain
	 * \param coord : coordinate (x,y)
	 * \return True if the coordinate lies into the terrain, false otherwise
	 */
	bool isInBounds( Double2 coord );

	/*!
	 * \fn int index2d( Int2 coord )
	 * \brief Return the index of the terrain corresponding to this coord (x,y)
	 * \param coord : coordinate (x,y)
	 * \return The index of the domain corresponding to this coord (x,y)
	 */
	int index2d( Int2 coord );

	/*!
	 * \fn std::vector< std::vector< int > > getParticleDeposition( Int2 sampleSize )
	 * \brief get a vector storing the number of particle deposited in each area of the terrain for each particle family
	 * \brief position in the int vector == family id
	 * \brief the array is a contiguous (x,y) array with row major ordering
	 * \param sampleSize : the discrete size (x,y) of the sample to return
	 * \return the vector containing particle deposition
	 */
	virtual std::vector< std::vector< int > > getParticleDeposition( Int2 sampleSize );

	/*!
	 * \fn std::vector< std::vector< int > > getParticleDeposition()
	 * \brief get a vector storing the number of particle deposited in each area of the terrain for each particle family
	 * \brief position in the int vector == family id
	 * \brief the array is a contiguous (x,y) array with row major ordering corresponding to the discrete size of the terrain
	 * \return the vector containing particle depositions
	 */
	virtual std::vector< std::vector< int > > getParticleDeposition();

	/*!
	 * \fn Int2 getDiscreteSize()
	 * \brief get the size of the terrain (x,y)
	 * \return the discrete size of the terrain
	 */
	Int2 getDiscreteSize();

	int getXSize();

	int getYSize();

	/*!
	 * \fn double getDx()
	 * \brief get the spatial integration step of the terrain
	 * \return dx
	 */
	double getDx();

	Double2 getPosition();

	/*!
	 * \fn int getNumberOfParticles()
	 * \brief get the total number of deposited particles
	 * \return the total number of deposited particle
	 */
	virtual int getNumberOfParticles();

	virtual int countTerrainParticles();

	/*!
	* return the number of particles stored for each class
	*/
	virtual std::vector< int > getStoredParticlesCounts();

};

}//namespace piaf

#endif /* GRIDTERRAIN_H_ */
