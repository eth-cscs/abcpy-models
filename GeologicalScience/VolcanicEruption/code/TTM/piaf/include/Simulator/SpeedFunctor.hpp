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


#ifndef SPEEDFUNCTOR_H_
#define SPEEDFUNCTOR_H_

#include <memory>
#include "SimulatorTypes.hpp"
#include "GenericParticleFamily.hpp"

#include "boost/shared_ptr.hpp"

//using namespace piaf;

namespace piaf{

/*!
 * \class SpeedFunctor
 * \brief Abstract class used to define advecting field in the domain
 */
class SpeedFunctor {
public:
	/*!
	 * \fn SpeedFunctor()
	 * \brief Constructor of SpeedFunctor
	 */
	SpeedFunctor( std::string name );

	/*!
	 * \fn virtual ~SpeedFunctor()
	 * \brief Destructor of SpeedFunctor (pure virtual method)
	 */
	virtual ~SpeedFunctor();

	/*!
	 * \fn virtual Double3 compute( Double3 pos, double t, double dt, std::shared_ptr< GenericParticleFamily > family )
	 * \brief compute the speed for a given site, time and particle family (pure virtual method), dt is here only for diffusion... try to find a better way to do
	 */
	virtual Double3 operator()( Double3 pos, double t, double dt , std::shared_ptr< GenericParticleFamily > family ) = 0;

	/*!
	 * \fn void setSpeedType( ESpeedRepresentation sr )
	 * \brief set the type of the speed (EULERIAN_STATIC, EULERIAN_DYNAMIC or LAGRANGIAN_DYNAMIC)
	 * \todo : type contraint ? validite check par defaut ? si non faire a la main
	 */
	void setSpeedType( ESpeedType sr );

	/*!
	 * \fn ESpeedRepresentation getSpeedType()
	 * \brief get the type of the speed 
	 * \return the speed type 
	 */
	ESpeedType getSpeedType();

	bool isDiffusion();

	void setIsDiffusion( bool id );

	std::string getName();

protected:
	/*!
	 * \brief type of the speed 
	 */
	ESpeedType speedType_;

	bool isDiffusion_;

	std::string name_;
};

} // namespace piaf

#endif /* SPEEDFUNCTOR_H_ */
