/*
TEphra TRAsport Simulator (TETRAS)
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


#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include <assert.h>
#include <iostream>

/*!
 * \class Atmosphere
 * \brief this class is used to compute characteristics of the atmosphere such as temperature,
 * pressure and air density
 */

/*!
 * \todo add wind models to the atmosphere class
 */

class Atmosphere {
public:

	/*!
	 * \fn Atmosphere( double bottomTropopauseHeight, double topTropopauseHeight, double T0, double P0 )
	 * \brief constructor of the class Atmosphere
	 * \param tropopauseHeight : height of the tropopause
	 * \param topSteadyTemperatureHeight : height of the top of the steady temperature zone (above the tropopause)
	 * \param T0 : temperature at height 0
	 * \param P0 : pression at height 0
	 */
	Atmosphere( double tropopauseHeight, double topSteadyTemperatureHeight, double T0, double P0 );

	virtual ~Atmosphere();

	/*!
	 * \fn double getTropopauseHeight()
	 * \brief return the height of the tropopause
	 * \return tropopause height
	 */
	double getTropopauseHeight();

	/*!
	 * \fn double getTopSteadyHeight()
	 * \brief return the height of the top of the steady temperature zone
	 * \return top steady temperature height
	 */
	double getTopSteadyHeight();

	/*!
	 * \fn double getT0()
	 * \brief return the temperature at altitude 0
	 * \return T0
	 */
	double getT0();

	/*!
	 * \fn double getP0()
	 * \brief return the atmospheric pressure at altitude 0
	 * \return P0
	 */
	double getP0();

	/*!
	 * \fn double getRa()
	 * \brief return the gas constant of the air
	 * \return Ra
	 */
	double getRa();

	/*!
	 * \fn double getT( double z )
	 * \brief return the temperature at given height
	 * \param z : height
	 * \return temperature
	 */
	double getT( double z );

	/*!
	 * \fn double getP( double z )
	 * \brief return the pressure at given height
	 * \param z : height
	 * \return pressure
	 */
	double getP( double z );

	/*!
	 * \fn double getAlpha( double T, double P )
	 * \brief return the ambiant atmospheric density alpha
	 * \param T : temperature
	 * \param P : pressure
	 * \return alpha
	 */
	double getAlpha( double T, double P );

private:

	/*! \brief the height of the tropopause */
	double tropopauseHeight_;

	/*! \brief the height of the end of the zone of steady temperature () */
	double topSteadyTemperatureHeight_;

	/*! \brief temperature at height 0 */
	double T0_;
	/*! \brief pressure at height 0 */
	double P0_;

	/*! \brief the gas constant for the air */
	const double Ra_;
	/*! \brief temperature gradient in the troposphere */
	const double omegaT_;
	/*! \brief temperature gradient in the stratosphere */
	const double omegaS_;
};

#endif /* ATMOSPHERE_H_ */
