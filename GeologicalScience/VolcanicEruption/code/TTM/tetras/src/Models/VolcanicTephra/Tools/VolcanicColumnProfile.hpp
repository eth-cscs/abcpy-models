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


#ifndef VOLCANICCOLUMNPROFILE_H_
#define VOLCANICCOLUMNPROFILE_H_

#include <vector>
#include <math.h>
#include "Atmosphere.hpp"
#include <boost/shared_ptr.hpp>

#include <iostream>

/*!
 * \class VolcanicColumnProfile
 * \brief this class is used to compute the speed profile and L profile of the volcanic column
 * the computation is made inside the constructor (with RK1-euler for now) and profiles are available
 * with methods getUCentralProfile() and getLProfile()
 * the profile is computed with dx = 1.0, so there is one value per meter from the altitude "ventHeight"
 */

/*!
 * \todo allow the usage of custom dx
 */

class VolcanicColumnProfile {
public:

	/*!
	 * \fn VolcanicColumnProfile( double U0, double L0, double n0, double theta0, double ventHeight, Atmosphere* atmosphere )
	 * \brief constructor of the classe VolcanicColumnProfile
	 * \param U0 : ejection velocity
	 * \param L0 : plume radius at vent
	 * \param n0 : gas mass fraction at vent
	 * \param theta0 : plume temperature at vent
	 * \param ventHeight : altitude of vent
	 * \param atmosphere : the atmosphere
	 */
	VolcanicColumnProfile( double U0, double L0, double n0, double theta0, double ventHeight, std::shared_ptr< Atmosphere > atmosphere );

	virtual ~VolcanicColumnProfile();

	/*!
	 * \fn std::vector<double> getUCentralProfile()
	 * \brief return the central velocity profile of the volcanic plume, there is one value per dx along the Z axis (default dx = 1), first value is at vent height
	 * \return central velocity profile of the plume along Z axis
	 */
	std::vector<double> getUCentralProfile();

	/*!
	 * \fn std::vector<double> getLProfile()
	 * \brief return the plume radius, there is one value per dx along the Z axis (default dx = 1), first value is at vent height
	 * \return plume raidus along Z axis
	 */
	std::vector<double> getLProfile();

	/*!
	 * \fn double getL( double z )
	 * \brief return the plume radius at the given height ( 0 == vent height )
	 * \param z : height
	 * \return plume radius at given height
	 */
	double getL( double z );

	/*!
	 * \fn double getU( double z )
	 * \brief return the central velocity of the plume at the given height ( 0 == vent height )
	 * \param z : height
	 * \return plume central velocity at given height
	 */
	double getU( double z );

	/*!
	 * \fn double getHt()
	 * \brief return the maximum height of the plume ( Ht )
	 * \return Ht
	 */
	double getHt();

	/*!
	 * \fn double getN0()
	 * \brief return the gas mass fraction at vent (n0)
	 * \return n0
	 */
	double getN0();

	/*!
	 * \fn double getTheta0()
	 * \brief return the plume temperature at vent (theta0)
	 * \return theta0
	 */
	double getTheta0();

private:
	/*! \brief the resulting velocity profile */
	std::vector<double> uCentralProfile_;
	/*! \brief the resulting radius profile */
	std::vector<double> lProfile_;

	/*! \brief parameters to be set by the user (via constructor) */
	/*! \brief ejection velocity at vent */
	const double U0_;
	/*! \brief plume radius at vent */
	const double L0_;
	/*! \brief pas mass fraction at vent */
	const double n0_;
	/*! \brief plume temperature at vent */
	const double theta0_;
	/*! \brief height of the vent */
	const double ventHeight_;

	/*! \brief predefined parameters (see the constructor implementation to see the values) */
	/*! \brief delta x (correspond to delta z)*/
	const double dx_;
	/*! \brief specific heat at constant pressure of the air */
	const double Ca_;
	/*! \brief specific heat at constant pressure of the volcanic gas emited at the vent (= Cm) */
	const double Cp0_;
	/*! \brief the entertainment constant */
	const double k_;
	/*! \brief the gas constant for the volcanic gas (= Rm) */
	const double Rg0_;
	/*! \brief the density of the solid pyroclasts */
	const double sigma_;
	/*! \brief the bulk density in the column */
	const double beta0_;
	/*! \brief specific gas constant of the air */
	const double Ra_;

	/*!
	 * \fn static double dUModelA( double U, double L, double alpha, double beta )
	 * \brief compute delta (U) / delta z (delta z = 1) for model A
	 * \param U : speed
	 * \param L : plume radius
	 * \param alpha : atmospheric density (air density)
	 * \param beta : density of plume
	 * \return deltaU / deltaZ (according to model A)
	 */
	static double dUModelA( double U, double L, double alpha, double beta );

	/*!
	 * \fn static double dbetaUL2ModelA( double U, double L, double alpha, double beta, double dU )
	 * \brief compute delta (beta*U*L*L) / delta z (delta z = 1) for model A
	 * \param U : speed
	 * \param L : plume radius
	 * \param alpha : atmospheric density (air density)
	 * \param beta : density of plume
	 * \param dU : deltaU
	 * \return delta (beta*U*L*L) / delta z (according to model A)
	 */
	static double dbetaUL2ModelA( double U, double L, double alpha, double beta, double dU );

	/*!
	 * \fn static double dCpThetaModelAB( double U, double L, double alpha, double beta, double dbetaUL2, double CpTheta, double Ca, double T )
	 * \brief compute delta (Cp*theta) /delta z (delta z = 1) for model A
	 * \param U : speed
	 * \param L : plume radius
	 * \param alpha : atmospheric density (air density)
	 * \param beta : density of plume
	 * \param dbetaUL2 : delta (beta*U*L*L)
	 * \param CpTheta : Cp * theta (Cp = specific heat of the plume at vent)
	 * \param Ca : Ca (specific heat of the aire)
	 * \param T : temperature
	 * \return delta (Cp*theta) /delta z (delta z = 1) (according to model A)
	 */
	static double dCpThetaModelAB( double U, double L, double alpha, double beta, double dbetaUL2, double CpTheta, double Ca, double T );


	/*!
	 * \fn static double dUModelB( double U, double L, double alpha, double beta, double dbetaUL2 )
	 * \brief compute delta (U) / delta z (delta z = 1) for model B
	 * \param U : speed
	 * \param L : plume radius
	 * \param alpha : atmospheric density (air density)
	 * \param beta : density of plume
	 * \param dbetaUL2 : delta (beta * U * L * L)
	 * \return deltaU / deltaZ (according to model B)
	 */
	static double dUModelB( double U, double L, double alpha, double beta, double dbetaUL2 );


	/*!
	 * \fn static double dbetaUL2ModelB( double U, double L, double alpha, double k )
	 * \brief compute delta (beta*U*L*L) / delta z (delta z = 1) for model B
	 * \param U : speed
	 * \param L : plume radius
	 * \param alpha : density of plume
	 * \param k : entrainment constant
	 */
	static double dbetaUL2ModelB( double U, double L, double alpha, double k );

};

#endif /* VOLCANICCOLUMNPROFILE_H_ */
