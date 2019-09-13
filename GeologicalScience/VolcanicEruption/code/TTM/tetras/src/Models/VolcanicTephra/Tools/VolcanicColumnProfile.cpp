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


#include "VolcanicColumnProfile.hpp"

VolcanicColumnProfile::VolcanicColumnProfile( double U0, double L0, double n0, double theta0, double ventHeight, std::shared_ptr< Atmosphere > atmosphere ):
U0_			( U0 ),
L0_			( L0 ),
n0_			( n0 ),
theta0_		( theta0 ),
ventHeight_	( ventHeight ),
dx_			( 1.0 ),
Ca_			( 998.0 ),
Cp0_		( 1617.0 ),
k_			( 0.09 ),
Rg0_		( 462.0 ),
sigma_		( 2500.0 ),
beta0_		( 1.0 / ( ( 1.0 - n0 ) * ( 1.0 / sigma_ ) + ( n0 * Rg0_ * theta0 ) / ( atmosphere->getP0() ) ) ),
Ra_			( 285.0 )
{

	double Tn, Pn, alphan, Up1, betaUL2p1, CpThetap1, np1, Rgp1, betap1, thetap1, Lp1, Cpp1, dU, dbetaUL2, dCpTheta;

	double Un 		= U0;
	double betaUL2n = beta0_ * U0 * L0 * L0;
	double CpThetan = Cp0_ * theta0;
	double Ln 		= L0;
	double Zn 		= ventHeight;
	double nn		= n0;
	double Rgn		= Rg0_;
	double betan 	= beta0_;
	double thetan 	= theta0;

	// store the profile
	uCentralProfile_	. push_back( U0 );
	lProfile_ 			. push_back( L0 );
	// we use a tolerance, computation with RK1
	while( Un > 1.0 ){
		Tn 		= atmosphere->getT( Zn );
		Pn 		= atmosphere->getP( Zn );
		alphan 	= atmosphere->getAlpha( Tn, Pn );

		dU 			= 0.0;
		dbetaUL2	= 0.0;
		dCpTheta 	= 0.0;

		// gas thrust region, apply model A
		if( alphan < betan ){

			dU 			= dUModelA( Un, Ln, alphan, betan );
			dbetaUL2 	= dbetaUL2ModelA( Un, Ln, alphan, betan, dU );
			dCpTheta 	= dCpThetaModelAB( Un, Ln, alphan, betan, dbetaUL2, CpThetan, Ca_, Tn );

		}

		// convective region, apply model B
		else{

			dbetaUL2 	= dbetaUL2ModelB( Un , Ln, alphan, k_ );
			dU 			= dUModelB( Un, Ln, alphan, betan, dbetaUL2 );
			dCpTheta 	= dCpThetaModelAB( Un, Ln, alphan, betan, dbetaUL2, CpThetan, Ca_, Tn );

		}

		dU 			*= dx_;
		dbetaUL2 	*= dx_;
		dCpTheta 	*= dx_;

		np1 		= 1.0 + ( n0 - 1.0 ) * ( ( L0 * L0 * U0 * beta0_ ) / ( Ln * Ln * Un *betan ) );
		Rgp1 		= Ra_ + ( Rg0_ - Ra_ ) * ( ( 1.0 - nn ) / nn ) * ( n0 / ( 1.0 - n0 ) );
		Cpp1 		= Ca_ + ( Cp0_ - Ca_ ) * ( ( 1.0 - nn ) / ( 1.0 - n0 ) );
		betap1 		= ( 1.0 / ( ( 1.0 - nn ) * ( 1.0 / sigma_ ) + ( nn * Rgn * thetan ) / Pn ) );

		Up1 		= Un + dU;
		betaUL2p1 	= betaUL2n + dbetaUL2;
		CpThetap1 	= CpThetan + dCpTheta;

		thetap1 	= CpThetap1 / Cpp1;
		Lp1 		= sqrt( betaUL2p1 / ( betap1 * Up1 ) );

		nn 			= np1;
		Rgn 		= Rgp1;
		betan 		= betap1;
		Un 			= Up1;
		betaUL2n 	= betaUL2p1;
		CpThetan 	= CpThetap1;
		thetan 		= thetap1;
		Ln 			= Lp1;

		Zn 			+= dx_;

		// store the profile
		uCentralProfile_	. push_back( Un );
		lProfile_			. push_back( Ln );

	}

}

VolcanicColumnProfile::~VolcanicColumnProfile() {}

double VolcanicColumnProfile::dUModelA( double U, double L, double alpha, double beta ){
	return - ( U / ( 8.0 * L ) ) * sqrt( alpha / beta ) + ( ( 9.81 * ( alpha - beta ) ) / ( beta * U ) );
}

double VolcanicColumnProfile::dbetaUL2ModelA( double U, double L, double alpha, double beta, double dU ){
	return ( ( 9.81 * ( alpha - beta ) * L * L ) / U ) - ( beta * L * L ) * dU;
}

double VolcanicColumnProfile::dCpThetaModelAB( double U, double L, double alpha, double beta, double dbetaUL2, double CpTheta, double Ca, double T ){
	return ( 1.0 / ( beta * U * L * L ) ) * ( ( ( Ca * T ) + ( ( U * U ) / 2.0 ) - CpTheta ) * dbetaUL2 - alpha * U * L * L * 9.81 );
}

double VolcanicColumnProfile::dUModelB( double U, double L, double alpha, double beta, double dbetaUL2 ){
	return ( 1.0 / ( beta * U * L * L ) ) * ( ( 9.81 * ( alpha - beta ) * L * L ) - U * dbetaUL2 );
}

double VolcanicColumnProfile::dbetaUL2ModelB( double U, double L, double alpha, double k ){
	return 2.0 * k * U * L * alpha;
}

std::vector<double> VolcanicColumnProfile::getUCentralProfile(){
	return uCentralProfile_;
}

std::vector<double> VolcanicColumnProfile::getLProfile(){
	return lProfile_;
}

double VolcanicColumnProfile::getHt(){
	return ventHeight_ + lProfile_.size();
}

double VolcanicColumnProfile::getN0(){
	return n0_;
}

double VolcanicColumnProfile::getTheta0(){
	return theta0_;
}

double VolcanicColumnProfile::getL( double z ){
	if( int( z ) >= int( lProfile_.size() ) || int( z ) < 0 ) return 0.0;
	return lProfile_[ int( z ) ];
}

double VolcanicColumnProfile::getU( double z ){
	if( int( z ) >= int( uCentralProfile_.size() ) || int( z ) < 0 ) return 0.0;
	return uCentralProfile_[ int( z ) ];
}
