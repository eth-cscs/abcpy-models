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


#include "EruptionParameters.hpp"

EruptionParameters::EruptionParameters( double ht, double hb, Double3 ventPosition, double duration, double eruptedMass,
	std::shared_ptr< VolcanicColumnProfile > columnProfile, double columnVerticalDiffusion, double columnHorizontalDiffusion,
	double eccentricity, Double2 focusSource,
	std::vector<std::shared_ptr<TephraParticleFamily> > particleFamilies,
	std::shared_ptr< Atmosphere >  atmosphere, double maxWindSpeed, double windDirection, double atmosphereVerticalDiffusion,
	double atmosphereHorizontalDiffusion, bool useWoodsModel, std::vector< std::string > velocitiesNames, std::string windModel, Int3 eruptionDate
):
ht_						( ht ),
hb_						( hb ),
ventPosition_			( ventPosition ),
duration_				( duration ),
eruptedMass_			( eruptedMass ),
columnProfile_			( columnProfile ),
columnVerticalDiffusion_		( columnVerticalDiffusion ),
columnHorizontalDiffusion_		( columnHorizontalDiffusion ),
eccentricity_			( eccentricity ),
focusSource_			( focusSource ),
particleFamilies_		( particleFamilies ),
atmosphere_				( atmosphere ),
maxWindSpeed_			( maxWindSpeed ),
windDirection_			( windDirection ),
atmosphereVerticalDiffusion_	( atmosphereVerticalDiffusion ),
atmosphereHorizontalDiffusion_	( atmosphereHorizontalDiffusion ),
useWoodsModel_ ( useWoodsModel ),
plume_ ( WimPlume( ventPosition ) ),
velocitiesNames_ ( velocitiesNames ),
windModel_ ( windModel ),
eruptionDate_ ( eruptionDate )
{
	if( !useWoodsModel ){
		duration_ = -1.0;
		eruptedMass_ = -1.0;
	}
}

EruptionParameters::~EruptionParameters() {}

void EruptionParameters::setPlume( WimPlume plume ){
	plume_ = plume;
}

void EruptionParameters::setWind( WindProfile wind ){
	wind_ = wind;
}

void EruptionParameters::setEruptionEvents( std::vector< std::tuple< double, double > > eruptionEvents ){
	eruptionEvents_ = eruptionEvents;
}


EruptionParameters EruptionParameters::buildWimEruption(  std::vector<std::shared_ptr<TephraParticleFamily> > &families,
	Double3 &ventPosition, double &columnVerticalDiffusion, double &columnHorizontalDiffusion, double &tropopause,
	double &stratosphere, double &T0, double &P0, double &atmosphereVerticalDiffusion,
	double &atmosphereHorizontalDiffusion, std::vector< std::tuple< double, double > > &eruptionEvents, WimPlume &plume,
	WindProfile &wind, double eccentricity, Double2 focusSource, std::vector< std::string > velocitiesNames, std::string windModel, Int3 eruptionDate
){

	std::shared_ptr< Atmosphere > atmosphere = std::make_shared< Atmosphere >( tropopause, stratosphere, T0, P0 );

	EruptionParameters ep(0.0, 0.0, ventPosition, 0.0, 0.0, NULL, columnVerticalDiffusion, columnHorizontalDiffusion,
		eccentricity, focusSource,  families, atmosphere, 0.0, 0.0, atmosphereVerticalDiffusion, atmosphereHorizontalDiffusion,
		false, velocitiesNames, windModel, eruptionDate
	);

	ep.setEruptionEvents( eruptionEvents );
	ep.setWind( wind );
	ep.setPlume( plume );

	return ep;

}

EruptionParameters EruptionParameters::buildWoodsEruption( double ht, double hb, Double3 ventPosition, double duration, double eruptedMass,
	double U0, double L0, double n0, double theta0, double columnVerticalDiffusion, double columnHorizontalDiffusion,
	double eccentricity, Double2 focusSource,
	std::vector<std::shared_ptr<TephraParticleFamily> > particleFamilies,
	double tropopauseHeight, double topSteadyTemperatureHeight, double T0, double P0, double maxWindSpeed, double windDirection,
	double atmosphereVerticalDiffusion, double atmosphereHorizontalDiffusion, std::vector< std::string > velocitiesNames, Int3 eruptionDate
)
{

	std::cout << "Before computation, U0 : " << U0 << ", L0 : " << L0 << std::endl;

	/*! \todo handle with exception */
	if( ht < 0.0 && ( U0 < 0.0 || L0 < 0.0 || n0 < 0.0 || theta0 < 0.0 ) ){
		std::cout << "unable to construct a volcanic column profile with given parameters (HT is undefined and one or more column parameter is undefined)" << std::endl;
		exit( 0 );
	}

	std::shared_ptr< Atmosphere > atmosphere = std::make_shared< Atmosphere >( tropopauseHeight, topSteadyTemperatureHeight, T0, P0 );

	double htcalc;
	double hbcalc;
	double U0calc;
	double L0calc;
	double n0calc;
	double theta0calc;
	hbcalc = hb;
	htcalc = ht;
	U0calc = ( U0 > 0.0 ) ? U0 : double( ( rand() % 391 ) + 10 ) ;
	L0calc = ( L0 > 0.0 ) ? L0 : double( ( rand() % 181 ) + 20 ) ;
	n0calc = ( n0 > 0.0 ) ? n0 : ( ( double( rand() ) / double( RAND_MAX ) ) * 0.04 ) + 0.01 ;
	theta0calc = ( theta0 > 0.0 ) ? theta0 : double( ( rand() % 201 ) + 1100 ) ;
	std::shared_ptr< VolcanicColumnProfile > columnProfile = std::make_shared< VolcanicColumnProfile >( U0calc, L0calc, n0calc, theta0calc, ventPosition.z_, atmosphere );

	if( ht < 0.0 ) {
		htcalc = columnProfile->getHt();
	}
	else{

		int cpt = 0;

		double U0avg = 0.0;
		double L0avg = 0.0;
		double n0avg = 0.0;
		double theta0avg = 0.0;

		// take 10 valid column configurations, keep the sum in order to compute the average if wanted

		bool columnOK = false;

		while( !columnOK ){

			for( int i = 0; i < 10; i++ ){

				U0calc = ( U0 > 0.0 ) ? U0 : double( ( rand() % 391 ) + 10 ) ;
				L0calc = ( L0 > 0.0 ) ? L0 : double( ( rand() % 181 ) + 20 ) ;
				n0calc = ( n0 > 0.0 ) ? n0 : ( ( double( rand() ) / double( RAND_MAX ) ) * 0.04 ) + 0.01 ;
				theta0calc = ( theta0 > 0.0 ) ? theta0 : double( ( rand() % 201 ) + 1100 ) ;

				columnProfile = std::make_shared< VolcanicColumnProfile >( U0calc, L0calc, n0calc, theta0calc, ventPosition.z_, atmosphere );

				while( abs( columnProfile->getHt() - htcalc ) > 100.0 ) {

					U0calc = ( U0 > 0.0 ) ? U0 : double( ( rand() % 391 ) + 10 ) ;
					L0calc = ( L0 > 0.0 ) ? L0 : double( ( rand() % 181 ) + 20 ) ;
					n0calc = ( n0 > 0.0 ) ? n0 : ( ( double( rand() ) / double( RAND_MAX ) ) * 0.04 ) + 0.01 ;
					theta0calc = ( theta0 > 0.0 ) ? theta0 : double( ( rand() % 201 ) + 1100 ) ;
					columnProfile = std::make_shared< VolcanicColumnProfile >( U0calc, L0calc, n0calc, theta0calc, ventPosition.z_, atmosphere );
					cpt++;

					/*! \todo handle with exception */
					if( cpt == 1000 ){
						std::cout << "unable to construct a volcanic column profile with given parameters (HT seems not compatible with other column parameters)" << std::endl;
						exit( 0 );
					}
				}
				//std::cout << "column generated : U0 = " << columnProfile->getU(0) << " - L0 = " << columnProfile->getL(0) << " - n0 = " << columnProfile->getN0() << " - theta0 = " << columnProfile->getTheta0() << " - Ht = " << columnProfile->getHt() << std::endl;
				U0avg += U0calc;
				L0avg += L0calc;
				n0avg += n0calc;
				theta0avg += theta0calc;
			}
			U0avg = U0avg / 10.0;
			L0avg = L0avg / 10.0;
			n0avg = n0avg / 10.0;
			theta0avg = theta0avg / 10.0;

			//columnProfile = std::make_shared< VolcanicColumnProfile >( U0avg, L0avg, n0avg, theta0avg, ventPosition.z_, atmosphere );

			//std::cout << "averaged column generated : U0 = " << columnProfile->getU(0) << " - L0 = " << columnProfile->getL(0) << " - n0 = " << columnProfile->getN0() << " - theta0 = " << columnProfile->getTheta0() << " - Ht = " << columnProfile->getHt() << std::endl;

			// don't use the average but the last value
			columnProfile = std::make_shared< VolcanicColumnProfile >( U0calc, L0calc, n0calc, theta0calc, ventPosition.z_, atmosphere );

			/*! \todo handle with exception */
			if( abs( columnProfile->getHt() - htcalc ) > 100.0 ){
				std::cout << "averaging parameters gives an invalid column... retrying" << std::endl;
				//exit( 0 );
			}



			else columnOK = true;

		}

	}

	std::cout << "Value of U0 used : " << columnProfile->getU(0) << ", value of L0 used : " << columnProfile->getL(0) << std::endl;

	// bonnadona 2003
	if( ht < 0.0 ) hbcalc = 0.7 * htcalc;

	EruptionParameters ep( htcalc, hbcalc, ventPosition, duration, eruptedMass,
		columnProfile, columnVerticalDiffusion, columnHorizontalDiffusion,
		eccentricity, focusSource,
		particleFamilies,
		atmosphere, maxWindSpeed, windDirection, atmosphereVerticalDiffusion, atmosphereHorizontalDiffusion, true,
		velocitiesNames, "", eruptionDate
	);

	std::vector< std::tuple< double, double > > eruptionEvents_;

	ep.eruptionEvents_.push_back( make_tuple( 0.0, duration ) );
	ep.mer_ = eruptedMass / duration;

	return ep;
}

double EruptionParameters::getHt(){
	if( !useWoodsModel_ ) throw std::runtime_error("no ht parameter when using Wim model");
	return ht_;
}

double EruptionParameters::getHb(){
	if( !useWoodsModel_ ) throw std::runtime_error("no hb parameter when using Wim model");
	return hb_;
}

Double3 EruptionParameters::getVentPosition(){
	return ventPosition_;
}

double EruptionParameters::getDuration(){
	if( !useWoodsModel_ ){
		if( duration_ > 0.0 ) return duration_;
		duration_ = 0.0;
		for( std::tuple< double, double > ev : eruptionEvents_ ) duration_ += get<1>(ev);
	}
	return duration_;
}

double EruptionParameters::getEndOfEruptionEvents(){
    if( !useWoodsModel_ ){
        double max = 0.0;
        for( std::tuple< double, double > ev : eruptionEvents_ ){
            if( max < get<0>(ev) + get<1>(ev) ){
                max = get<0>(ev) + get<1>(ev);
            }
        }
        return max;
    }
    return duration_;
}

double EruptionParameters::getEruptedMass(){
	if( !useWoodsModel_ ) {
		if( eruptedMass_ > 0.0 ) return eruptedMass_;
		eruptedMass_ = getPlume().getMs( 0.0 )*getDuration();
	}
	return eruptedMass_;
}

std::shared_ptr< VolcanicColumnProfile > EruptionParameters::getColumnProfile(){
	if( !useWoodsModel_ ) throw std::runtime_error("no column profile parameter when using Wim model");
	return columnProfile_;
}

double EruptionParameters::getColumnVerticalDiffusion(){
	return columnVerticalDiffusion_;
}

double EruptionParameters::getColumnHorizontalDiffusion(){
	return columnHorizontalDiffusion_;
}

double EruptionParameters::getEccentricity(){
	return eccentricity_;
}

Double2 EruptionParameters::getFocusSource(){
	return focusSource_;
}

std::vector<std::shared_ptr<TephraParticleFamily> > EruptionParameters::getParticleFamilies(){
	return particleFamilies_;
}

std::shared_ptr< Atmosphere > EruptionParameters::getAtmosphere(){
	return atmosphere_;
}

double EruptionParameters::getMaxWindSpeed(){
	if( !useWoodsModel_ ) throw std::runtime_error("no maximum wind speed parameter when using Wim model");
	return maxWindSpeed_;
}

double EruptionParameters::getWindDirection(){
	if( !useWoodsModel_ ) throw std::runtime_error("no wind direction parameter when using Wim model");
	return windDirection_;
}

double EruptionParameters::getAtmosphereVerticalDiffusion(){
	return atmosphereVerticalDiffusion_;
}

double EruptionParameters::getAtmosphereHorizontalDiffusion(){
	return atmosphereHorizontalDiffusion_;
}

double EruptionParameters::getParticlePerIteration(int familyId, double dt){
	//if( !useWoodsModel_ ) throw std::runtime_error("unable to use getParticlePerIteration when using Wim model");
	return ( ( particleFamilies_[ familyId ]->originalWeightPercent_ / 100.0) * getEruptedMass() ) / ( particleFamilies_[ familyId ]->mass_ ) / ( getDuration() / dt);
}

double EruptionParameters::getFamilyTotalMass(int familyId){
	//if( !useWoodsModel_ ) throw std::runtime_error("unable to use getFamilyTotalMass when using Wim model");
	return ( particleFamilies_[ familyId ]->originalWeightPercent_ / 100.0) * getEruptedMass();
}

double EruptionParameters::getFamilyTotalParticles(int familyId){
	//if( !useWoodsModel_ ) throw std::runtime_error("unable to use getFamilyTotalParticles when using Wim model");
	return ( ( particleFamilies_[ familyId ]->originalWeightPercent_ / 100.0) * getEruptedMass() ) / ( particleFamilies_[ familyId ]->mass_ );
}

WimPlume EruptionParameters::getPlume(){
	if( useWoodsModel_ ) throw std::runtime_error("no Wim plume profile when using woods model");
	return plume_;
}

WindProfile EruptionParameters::getWind(){
	if( useWoodsModel_ ) throw std::runtime_error("no wind profile when using woods model");
	return wind_;
}

std::vector< std::tuple< double, double > > EruptionParameters::getEruptionEvents(){
	if( useWoodsModel_ ) throw std::runtime_error("no eruption events when using woods model");
	return eruptionEvents_;
}

bool EruptionParameters::useWoodsModel(){
	return useWoodsModel_;
}

double EruptionParameters::getMer( double t ){

	//std::cout << " t : " << t << std::endl;
	//std::cout << "eruption events size : " << eruptionEvents_.size() << std::endl;

	if( eruptionEvents_.size() == 0 || ( t > get<0>( eruptionEvents_.back() ) + get<1>( eruptionEvents_.back() ) ) ) return -1.0;

	for( auto ee : eruptionEvents_ ){
		if( t >= get<0>( ee ) &&  t <= get<0>( ee )+get<1>( ee ) ) {
			if( useWoodsModel_ ) return mer_;
			return plume_.getMs( t );
		}
	}

	return 0.0;
}

std::vector< std::string > EruptionParameters::getVelocitiesNames(){
	return velocitiesNames_;
}
