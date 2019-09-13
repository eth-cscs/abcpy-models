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


#ifndef ERUPTIONPARAMETERS_H_
#define ERUPTIONPARAMETERS_H_

#include "../../Models/VolcanicTephra.hpp"
#include "boost/make_shared.hpp"

#include <tuple>

class EruptionParameters {
private:

	EruptionParameters( double ht, double hb, Double3 ventPosition, double duration, double eruptedMass,
		std::shared_ptr< VolcanicColumnProfile > columnProfile, double columnVerticalDiffusion, double columnHorizontalDiffusion,
		double eccentricity, Double2 focusSource,
		std::vector<std::shared_ptr<TephraParticleFamily> > particleFamilies,
		std::shared_ptr< Atmosphere > atmosphere, double maxWindSpeed, double windDirection, double atmosphereVerticalDiffusion, double atmosphereHorizontalDiffusion,
		bool useWoodsModel, std::vector< std::string > velocitiesNames, std::string windModel, Int3 eruptionDate
	);


	void setPlume( WimPlume plume );
	void setWind( WindProfile wind );
	void setEruptionEvents( std::vector< std::tuple< double, double > > eruptionEvents );

	// general parameters
	const double ht_;
	const double hb_;
	const Double3 ventPosition_;


	double duration_;
	double eruptedMass_;

	// volcanic column
	//const double U0_;
	//const double L0_;
	//const double n0_;
	//const double theta0_;
	const std::shared_ptr< VolcanicColumnProfile > columnProfile_;
	const double columnVerticalDiffusion_;
	const double columnHorizontalDiffusion_;

	// umbrella cloud
	const double eccentricity_;
	const Double2 focusSource_;

	// particle families
	const std::vector<std::shared_ptr<TephraParticleFamily> > particleFamilies_;

	// atmosphere
	//const double tropopauseHeight_;
	//const double thicknessSteadyZone_;
	//const double temperatureASL_;
	//const double pressureASL_;
	const std::shared_ptr< Atmosphere > atmosphere_;
	const double maxWindSpeed_;
	const double windDirection_;
	const double atmosphereVerticalDiffusion_;
	const double atmosphereHorizontalDiffusion_;

	const bool useWoodsModel_;


	WimPlume plume_;
	WindProfile wind_;
	std::vector< std::tuple< double, double > > eruptionEvents_;
	double mer_;

	std::vector< std::string > velocitiesNames_;

	std::string windModel_;

	Int3 eruptionDate_;

public:

	static EruptionParameters buildWoodsEruption( double ht, double hb, Double3 ventPosition, double duration,
		double eruptedMass, double U0, double L0, double n0, double theta0, double columnVerticalDiffusion,
		double columnHorizontalDiffusion, double eccentricity, Double2 focusSource,
		std::vector<std::shared_ptr<TephraParticleFamily> > particleFamilies, double tropopauseHeight,
		double topSteadyTemperatureHeight, double T0, double P0, double maxWindSpeed, double windDirection,
		double atmosphereVerticalDiffusion, double atmosphereHorizontalDiffusion, std::vector< std::string > velocitiesNames, Int3 eruptionDate
	);

	static EruptionParameters buildWimEruption( std::vector<std::shared_ptr<TephraParticleFamily> > &families,
		Double3 &ventPosition, double &columnVerticalDiffusion, double &columnHorizontalDiffusion, double &tropopause,
		double &stratosphere, double &T0, double &P0, double &atmosphereVerticalDiffusion,
		double &atmosphereHorizontalDiffusion, std::vector< std::tuple< double, double > > &eruptionEvents, WimPlume &plume,
		WindProfile &wind, double eccentricity, Double2 focusSource,std::vector< std::string > velocitiesNames, std::string windModel, Int3 eruptionDate
	);

	virtual ~EruptionParameters();

	double getHt();

	double getHb();

	Double3 getVentPosition();

	double getDuration();

    double getEndOfEruptionEvents();

	double getEruptedMass();

	std::shared_ptr< VolcanicColumnProfile > getColumnProfile();

	double getColumnVerticalDiffusion();

	double getColumnHorizontalDiffusion();

	double getEccentricity();

	Double2 getFocusSource();

	std::vector< std::shared_ptr< TephraParticleFamily> > getParticleFamilies();

	std::shared_ptr< Atmosphere > getAtmosphere();

	double getMaxWindSpeed();

	double getWindDirection();

	double getAtmosphereVerticalDiffusion();

	double getAtmosphereHorizontalDiffusion();

	double getParticlePerIteration(int familyId, double dt);

	double getFamilyTotalMass(int familyId);

	double getFamilyTotalParticles(int familyId);

	WimPlume getPlume();

	WindProfile getWind();

	std::vector< std::tuple< double, double > > getEruptionEvents();

	bool useWoodsModel();

	double getMer( double time );

	std::vector< std::string > getVelocitiesNames();

};

#endif /* ERUPTION_H_ */
