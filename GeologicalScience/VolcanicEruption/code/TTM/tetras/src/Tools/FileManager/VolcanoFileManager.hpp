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


#ifndef VOLCANOFILEMANAGER_H_
#define VOLCANOFILEMANAGER_H_

#include "../../Models/VolcanicTephra/VolcanicTephraParticleFamilies.hpp"
#include "../EruptionParameters/EruptionParameters.hpp"

#include "piaf.hpp"

using namespace piaf;

class VolcanoFileManager : public FileManager {
public:
	VolcanoFileManager();
	virtual ~VolcanoFileManager();

    virtual void write( string filename, GridTerrain* terrain, std::vector<double> factors, EruptionParameters parameters, std::vector< int > injectedParticles, double dx, double dt, ESimulatorType simType, double eruptionTotalDuration, std::vector<std::string> vNames );

	EruptionParameters loadEruptionParameters( string filename, double U0, double L0 );

private:

	template <class T>
	std::vector<T> readVector( std::vector<std::string> tokens );

	void readPlumeData( WimPlume& plume, string filename );

	void readWindData( WindProfile& wind, string filename );

	std::vector<std::shared_ptr<TephraParticleFamily> > interpolateParticleFamilies( std::vector<std::shared_ptr<TephraParticleFamily> > families, uint n );

	// read the parameters of the eruption
	// particles families must be loaded in increasing order of size
	void loadCommonEruptionParameters( string filename, std::vector<std::shared_ptr<TephraParticleFamily> > &families,
		Double3 &ventPosition, double &columnVerticalDiffusion, double &columnHorizontalDiffusion, double &tropopause,
		double &stratosphere, double &T0, double &P0, double &atmosphereVerticalDiffusion,
		double &atmosphereHorizontalDiffusion, std::vector< std::string > &velocitiesNames, std::string &windModel, Int3 eruptionDate );

	EruptionParameters loadWoodsEruptionParameters( string filename, double U0p, double L0p );

	EruptionParameters loadWimEruptionParameters( string filename );

};

#endif /* VOLCANOFILEMANAGER_H_ */
