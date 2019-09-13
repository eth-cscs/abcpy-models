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


#include "VolcanoMPIFileManager.hpp"

VolcanoMPIFileManager::VolcanoMPIFileManager(MPITopology* topology):
MPIFileManager( topology )
 {}

VolcanoMPIFileManager::~VolcanoMPIFileManager() {}

void VolcanoMPIFileManager::write( string filename, MPIGridTerrain* terrain, std::vector<double> factors, EruptionParameters parameters, std::vector< int > injectedParticles, double dx, double dt, ESimulatorType simType, double eruptionTotalDuration, std::vector<std::string> vNames ){
    if( terrain->getTopology()->getRank() == 0 ) VolcanoFileManager::write( filename, terrain, factors, parameters, injectedParticles, dx, dt, simType, eruptionTotalDuration, vNames );
	// else, still perform the collective operation
	else terrain->getParticleDeposition();
}
