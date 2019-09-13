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


#ifndef VOLCANOMPIFILEMANAGER_H_
#define VOLCANOMPIFILEMANAGER_H_

#include "VolcanoFileManager.hpp"
#include "../EruptionParameters/EruptionParameters.hpp"

#include "piaf.hpp"

using namespace piaf;

class VolcanoMPIFileManager : public MPIFileManager, public VolcanoFileManager {
public:
	VolcanoMPIFileManager(MPITopology* topology);
	virtual ~VolcanoMPIFileManager();

    virtual void write(string filename, MPIGridTerrain* terrain, std::vector<double> factors, EruptionParameters parameters, std::vector< int > injectedParticles, double dx, double dt, ESimulatorType simType, double eruptionTotalDuration, std::vector<std::string> vNames);

};

#endif /* VOLCANOMPIFILEMANAGER_H_ */
