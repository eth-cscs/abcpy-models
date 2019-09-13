/*
TEphra TRAnsport Simulator (tetras)
Copyright (C) 2015  University of Geneva, Switzerland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef TOOLSTYPES_H_
#define TOOLSTYPES_H_

#include <string>
//#include "../Simulator/SimulatorTypes.hpp"
#include "../Models/VolcanicTephra/Tools/Atmosphere.hpp"

#include "piaf.hpp"

using namespace piaf;

struct params{
	ESimulatorType simulatorType_;
	std::string eruptionFile_;
	std::string outputFile_;
	std::string version_;
#ifdef DISPLAY
	bool display_;
	bool video_;
#endif
#ifdef CLI
	bool isInteractive_;
#endif
	int npart_;
	double dx_;
	double dt_;
	bool ddt_;
	double ddtValue_;
	double dxT_;
	//double simTime_;
	//double d_;
	double domainSizeX_;
	double domainSizeY_;
	double domainSizeZ_;
	//Atmosphere *atmosphere_;
	bool test_;
	int blockCyclic_;
	int simpleDom_;
	int seed_;
	bool id_;
	Double3 origin_;
	bool track_;
	bool nodif_;
	double trackinterval_;
	double trackwriteinterval_;
    std::string trackFileBaseName_;
    std::string trackFileType_;
	bool particlescsv_;
	bool injectAtCrater_;
	double timeOut_;

#ifdef MUSCLE
    bool isMuscle_;
#endif
};

#endif /* TOOLSTYPES_H_ */
