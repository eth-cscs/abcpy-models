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


#ifndef FILEMANAGER_H_
#define FILEMANAGER_H_

#include <string>
#include <fstream>
#include <typeinfo>
#include <tuple>

#include <boost/algorithm/string.hpp>
#include <boost/make_shared.hpp>

#include "hdf5.h"
#include "hdf5_hl.h"

#include "../../Simulator/AbstractEulerianDomain.hpp"
#include "../../Simulator/Domain.hpp"

#include "../../Simulator/VelocityField.hpp"

namespace piaf{

class FileManager {

public:
	FileManager();
	virtual ~FileManager();

	//virtual void write( std::string filename, AbstractEulerianDomain *dom, bool sumParticles );

    virtual void write( std::string filename, std::vector< std::vector<double> > & domainParticlesNumber, Int3 dim );

    virtual void write( std::string filename, AbstractEulerianDomain *dom );

	virtual void write( std::string filename, GridTerrain *terrain );

    virtual void writeHdf5( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles );

    virtual void writeVtk( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles );

    virtual void writeCsv( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles );

    virtual void writeParticlesCsv( std::string filename,  std::vector< Particle >& particles );

    virtual void writeVelocityField( std::string filename, VelocityField& vc, double t, double dt );

private:

    int vtkTrackFileCounter_ = 0;
    int hdf5TrackFileCounter_ = 0;
    int csvTrackFileCounter_ = 0;

};

} // namespace piaf

#endif /* FILEMANAGER_H_ */
