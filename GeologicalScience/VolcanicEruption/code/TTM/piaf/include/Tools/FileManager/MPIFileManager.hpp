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

#ifndef MPIFILEMANAGER_H_
#define MPIFILEMANAGER_H_

#include "FileManager.hpp"
#include "../../Simulator/MPI/MPI3DBlockCyclicDomain.hpp"
#include "../../Simulator/MPI/MPIGridTerrain.hpp"


namespace piaf{

class MPIFileManager : public FileManager {
public:
    MPIFileManager( MPITopology* topology );
    virtual ~MPIFileManager();

    //virtual void write( std::string filename, AbstractEulerianDomain *dom, bool sumParticles );

    virtual void write( std::string filename, std::vector< std::vector<double> > & domainParticlesNumber, Int3 dim );

    virtual void write( std::string filename, AbstractEulerianDomain *dom );

    virtual void write( std::string filename, MPIGridTerrain *terrain );

    virtual void writeHdf5( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles );

    virtual void writeVtk( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles );

    virtual void writeCsv( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles );

    virtual void writeVelocityField( std::string filename, VelocityField& vc, double t, double dt );

private:

    MPITopology* topology_;

};

} // namespace piaf
#endif /* MPIFILEMANAGER_H_ */
