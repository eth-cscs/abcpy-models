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


#include "../../include/Tools/FileManager/MPIFileManager.hpp"

namespace piaf{

MPIFileManager::MPIFileManager(MPITopology* topology):
    topology_(topology)
{}

MPIFileManager::~MPIFileManager() {}

/*
    void MPIFileManager::write(std::string filename, AbstractEulerianDomain *dom, bool sumParticles){
        //std::cout << "before gather rank "<< topology_->getRank() << std::endl;
        ( ( MPI3DBlockCyclicDomain* )dom )->gather();
        //std::cout << "after gather rank "<< topology_->getRank() << std::endl;
        if( ( ( MPI3DBlockCyclicDomain* )dom )->getTopology()->getRank() == 0 ){
            Domain gatheredDomain = ( ( MPI3DBlockCyclicDomain* )dom )->getGatheredDomain();
            FileManager::write( filename, &gatheredDomain, sumParticles );
        }
    }
    */

void MPIFileManager::write( std::string filename, std::vector< std::vector<double> > & domainParticlesNumber, Int3 dim ){
    if( topology_->getRank() == 0 ) FileManager::write( filename, domainParticlesNumber, dim );
}


void MPIFileManager::write( std::string filename, AbstractEulerianDomain *dom ){
    try{
        ((MPIAbstractEulerianDomain*)dom)->gather();
    }catch(exception& e){
        std::cout << "Exception : " << e.what() << std::endl;
    }
    if( topology_->getRank() == 0 ) {
        auto gatheredDomain = ((MPIAbstractEulerianDomain*)dom)->getGatheredDomain();
        FileManager::write( filename, &gatheredDomain );
    }
}


void MPIFileManager::write( std::string filename, MPIGridTerrain *terrain ){
    if( terrain->getTopology()->getRank() == 0 ) FileManager::write( filename, terrain );
    // else, still perform the collective operation
    else terrain->getParticleDeposition();
}


void MPIFileManager::writeHdf5( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles ){
    if( topology_->getRank() == 0 ) FileManager::writeHdf5( filename, trackedParticles );
}

void MPIFileManager::writeVtk( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles ){
    if( topology_->getRank() == 0 ) FileManager::writeVtk( filename, trackedParticles );
}

void MPIFileManager::writeCsv( std::string filename, const std::tuple< std::vector< double >, std::vector< std::vector< Particle > > >& trackedParticles ){
    if( topology_->getRank() == 0 ) FileManager::writeCsv( filename, trackedParticles );
}

void MPIFileManager::writeVelocityField( std::string filename, VelocityField& vc, double t, double dt ){
    if( topology_->getRank() == 0 ) FileManager::writeVelocityField(filename, vc, t, dt);
}

} // namespace piaf
