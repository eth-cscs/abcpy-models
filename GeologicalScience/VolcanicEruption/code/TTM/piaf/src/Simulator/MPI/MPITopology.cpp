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

#include "../../../include/Simulator/MPI/MPITopology.hpp"

namespace piaf{

MPITopology::MPITopology(int argc, char** argv,  shared_ptr<MpiManager>mpiManager)
{
    this->mpiManager=mpiManager;
}

MPITopology::~MPITopology() {}

void MPITopology::init(){

}

int MPITopology::getSize(){
    return this->mpiManager->getSize();
}

int MPITopology::getRank(){
    return this->mpiManager->getRank();
}


/*MpiManager * MPITopology::getMpiManagerPtr() const{
    return (this->mpiManager);
}*/

MpiManager & MPITopology::getMpiManager() const{
    return *(this->mpiManager.get());
}



void MPITopology::setMpiManager(shared_ptr<MpiManager> newmpiManager){

    //if(this->mpiManager)
      //  delete this->mpiManager;
    this->mpiManager= newmpiManager;
}

} // namespace piaf
