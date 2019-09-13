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

#include "Tools/Simulation/Simulation.hpp"

#ifdef MUSCLE
#include "Tools/Simulation/distributedsimulation.h"
#endif

#include "parseargs.hpp"


int main( int argc, char **argv ) {


#ifdef MUSCLE
    FactoryDistributedContinent * factory;
    DistributedContinent * sim;
    factory = new FactorMuscleDistributedContinent();
    sim = factory->newInstance(argc, argv);
#else
    FactoryLocalContinent * factory;
    LocalContinent * sim;
    factory=new FactoryLocalContinent();
    sim = factory->newInstance(argc, argv, MPI_COMM_WORLD, -1.0, -1.0);
#endif
    delete factory;
    sim->mainLoop();

    //MpiManager * mpiManager= & (sim->getMpiTopology()->getMpiManager());
    delete sim;
    //if(mpiManager)
     //   delete mpiManager;

}
