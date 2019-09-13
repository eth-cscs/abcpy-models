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


#ifdef CLI
#include "Tools/ConsoleInterface/ConsoleInterface.hpp"
#endif

#include "Tools/Simulation/Simulation.hpp"

#ifdef MUSCLE
    #include "Tools/Simulation/distributedsimulation.h"
#endif

#include "Tools/Test/test.cpp"

#include "parseargs.hpp"

void model( MPI_Comm communicator, double U0, double L0, int seed, char* outputfile, int len, int* success ){
    std::vector<std::string> arguments = {"./tetras", "-s", "EXACT", "-x", "60000", "-y", "60000", "-z", "30000", "-ox", "-30000", "-oy", "-30000", "-oz", "0", "-dx", "500", "-dt", "0.5", "-ddt", "-dx_t", "500", "-f", "/users/duttar/abcvolcano/eruption-data/eruptions/pululagua/pululagua2450BP-constantdensity-abcpy.csv", "-npart", "100000", "-blockcyclic", "5", "-o", std::string(outputfile), "-id", "-timeout", "20", "-determinist", std::to_string(seed)};

    //std::vector<std::string> arguments = {"./tetras", "-s", "EXACT", "-x", "60000", "-y", "60000", "-z", "30000", "-ox", "-30000", "-oy", "-30000", "-oz", "0", "-dx", "500", "-dt", "0.5", "-ddt", "-dx_t", "500", "-f", "/users/duttar/abcvolcano/eruption-data/eruptions/pululagua/pululagua2450BP-abcpy.csv", "-npart", "100000", "-blockcyclic", "5", "-o", std::string(outputfile), "-id", "-timeout", "20", "-determinist", std::to_string(seed)};

    std::vector<char*> argv;
    for (const auto& arg : arguments)
        argv.push_back((char*)arg.data());
    argv.push_back(nullptr);

    FactoryLocalPlume * factory;
    LocalPlume * sim;
    factory=new FactoryLocalPlume();

    sim = factory->newInstance(argv.size()-1, argv.data(), communicator, U0, L0);
    delete factory;
    
    *success = 1;
    
    try{
        sim->mainLoop();
    }
    catch(std::string const& e){
        std::cout << "Error during computation : " << e << std::endl;
        *success = 0;
    }
    
    delete sim;
}

/*
int main(int argc, char **argv){

  MPI_Init(NULL, NULL);

  MPI_Comm communicator = MPI_COMM_WORLD;
  double U0 = stod(argv[1]);
  double L0 = stod(argv[2]);
  int seed = atoi(argv[3]);
  char* outputfile = argv[4];

    std::vector<std::string> arguments = {"./tetras", "-s", "EXACT", "-x", "60000", "-y", "60000", "-z", "30000", "-ox", "-30000", "-oy", "-30000", "-oz", "0", "-dx", "500", "-dt", "0.5", "-ddt", "-dx_t", "500", "-f", "/users/duttar/abcvolcano/eruption-data/eruptions/pululagua/pululagua2450BP-constantdensity-abcpy.csv", "-npart", "100000", "-blockcyclic", "5", "-o", std::string(outputfile), "-id", "-timeout", "20", "-determinist", std::to_string(seed)};

    std::vector<char*> newargv;
    for (const auto& arg : arguments)
        newargv.push_back((char*)arg.data());
    newargv.push_back(nullptr);

    FactoryLocalPlume * factory;
    LocalPlume * sim;
    factory=new FactoryLocalPlume();

    sim = factory->newInstance(newargv.size()-1, newargv.data(), communicator, U0, L0);
    delete factory;

    try{
        sim->mainLoop();
    }
    catch(std::string const& e){
        std::cout << "Error during computation : " << e << std::endl;
    }
    
    delete sim;


  MPI_Finalize();

  return 0;

}
*/


int main( int argc, char **argv ) {

#ifdef MUSCLE
    FactoryDistributedPlume * factory;
    DistributedPlume * sim;
    factory = new FactorMuscleDistributedPlume();
    sim = factory->newInstance(argc, argv, MPI_COMM_WORLD, -1.0, -1.0);
#else
    FactoryLocalPlume * factory;
    LocalPlume * sim;
    factory=new FactoryLocalPlume();
    sim = factory->newInstance(argc, argv, MPI_COMM_WORLD, -1.0, -1.0);
#endif
    delete factory;

    sim->mainLoop();

    //MpiManager * mpiManager= & (sim->getMpiTopology()->getMpiManager());
    delete sim;
    //if(mpiManager)
       // delete mpiManager;

}

