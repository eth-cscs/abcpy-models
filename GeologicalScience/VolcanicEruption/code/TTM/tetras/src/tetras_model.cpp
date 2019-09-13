#include "mpi.h"
#include <iostream>

int model( double* results, unsigned int rsize, MPI_Comm communicator, double U0, double L0 ){
    
    std::cout << "Hello from model" << std::endl;
    
    // FactoryLocalPlume * factory;
    // LocalPlume * sim;
    // factory=new FactoryLocalPlume();
    // sim = factory->newInstance(argc, argv, communicator, U0, L0);
    // delete factory;
    // sim->mainLoop();
    // delete sim;
    return 0;
}