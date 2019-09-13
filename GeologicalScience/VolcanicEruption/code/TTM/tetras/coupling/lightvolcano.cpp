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

@author: Mohamed Ben Belgacem
*/


#include "../src/Tools/Simulation/Simulation.hpp"

#include "../src/parseargs.hpp"
#include "../src/Tools/communication/mapper.h"
#include "musclehpc/conduits/mmsf.h"
#include "../../tetras/src/Tools/Simulation/atmosphere.h"
#include "memory.h"

using namespace std;

int main3( int argc, char **argv ) {

    //-------- MMSF_MPI handler --------------------
    // do it as first statement
    MMSF_MPI * mmsf = new MMSF_MPI(argc, argv, MPI_COMM_WORLD);

    //--------- create kernels ---------------------
    shared_ptr<FactoryDistributedAtmospheredPlume> fplume (make_shared< FactoryLightDistributedAtmospheredPlume>());
    shared_ptr<FactoryDistributedContinent> fcontinent (make_shared< FactoryLightDistributedContinent>());

    LightDistributedAtmospheredPlume * plume =
                 dynamic_cast<LightDistributedAtmospheredPlume *> (fplume->newInstance(argc, argv, MPI_COMM_WORLD, -1.0, -1.0));
    LightDistributedContinent * continent =
                 dynamic_cast<LightDistributedContinent *> (fcontinent->newInstance(argc, argv, MPI_COMM_WORLD, -1.0, -1.0));
    LightTransportInterpolation * interpoler = new LightTransportInterpolation();
    LightDistributedAtmosphereSubModel *atm = new LightDistributedAtmosphereSubModel(1.0, 2.0);

    //--------- register kernels -------------------
    mmsf->addKernel(plume, "S1");
    mmsf->addKernel(continent, "S2");
    mmsf->addKernel(interpoler, "G");
    mmsf->addKernel(atm, "atm");

    //--------- setup entrance connections ----------
    mmsf->loadCxACoupling();

    //--------- start simulation -------------------
    mmsf->compute();

    //--------- free kernels then the mmsf ---------
    delete plume;
    delete continent;
    delete interpoler;
    delete atm;

    //--------- do it as  last statement -----------
    delete mmsf;
    return 0;

}


int main( int argc, char **argv ) {

    //-------- MMSF_MPI handler --------------------
    // do it as first statement
    MMSF_MPI * mmsf = new MMSF_MPI(argc, argv, MPI_COMM_WORLD);

    //--------- create kernels ---------------------
    shared_ptr<FactoryDistributedPlume> fplume (make_shared< FactoryLightDistributedPlume>());
    shared_ptr<FactoryDistributedContinent> fcontinent (make_shared< FactoryLightDistributedContinent>());

    std::cout << "Create plume, continent and interpoler" << std::endl;

    LightDistributedPlume * plume =
                 dynamic_cast<LightDistributedPlume *> (fplume->newInstance(argc, argv, MPI_COMM_WORLD, -1.0, -1.0));
    LightDistributedContinent * continent =
                 dynamic_cast<LightDistributedContinent *> (fcontinent->newInstance(argc, argv, MPI_COMM_WORLD, -1.0, -1.0));
    LightTransportInterpolation * interpoler = new LightTransportInterpolation();

    std::cout << "Register kernels" << std::endl;

    //--------- register kernels -------------------
    mmsf->addKernel(plume, "S1");
    mmsf->addKernel(continent, "S2");
    mmsf->addKernel(interpoler, "G");

    //--------- setup entrance connections ----------
    /*mmsf->connect<char>("S1.f_out")->To("G.f1_in")->with_MTM_MPI("conduit1");
    mmsf->connect<char>("S2.f_out")->To("G.f2_in")->with_MTM_MPI("conduit2");
    mmsf->connect<char>("G.f1_out")->To("S1.f_in")->with_MTM_MPI("conduit3");
    mmsf->connect<char>("G.f2_out")->To("S2.f_in")->with_MTM_MPI("conduit4");*/

    std::cout << "Setup connections" << std::endl;

    //--------- setup entrance connections ----------
    mmsf->loadCxACoupling();

    std::cout << "Compute" << std::endl;

    //--------- start simulation -------------------
    mmsf->compute();

    //--------- free kernels then the mmsf ---------
    // delete plume;
    // delete continent;
    // delete interpoler;

    std::cout << "Done" << std::endl;

    //--------- do it as  last statement -----------
    delete mmsf;
    return 0;
}

int main2( int argc, char **argv ) {

    //-------- MMSF_MPI handler --------------------
    // do it as first statement
    MMSF_MPI * mmsf = new MMSF_MPI(argc, argv, MPI_COMM_WORLD);

    //--------- create kernels ---------------------
    LightDistributedAtmosphereRequester * requester = new LightDistributedAtmosphereRequester(1.0, 2.0);
    LightDistributedAtmosphereSubModel *responser = new LightDistributedAtmosphereSubModel(1.0, 2.0);

    //--------- register kernels -------------------
    mmsf->addKernel(requester, "S1", 3);
    mmsf->addKernel(responser, "S2", 4);

    //--------- setup entrance connections ----------
    /*mmsf->connect<char>("S1.f_out")->To("S2.f_in")->with_MTM_MPI("conduit1");
    mmsf->connect<char>("S2.f_out")->To("S1.f_in")->with_MTM_MPI("conduit2");*/

    //--------- setup entrance connections ----------
    mmsf->loadCxACoupling();

    //--------- start simulation -------------------
    mmsf->compute();

    //--------- free kernels then the mmsf ---------
    delete requester;
    delete responser;

    //--------- do it as  last statement -----------
    delete mmsf;
    return 0;
}
