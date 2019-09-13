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


#ifndef TESTCPP
#define TESTCPP

#ifdef DISPLAY
#include "../Display2D/Display2D.hpp"
#endif

#include "../ToolsTypes.hpp"
#include "../../Models/VolcanicTephra.hpp"
#include "../FileManager/VolcanoFileManager.hpp"

#include "piaf.hpp"

#include "boost/scoped_ptr.hpp"

class TestSpeed: public SpeedFunctor{
public:
	Double3 compute( Double3 pos, double t, std::shared_ptr<GenericParticleFamily> family ){
		Double3 res = { 0.0, 0.0, -100.0 };
		return res;
	}
};

inline void test( params p, MPITopology* topology ){
/*
	// General initializations
	std::shared_ptr< TestSpeed > testSpeed( new TestSpeed() );
	testSpeed->setSpeedType( LAGRANGIAN_DYNAMIC );
	//Double3 ventPosition = { p.domainSizeX_ / 2.0, p.domainSizeY_ / 2.0, 500.0 };
	Double3 ventPosition = { p.domainSizeX_ / 2.0, p.domainSizeY_ / 2.0, p.domainSizeZ_ / 2.0 };
	FileManager fm;

	double D = 0.01;
	std::shared_ptr<DiffusionSpeed> diffusion(new DiffusionSpeed());
	//diffusion->v_ = p.dx_/p.dt_;
	//diffusion->r_ = (D * 3.14159265 * sqrt(p.dt_));
	diffusion->urAtm_ = sqrt( ( 2.0 * D ) / p.dt_ );
	diffusion->setSpeedType(LAGRANGIAN_DYNAMIC);

	std::shared_ptr< Domain > dom;// Domain( { 0.0, 0.0, 0.0 }, { p.domainSizeX_ , p.domainSizeY_ , p.domainSizeZ_ } , p.dx_ );

	// in order to make the member function "equals" callable by ervery process ( but the result is relevant only for process number 0 )
	dom = std::shared_ptr< Domain >( new Domain( { 0.0, 0.0, 0.0 }, { p.domainSizeX_ , p.domainSizeY_ , p.domainSizeZ_ } , p.dx_ ) );

	std::shared_ptr< AbstractSimulator > sim;
	std::shared_ptr< ParticleRepository > repository;
	std::shared_ptr< GridTerrain > terrain;

	int cpt = 0;

	int totalPart = 1000000;
	int iterations = 1000;
*/
	/*if( topology->getRank() == 0 ){
		// Sequential simulation initializations
		Eruption eruption = fm.loadEruption( p.eruptionFile_, dom.get() );
		dom->addSpeed( testSpeed );
		dom->addSpeed( diffusion );
		repository = std::shared_ptr< ParticleRepository >( new ParticleRepository( dom->getParticleFamilies() ) );
		terrain = std::shared_ptr< GridTerrain >( new GridTerrain( { p.domainSizeX_, p.domainSizeY_ }, p.dx_, dom->getParticleFamilies() ) );
		repository->setActive( true );
		terrain->setActive( true );
		sim = std::shared_ptr< AbstractSimulator >( new ExactSimulator( dom.get(), p.dt_, terrain.get(), repository.get() ) );

		// put some particles int the domain
		//dom->insertParticles( ventPosition, 1000000, 0 );

		//Display2D disp(dom.get(), NULL, (p.domainSizeY_/2.0)/p.dx_, false, 30);
		//disp.show();

		//for( int i = 0; i < iterations; i++ ){

		while( true ){
			if( cpt < iterations )dom->insertParticles( ventPosition, totalPart / iterations, 0 );
			std::cout << cpt << std::endl;
			if( cpt % 10 == 0 ){
				if( ! dom->containsParticles() ) break;
			}
			sim->step();
			//dom->computeMaxNumPart();
			//disp.show();
			cpt++;
		}
		// just one class and we don't care about scaling to correct mass for the test
		std::vector<double> factors;
		factors.push_back( 1.0 );
		fm.write("test_output_seq.h5",terrain.get(), factors);
	}*/

/*
	// MPI simulation initializations
	MPI3DBlockCyclicDomain mpiDom = MPI3DBlockCyclicDomain( { 0.0 , 0.0 , 0.0 } , { p.domainSizeX_ , p.domainSizeY_ , p.domainSizeZ_ } , p.dx_ , topology );
	EruptionParameters mpiEruption = fm.loadEruptionParameters( p.eruptionFile_ );
	mpiDom.addSpeed( testSpeed );
	mpiDom.addSpeed( diffusion );
	std::shared_ptr< ParticleRepository > mpiRepository = std::shared_ptr< ParticleRepository >( new ParticleRepository( mpiDom.getParticleFamilies() ) );
	std::shared_ptr< GridTerrain > mpiTerrain = std::shared_ptr< GridTerrain >( new GridTerrain( { 0.0, 0.0 }, { p.domainSizeX_, p.domainSizeY_ }, p.dx_, mpiDom.getParticleFamilies() ) );
	mpiRepository->setActive( true );
	mpiTerrain->setActive( true );
	std::shared_ptr< AbstractSimulator > mpiSim = std::shared_ptr< AbstractSimulator >( new MPIExactSimulator( &mpiDom, p.dt_, mpiTerrain.get(), mpiRepository.get() ) );

	// put some particles int the domain
	//mpiDom.insertParticles( ventPosition, totalPart, 0 );

	cpt = 0;
	while( true ){
	//for( int i = 0; i < iterations; i++ ){
		if( cpt < iterations ) mpiDom.insertParticles( ventPosition, totalPart / iterations , 0 );
		if( topology->getRank() == 0 )std::cout << cpt << std::endl;
		if( cpt % 10 == 0 ){
			if( ! mpiDom.containsParticles() ) {
				//std::cout << "le domaine ne contient plus de particule" << std::endl;
				break;
			}
		}
		mpiSim->step();
		//dom->computeMaxNumPart();
		//disp.show();
		cpt++;
	}


	std::vector<double> factors;
	factors.push_back( 1.0 );
	//fm.write("test_output_mpi.h5",mpiTerrain.get(), factors);
*/
	/*bool testSucceeded = true;
	//for( int i = 0; i < 100; i++ )mpiSim->step();

	for( int i = 0; i < 100; i++ ){

		std::cout << i << std::endl;


		if( topology->getRank() == 0 ) sim->step();

		//if( topology->getRank() == 0 ) dom->test();

		mpiSim->step();

		if( topology->getRank() == 0 ) dom->test();

		//if(dom->equals( &mpiDom )){}

		if( ! dom->equals( &mpiDom ) ) {
			if( topology->getRank() == 0 && testSucceeded ) {
				std::cout << "dom != mpiDom after " << i+1 << " iterations" << std::endl;
				testSucceeded = false;
			}
			//return;
		}
	}

	if( topology->getRank() == 0 && testSucceeded ) std::cout << "test succeeded" << std::endl;
	else if ( topology->getRank() == 0 ) std::cout << "test failed" << std::endl;

	 */

}





#endif
