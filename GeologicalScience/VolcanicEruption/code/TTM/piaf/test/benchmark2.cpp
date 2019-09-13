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

#include "piaf.hpp"
#include <numeric>

struct params{
  double time_;
  std::string filename_;
  bool ddt_;
};

void printUsage(){
  std::cout<<"Usage: mpirun -np [n] -t [TIME] ..."<<std::endl;
  std::cout<<" -h, --help         help, show this message"<<std::endl;
  std::cout<<" -o [filename]      output filename"<<std::endl;
  std::cout<<" -ddt               enable dynamic delta t adjustement"<<std::endl;
}

void parseArgs(int argc, char **argv, params &p){
  p.time_ = -1;
  p.filename_ = "\0";
  p.ddt_ = false;

  for (int i = 1; i < argc; i++) {

    if((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"--help") == 0)){
      printUsage();
      throw -1;
    }

    if(strcmp(argv[i],"-t") == 0){
      try{
        std::string text = argv[i+1];
        p.time_ = boost::lexical_cast< double >( text );
      }
      catch(...){
        std::cout<<"Invalid value for parameter t"<<std::endl<<std::endl;
        printUsage();
        throw -1;
      }
      i++;
    }

    if(strcmp(argv[i],"-o") == 0){
      std::string text = argv[i+1];
      p.filename_ = text;
      i++;
    }

    if(strcmp(argv[i],"-ddt") == 0){
      p.ddt_ = true;
    }

  }

  if(p.time_ < 0){
    std::cout<<"please choose a valid value for parameter t"<<std::endl; printUsage(); throw -1;
  }

}

int main( int argc, char **argv ) {

  params p;
  parseArgs( argc, argv, p );

  // parameters
  const double dx = 10.0;
  const double dt = 1.0;
  const double duration = 10.0; // duration of the "eruption"
  const piaf::Double3 speed = { 0.0, 0.0, -1.0 };
  const double simTime = p.time_;
  const piaf::Double3 position = { 0.0, 0.0, 0.0 };
  const piaf::Double3 size = { 1000.0, 1000.0, 1000.0 };
  const bool ddt = p.ddt_;
  const double ddtValue = 0.1;
  std::string filename = p.filename_;
  int nFam = 10;
  double diffusion = 5.0;
  piaf::Double3 ventPos = { 500.0, 500.0, 900.0 };
  std::vector<int> partNumbers;
  //for(int i = 0; i<nFam; i++)partNumbers.push_back( i*1000 );
  for(int i = 0; i<nFam; i++)partNumbers.push_back( 100000 );
  std::vector<int> insertedPart( nFam, 0 );

  boost::posix_time::ptime timeStart;
  boost::posix_time::ptime timeEnd;
  boost::posix_time::time_duration timeDif;

  timeStart = boost::posix_time::ptime(boost::posix_time::microsec_clock::local_time());

  // create topology, speed, domain, simulator and particle families (mendatory data structures for a simulation)
  piaf::MPITopology* topology = new piaf::MPI3DGridTopology( argc, argv );
  std::shared_ptr< piaf::ConstantSpeed > constantSpeed( new piaf::ConstantSpeed( speed ) );
  std::shared_ptr< piaf::Diffusion > diffusionSpeed( new piaf::Diffusion( diffusion ) );
  piaf::MPIAbstractDomain* domain = new piaf::MPI3DBlockCyclicDomain( position , size , dx , (piaf::MPI3DGridTopology*)topology );
  std::vector< std::shared_ptr< piaf::GenericParticleFamily > > families;

  for( int i=0; i<nFam; i++ ) families.push_back( std::make_shared<piaf::GenericParticleFamily>( piaf::GenericParticleFamily( i ) ) );

  domain->addSpeed( constantSpeed );
  domain->addSpeed( diffusionSpeed );
  for( int i=0; i<nFam; i++ ) domain->addParticleFamily( families[ i ] );

  piaf::MPIGridTerrain* terrain = NULL;
  piaf::MPIParticleRepository* repository = NULL;

  if( filename != "\0"){
    terrain = new piaf::MPIGridTerrain( { 0.0 , 0.0 }, { size.x_ , size.y_ }, dx, domain->getParticleFamilies(), (piaf::MPI3DGridTopology*)topology );
    terrain->setActive( true );
    repository = new piaf::MPIParticleRepository( domain->getParticleFamilies(),(piaf::MPI3DGridTopology*)topology );
    repository->setActive( true );
  }

  piaf::MPIAbstractSimulator* simulator = new piaf::MPIExactSimulator( domain, 0.0, dt, terrain, repository, ddt );

  if( ddt ) simulator->setDdtValue( ddtValue );

  timeEnd = boost::posix_time::ptime(boost::posix_time::microsec_clock::local_time());
  timeDif = boost::posix_time::time_duration(timeEnd - timeStart);

  std::cout<<"initialization time : " << timeDif << std::endl;

  timeStart = boost::posix_time::ptime(boost::posix_time::microsec_clock::local_time());

  int incr, rest;
  int it = 0;

  while( simulator->getTime() < p.time_ ){
    for(int i=0; i<partNumbers.size(); i++){
      incr = partNumbers[i] / ( duration/simulator->getDt() );
      rest = partNumbers[i] - insertedPart[i];
      if(incr == 0 || incr > rest) incr = rest;
      insertedPart[i] += incr;
      domain->insertParticles( ventPos, incr,i );
    }
    simulator->step();
    it++;
  }

  timeEnd = boost::posix_time::ptime(boost::posix_time::microsec_clock::local_time());
  timeDif = boost::posix_time::time_duration(timeEnd - timeStart);

  if( topology->getRank() == 0 ){
    std::cout << it << " iterations completed - computation time : " << timeDif << " ( " << timeDif.total_seconds() << " [s] ) with " << topology->getSize() << " processors" << ", simTime : " << simTime << std::endl << ", total injected particles : "<< std::accumulate(insertedPart.begin(), insertedPart.end(), 0) << ", domain size : (" << size.x_ << ", " << size.y_ << ", " << size.z_ << "), dx : " << dx << ", domain discrete size : (" << domain->getSize().x_ << ", " << domain->getSize().y_ << ", " << domain->getSize().z_ << ")" << std::endl ;
  }

  if( filename != "\0"){
    piaf::MPIFileManager fm(topology);
    fm.write( filename, terrain );
  }

  delete terrain;
  delete repository;
  delete simulator;
  delete topology;

  return 0;

}
