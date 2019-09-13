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


#include "../../../include/Simulator/MPI/MPIAbstractEulerianDomain.hpp"

#include <exception>

namespace piaf{

MPIAbstractEulerianDomain::MPIAbstractEulerianDomain( Double3 pos, Double3 size, double dx , MPITopology* topology ):
    AbstractEulerianDomain( pos, size, dx ),
    topology_ (topology)
{
    // TODO Auto-generated constructor stub

}

MPIAbstractEulerianDomain::MPIAbstractEulerianDomain() {

}

MPIAbstractEulerianDomain::~MPIAbstractEulerianDomain() {
    // TODO Auto-generated destructor stub
}

std::vector< std::string > MPIAbstractEulerianDomain::getSpeedsName(){
    std::vector< std::string > res;
    if( domainBlocks_.size() > 0 ) res = domainBlocks_[0].getSpeedsName();
    return res;
}

std::vector< piaf::Particle > MPIAbstractEulerianDomain::getParticlesInBox( piaf::Double3 p1, piaf::Double3 p2 ){
    std::vector< piaf::Particle > resLoc;
    std::vector< piaf::Particle > res;

    for(Domain& d : domainBlocks_){
        std::vector< piaf::Particle > box = d.getParticlesInBox( p1, p2 );
        resLoc.insert( resLoc.end(), box.begin(), box.end() );
    }

    if(topology_->getSize() == 1) return resLoc;

    std::vector<int> sizes;
    if( topology_->getRank() == 0 ) sizes.resize(topology_->getSize());
    int resLocSize = resLoc.size();
 
    MpiManager & mpiManager= this->topology_->getMpiManager();
    mpiManager.gather<int>(&resLocSize, 1, sizes.data(), 1, 0);
    int globalSize = 0;
    for(int i : sizes) globalSize += i;
    std::vector<int> displs(1,0);
    if( topology_->getRank() == 0 ) for(uint i=0; i<sizes.size()-1; i++) displs.push_back(displs.at(i)+sizes.at(i)*sizeof(Particle));

    res.resize(globalSize);

    for(auto& i : sizes) i *= sizeof(Particle);

    mpiManager.gatherV<char>((char*)(resLoc.data()), resLoc.size()*sizeof(Particle), (char*)(res.data()), sizes.data(), displs.data(), 0);

    return res;
}

Int3 MPIAbstractEulerianDomain::getGlobalSize(){
    return size_;
}

Int3 MPIAbstractEulerianDomain::getBlockSize(){
    return blockSize_;
}

Int3 MPIAbstractEulerianDomain::getNumBlockLocal(){
    return numBlockLocal_;
}

MPITopology* MPIAbstractEulerianDomain::getTopology(){
    return topology_;
}

bool MPIAbstractEulerianDomain::isInBounds( Int3 coord ){
    if( coord.x_ < size_.x_ && coord.y_ < size_.y_ && coord.z_ < size_.z_ && coord.x_ >= 0 && coord.y_ >= 0 && coord.z_ >= 0) return true;
    else return false;
}

void MPIAbstractEulerianDomain::addSpeed( std::shared_ptr< SpeedFunctor > newSpeed ){
    for( int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ ; i++ ) domainBlocks_[ i ].addSpeed( newSpeed );
}

void MPIAbstractEulerianDomain::removeSpeed( std::string name ){
    for( int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ ; i++ ) domainBlocks_[ i ].removeSpeed( name );
}

void MPIAbstractEulerianDomain::addParticleFamily( std::shared_ptr< GenericParticleFamily > newFamily ){
    for( int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ ; i++ ) domainBlocks_[ i ].addParticleFamily( newFamily );
}

int MPIAbstractEulerianDomain::getNextFamilyId(){
    return domainBlocks_[ 0 ].getNextFamilyId();
}

std::vector< std::shared_ptr< GenericParticleFamily> > MPIAbstractEulerianDomain::getParticleFamilies(){
    return domainBlocks_[0].getParticleFamilies();
}

void MPIAbstractEulerianDomain::setTerrain( GridTerrain *terrain ){
    for( int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++ ) domainBlocks_[ i ].setTerrain( terrain );
}

bool MPIAbstractEulerianDomain::containsParticles(){
    // int containsParticlesInt;
    //MPI_Status status;
    int i, containsParticlesAll;
    containsParticles_ = 0;
    for( i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++){
        if( domainBlocks_[ i ].containsParticles() ) containsParticles_ = 1;
    }

    MpiManager & mpiManager= this->topology_->getMpiManager();
    std::vector<int> sendRecvVal;
    sendRecvVal.push_back(containsParticles_);
    mpiManager.allReduceVect<int>(sendRecvVal, MPI_MAX);
    containsParticlesAll=sendRecvVal.at(0);
    sendRecvVal.clear();

    if( containsParticlesAll == 1 ) return true;
    return false;
}


// does not work, fix it

void MPIAbstractEulerianDomain::gather(){

    // gather all domain blocks
    std::vector< std::vector< Domain > > gatheredDomainBlocks;
    //MPI_Request request;
    // serialize data to send
    std::ostringstream oss;
    {
        boost::archive::binary_oarchive oa( oss );
        oa << BOOST_SERIALIZATION_NVP( domainBlocks_ );
    }
    int sendSize = oss.str().size();
    // gather send sizes
    // mohamed
    MpiManager & mpiManager = this->topology_->getMpiManager();
    std::vector< int > recvSizes;
    recvSizes.resize( topology_->getSize() );
    mpiManager.gather<int>(&sendSize, 1, &recvSizes[ 0 ], 1, mpiManager.bossId());
    // everybody sends its data to 0
    int tag_send=0; int tag_rec=0;
    int dest=0;

    if(mpiManager.getRank() != 0){
        //std::cout << "p" << mpiManager.getRank() << " send of size " << sendSize << std::endl;
        //mpiManager.send<char>((char*)oss.str().c_str(), sendSize, dest,&request, tag_send);
        mpiManager.send<char>((char*)oss.str().c_str(), sendSize, dest, tag_send);
    }

    // 0 receives all data
    if( topology_->getRank() == 0 ){
        for( int n = 0; n < topology_->getSize(); n++ ){

            if(n==0){
                //std::cout << "p" << mpiManager.getRank() << " getting domain blocks locally" << std::endl;
                gatheredDomainBlocks.push_back( domainBlocks_ );
                //std::cout << "p" << mpiManager.getRank() << " done getting domain blocks locally" << std::endl;
            } else {
                std::vector< Domain > domainBlocks;
                std::vector< char > incomingBuffer( recvSizes[ n ] );
                //std::cout << "p" << mpiManager.getRank() << " will now receive data from p" << n << std::endl;
                mpiManager.receive<char>(&incomingBuffer[ 0 ], recvSizes[ n ], n, tag_rec);
                //std::cout << "p" << mpiManager.getRank() << " finished receive from p" << n << " of size " << recvSizes[ n ] << std::endl;
                // deserialize and put data into gatheredDomainBlock
                try{
                    //std::cout << "p" << mpiManager.getRank() << " tries to deserialize data  ..." << std::endl;
                    std::istringstream iss( std::string( &incomingBuffer[ 0 ], incomingBuffer.size() ) );
                    boost::archive::binary_iarchive ia(iss);
                    ia >> BOOST_SERIALIZATION_NVP( domainBlocks );
                    gatheredDomainBlocks.push_back( domainBlocks );
                    //std::cout << "p" << mpiManager.getRank() << " done deserialization " << std::endl;
                }
                catch(std::exception& e){
                    std::cout << "p" << mpiManager.getRank() << " there was an exception during deserialization" << std::endl;
                }
            }
        }
    }

    // if rank == 0, put all particles in the gatheredDomain_
    if( topology_->getRank() == 0 ){
        gatheredDomain_ = Domain( position_, { size_.x_ * dx_, size_.y_ * dx_, size_.z_ * dx_ } , dx_ );
        for( auto& blocks : gatheredDomainBlocks){
            for( auto& dom : blocks ){
                for( int z = 0; z < dom.getZSize(); z++ ){
                    for( int y = 0; y < dom.getYSize(); y++ ){
                        for( int x = 0; x < dom.getXSize(); x++ ){
                            for( auto& part : dom.domain_[ dom.index3d( { x, y, z } ) ] ){
                                part.displacement_.x_ += dom.position_.x_ + x * dx_;
                                part.displacement_.y_ += dom.position_.y_ + y * dx_;
                                part.displacement_.z_ += dom.position_.z_ + z * dx_;
                                gatheredDomain_.putParticle( part );
                            }
                        }
                    }
                }
            }
        }
        for(auto pf : getParticleFamilies()){
            gatheredDomain_.addParticleFamily(pf);
        }
    }

}


Domain MPIAbstractEulerianDomain::getGatheredDomain(){
    return gatheredDomain_;
}

void MPIAbstractEulerianDomain::swapBuffer(){
    for(int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++) domainBlocks_[ i ].swapBuffer();
}

void MPIAbstractEulerianDomain::computeSpeeds( double t, double dt ){
    for(int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++) domainBlocks_[ i ].computeSpeeds( t, dt );
}

void MPIAbstractEulerianDomain::computeSpeedsOnBorder( double t, double dt ){
    for(int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++) domainBlocks_[ i ].computeSpeedsOnBorder( t, dt );
}

void MPIAbstractEulerianDomain::computeSpeedsExceptBorder( double t, double dt ){
    for(int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++) domainBlocks_[ i ].computeSpeedsExceptBorder( t, dt );
}

void MPIAbstractEulerianDomain::computeStaticSpeeds(){
    for(int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++) domainBlocks_[ i ].computeStaticSpeeds();
}


int MPIAbstractEulerianDomain::countDomainParticles(){
    int cpt = 0;
    int totalCpt;

    MpiManager & mpiManager = this->topology_->getMpiManager();

    // std::cout << "--- MPIAbstractEulerianDomain::countDomainParticles() ---" << std::endl;
    // std::cout << "number of domain blocks : " << domainBlocks_.size() << std::endl;

    int blockCounter = 0;
    int incr;

    for(auto &block: domainBlocks_){
        incr = block.countDomainParticles();
        // if(incr > 0){
        //     std::cout << incr << "particles found in block " << blockCounter << std::endl;
        // }
        cpt += incr;
        blockCounter++;
    }

    mpiManager.reduce<int>(cpt, totalCpt, MPI_SUM, mpiManager.bossId());

    return totalCpt;
}

std::vector< double > MPIAbstractEulerianDomain::getNumberOfparticles(){

    std::vector< double > resTmp(getParticleFamilies().size(), 0.0);

    for(auto &block: domainBlocks_){
        std::transform( resTmp.begin(), resTmp.end(), block.getNumberOfparticles().begin(),
                        resTmp.begin(), std::plus<double>() );
    }


    MpiManager & mpiManager = this->topology_->getMpiManager();

    std::vector< double > res(getParticleFamilies().size(), 0.0);
    mpiManager.reduceVect<double>( resTmp, res, MPI_SUM, mpiManager.bossId() );


    return res;
}

std::shared_ptr< std::vector< Particle > > MPIAbstractEulerianDomain::getParticles(){

    MpiManager & mpiManager = this->topology_->getMpiManager();



    std::vector< Particle > resLoc;

    for( auto &block: domainBlocks_ ){
        auto tmp = block.getParticles();
        //resLoc.insert(resLoc.end(), tmp->begin(), tmp->end());
        for( auto p : *tmp ) resLoc.push_back( p );
    }

    //std::shared_ptr <std::vector< Particle >> res;
    //std::vector< Particle > res;
    
    auto res = new std::vector< Particle >();

    std::vector< int > sizes(mpiManager.getSize(), 0);
    std::vector< int > displacements(mpiManager.getSize(), 0);
    int sizeLoc = resLoc.size()*sizeof(Particle);

    mpiManager.gather<int>(&sizeLoc, 1, sizes.data(), 1, mpiManager.bossId());



    for(uint i=0; i<displacements.size(); i++){
        if(i>0){
            displacements[i]=displacements[i-1]+sizes[i-1];
        }
    }

    int totalSize = 0;
    for(auto v : sizes) totalSize += v;


    //if( mpiManager.isMainProcessor() ) res = std::make_shared< std::vector< Particle > >( std::vector< Particle >( ) );
    //res->resize(totalSize/sizeof(Particle));
    if( mpiManager.isMainProcessor() ){
        //std::cout << "before resize res on p0, will resize for " << totalSize/sizeof(Particle) << " particles, which means " << totalSize << " bytes" << std::endl;
        //res.resize( totalSize/sizeof(Particle) );
        //res->resize( totalSize/sizeof(Particle) );
        delete res;
        res = new std::vector< Particle >( totalSize/sizeof(Particle) );
        //std::cout << "after resize res on p0" << std::endl;
    }

    //std::cout << "res size : " << res.size() << std::endl;

    //std::cout << "p" << mpiManager.getRank() << " -  1 getParticles()" << std::endl;

    mpiManager.barrier();
    //mpiManager.gatherV< char >((char*)resLoc.data(), sizeLoc, (char*)res.data(), sizes.data(), displacements.data(), mpiManager.bossId());

    mpiManager.gatherV< char >((char*)resLoc.data(), sizeLoc, (char*)res->data(), sizes.data(), displacements.data(), mpiManager.bossId());

    //mpiManager.gatherV< char >((char*)resLoc.data(), sizeLoc, (char*)res.get(), sizes.data(), displacements.data(), mpiManager.bossId());

    //std::cout << "p" << mpiManager.getRank() << " -  2 getParticles()" << std::endl;

    //return std::make_shared< std::vector< Particle > >( res );

    std::shared_ptr< std::vector< Particle > > resShared(res);
    return resShared;
    
    //return res;
}

std::vector< std::vector<double> > MPIAbstractEulerianDomain::getDomainParticlesNumber(){
    std::vector< std::vector<double> > res;
    gather();
    if(topology_->getRank() == 0) {
        res = getGatheredDomain().getDomainParticlesNumber();
    }
    return res;
}

void MPIAbstractEulerianDomain::insertParticlesRandomPosition ( int incr, int familyId, double scaling ){

    // first broadcast the seed
    auto seed = time(NULL);
    MpiManager & mpiManager = this->topology_->getMpiManager();
    mpiManager.bCast(&seed, 1, mpiManager.bossId());
    srand(seed);

    // then insert at random positions
    double xMin = position_.x_;
    double xMax = position_.x_ + size_.x_*dx_;
    double xInterval = xMax - xMin;

    // removing 10% on each side... avoid insertion bug... fix this...
    xMin += 0.1*xInterval;
    xMax -= 0.1*xInterval;
    xInterval *= 0.8;

    double yMin = position_.y_;
    double yMax = position_.y_ + size_.y_*dx_;
    double yInterval = yMax - yMin;

    yMin += 0.1*yInterval;
    yMax -= 0.1*yInterval;
    yInterval *= 0.8;

    double zMin = position_.z_;
    double zMax = position_.z_ + size_.z_*dx_;
    double zInterval = zMax - zMin;

    zMin += 0.1*zInterval;
    zMax -= 0.1*zInterval;
    zInterval *= 0.8;

    for(int i=0; i<incr; i++){
        //std::cout << "xMin = " << xMin << ", rand() = " << rand() << ", RAND_MAX = " << RAND_MAX << ", rand()/RAND_MAX = " << double(rand())/double(RAND_MAX) << ", xInterval = " << xInterval << std::endl;
        Double3 coord = { xMin + double(rand())/double(RAND_MAX) *  xInterval, yMin + double(rand())/double(RAND_MAX) *  yInterval, zMin + double(rand())/double(RAND_MAX) *  zInterval };
        //std::cout << "Inserting at : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
        insertParticles( coord, 1, familyId, scaling );
    }
}

void MPIAbstractEulerianDomain::insertParticlesRandomPositionSphere( int incr, int familyId, double scaling, Double3 sphereCenter, double sphereRadius ){

    // first broadcast the seed
    auto seed = time(NULL);
    MpiManager & mpiManager = this->topology_->getMpiManager();
    mpiManager.bCast(&seed, 1, mpiManager.bossId());
    //srand(seed);

    // then insert at random positions

    double xMin = position_.x_;
    double xMax = position_.x_ + size_.x_*dx_;
    double xInterval = xMax - xMin;

    // removing 10% on each side... avoid insertion bug... fix this...
    xMin += 0.1*xInterval;
    xMax -= 0.1*xInterval;
    xInterval *= 0.8;

    double yMin = position_.y_;
    double yMax = position_.y_ + size_.y_*dx_;
    double yInterval = yMax - yMin;

    yMin += 0.1*yInterval;
    yMax -= 0.1*yInterval;
    yInterval *= 0.8;

    double zMin = position_.z_;
    double zMax = position_.z_ + size_.z_*dx_;
    double zInterval = zMax - zMin;

    zMin += 0.1*zInterval;
    zMax -= 0.1*zInterval;
    zInterval *= 0.8;

    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> xDistribution(xMin,xMax);
    std::uniform_real_distribution<double> yDistribution(yMin,yMax);
    std::uniform_real_distribution<double> zDistribution(zMin,zMax);

    for(int i=0; i<incr; i++){
        // pick random position until it lies within the sphere
        //Double3 coord = { xMin + double(rand())/double(RAND_MAX) *  xInterval, yMin + double(rand())/double(RAND_MAX) *  yInterval, zMin + double(rand())/double(RAND_MAX) *  zInterval };
        Double3 coord = { xDistribution(generator), yDistribution(generator), zDistribution(generator) };
        //        std::cout << "coord picked : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
        //        std::cout << "sphere center : " << sphereCenter.x_ << ", " << sphereCenter.y_ << ", " << sphereCenter.z_ << std::endl;
        //        std::cout << "sphere radius : " << sphereRadius << std::endl;
        while( !coordInSphere(coord, sphereCenter, sphereRadius) ){
            //coord = { xMin + double(rand())/double(RAND_MAX) *  xInterval, yMin + double(rand())/double(RAND_MAX) *  yInterval, zMin + double(rand())/double(RAND_MAX) *  zInterval };
            coord = { xDistribution(generator), yDistribution(generator), zDistribution(generator) };
            //            std::cout << "coord picked : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
            //            std::cout << "sphere center : " << sphereCenter.x_ << ", " << sphereCenter.y_ << ", " << sphereCenter.z_ << std::endl;
            //            std::cout << "sphere radius : " << sphereRadius << std::endl;
        }
        //        std::cout << "insetring at coord : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
        insertParticles( coord, 1, familyId, scaling );
    }
}

void MPIAbstractEulerianDomain::insertParticlesRandomPositionSphereGaussian( int incr, int familyId, double scaling, Double3 sphereCenter, double sphereRadius, double stdev ){

    // first broadcast the seed
    auto seed = time(NULL);
    MpiManager & mpiManager = this->topology_->getMpiManager();
    mpiManager.bCast(&seed, 1, mpiManager.bossId());
    //srand(seed);

    // then insert at random positions

    double xMin = position_.x_;
    double xMax = position_.x_ + size_.x_*dx_;
    double xInterval = xMax - xMin;

    // removing 10% on each side... avoid insertion bug... fix this...
    xMin += 0.1*xInterval;
    xMax -= 0.1*xInterval;
    xInterval *= 0.8;

    double yMin = position_.y_;
    double yMax = position_.y_ + size_.y_*dx_;
    double yInterval = yMax - yMin;

    yMin += 0.1*yInterval;
    yMax -= 0.1*yInterval;
    yInterval *= 0.8;

    double zMin = position_.z_;
    double zMax = position_.z_ + size_.z_*dx_;
    double zInterval = zMax - zMin;

    zMin += 0.1*zInterval;
    zMax -= 0.1*zInterval;
    zInterval *= 0.8;

    std::default_random_engine generator;
    std::normal_distribution<double> xDistribution(xMin+xInterval/2.0,stdev);
    std::normal_distribution<double> yDistribution(yMin+yInterval/2.0,stdev);
    std::normal_distribution<double> zDistribution(zMin+zInterval/2.0,stdev);

    for(int i=0; i<incr; i++){
        // pick random position until it lies within the sphere
        //Double3 coord = { xMin + double(rand())/double(RAND_MAX) *  xInterval, yMin + double(rand())/double(RAND_MAX) *  yInterval, zMin + double(rand())/double(RAND_MAX) *  zInterval };
        Double3 coord = { xDistribution(generator), yDistribution(generator), zDistribution(generator) };
        //        std::cout << "coord picked : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
        //        std::cout << "sphere center : " << sphereCenter.x_ << ", " << sphereCenter.y_ << ", " << sphereCenter.z_ << std::endl;
        //        std::cout << "sphere radius : " << sphereRadius << std::endl;
        while( !coordInSphere(coord, sphereCenter, sphereRadius) ){
            //coord = { xMin + double(rand())/double(RAND_MAX) *  xInterval, yMin + double(rand())/double(RAND_MAX) *  yInterval, zMin + double(rand())/double(RAND_MAX) *  zInterval };
            coord = { xDistribution(generator), yDistribution(generator), zDistribution(generator) };
            //            std::cout << "coord picked : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
            //            std::cout << "sphere center : " << sphereCenter.x_ << ", " << sphereCenter.y_ << ", " << sphereCenter.z_ << std::endl;
            //            std::cout << "sphere radius : " << sphereRadius << std::endl;
        }
        //        std::cout << "insetring at coord : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
        insertParticles( coord, 1, familyId, scaling );
    }
}

} // namespace piaf
