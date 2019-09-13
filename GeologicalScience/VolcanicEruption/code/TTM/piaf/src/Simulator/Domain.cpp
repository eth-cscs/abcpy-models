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


#include "../../include/Simulator/Domain.hpp"

namespace piaf {

/*
    * Implementation of constructors and destructor
    */

Domain::Domain( Double3 pos, Double3 size, double dx ) :
    AbstractEulerianDomain  ( pos, size, dx ),
    domain_    ( std::vector< std::vector< Particle > >( size_.x_ * size_.y_ * size_.z_ ) ),
    buffer_    ( std::vector< std::vector< Particle > >( size_.x_ * size_.y_ * size_.z_ ) ),
    eulerianStaticSpeeds_ ( std::vector< std::vector< Double3 > >( size_.x_ * size_.y_ * size_.z_ ) ),
    eulerianDynamicSpeeds_ ( std::vector< std::vector< Double3 > >( size_.x_ * size_.y_ * size_.z_ ) ),
    isInGroundContact_ ( std::vector< bool >( size_.x_ * size_.y_ * size_.z_ ) )
{
    maxNParticle_ = 0;
    for( int i = 0; i < size_.x_ * size_.y_ * size_.z_; i++ ) isInGroundContact_[ i ] = false;
}

Domain::Domain() :
    AbstractEulerianDomain(),
    maxNParticle_  ( 0 )

{
}

Domain::~Domain(){
}

void Domain::test(){
    /*bool partFound = false;
        for( std::vector< Particle >& site : domain_ ){
        if( site.size() > 0 ) partFound = true;
    }
    if( partFound ) std::cout << "il y a des particules" << std::endl;
    else std::cout << "pas trouve de particule" << std::endl;*/

    for( Uint i = 0; i < domain_.size(); i++ ) {
        if( domain_[ i ].size() > 0 ) std::cout << "particule trouvee au site n° " << i << std::endl;
    }
}

std::vector< std::vector< Particle > > Domain::getDomainParticles(){
    // return a copy of the current domain
    return domain_;
}

void Domain::computeMaxNumPart(){
    maxNParticle_ = 0;
    for( int i = 0; i < size_.x_ * size_.y_ * size_.z_; i++ ) if( domain_[ i ].size() > Uint( maxNParticle_ ) ) maxNParticle_ = domain_[ i ].size();
}

void Domain::init( Double3 pos, Double3 size, double dx ){
    dx_ = dx;
    size_ = { int( size.x_ / dx_ ), int( size.y_ / dx_ ), int( size.z_ / dx_ ) };
    position_ = pos;
    domain_ = std::vector< std::vector< Particle > >( size_.x_ * size_.y_ * size_.z_ );
    buffer_ = std::vector< std::vector< Particle > >( size_.x_ * size_.y_ * size_.z_ );
    eulerianStaticSpeeds_ = std::vector< std::vector< Double3 > >( size_.x_ * size_.y_ * size_.z_ );
    eulerianDynamicSpeeds_ = std::vector< std::vector< Double3 > >( size_.x_ * size_.y_ * size_.z_ );
    isInGroundContact_ = std::vector< bool >( size_.x_ * size_.y_ * size_.z_ );
    maxNParticle_ = 0;
    for( int i = 0; i < size_.x_ * size_.y_ * size_.z_; i++ ) isInGroundContact_[ i ] = false;
}

void Domain::addSpeed( std::shared_ptr< SpeedFunctor > newSpeed ){
    if( newSpeed->getSpeedType() == LAGRANGIAN_DYNAMIC ) lagrangianDynamicSpeedFunctors_.push_back( newSpeed );
    if( newSpeed->getSpeedType() == EULERIAN_DYNAMIC ) eulerianDynamicSpeedFunctors_.push_back( newSpeed );
    if( newSpeed->getSpeedType() == EULERIAN_STATIC ) {
        AbstractEulerianDomain::addSpeed( newSpeed );
        eulerianStaticSpeedFunctors_.push_back( newSpeed );
    }
}

void Domain::removeSpeed( std::string name ){
    AbstractEulerianDomain::removeSpeed( name );
    for(uint i=0; i<lagrangianDynamicSpeedFunctors_.size(); i++){
        if( lagrangianDynamicSpeedFunctors_[i]->getName().compare(name) == 0 )
            lagrangianDynamicSpeedFunctors_.erase(lagrangianDynamicSpeedFunctors_.begin()+i);
    }
    for(uint i=0; i<eulerianDynamicSpeedFunctors_.size(); i++){
        if( eulerianDynamicSpeedFunctors_[i]->getName().compare(name) == 0 )
            eulerianDynamicSpeedFunctors_.erase(eulerianDynamicSpeedFunctors_.begin()+i);
    }
    for(uint i=0; i<eulerianStaticSpeedFunctors_.size(); i++){
        if( eulerianStaticSpeedFunctors_[i]->getName().compare(name) == 0 )
            eulerianStaticSpeedFunctors_.erase(eulerianStaticSpeedFunctors_.begin()+i);
    }
}

void Domain::addParticleFamily( std::shared_ptr< GenericParticleFamily > newFamily ){
    AbstractEulerianDomain::addParticleFamily( newFamily );
    particleFamilies_.push_back( newFamily );
    // add new speeds for this family
    for( int z = 0; z < size_.z_; z++ ) {
        for( int y = 0; y < size_.y_; y++ ) {
            for( int x = 0; x < size_.x_; x++ ) {
                /*! \todo : for gcc 4.4 compatibility */
                const Double3 elem = { 0.0, 0.0, 0.0 };
                eulerianStaticSpeeds_[ index3d( { x, y, z } ) ].push_back( elem );
                eulerianDynamicSpeeds_[ index3d( { x, y, z } ) ].push_back( elem );
            }
        }
    }
    //this->computeStaticSpeeds();
    if( particleFamilies_.size() - 1 != Uint( newFamily->familyId_ ) ) {
        std::cout << "new family id does not match with index" << std::endl;
        exit( 0 );
    }
}

int Domain::getNextFamilyId(){
    return particleFamilies_.size();
}

int Domain::getNParticles(){
    int cpt;
    cpt = 0;
    for( int z = 0; z < getZSize(); z++ ) {
        for( int y = 0; y < getYSize(); y++ ) {
            for( int x = 0; x < getXSize(); x++ ) {
                cpt += domain_[ index3d( { x, y, z } ) ].size();
            }
        }
    }
    return cpt;
}

bool Domain::containsParticles( ){
    if( containsParticles_ == 1 ) {
        for( int z = 0; z < getZSize(); z++ ) {
            for( int y = 0; y < getYSize(); y++ ) {
                for( int x = 0; x < getXSize(); x++ ) {
                    if( domain_[ index3d( { x, y, z } ) ].size() > 0) return true;
                }
            }
        }
    }
    else return false;
    containsParticles_ = 0;
    return false;
}

bool Domain::isInGroundContactAt( Int3 coord ){
    return isInGroundContact_[ index3d( coord ) ];
}

int Domain::getNParticleAt( Int3 coord ){
    if( isInBounds( coord ) ) {
        return domain_[ index3d( coord ) ].size();
    }
    else {
        std::cout << "array index out of bound in Domain::getNParticleAt( Int3 coord )" << std::endl;
        exit( 0 );
    }
}

int Domain::getNParticleAt( Int3 coord, int familyId ){
    if( isInBounds( coord ) ) {
        int cpt = 0;
        for( Uint i = 0; i < domain_[ index3d( coord ) ].size(); i++ ) {
            if( domain_[ index3d( coord ) ].at( i ).familyId_ == familyId ) cpt++;
        }
        return cpt;
    }
    else {
        std::cout << "array index out of bound in Domain::getNParticleAt( Int3 coord, int familyId )" << std::endl;
        exit( 0 );
    }
}

bool Domain::isInBounds( Int3 coord ){
    if( coord.x_ < size_.x_ && coord.y_ < size_.y_ && coord.z_ < size_.z_ && coord.x_ >= 0 && coord.y_ >= 0 && coord.z_ >= 0){
        //std::cout << "coordinate : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << " is in bound" << std::endl;
        //std::cout << "size of the domain : " << size_.x_ << ", " << size_.y_ << ", " << size_.z_ << std::endl;
        return true;
    }
    else{
        // std::cout << "coordinate : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << " is out of bound" << std::endl;
        // std::cout << "size of the domain : " << size_.x_ << ", " << size_.y_ << ", " << size_.z_ << std::endl;
        return false;
    }
}

int Domain::index3d( Int3 coord ){
    return coord.x_ + coord.y_ * size_.x_ + coord.z_ * size_.x_ * size_.y_;
}

std::vector< std::shared_ptr< GenericParticleFamily> > Domain::getParticleFamilies(){
    return particleFamilies_;
}

void Domain::insertParticlesRandomPosition ( int incr, int familyId, double scaling ){

    double xMin = position_.x_;
    double xMax = position_.x_ + size_.x_*dx_;
    double xInterval = xMax - xMin;

    double yMin = position_.y_;
    double yMax = position_.y_ + size_.y_*dx_;
    double yInterval = yMax - yMin;

    double zMin = position_.z_;
    double zMax = position_.z_ + size_.z_*dx_;
    double zInterval = zMax - zMin;

    for(int i=0; i<incr; i++){
        Double3 coord = { xMin + rand()/RAND_MAX *  xInterval, yMin + rand()/RAND_MAX *  yInterval, zMin + rand()/RAND_MAX *  zInterval };
        insertParticles( coord, 1, familyId, scaling );
    }
}

void Domain::insertParticlesRandomPositionSphere( int incr, int familyId, double scaling, Double3 sphereCenter, double sphereRadius ){
    double xMin = position_.x_;
    double xMax = position_.x_ * size_.x_*dx_;
    double xInterval = xMax - xMin;

    double yMin = position_.y_;
    double yMax = position_.y_ * size_.y_*dx_;
    double yInterval = yMax - yMin;

    double zMin = position_.z_;
    double zMax = position_.z_ * size_.z_*dx_;
    double zInterval = zMax - zMin;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> xDistribution(xMin,xMax);
    std::uniform_real_distribution<double> yDistribution(yMin,yMax);
    std::uniform_real_distribution<double> zDistribution(zMin,zMax);

    for(int i=0; i<incr; i++){
        // pick random position until it lies within the sphere
        //Double3 coord = { xMin + rand()/RAND_MAX *  xInterval, yMin + rand()/RAND_MAX *  yInterval, zMin + rand()/RAND_MAX *  zInterval };
        Double3 coord = { xDistribution(generator), yDistribution(generator), zDistribution(generator) };
        while( !coordInSphere(coord, sphereCenter, sphereRadius) ){
            coord = { xMin + rand()/RAND_MAX *  xInterval, yMin + rand()/RAND_MAX *  yInterval, zMin + rand()/RAND_MAX *  zInterval };
        }
        insertParticles( coord, 1, familyId, scaling );
    }
}

void Domain::insertParticlesRandomPositionSphereGaussian( int incr, int familyId, double scaling, Double3 sphereCenter, double sphereRadius, double stdev ){
    double xMin = position_.x_;
    double xMax = position_.x_ * size_.x_*dx_;
    double xInterval = xMax - xMin;

    double yMin = position_.y_;
    double yMax = position_.y_ * size_.y_*dx_;
    double yInterval = yMax - yMin;

    double zMin = position_.z_;
    double zMax = position_.z_ * size_.z_*dx_;
    double zInterval = zMax - zMin;

    std::default_random_engine generator;
    std::normal_distribution<double> xDistribution(xMin+xInterval/2.0,stdev);
    std::normal_distribution<double> yDistribution(yMin+yInterval/2.0,stdev);
    std::normal_distribution<double> zDistribution(zMin+zInterval/2.0,stdev);

    for(int i=0; i<incr; i++){
        // pick random position until it lies within the sphere
        //Double3 coord = { xMin + rand()/RAND_MAX *  xInterval, yMin + rand()/RAND_MAX *  yInterval, zMin + rand()/RAND_MAX *  zInterval };
        Double3 coord = { xDistribution(generator), yDistribution(generator), zDistribution(generator) };
        while( !coordInSphere(coord, sphereCenter, sphereRadius) ){
            coord = { xMin + rand()/RAND_MAX *  xInterval, yMin + rand()/RAND_MAX *  yInterval, zMin + rand()/RAND_MAX *  zInterval };
        }
        insertParticles( coord, 1, familyId, scaling );
    }
}

// result for computing coordu is rounded in order to keep the particle linked to its closest possible site ( like in the simulator functions )
// void Domain::insertParticles( Double3 coord, int incr, int familyId ){
//     insertParticles( coord, incr, familyId, 1.0 );
// }

void Domain::insertParticles( Double3 coord, int incr, int familyId, double scaling ){
    //std::cout << "insertParticles class Domain" << std::endl;

    containsParticles_ = 1;

    //Int3 coordu = { int( coord.x_ / dx_ - position_.x_ ) , int( coord.y_ / dx_ - position_.y_ ) , int( coord.z_ / dx_  - position_.z_ ) };
    //Int3 coordu = { int( ( coord.x_ - position_.x_ ) / dx_ ) , int( ( coord.y_ - position_.y_ ) / dx_  ) , int( ( coord.z_ - position_.z_ ) / dx_ ) };
    Int3 coordu = { int( roundhalfup( ( coord.x_ - position_.x_ ) / dx_ ) ), int( roundhalfup( ( coord.y_ - position_.y_ ) / dx_  ) ), int( roundhalfup( ( coord.z_ - position_.z_ ) / dx_ ) ) };
    //Int3 coordu = { int( ( ( coord.x_ - position_.x_ ) / dx_ ) ) , int( ( ( coord.y_ - position_.y_ ) / dx_  ) ) , int( ( ( coord.z_ - position_.z_ ) / dx_ ) ) };

    if( isInBounds( coordu ) ) {
        // if incr>0, add incr new particles to domain[index3d(coord)]
        if( incr > 0 ) {
            //std::cout << "inserting " << incr << " particles at position " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
            for( int i = 0; i < incr; i++ ) {
                //domain_[index3d(coord)].push_back({familyId,{0.0,0.0,0.0}});
                //domain_[ index3d( coordu ) ].push_back( Particle( familyId , { 0.0 , 0.0 , 0.0 } , { 0.0 , 0.0 , 0.0 } ) );
                domain_[ index3d( coordu ) ].push_back( Particle( familyId, { coord.x_ - ( coordu.x_*dx_ + position_.x_ ), coord.y_ - ( coordu.y_*dx_ + position_.y_ ), coord.z_ - ( coordu.z_*dx_ + position_.z_ ) }, { 0.0, 0.0, 0.0 }, scaling ) );
            }
        }
        // if incr<0, remove incr Particle from domain[index3d(coord)] with familyId_ == familyId
        // if there is not enough Particle with familyId_ == familyId, remove all these particles
        if( incr < 0 ) {
            for( int i = 0; i > incr; i-- ) {
                //std::cout << "removing " << incr << " particles at position " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
                //for(std::vector<Particle>::iterator j=domain_[index3d(coord)].begin();j!=domain_[index3d(coord)].end();++j){
                for( std::vector< Particle >::iterator j = domain_[ index3d( coordu ) ].begin(); j != domain_[ index3d( coordu ) ].end(); ++j ) {
                    if( j->familyId_ == familyId ) {
                        domain_[ index3d( coordu ) ].erase( j );
                        //delete ( *j );
                    }
                    if( j == domain_[ index3d( coordu ) ].end() ) i = incr-1;
                }
            }
        }
    }
    else {
        std::cout << "array index out of bound in Domain::insertParticles( Double3 coord, int incr, int familyId )" << std::endl;
        std::cout << "real coordinate : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
        std::cout << "real position of the domain : " << position_.x_ << ", " << position_.y_ << ", " << position_.z_ << std::endl;
        std::cout << "real size of the domain : " << dx_*size_.x_ << ", " << dx_*size_.y_ << ", " << dx_*size_.z_ << std::endl;
        std::cout << "coordinate : " << coordu.x_ << ", " << coordu.y_ << ", " << coordu.z_ << std::endl;
        std::cout << "size of the domain : " << size_.x_ << ", " << size_.y_ << ", " << size_.z_ << std::endl;
        exit(0);
    }
}

int Domain::getMaxNParticles(){
    return maxNParticle_;
}

Double3 Domain::getSpeed( Int3 coord, int familyId ){
    return eulerianStaticSpeeds_[ index3d( coord ) ].at( familyId ) + eulerianDynamicSpeeds_[ index3d( coord ) ].at( familyId );
}


void Domain::insertParticlesInBuffer( Int3 coord, int incr, int familyId ){
    containsParticles_ = 1;
    if( isInBounds( coord ) ) {
        // if incr>0, add incr new particles to buffer[index3d(coord)]
        if( incr > 0 ) {
            for( int i = 0; i < incr; i++ ) {
                //buffer_[index3d(coord)].push_back({familyId,{0.0,0.0,0.0}});
                buffer_[ index3d( coord ) ].push_back( Particle( familyId, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } ) );
            }
        }
        // if incr<0, remove incr Particle from buffer[index3d(coord)] with familyId_ == familyId
        // if there is not enough Particle with familyId_ == familyId, remove all these particles
        if( incr < 0 ) {
            for( int i = 0; i > incr; i-- ) {
                for( std::vector< Particle >::iterator j = buffer_[ index3d( coord ) ].begin(); j != buffer_[ index3d( coord ) ].end(); ++j ) {
                    if( j->familyId_ == familyId ) {
                        buffer_[ index3d( coord ) ].erase( j );
                        //delete ( *j );
                    }
                    if( j == buffer_[ index3d( coord ) ].end() ) i = incr - 1;
                }
            }
        }
    }
    else {
        std::cout << "array index out of bound in Domain::insertParticlesInBuffer( Int3 coord, int incr, int familyId )" << std::endl;
        exit( 0 );
    }
}

void Domain::setAt( Int3 coord, int val, int familyId ){
    containsParticles_ = 1;
    if( isInBounds( coord ) ) {
        int numPartThisFamily = 0;
        for( std::vector< Particle >::iterator j = domain_[ index3d( coord ) ].begin(); j != domain_[ index3d( coord ) ].end(); ++j ) {
            if( j->familyId_ == familyId ) numPartThisFamily++;
        }
        // if there are too many particles, delete some
        if( numPartThisFamily > val ) {
            for( int i = numPartThisFamily; i > val; i-- ) {
                for( std::vector< Particle >::iterator j = domain_[ index3d( coord ) ].begin(); j != domain_[ index3d( coord ) ].end(); ++j ) {
                    if( j->familyId_ == familyId ) {
                        domain_[ index3d( coord ) ].erase( j );
                        break;
                        //delete ( *j );
                    }
                }
            }
        }
        // if there are too few particles, add some
        if( numPartThisFamily < val ) {
            for( int i = numPartThisFamily; i < val; i++ ) {
                //domain_[index3d(coord)].push_back({familyId,{0.0,0.0,0.0}});
                domain_[ index3d( coord ) ].push_back( Particle( familyId, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } ) );
            }
        }
    }
    else {
        std::cout << "array index out of bound in Domain::setAt( Int3 coord, int val, int familyId )" << std::endl;
        exit( 0 );
    }
}

void Domain::eraseAt( Int3 coord ){
    if( isInBounds( coord ) ) {

        // clear and force memory reallocation
        std::vector< Particle >().swap( domain_[ index3d( coord ) ] );

        //domain_[ index3d( coord ) ].clear();
    }
    else {
        std::cout << "array index out of bound in Domain::eraseAt( Int3 coord )" << std::endl;
        exit( 0 );
    }
}

void Domain::eraseBufferAt( Int3 coord ){
    if( isInBounds( coord ) ) {
        //domain_[index3d(coord)].clear();

        // clear and force memory reallocation
        std::vector< Particle >().swap( buffer_[ index3d( coord ) ] );

        //buffer_[ index3d( coord ) ].clear();
    }
    else {
        std::cout << "array index out of bound in Domain::eraseBufferAt( Int3 coord )" << std::endl;
        exit( 0 );
    }
}


void Domain::swapBuffer(){
    swap( domain_, buffer_ );
}

void Domain::computeStaticSpeeds(){
    double index3d;
    Double3 speed;
    // for each site
    for( int z = 0; z < size_.z_; z++ ) {
        for( int y = 0; y < size_.y_; y++ ) {
            for( int x = 0; x < size_.x_; x++ ) {
                index3d = this->index3d( { x, y, z } );
                // for each Particle family
                for( std::vector< std::shared_ptr< GenericParticleFamily > >::iterator fam = particleFamilies_.begin(); fam != particleFamilies_.end(); ++fam ) {
                    eulerianStaticSpeeds_[ index3d ].at( ( *fam )->familyId_ ) = { 0.0, 0.0, 0.0 };
                    for( auto& speedFunct : eulerianStaticSpeedFunctors_ ) {
                        //speed = speedFunct->compute( { x * dx_ + position_.x_ , y * dx_ + position_.y_ , z * dx_ + position_.z_ } ,0.0 , 0.0 , *fam );
                        speed = ( *speedFunct )( { x * dx_ + position_.x_, y * dx_ + position_.y_, z * dx_ + position_.z_ },0.0, 0.0, *fam );
                        eulerianStaticSpeeds_[ index3d ].at( ( *fam )->familyId_ ).x_ += speed.x_;
                        eulerianStaticSpeeds_[ index3d ].at( ( *fam )->familyId_ ).y_ += speed.y_;
                        eulerianStaticSpeeds_[ index3d ].at( ( *fam )->familyId_ ).z_ += speed.z_;
                    }
                }
            }
        }
    }
}

void Domain::computeSpeedsOnBorder( double t, double dt ){
    double index3d;
    Double3 speed;
    // for each site on the border
    int incrZ = ( size_.z_ - 1 > 0 ) ? size_.z_ -1 : 1;
    int incrY = ( size_.y_ - 1 > 0 ) ? size_.y_ -1 : 1;
    int incrX = ( size_.x_ - 1 > 0 ) ? size_.x_ -1 : 1;

    for( int x = 0; x < size_.x_; x += incrX ) {
        for( int y = 0; y < size_.y_; y++ ) {
            for( int z = 0; z < size_.z_; z++ ) {
                index3d = this->index3d( { x, y, z } );
                if( domain_[ index3d ].size() > 0 ) {
                    for( auto& part : domain_[ index3d ] ) {
                        part.lagrangianSpeed_ = { 0.0, 0.0, 0.0 };
                    }
                    for( std::vector< std::shared_ptr< GenericParticleFamily > >::iterator fam = particleFamilies_.begin(); fam != particleFamilies_.end(); ++fam ) {
                        eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ) = { 0.0, 0.0, 0.0 }; //eulerianStaticSpeeds_[ index3d ].at( ( *fam )->familyId_ );
                        for( Uint s = 0; s < eulerianDynamicSpeedFunctors_.size(); s++ ) {
                            //speed = eulerianDynamicSpeedFunctors_.at( s )->compute( { x * dx_ + position_.x_ , y * dx_ + position_.y_ , z * dx_ + position_.z_ } , t, dt , *fam );
                            speed = ( *eulerianDynamicSpeedFunctors_.at ( s ) )( { x * dx_ + position_.x_, y * dx_ + position_.y_, z * dx_ + position_.z_ }, t, dt, *fam );
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).x_ += speed.x_;
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).y_ += speed.y_;
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).z_ += speed.z_;
                        }
                    }
                    for( Uint s = 0; s < lagrangianDynamicSpeedFunctors_.size(); s++ ) {
                        for( auto& part : domain_[ index3d ] ) {
                            //speed = lagrangianDynamicSpeedFunctors_.at( s )->compute( { x * dx_ + position_.x_ + part.displacement_.x_ , y * dx_ + position_.y_ + part.displacement_.y_, z * dx_ + position_.z_ + part.displacement_.z_ } , t, dt , particleFamilies_[ part.familyId_ ] );
                            speed = ( *lagrangianDynamicSpeedFunctors_.at ( s ) )( { x * dx_ + position_.x_ + part.displacement_.x_, y * dx_ + position_.y_ + part.displacement_.y_, z * dx_ + position_.z_ + part.displacement_.z_ }, t, dt, particleFamilies_[ part.familyId_ ] );
                            part.lagrangianSpeed_.x_ += speed.x_;
                            part.lagrangianSpeed_.y_ += speed.y_;
                            part.lagrangianSpeed_.z_ += speed.z_;
                        }
                    }
                }
            }
        }
    }

    for( int y = 0; y < size_.y_; y += incrY ) {
        for( int x = 1; x < size_.x_ - 1; x++ ) {
            for( int z = 0; z < size_.z_; z++ ) {
                index3d = this->index3d( { x, y, z } );
                if( domain_[ index3d ].size() > 0 ) {
                    for( auto& part : domain_[ index3d ] ) {
                        part.lagrangianSpeed_ = { 0.0, 0.0, 0.0 };
                    }
                    for( std::vector< std::shared_ptr< GenericParticleFamily > >::iterator fam = particleFamilies_.begin(); fam != particleFamilies_.end(); ++fam ) {
                        eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ) = { 0.0, 0.0, 0.0 }; //eulerianStaticSpeeds_[ index3d ].at( ( *fam )->familyId_ );
                        for( Uint s = 0; s < eulerianDynamicSpeedFunctors_.size(); s++ ) {
                            //speed = eulerianDynamicSpeedFunctors_.at( s )->compute( { x * dx_ + position_.x_ , y * dx_ + position_.y_ , z * dx_ + position_.z_ } , t, dt , *fam );
                            speed = ( *eulerianDynamicSpeedFunctors_.at ( s ) )( { x * dx_ + position_.x_, y * dx_ + position_.y_, z * dx_ + position_.z_ }, t, dt, *fam );
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).x_ += speed.x_;
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).y_ += speed.y_;
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).z_ += speed.z_;
                        }
                    }
                    for( Uint s = 0; s < lagrangianDynamicSpeedFunctors_.size(); s++ ) {
                        for( auto& part : domain_[ index3d ] ) {
                            //speed = lagrangianDynamicSpeedFunctors_.at( s )->compute( { x * dx_ + position_.x_ + part.displacement_.x_ , y * dx_ + position_.y_ + part.displacement_.y_, z * dx_ + position_.z_ + part.displacement_.z_ } , t, dt , particleFamilies_[ part.familyId_ ] );
                            speed = ( *lagrangianDynamicSpeedFunctors_.at ( s ) )( { x * dx_ + position_.x_ + part.displacement_.x_, y * dx_ + position_.y_ + part.displacement_.y_, z * dx_ + position_.z_ + part.displacement_.z_ }, t, dt, particleFamilies_[ part.familyId_ ] );
                            part.lagrangianSpeed_.x_ += speed.x_;
                            part.lagrangianSpeed_.y_ += speed.y_;
                            part.lagrangianSpeed_.z_ += speed.z_;
                        }
                    }
                }
            }
        }
    }

    for( int z = 0; z < size_.z_; z += incrZ ) {
        for( int x = 1; x < size_.x_ - 1; x++ ) {
            for( int y = 1; y < size_.y_ - 1; y++ ) {
                index3d = this->index3d( { x, y, z } );
                if( domain_[ index3d ].size() > 0 ) {
                    for( auto& part : domain_[ index3d ] ) {
                        part.lagrangianSpeed_ = { 0.0, 0.0, 0.0 };
                    }
                    for( std::vector< std::shared_ptr< GenericParticleFamily > >::iterator fam = particleFamilies_.begin(); fam != particleFamilies_.end(); ++fam ) {
                        eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ) = { 0.0, 0.0, 0.0 }; //eulerianStaticSpeeds_[ index3d ].at( ( *fam )->familyId_ );
                        for( Uint s = 0; s < eulerianDynamicSpeedFunctors_.size(); s++ ) {
                            //speed = eulerianDynamicSpeedFunctors_.at( s )->compute( { x * dx_ + position_.x_ , y * dx_ + position_.y_ , z * dx_ + position_.z_ } , t, dt , *fam );
                            speed = ( *eulerianDynamicSpeedFunctors_.at ( s ) )( { x * dx_ + position_.x_, y * dx_ + position_.y_, z * dx_ + position_.z_ }, t, dt, *fam );
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).x_ += speed.x_;
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).y_ += speed.y_;
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).z_ += speed.z_;
                        }
                    }
                    for( Uint s = 0; s < lagrangianDynamicSpeedFunctors_.size(); s++ ) {
                        for( auto& part : domain_[ index3d ] ) {
                            //speed = lagrangianDynamicSpeedFunctors_.at( s )->compute( { x * dx_ + position_.x_ + part.displacement_.x_ , y * dx_ + position_.y_ + part.displacement_.y_, z * dx_ + position_.z_ + part.displacement_.z_ } , t, dt , particleFamilies_[ part.familyId_ ] );
                            speed = ( *lagrangianDynamicSpeedFunctors_.at ( s ) )( { x * dx_ + position_.x_ + part.displacement_.x_, y * dx_ + position_.y_ + part.displacement_.y_, z * dx_ + position_.z_ + part.displacement_.z_ }, t, dt, particleFamilies_[ part.familyId_ ] );
                            part.lagrangianSpeed_.x_ += speed.x_;
                            part.lagrangianSpeed_.y_ += speed.y_;
                            part.lagrangianSpeed_.z_ += speed.z_;
                        }
                    }
                }
            }
        }
    }

}

void Domain::computeSpeedsExceptBorder( double t, double dt ){

    AbstractEulerianDomain::computeSpeeds( t, dt );

    double index3d;
    Double3 speed;
    for( int z = 1; z < size_.z_ - 1; z++ ) {
        for( int y = 1; y < size_.y_ - 1; y++ ) {
            for( int x = 1; x < size_.x_ - 1; x++ ) {
                index3d = this->index3d( { x, y, z } );
                if( domain_[ index3d ].size() > 0 ) {
                    for( auto& part : domain_[ index3d ] ) {
                        part.lagrangianSpeed_ = { 0.0, 0.0, 0.0 };
                    }
                    for( std::vector< std::shared_ptr< GenericParticleFamily > >::iterator fam = particleFamilies_.begin(); fam != particleFamilies_.end(); ++fam ) {
                        eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ) = { 0.0, 0.0, 0.0 }; //eulerianStaticSpeeds_[ index3d ].at( ( *fam )->familyId_ );
                        for( Uint s = 0; s < eulerianDynamicSpeedFunctors_.size(); s++ ) {
                            //speed = eulerianDynamicSpeedFunctors_.at( s )->compute( { x * dx_ + position_.x_ , y * dx_ + position_.y_ , z * dx_ + position_.z_ } , t, dt , *fam );
                            speed = ( *eulerianDynamicSpeedFunctors_.at ( s ) )( { x * dx_ + position_.x_, y * dx_ + position_.y_, z * dx_ + position_.z_ }, t, dt, *fam );
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).x_ += speed.x_;
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).y_ += speed.y_;
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).z_ += speed.z_;
                        }
                    }
                    for( Uint s = 0; s < lagrangianDynamicSpeedFunctors_.size(); s++ ) {
                        for( auto& part :
                             domain_[ index3d ] ) {
                            //speed = lagrangianDynamicSpeedFunctors_.at( s )->compute( { x * dx_ + position_.x_ + part.displacement_.x_ , y * dx_ + position_.y_ + part.displacement_.y_, z * dx_ + position_.z_ + part.displacement_.z_ } , t, dt , particleFamilies_[ part.familyId_ ] );
                            speed = ( *lagrangianDynamicSpeedFunctors_.at ( s ) )( { x * dx_ + position_.x_ + part.displacement_.x_, y * dx_ + position_.y_ + part.displacement_.y_, z * dx_ + position_.z_ + part.displacement_.z_ }, t, dt, particleFamilies_[ part.familyId_ ] );
                            part.lagrangianSpeed_.x_ += speed.x_;
                            part.lagrangianSpeed_.y_ += speed.y_;
                            part.lagrangianSpeed_.z_ += speed.z_;
                        }
                    }
                }
            }
        }
    }
}

void Domain::computeSpeeds( double t, double dt ){



    AbstractEulerianDomain::computeSpeeds( t, dt );

    //boost::posix_time::ptime    timeStart;
    //boost::posix_time::ptime    timeEnd;
    //boost::posix_time::time_duration  duration;

    // set the initial speed to static speed
    double index3d;
    Double3 speed; // = {0.0,0.0,0.0};
    // for each site
    for( int z = 0; z < size_.z_; z++ ) {
        for( int y = 0; y < size_.y_; y++ ) {
            for( int x = 0; x < size_.x_; x++ ) {
                index3d = this->index3d( { x, y, z } );
                // if there is no particle in this site, no need to compute

                if( domain_[ index3d ].size() > 0 ) {
                    // set the speed to 0 for each Particle in the site
                    for( auto& part : domain_[ index3d ] ) {
                        part.lagrangianSpeed_ = { 0.0, 0.0, 0.0 };
                        part.diffusionSpeed_ = { 0.0, 0.0, 0.0 };
                    }

                    // for each Particle family
                    for( std::vector< std::shared_ptr< GenericParticleFamily > >::iterator fam = particleFamilies_.begin(); fam != particleFamilies_.end(); ++fam ) {
                        // set the initial speed to 0
                        eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ) = { 0.0, 0.0, 0.0 }; //eulerianStaticSpeeds_[ index3d ].at( ( *fam )->familyId_ );
                        // add all speeds to this site
                        for( Uint s = 0; s < eulerianDynamicSpeedFunctors_.size(); s++ ) {
                            // if the speed is eulerian, apply it to the sites speed
                            speed = ( *eulerianDynamicSpeedFunctors_.at ( s ) )( { x * dx_ + position_.x_, y * dx_ + position_.y_, z * dx_ + position_.z_ }, t, dt, *fam );
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).x_ += speed.x_;
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).y_ += speed.y_;
                            eulerianDynamicSpeeds_[ index3d ].at( ( *fam )->familyId_ ).z_ += speed.z_;
                        }
                    }
                    for( Uint s = 0; s < lagrangianDynamicSpeedFunctors_.size(); s++ ) {
                        for( auto& part : domain_[ index3d ] ) {
                            std::shared_ptr< SpeedFunctor > speedfunctor = lagrangianDynamicSpeedFunctors_.at( s );
                            speed = ( *speedfunctor )( { x * dx_ + position_.x_ + part.displacement_.x_, y * dx_ + position_.y_ + part.displacement_.y_, z * dx_ + position_.z_ + part.displacement_.z_ }, t, dt, particleFamilies_[ part.familyId_ ] );

                            if( (*lagrangianDynamicSpeedFunctors_.at( s )).isDiffusion() ) {
                                part.diffusionSpeed_.x_ += speed.x_;
                                part.diffusionSpeed_.y_ += speed.y_;
                                part.diffusionSpeed_.z_ += speed.z_;
                            }

                            part.lagrangianSpeed_.x_ += speed.x_;
                            part.lagrangianSpeed_.y_ += speed.y_;
                            part.lagrangianSpeed_.z_ += speed.z_;

                        }
                    }
                }
            }
        }
    }



}

void Domain::putParticleInBuffer( Particle part, Int3 coord ){
    containsParticles_ = 1;
    if( isInBounds( coord ) ) {
        buffer_[ index3d( coord ) ].push_back( part );
    }
    else {
        std::cout << "array index out of bound in Domain::putParticleInBuffer( Particle part, Int3 coord )" << std::endl;
        exit( 0 );
    }
}

void Domain::putParticleInBuffer( Particle part ){
    //Int3 coord = { int( floor( ( part.displacement_.x_ - position_.x_ ) / dx_ ) ), int( floor( ( part.displacement_.y_ - position_.y_ ) / dx_ ) ), int( floor( ( part.displacement_.z_ - position_.z_ ) / dx_ ) ) };
    Int3 coord = { int( round( ( part.displacement_.x_ - position_.x_ ) / dx_ ) ), int( round( ( part.displacement_.y_ - position_.y_ ) / dx_ ) ), int( round( ( part.displacement_.z_ - position_.z_ ) / dx_ ) ) };
    //Int3 coord = { int( ( ( part.displacement_.x_ - position_.x_ ) / dx_ ) ), int( ( ( part.displacement_.y_ - position_.y_ ) / dx_ ) ), int( ( ( part.displacement_.z_ - position_.z_ ) / dx_ ) ) };
    containsParticles_ = 1;
    part.displacement_.x_ -= coord.x_ * dx_ + position_.x_;
    part.displacement_.y_ -= coord.y_ * dx_ + position_.y_;
    part.displacement_.z_ -= coord.z_ * dx_ + position_.z_;
    if( isInBounds( coord ) ) {
        buffer_[ index3d( coord ) ].push_back( part );
    }
    else{
        std::cout << "array index out of bound in Domain::putParticleInBuffer" << std::endl;
        exit( 0 );
    }
}

void Domain::putParticle( Particle part, Int3 coord ){
    containsParticles_ = 1;
    if( isInBounds( coord ) ) {
        domain_[ index3d( coord ) ].push_back( part );
    }
    else {
        std::cout << "array index out of boundin Domain::putParticle( Particle part, Int3 coord )" << std::endl;
        exit( 0 );
    }
}

void Domain::putParticle( Particle part ){

    /*std::cout << "placement d'une particule, position : " << part.displacement_.x_ << ", " << part.displacement_.y_ << ", " << part.displacement_.z_ << std::endl;
        std::cout << "position du domaine : " << position_.x_ << ", " << position_.y_ << ", " << position_.z_ << std::endl;
        std::cout << "taille du domaine : " << size_.x_ * dx_ << ", " << size_.y_ * dx_ << ", " << size_.z_ * dx_ << std::endl;
        std::cout << "taille du domaine (discret) : " << size_.x_ << ", " << size_.y_ << ", " << size_.z_ << std::endl;*/

    //Int3 coord = { int( floor( ( part.displacement_.x_ - position_.x_ ) / dx_ ) ), int( floor( ( part.displacement_.y_ - position_.y_ ) / dx_ ) ), int( floor( ( part.displacement_.z_ - position_.z_ ) / dx_ ) ) };
    //Int3 coord = { int( ( part.displacement_.x_ - position_.x_ ) / dx_ ), int( ( part.displacement_.y_ - position_.y_ ) / dx_ ), int( ( part.displacement_.z_ - position_.z_ ) / dx_ ) };
    //Int3 coord = { int( round( ( part.displacement_.x_ - position_.x_ ) / dx_ ) ), int( round( ( part.displacement_.y_ - position_.y_ ) / dx_ ) ), int( round( ( part.displacement_.z_ - position_.z_ ) / dx_ ) ) };
    //Int3 coord = { int( ( ( part.displacement_.x_ - position_.x_ ) / dx_ ) ), int( ( ( part.displacement_.y_ - position_.y_ ) / dx_ ) ), int( ( ( part.displacement_.z_ - position_.z_ ) / dx_ ) ) };
    //Int3 coord = { int( roundhalftowardzero( ( part.displacement_.x_ - position_.x_ ) / dx_ ) ), int( roundhalftowardzero( ( part.displacement_.y_ - position_.y_ ) / dx_ ) ), int( roundhalftowardzero( ( part.displacement_.z_ - position_.z_ ) / dx_ ) ) };
    Int3 coord = { int( roundhalfup( ( part.displacement_.x_ - position_.x_ ) / dx_ ) ), int( roundhalfup( ( part.displacement_.y_ - position_.y_ ) / dx_ ) ), int( roundhalfup( ( part.displacement_.z_ - position_.z_ ) / dx_ ) ) };
    containsParticles_ = 1;
    part.displacement_.x_ -= coord.x_ * dx_ + position_.x_;
    part.displacement_.y_ -= coord.y_ * dx_ + position_.y_;
    part.displacement_.z_ -= coord.z_ * dx_ + position_.z_;

    //std::cout << "coord ou va etre place la particule : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
    //std::cout << "nouveau deplacement de la particule : " << part.displacement_.x_ << ", " << part.displacement_.y_ << ", " << part.displacement_.z_ << std::endl;

    if( isInBounds( coord ) ) {
        domain_[ index3d( coord ) ].push_back( part );
    }
    else{
        std::cout << "array index out of bound in Domain::putParticle( Particle part )" << std::endl;
        exit( 0 );
    }
}

void Domain::setTerrain( GridTerrain *terrain ){
    for( int z = 0; z < size_.z_; z++ ) {
        for( int y = 0; y < size_.y_; y++ ) {
            for( int x = 0; x < size_.x_; x++ ) {
                // check if one of the 4 lowest points of the domain is in the ground
                if( terrain->isInGround( { x * dx_ + position_.x_, y * dx_ + position_.y_, z * dx_ + position_.z_ } ) ||
                        terrain->isInGround( { ( x + 1 ) * dx_ + position_.x_, y * dx_ + position_.y_, z * dx_ + position_.z_ } ) ||
                        terrain->isInGround( { x * dx_ + position_.x_, ( y + 1 ) * dx_ + position_.y_, z * dx_ + position_.z_ } ) ||
                        terrain->isInGround( { ( x + 1 ) * dx_ + position_.x_, ( y + 1 ) * dx_ + position_.y_, z * dx_ + position_.z_ } ) ) {
                    isInGroundContact_[ index3d( { x, y, z } ) ] = true;
                }
                else isInGroundContact_[ index3d( { x, y, z } ) ] = false;
            }
        }
    }
}

int Domain::countDomainParticles(){

    //std::cout << "--- Domain::countDomainParticles() ---" << std::endl;

    int cpt = 0;
    int incr;
    int siteCount = 0;
    for(auto &site : domain_) {
        incr = site.size();
        cpt += incr;
        // if(incr > 0){
        //     std::cout << "found " << incr << " particles in site " << siteCount << std::endl;
        // }
        siteCount++;
    }

    //std::cout << "--- --- ---" << std::endl;

    return cpt;
}


std::vector< piaf::Particle > Domain::getParticlesInBox( piaf::Double3 p1, piaf::Double3 p2 ){
    std::vector< piaf::Particle > res;
    double startX = std::min(p1.x_, p2.x_); double startY = std::min(p1.y_, p2.y_); double startZ = std::min(p1.z_, p2.z_);
    double endX = std::max(p1.x_, p2.x_); double endY = std::max(p1.y_, p2.y_); double endZ = std::max(p1.z_, p2.z_);

    for( double z = startZ; z < endZ; z += dx_ ) {
        for( double y = startY; y < endY; y += dx_ ) {
            for( double x = startX; x < endX; x += dx_ ) {
                Int3 coord = { int( roundhalfup( ( x - position_.x_ ) / dx_ ) ),
                               int( roundhalfup( ( y - position_.y_ ) / dx_ ) ),
                               int( roundhalfup( ( z - position_.z_ ) / dx_ ) ) };
                if( isInBounds( coord ) ) {
                    // res.insert( res.end(), domain_[ index3d( coord ) ].begin(), domain_[ index3d( coord ) ].end() );
                    for( auto part : domain_[ index3d( coord ) ] ) {
                        part.displacement_.x_ += position_.x_ + coord.x_ * dx_;
                        part.displacement_.y_ += position_.y_ + coord.y_ * dx_;
                        part.displacement_.z_ += position_.z_ + coord.z_ * dx_;
                        if( isInside( part.displacement_ , p1, p2 ) ) res.push_back( part );
                    }
                }
            }
        }
    }

    return res;
}

std::vector < std::string > Domain::getSpeedsName(){
    std::vector < std::string > res;

    for( auto s : eulerianDynamicSpeedFunctors_ ) res.push_back( s->getName() );
    for( auto s : eulerianStaticSpeedFunctors_ ) res.push_back( s->getName() );
    for( auto s : lagrangianDynamicSpeedFunctors_ ) res.push_back( s->getName() );

    return res;
}

std::vector< double > Domain::getNumberOfparticles(){
    std::vector< double > res(particleFamilies_.size(), 0.0);

    for(auto& v : domain_){
        for(auto& p : v){
            res[p.familyId_] += p.scaling_;
        }
    }

    return res;
}

std::shared_ptr<std::vector< Particle >> Domain::getParticles(){
    std::shared_ptr <std::vector< Particle >> res = std::make_shared< std::vector< Particle > >(std::vector< Particle >());
    for(auto& v : domain_){
        for(auto& p : v){
            (*res).push_back(p);
        }
    }
    return res;
}

std::vector< std::vector< double > > Domain::getDomainParticlesNumber(){
    std::vector< std::vector<double> > res = std::vector< std::vector<double> >(particleFamilies_.size());
    for( auto& d : res ){
        d = std::vector< double >(domain_.size());
    }
    int i = 0;
    int lastParticleFamily = 0;
    for( int x = 0; x < size_.x_; x++ ) {
        for( int y = 0; y < size_.y_; y++ ) {
            for( int z = 0; z < size_.z_; z++ ) {
                for( auto& part : domain_[index3d({x, y, z})] ) {
                     res[part.familyId_][i] += part.scaling_;
                     lastParticleFamily = std::max(part.familyId_, lastParticleFamily);
                }
                i++;
            }
        }
    }
    //std::cout << "lastParticleFamily : " << lastParticleFamily << std::endl;
    res.erase(res.begin()+lastParticleFamily+1, res.end());
    return res;
}

}    //namespace piaf
