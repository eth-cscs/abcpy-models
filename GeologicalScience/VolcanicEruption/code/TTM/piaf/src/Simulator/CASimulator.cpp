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


#include "../../include/Simulator/CASimulator.hpp"

namespace piaf{

/*
 * Implementation of constructor and destructor
 */

CASimulator::CASimulator( AbstractEulerianDomain *d, double initialTime, double dt, GridTerrain *terrain, ParticleRepository *repository, bool ddt ):
		AbstractSimulator	( d, initialTime, dt, terrain, repository, ddt )
{}

CASimulator::~CASimulator() {}

double CASimulator::computeProba( double speed, double dx ){
	assert( abs( speed * dt_ / dx ) <= 1.0 );
	return speed * dt_ / dx;
}

/*
 * Implemenation of the method step()
 */

void CASimulator::step( ){

	Domain* domain = static_cast< Domain* >( domain_ );

	// update advecting field
	domain->computeSpeeds( currentTime_, dt_ );

	int neighborX;
	int neighborY;
	int neighborZ;

	bool moveX;
	bool moveMX;
	bool moveY;
	bool moveMY;
	bool moveZ;
	bool moveMZ;

	Double3 proba;
	Double3 speed;

	Double3 domainPosition = domain->getPosition();

	// for each domain site
	for( int z = 0 ; z < domain->getZSize() ; z++ ){
		for( int y = 0 ; y < domain->getYSize() ; y++ ){
			for( int x = 0 ; x < domain->getXSize() ; x++ ){

				if( ! domain->domain_[ domain->index3d( { x , y , z } ) ].empty() ){

					if( domain->isInGroundContact_ [ domain->index3d( { x , y , z } ) ] ) {
						/*
						 * \todo factorize CASimulator and ExactSimulator... deposit process can be a (inline) function
						 */
						//for( std::vector<Particle*>::iterator part = domain->domain_[ domain->index3d( { x , y , z } ) ].begin() ; part != domain->domain_[ domain->index3d( { x , y , z } ) ].end() ; part++ ){
						//for( Particle& part :  domain->domain_[ domain->index3d( { x , y , z } ) ] ){
						BOOST_FOREACH( Particle& part, domain->domain_[ domain->index3d( { x , y , z } ) ] ){
							// if the particle is in a site that is in contact with the ground, deposit the particle
							/*! \todo : give the right speed to the particle (EULERIAN + LAGRANGIAN->LAGRANGIAN) */
							//if( terrain_ != NULL ){
							//std::cout<<"particle deposit : "<<x<<";"<<y<<";"<<z<<std::endl;
							// the displacement of the particle becomes its absolute position and its speed becomes its total speed
							// mixed eulerian/lagrangian representation->lagrangian
							part.displacement_.x_ += x * domain->dx_ + domainPosition.x_;
							part.displacement_.y_ += y * domain->dx_ + domainPosition.y_;
							part.displacement_.z_ += z * domain->dx_ + domainPosition.z_;
							// particle speed (lagrangian) = eulerian + lagrangian
							speed = domain->getSpeed( { x , y , z } , part.familyId_ );
							//part.lagrangianSpeed_.x_ += speed.x_;
							//part.lagrangianSpeed_.y_ += speed.y_;
							//part.lagrangianSpeed_.z_ += speed.z_;
							part.lagrangianSpeed_ = speed;
							// deposit the particle
							terrain_->depositParticle( part );
							//delete (*part);
							//}
						}
					}
					else{
						// for each Particle in the site
						//for( std::vector<Particle*>::iterator part = domain->domain_[ domain->index3d( { x , y , z } ) ].begin() ; part != domain->domain_[ domain->index3d( { x , y , z } ) ].end() ; part++ ){
						//for( Particle& part : domain->domain_[ domain->index3d( { x , y , z } ) ] ){
						BOOST_FOREACH( Particle& part, domain->domain_[ domain->index3d( { x , y , z } ) ] ){
							// else, move the particle
							/*gMoveX = (double(rand())/double(RAND_MAX))<domain->ProbaX({x,y,z},part->familyId_,dt_);
					gMoveMX = (double(rand())/double(RAND_MAX))<domain->ProbaMX({x,y,z},part->familyId_,dt_);
					gMoveY = (double(rand())/double(RAND_MAX))<domain->ProbaY({x,y,z},part->familyId_,dt_);
					gMoveMY = (double(rand())/double(RAND_MAX))<domain->ProbaMY({x,y,z},part->familyId_,dt_);
					gMoveZ = (double(rand())/double(RAND_MAX))<domain->ProbaZ({x,y,z},part->familyId_,dt_);
					gMoveMZ = (double(rand())/double(RAND_MAX))<domain->ProbaMZ({x,y,z},part->familyId_,dt_);*/

							// We use lagrangian speed of the particle to compute probability !
							proba.x_ = computeProba( domain->getSpeed( { x , y , z } , part.familyId_ ).x_ + part.lagrangianSpeed_.x_ , domain->getDx() );
							proba.y_ = computeProba( domain->getSpeed( { x , y , z } , part.familyId_ ).y_ + part.lagrangianSpeed_.y_ , domain->getDx() );
							proba.z_ = computeProba( domain->getSpeed( { x , y , z } , part.familyId_ ).z_ + part.lagrangianSpeed_.z_ , domain->getDx() );


							moveX = moveMX = moveY = moveMY = moveZ = moveMZ = false;

							if( proba.x_ < 0.0 && ( double( rand() ) / double( RAND_MAX ) ) < std::abs( proba.x_ ) ) moveMX = true;
							else if( ( double( rand() ) / double( RAND_MAX ) ) < proba.x_ ) moveX = true;

							if( proba.y_ < 0.0 && ( double( rand() ) / double( RAND_MAX ) ) < std::abs( proba.y_ ) ) moveMY = true;
							else if( ( double( rand() ) / double( RAND_MAX ) ) < proba.y_ ) moveY = true;

							if( proba.z_ < 0.0 && ( double( rand() ) / double( RAND_MAX ) ) < std::abs( proba.z_ ) ) moveMZ = true;
							else if( ( double( rand() ) / double( RAND_MAX ) ) < proba.z_ ) moveZ = true;


							// choose the neighbor to move the Particle, according to the chosen movement
							neighborX = x;
							neighborY = y;
							neighborZ = z;
							if( moveX ) 	neighborX++;
							if( moveMX ) 	neighborX--;
							if( moveY ) 	neighborY++;
							if( moveMY ) 	neighborY--;
							if( moveZ ) 	neighborZ++;
							if( moveMZ ) 	neighborZ--;
							// move the Particle (in the buffer)
							if( domain->isInBounds( { neighborX , neighborY , neighborZ } ) ) {
								//domain->IncrBufferAt({gNeighborX, gNeighborY,gNeighborZ}, 1,part->familyId_);
								domain->putParticleInBuffer( part , { neighborX , neighborY , neighborZ } );
								//buffer_->putParticle( *part , { neighborX , neighborY , neighborZ } );
							}
							// if the particle go out of the buffer, place it in the repository
							else {
								repository_->depositParticle( part );
								//delete (*part);
							}

							//}
						}
					}
					// reinitialize the domain for the next step
					//domain->setAt({i,j,0},0);
					domain->eraseAt( { x , y , z } );
				}
			}
		}
	}
	// swap domain and buffer and update time
	domain->swapBuffer();
	//std::swap( domain, buffer_ );
	currentTime_ += dt_;
}

} // namespace piaf
