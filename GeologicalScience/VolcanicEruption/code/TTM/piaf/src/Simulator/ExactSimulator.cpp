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


#include "../../include/Simulator/ExactSimulator.hpp"

namespace piaf{

ExactSimulator::ExactSimulator(AbstractEulerianDomain *d, double initialTime, double dt, GridTerrain *terrain, ParticleRepository *repository, bool ddt ):
		AbstractSimulator	( d, initialTime, dt, terrain, repository, ddt )
//diffusionSpeed_(d->getDx()/dt)
{}


ExactSimulator::~ExactSimulator() {}

void ExactSimulator::step( ){

	bool iterFirst = true;
	static int iterCpt = 0;

	// because domain_ is an AbstractDomain
	Domain* domain = static_cast< Domain* >( domain_ );


	STARTITER:

	if( ( iterCpt % 100 == 0 ) && iterFirst && ddt_ ) dt_ += ddtValue_;


	std::vector< Particle > terrainBuffer;
	std::vector< Particle > repositoryBuffer;

#ifdef TRACE
	static int iterStart;
	static clock_t clockCpt, clockStart;
	if( iterCpt % 10 == 0 ){
		clockCpt = 0;
		iterStart = iterCpt;
	}
	clockStart = clock();
#endif

	/*! \todo : une grosse partie du code de cette fonction est semblable au code de MPIExactSimulator::step() ... essayer de regrouper */



	// update advecting field (eulerian speeds) and Particle velocities (lagrangian speeds)
	domain->computeSpeeds( currentTime_, dt_ );

	Double3 			speed;
	Double3 			displacement;
	Int3 			gridDisplacement = { 0 , 0 , 0 };
	Int3 			coord;
	//Particle* 		particle;
	const double 	dx = domain->getDx();
	int				index3d;
	Double3			domainPosition = domain->getPosition();
	Particle newPart;

	if( ddt_ ){
		for( int z = 0 ; z < domain->getZSize() ; z++ ){ for( int y = 0 ; y < domain->getYSize() ; y++ ){ for( int x = 0 ; x < domain->getXSize() ; x++ ){
			domain->eraseBufferAt( { x, y, z } );
		} } }
	}

	// for each domain site
	for( int z = 0 ; z < domain->getZSize() ; z++ ){ for( int y = 0 ; y < domain->getYSize() ; y++ ){ for( int x = 0 ; x < domain->getXSize() ; x++ ){

		// if the site contains no particles, there is nothing to do
		/*!
		 * \todo there must be a most optimized way to do that !!!
		 */
		index3d = domain->index3d( { x , y , z } );


		//if( domain->domain_[ index3d ].size() > 0 ){
		if( ! domain->domain_[ index3d ].empty() ){
			if( domain->isInGroundContact_[ index3d ] ){
				//for(std::vector< Particle* >::iterator part = domain->domain_[ index3d ].begin() ; part != domain->domain_[ index3d ].end() ; part++ ){
				//for( Particle& part : domain->domain_[ index3d ] ){
				BOOST_FOREACH( Particle& part , domain->domain_[ index3d ] ){
					newPart = part;
					// if the particle is in a site that is in contact with the ground, deposit the particle
					//if( terrain_ != NULL ){
					//std::cout<<"particle deposit : "<<x<<";"<<y<<";"<<z<<std::endl;
					// the displacement of the particle becomes its absolute position and its speed becomes its total speed
					// mixed eulerian/lagrangian representation -> lagrangian
					newPart.displacement_.x_ += x * dx + domainPosition.x_;
					newPart.displacement_.y_ += y * dx + domainPosition.y_;
					newPart.displacement_.z_ += z * dx + domainPosition.z_;
					// particle speed (lagrangian) = eulerian + lagrangian
					speed = domain->getSpeed( { x , y , z } , newPart.familyId_ );
					//newPart.lagrangianSpeed_.x_ += speed.x_;
					//newPart.lagrangianSpeed_.y_ += speed.y_;
					//newPart.lagrangianSpeed_.z_ += speed.z_;
					newPart.lagrangianSpeed_ = speed;
					// deposit the particle
					//terrain_->depositParticle( part );
					terrainBuffer.push_back( newPart );
					//delete (*part);
					//}
				}
			}
			// else move the particle
			// for each Particle in the site
			else {
				//for( std::vector< Particle* >::iterator part = domain->domain_[ index3d ].begin() ; part != domain->domain_[ index3d ].end() ; part++ ){
				//for( Particle& part : domain->domain_[ index3d ] ){
				BOOST_FOREACH( Particle& part , domain->domain_[ index3d ] ){
					newPart = part;
					//particle = ( *part );
					speed = domain->getSpeed( { x , y , z } , newPart.familyId_ );
					//std::cout<<gSpeed.x_<<";"<<gSpeed.y_<<";"<<gSpeed.z_<<std::endl;

					// compute new gDisplacement of the Particle (taking into account eulerian and lagrangian gSpeed)
					//displacement = speed*dt_ + (*part)->lagrangianSpeed_*dt_ + (*part)->displacement_;
					//displacement = speed*dt_ + part.lagrangianSpeed_*dt_ + part.displacement_;
					displacement.x_ = speed.x_ * dt_ + newPart.lagrangianSpeed_.x_ * dt_ + newPart.displacement_.x_;
					displacement.y_ = speed.y_ * dt_ + newPart.lagrangianSpeed_.y_ * dt_ + newPart.displacement_.y_;
					displacement.z_ = speed.z_ * dt_ + newPart.lagrangianSpeed_.z_ * dt_ + newPart.displacement_.z_;

					// compute grid gDisplacement
					gridDisplacement.x_ = roundhalfup( displacement.x_ / dx );
					gridDisplacement.y_ = roundhalfup( displacement.y_ / dx );
					gridDisplacement.z_ = roundhalfup( displacement.z_ / dx );
					/*gridDisplacement.x_ = displacement.x_ / dx;
					gridDisplacement.y_ = displacement.y_ / dx;
					gridDisplacement.z_ = displacement.z_ / dx;*/

					//std::cout << "ExactSimulator a calcule comme gridDisplacement : " << gridDisplacement.x_ << ", " << gridDisplacement.y_ << ", " << gridDisplacement.z_ << " a t = " << currentTime_<< std::endl;

					//assert( gridDisplacement.x_ <= 1 && gridDisplacement.y_ <= 1 && gridDisplacement.z_ <= 1 && gridDisplacement.x_ >= -1 && gridDisplacement.y_ >= -1 && gridDisplacement.z_ >= -1 );
					if( gridDisplacement.x_ > 1 || gridDisplacement.y_ > 1 || gridDisplacement.z_ > 1 || gridDisplacement.x_ < -1 || gridDisplacement.y_ < -1 || gridDisplacement.z_ < -1 ){
						if( ddt_ ){
							iterFirst = false;
							dt_ -= ddtValue_;
							/*! \todo don't use a goto... */
							goto STARTITER;
						}
						else{
							/*! \todo handle with and exception */
							std::cout << "gridDisplacement > 1" << std::endl;
							exit(0);
						}
					}

					// compute new gDisplacement of the Particle
					displacement.x_ = displacement.x_ - gridDisplacement.x_ * dx;
					displacement.y_ = displacement.y_ - gridDisplacement.y_ * dx;
					displacement.z_ = displacement.z_ - gridDisplacement.z_ * dx;
					//(*part)->displacement_ = displacement;
					//part.displacement_ = displacement;
					newPart.displacement_.x_ = displacement.x_;
					newPart.displacement_.y_ = displacement.y_;
					newPart.displacement_.z_ = displacement.z_;

					// placing the Particle in the next domain
					//coord = {x+gridDisplacement.x_, y+gridDisplacement.y_, z+gridDisplacement.z_};
					coord.x_ = x + gridDisplacement.x_;
					coord.y_ = y + gridDisplacement.y_;
					coord.z_ = z + gridDisplacement.z_;

					if( domain->isInBounds( coord ) ) {

						//std::cout << "exact : une particule est placee dans le buffer a la coordonnee "<< domain->index3d( coord ) << " et elle venait de la coordonee " << domain->index3d( {x,y,z} ) << std::endl;

						domain->putParticleInBuffer( newPart, coord );
						//buffer_->putParticle( particle, coord );
					}

					// if the particle go out of the buffer, place it in the repository
					else{
						if( repository_ != NULL ) {
							repositoryBuffer.push_back( newPart );
							//repository_->depositParticle( part );
						}
						//delete (*part);
					}
				}
			}
			// reinitialize the domain for the next step
			if( ! ddt_ ) domain->eraseAt( { x , y , z } );
		}
	} } }

	BOOST_FOREACH( Particle& part, terrainBuffer ) terrain_->depositParticle( part );
	BOOST_FOREACH( Particle& part, repositoryBuffer ) repository_->depositParticle( part );
	// swap domain and buffer and update time
	domain->swapBuffer();

#ifdef TRACE
	clockCpt += clock() - clockStart;
	if( iterCpt % 10 == 9 ){
		trace("computeClockCount, iterPerSample : 10, iterStart : %1%, iterEnd : %2%, clockCpt : %3%, process : %4%", "0") % iterStart % iterCpt % clockCpt % 0;
	}
#endif

	//std::swap( domain, buffer_ );
	currentTime_ += dt_;

	iterCpt++;


	//terrainBuffer.clear();
	//repositoryBuffer.clear();
}

}//namespace piaf
