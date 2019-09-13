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


#include "../../../include/Simulator/MPI/MPIExactSimulator.hpp"

namespace piaf{

    MPIExactSimulator::MPIExactSimulator( AbstractEulerianDomain *d, double initialTime, double dt, GridTerrain *terrain, MPIParticleRepository *repository, bool ddt, MpiManager & mpiManager):
    MPIAbstractSimulator	( d, initialTime, dt, terrain, repository, ddt, mpiManager  )
	{}

		MPIExactSimulator::~MPIExactSimulator() {}

		/*! \todo : une grosse partie du code de cette fonction est semblable au code de ExactSimulator::step() ... essayer de regrouper */
		/*! \todo : remplacer structure goto par fonction tryMove */
		void MPIExactSimulator::step( ){

			//if( partTracker_ != NULL ) ((MPIParticleTracker*)partTracker_)->addTime( currentTime_ );

			MPIAbstractEulerianDomain* domain = static_cast< MPIAbstractEulerianDomain* >( domain_ );
			bool iterFirst = true;
			int iterOk;
			static int iterCpt = 0;
			#ifdef TRACE
			static int iterStart;
			static clock_t clockCpt, clockStart;
			if( iterCpt % 10 == 0 ){
				clockCpt = 0;
				iterStart = iterCpt;
			}
			clockStart = clock();
			#endif
			STARTITER:

			if(dt_ <= 0) {
				std::cout << "Too many tries adjusting dt, reached value <= 0" << std::endl;
				throw std::string("Too many tries adjusting dt, reached value <= 0");
			}

			iterOk = 1;
			if( ( iterCpt % 100 == 0 ) && iterFirst && ddt_ ) dt_ += ddtValue_;
			std::vector< Particle > terrainBuffer;
			std::vector< Particle > repositoryBuffer;

			// update advecting field (eulerian speeds) and Particle velocities (lagrangian speeds)
			domain->computeSpeeds( currentTime_, dt_ );

			Domain* domainBlock;
			Double3 		speed;
			Double3 		displacement;
			Int3 			gridDisplacement = { 0 , 0 , 0 };
			Int3 			coord;
			const double 	dx = domain->getDx();
			int			index3d;
			Double3			domainBlockPosition;
			Particle newPart;
			// for each domain block
			for( int zBlock = 0; zBlock <  domain->getNumBlockLocal().z_; zBlock++ ) {
				for( int yBlock = 0; yBlock < domain->getNumBlockLocal().y_; yBlock++ ) {
					for( int xBlock = 0; xBlock < domain->getNumBlockLocal().x_; xBlock++ ) {
						domainBlock = domain->getDomainBlockAtLocalIndex3d( { xBlock , yBlock , zBlock } );
						domainBlockPosition = domainBlock->getPosition();
						if( ddt_ ){
							for( int z = 0 ; z < domainBlock->getZSize() ; z++ ){
								for( int y = 0 ; y < domainBlock->getYSize() ; y++ ){
									for( int x = 0 ; x < domainBlock->getXSize() ; x++ ){
										domainBlock->eraseBufferAt( { x, y, z } );
									}
								}
							}
						}
						// for each domain site of the domain block
						for( int z = 0 ; z < domainBlock->getZSize() ; z++ ){
							for( int y = 0 ; y < domainBlock->getYSize() ; y++ ){
								for( int x = 0 ; x < domainBlock->getXSize() ; x++ ){
									index3d = domainBlock->index3d( { x, y, z } );
									if( ! domainBlock->domain_[ index3d ].empty() ){
										if( domainBlock->isInGroundContact_[ index3d ] ){
											BOOST_FOREACH( Particle& part , domainBlock->domain_[ index3d ] ){
												newPart = part;
												// if the particle is in a site that is in contact with the ground, deposit the particle
												// the displacement of the particle becomes its absolute position and its speed becomes its total speed
												// mixed eulerian/lagrangian representation -> lagrangian
												newPart.displacement_.x_ += x * dx + domainBlockPosition.x_;
												newPart.displacement_.y_ += y * dx + domainBlockPosition.y_;
												newPart.displacement_.z_ += z * dx + domainBlockPosition.z_;
												// eulerian + lagrangian speed (problem with diffusion ?)
												speed = domainBlock->getSpeed( { x , y , z } , newPart.familyId_ );
												newPart.lagrangianSpeed_ = newPart.lagrangianSpeed_ + speed;
												//newPart.lagrangianSpeed_ = speed;
												// deposit the particle
												terrainBuffer.push_back( newPart );
												if( partTracker_ != NULL ) {
													//if( part.familyId_ == 8 ) std::cout << part.displacement_.x_+x * dx + domainBlockPosition.x_ << "," << part.displacement_.y_+y * dx + domainBlockPosition.y_ << "," << part.displacement_.z_+z * dx + domainBlockPosition.z_ << std::endl;
													//partTracker_->addParticle( newPart, currentTime_ );
													trackBuffer_.push_back( newPart );
												}
											}
										}
										else{
											BOOST_FOREACH( Particle& part , domainBlock->domain_[ index3d ] ){


												if( partTracker_ != NULL ) {
													//if( part.familyId_ == 8 ) std::cout << part.displacement_.x_+x * dx + domainBlockPosition.x_ << "," << part.displacement_.y_+y * dx + domainBlockPosition.y_ << "," << part.displacement_.z_+z * dx + domainBlockPosition.z_ << std::endl;
													Particle newPart = part;
													newPart.displacement_.x_ += x * dx + domainBlockPosition.x_;
													newPart.displacement_.y_ += y * dx + domainBlockPosition.y_;
													newPart.displacement_.z_ += z * dx + domainBlockPosition.z_;
                          							speed = domainBlock->getSpeed( { x , y , z } , newPart.familyId_ );
  													newPart.lagrangianSpeed_ = newPart.lagrangianSpeed_ + speed;
													//partTracker_->addParticle( newPart, currentTime_ );
													trackBuffer_.push_back( newPart );
												}

												newPart = part;
												speed = domainBlock->getSpeed( { x , y , z } , newPart.familyId_ );
												// compute new gDisplacement of the Particle (taking into account eulerian and lagrangian gSpeed)
												displacement.x_ = speed.x_ * dt_ + newPart.lagrangianSpeed_.x_ * dt_ + newPart.displacement_.x_;
												displacement.y_ = speed.y_ * dt_ + newPart.lagrangianSpeed_.y_ * dt_ + newPart.displacement_.y_;
												displacement.z_ = speed.z_ * dt_ + newPart.lagrangianSpeed_.z_ * dt_ + newPart.displacement_.z_;
												// compute grid gDisplacement
												gridDisplacement.x_ = roundhalfup( displacement.x_ / dx );
												gridDisplacement.y_ = roundhalfup( displacement.y_ / dx );
												gridDisplacement.z_ = roundhalfup( displacement.z_ / dx );
												if( gridDisplacement.x_ > 1 || gridDisplacement.y_ > 1 || gridDisplacement.z_ > 1 || gridDisplacement.x_ < -1 || gridDisplacement.y_ < -1 || gridDisplacement.z_ < -1 ){
													if( ddt_ ){
														/*! \todo don't use a goto, wrap into a function instead */
														iterOk = 0;
														goto ENDITER;
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
												newPart.displacement_.x_ = displacement.x_;
												newPart.displacement_.y_ = displacement.y_;
												newPart.displacement_.z_ = displacement.z_;
												// placing the Particle in the next domain
												coord.x_ = x + gridDisplacement.x_;
												coord.y_ = y + gridDisplacement.y_;
												coord.z_ = z + gridDisplacement.z_;
												if( domainBlock->isInBounds( coord ) ) {
													domainBlock->putParticleInBuffer( newPart, coord );
												}
												// if the particle leaves the domain block ...
												else{
													// if the particle stays in the global domain, put it in the corresponding packet to send
													if( domain->isInBounds( { coord.x_ + int( round( ( domainBlock->position_.x_ - domain->position_.x_ ) / dx ) ), coord.y_ + int( round( ( domainBlock->position_.y_ - domain->position_.y_ ) / dx ) ), coord.z_ + int( round( ( domainBlock->position_.z_ - domain->position_.z_ ) / dx ) )  } ) ){
														// the displacement of the particle becomes its absolute position
														// but we don't care about its speed because the particle will be handled in the domain by the next process
														// mixed eulerian/lagrangian representation -> lagrangian
														newPart.displacement_.x_ += coord.x_ * dx + domainBlockPosition.x_;
														newPart.displacement_.y_ += coord.y_ * dx + domainBlockPosition.y_;
														newPart.displacement_.z_ += coord.z_ * dx + domainBlockPosition.z_;
														domain->putParticleInPacket( newPart, domainBlock );
													}
													// if the particle leaves the global domain
													else{
														if( repository_ != NULL ) {
															// the displacement of the particle becomes its absolute position
															// mixed eulerian/lagrangian representation -> lagrangian
															newPart.displacement_.x_ += coord.x_ * dx + domainBlockPosition.x_;
															newPart.displacement_.y_ += coord.y_ * dx + domainBlockPosition.y_;
															newPart.displacement_.z_ += coord.z_ * dx + domainBlockPosition.z_;
															// eulerian + lagrangian speed (problem with diffusion ?)
															speed = domainBlock->getSpeed( { x , y , z } , newPart.familyId_ );
															newPart.lagrangianSpeed_ = newPart.lagrangianSpeed_ + speed;
															//newPart.lagrangianSpeed_ = speed;
															repositoryBuffer.push_back( newPart );
														}
													}
												}
											}
										}
									}
									// reinitialize the domain for the next step
									if( ! ddt_ ) domainBlock->eraseAt( { x , y , z } );

								}
							}
						}
					}
				}
			}
			ENDITER:
			if( ddt_ ){
				int iterOkAll;
				#ifdef TRACE
				clockCpt += clock() - clockStart;
				#endif
                std::vector<int> sendRecvVal;
                sendRecvVal.push_back(iterOk);
                mpiManager.allReduceVect<int>(sendRecvVal, MPI_MIN);
                iterOkAll=sendRecvVal.at(0);
                sendRecvVal.clear();

				#ifdef TRACE
				clockStart = clock();
				#endif
				if( iterOkAll == 0 ){
					// clear particle packets and buffers
					trackBuffer_.clear();
					terrainBuffer.clear();
					repositoryBuffer.clear();
					BOOST_FOREACH( std::vector< Particle >& packet, domain->particlePackets_) packet.clear();
					dt_ -= ddtValue_;
					iterFirst = false;
					goto STARTITER;
				}
			}

			// update time
			currentTime_ += dt_;
			iterCpt++;

			if( partTracker_ != NULL ){
                if( (currentTime_ - partTracker_->getLastTrackTime()) >= partTracker_->getInterval() || relativeDif(currentTime_ - partTracker_->getLastTrackTime(), partTracker_->getInterval()) < 1e-6 ){
					((MPIParticleTracker*)partTracker_)->addTime( currentTime_ );
					for( auto& part : trackBuffer_ ){
						partTracker_->addParticle( part, currentTime_ );
					}
				}
			}
			trackBuffer_.clear();

			BOOST_FOREACH( Particle& part, terrainBuffer ){
				//if(part.familyId_ == 1){
				//	std::cout << "deposit particle : " << part.displacement_.x_ << "," << part.displacement_.y_ << "," << part.displacement_.z_ << " ; " << part.lagrangianSpeed_.x_ << "," << part.lagrangianSpeed_.y_ << "," << part.lagrangianSpeed_.z_ << std::endl;
				//}
				terrain_->depositParticle( part );
			}
			BOOST_FOREACH( Particle& part, repositoryBuffer ){
				repository_->depositParticle( part );
				repository_->depositParticle( BoundaryParticle( part, currentTime_ ) );
			}
			// swap buffers
			domain->swapBuffer();
			#ifdef TRACE
			clockCpt += clock() - clockStart;
			std::ostringstream oss;
			oss << domain->getTopology()->getRank();
			if( iterCpt % 10 == 9 ){
				trace("computeClockCount, iterPerSample : 10, iterStart : %1%, iterEnd : %2%, clockCpt : %3%, process : %4%", oss.str().c_str() ) % iterStart % iterCpt % clockCpt % domain->getTopology()->getRank();
			}
			#endif

			// exchange particles
			domain->exchangeParticles();
		}

	} // namespace piaf
