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

#include "../../include/Simulator/GridTerrain.hpp"

namespace piaf{

	GridTerrain::GridTerrain( Double2 pos, Double2 size, double dx, std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies ):
	ParticleRepository						( particleFamilies ),
	dx_										( dx ),
	discreteSize_							( { int( size.x_ / dx_ ) , int( size.y_ / dx_ ) } ),
	position_								( pos ),
	terrain_( new double[ discreteSize_.x_ * discreteSize_.y_ ] )
	{
		for( int y = 0 ; y < discreteSize_.y_ ; y++ ){
			for( int x = 0 ; x < discreteSize_.x_ ; x++ ){
				terrain_[ index2d( { x , y } ) ] = 0.1 * dx_ ;
			}
		}
	}

	GridTerrain::~GridTerrain() {
		delete[] terrain_;
	}

	bool GridTerrain::isInBounds( Int2 coord ){
		if( coord.x_ < discreteSize_.x_ && coord.y_ < discreteSize_.y_ ) return true;
		else return false;
	}

	bool GridTerrain::isInBounds( Double2 coord ){
		if( coord.x_ - position_.x_ < ( discreteSize_.x_ * dx_ ) && coord.y_ - position_.y_ < ( discreteSize_.y_ * dx_ ) && coord.x_ - position_.x_ >= 0.0 && coord.y_ - position_.y_ >= 0.0 ) return true;
		else return false;
	}

	int GridTerrain::index2d( Int2 coord ){
		return coord.x_ + coord.y_ * discreteSize_.x_;
	}

	void GridTerrain::depositParticle( Particle part ){
		if( isActive_ ){

			part.lagrangianSpeed_ = part.lagrangianSpeed_ - part.diffusionSpeed_;

			bool particleDeposited;
			bool particleOutOfTerrain;

			particleDeposited = false;
			particleOutOfTerrain = false;

			// compute site terrain border
			double yDestUp;
			double yDestInf;
			double xDestRight;
			double xDestLeft;

			double yMoveUp;
			double yMoveDown;
			double xMoveRight;
			double xMoveLeft;

			double xFact;
			double yFact;
			double zFact;

			Double3 nextPos;

			Int2 currentPointCoord;

			particleOutOfTerrain = !isInBounds( Double2( { part.displacement_.x_ , part.displacement_.y_ } ) );

			while( !particleDeposited && !particleOutOfTerrain ){

				//std::cout<<"loop"<<std::endl;

				currentPointCoord = { int( floor( ( part.displacement_.x_ - position_.x_ ) / dx_ ) ) , int( floor( ( part.displacement_.y_ - position_.y_ ) / dx_ ) ) };

				//assert( index2d( currentPointCoord ) >= 0 && index2d( currentPointCoord ) < discreteSize_.x_ * discreteSize_.y_ );

				xDestLeft = floor( ( part.displacement_.x_ ) / ( dx_ / 2.0 ) ) * ( dx_ / 2.0 );
				if( floor( xDestLeft / dx_ ) == ( xDestLeft / dx_ ) ) xDestLeft -= dx_ / 2.0;
				xDestRight = xDestLeft + dx_;

				yDestInf = floor( ( part.displacement_.y_ ) / ( dx_ / 2.0 ) ) * ( dx_ / 2.0 );
				if( floor( yDestInf / dx_ ) == ( yDestInf / dx_ ) ) yDestInf -= dx_ / 2.0;
				yDestUp = yDestInf + dx_;

				// if the particle is exactly in a border, force the destination to the next border
				if( yDestUp 	== ( part.displacement_.y_ ) ) yDestUp += dx_;
				if( yDestInf 	== ( part.displacement_.y_ ) ) yDestInf -= dx_;
				if( xDestRight 	== ( part.displacement_.x_ ) ) xDestRight += dx_;
				if( xDestLeft 	== ( part.displacement_.x_ ) ) xDestLeft -= dx_;

				// if the particle is immobile ( along x and y axis )
				if( part.lagrangianSpeed_.x_ == 0.0 && part.lagrangianSpeed_.y_ == 0.0 ){
					part.displacement_.z_ = terrain_[ index2d( currentPointCoord ) ];
					particleDeposited = true;
					// std::cout<<"particle deposited : no x and y speed"<<std::endl;
				}
				//  if the particle is already under the terrain ( the particle colide a "vertical side" )
				else if( part.displacement_.z_ <= terrain_[ index2d( currentPointCoord ) ] ){
					particleDeposited = true;
					// std::cout<<"particle deposited : vertical side collision"<<std::endl;
				}
				//  else move the particle to the next "terrain border" ( the particle is currently above the ground )
				else{
					//  if the particle has only a velocity along x axis
					if( part.lagrangianSpeed_.x_ != 0.0 && part.lagrangianSpeed_.y_ == 0.0 ){
						xMoveRight = xDestRight - part.displacement_.x_;
						xMoveLeft =  xDestLeft - part.displacement_.x_;
						if( xMoveRight / part.lagrangianSpeed_.x_ > 0.0 ) xFact = xMoveRight / part.lagrangianSpeed_.x_;
						else if( xMoveLeft / part.lagrangianSpeed_.x_ > 0.0 ) xFact = xMoveLeft / part.lagrangianSpeed_.x_;
						else{
							std::cout<<"there is a problem !!! case 1"<<std::endl;
							std::cout<<"dx : "<<dx_<<std::endl;
							std::cout<<"particle x position : "<<part.displacement_.x_<<std::endl;
							std::cout<<"particle speed x : "<<part.lagrangianSpeed_.x_<<std::endl;
							std::cout<<"dest right : "<<xDestRight<<std::endl;
							std::cout<<"dest left : "<<xDestLeft<<std::endl;
							std::cout<<"move right : "<<xMoveRight<<std::endl;
							std::cout<<"move left : "<<xMoveLeft<<std::endl;
							exit( 0 );
						}
						//  to avoid selection of gYFact in the test abs( gXFact ) < abs( gYFact )
						yFact = xFact + 1.0;
					}
					//  if the particle has only a velocity along y axis
					else if( part.lagrangianSpeed_.x_ == 0.0 && part.lagrangianSpeed_.y_ != 0.0 ){
						yMoveUp =  yDestUp - part.displacement_.y_;
						yMoveDown = yDestInf - part.displacement_.y_ ;
						if( yMoveUp / part.lagrangianSpeed_.y_ > 0.0 ) yFact = yMoveUp / part.lagrangianSpeed_.y_;
						else if( yMoveDown / part.lagrangianSpeed_.y_ > 0.0 ) yFact = yMoveDown / part.lagrangianSpeed_.y_;
						else{std::cout<<"there is a problem !!! case 2"<<std::endl;exit( 0 );}
						//  to avoid selection of gXFact in the test abs( gXFact ) < abs( gYFact )
						xFact = yFact + 1.0;
					}
					//  if the particle has a velocity along x and y axis
					else{
						xMoveRight = xDestRight - part.displacement_.x_;
						xMoveLeft =  xDestLeft - part.displacement_.x_;
						if( xMoveRight / part.lagrangianSpeed_.x_ > 0.0 ) xFact = xMoveRight / part.lagrangianSpeed_.x_;
						else if( xMoveLeft / part.lagrangianSpeed_.x_ > 0.0 ) xFact = xMoveLeft / part.lagrangianSpeed_.x_;
						else{std::cout<<"there is a problem !!! case 3"<<std::endl;exit( 0 );}

						yMoveUp =  yDestUp - part.displacement_.y_;
						yMoveDown =  yDestInf - part.displacement_.y_;
						if( yMoveUp / part.lagrangianSpeed_.y_ > 0.0 ) yFact = yMoveUp / part.lagrangianSpeed_.y_;
						else if( yMoveDown / part.lagrangianSpeed_.y_ > 0.0 ) yFact = yMoveDown / part.lagrangianSpeed_.y_;
						else{std::cout<<"there is a problem !!! case 4"<<std::endl;exit( 0 );}
					}

					//  we take the smallest movement ( to the next border, ie along x or y axis ) and perform the movement
					if( abs( xFact ) < abs( yFact ) ){
						nextPos = part.displacement_ + xFact  *  part.lagrangianSpeed_;
					}
					else{
						nextPos = part.displacement_ + yFact  *  part.lagrangianSpeed_;
					}

					//  if the particle is under the current cell, deposit the particle in the right place in the cell
					if( nextPos.z_ < terrain_[ index2d( currentPointCoord ) ] ){
						zFact = ( terrain_[ index2d( currentPointCoord ) ] - part.displacement_.z_ ) / part.lagrangianSpeed_.z_;
						part.displacement_ = part.displacement_ + zFact * part.lagrangianSpeed_;
						// particle collided a "horizontal side"
						particleDeposited = true;
					}
					//  else, continue ( the particle could be over the next cell, or has colided a "vertical side" )
					else{
						part.displacement_ = nextPos;
					}
				}

				particleOutOfTerrain = !isInBounds( Double2( { part.displacement_.x_, part.displacement_.y_ } ) );

				//currentPointCoord = { int( floor( part.displacement_.x_ / dx_ ) ) , int( floor( part.displacement_.y_ / dx_ ) ) };
				//particleOutOfTerrain = ( index2d( currentPointCoord ) < 0 ) || ( index2d( currentPointCoord ) >= discreteSize_.x_ * discreteSize_.y_ );

			}
			//  add the particle to the list of deposited particles...
			if( !particleOutOfTerrain ){
				storedParticles_.push_back( part );
				//if( part.familyId_ == 1 ){
				//	std::cout << part.displacement_.x_ << "," << part.displacement_.y_ << "," << part.displacement_.z_ << std::endl;
				//}
			}
			else garbage_.push_back( part );
		}
	}


	bool GridTerrain::isInGround( Double3 coord ){
		//  if the given point is under the nearest point of the terrain, the given point is on the ground

		double distSupLeft;
		double distSupRight;
		double distInfLeft;
		double distInfRight;

		Int2 coordSupLeft;
		Int2 coordSupRight;
		Int2 coordInfLeft;
		Int2 coordInfRight;

		Int2 coordClosest;

		//  compute the 4 nearest points of the terrain ( remember the terrain points are "in the middle" of the sites )
		coordInfLeft = { int( floor( ( coord.x_ - position_.x_ - ( dx_ / 2.0 ) ) / dx_ ) ) , int( floor( ( coord.y_ - position_.y_ - ( dx_ - 2.0 ) ) / dx_ ) ) };
		coordInfRight = { coordInfLeft.x_ + 1 , coordInfLeft.y_ };
		coordSupLeft = { coordInfLeft.x_ , coordInfLeft.y_ + 1 };
		coordSupRight = { coordInfLeft.x_ + 1 , coordInfLeft.y_ + 1 };

		//  compute distance between coord and points of the terrain
		distInfLeft = sqrt( ( coord.x_ - position_.x_ - coordInfLeft.x_ * dx_ + dx_ / 2.0 ) * ( coord.x_ - position_.x_ - coordInfLeft.x_ * dx_ + dx_ / 2.0 ) + ( coord.y_ - position_.y_ - coordInfLeft.y_ * dx_ + dx_ / 2.0 ) * ( coord.y_ - position_.y_ - coordInfLeft.y_ * dx_ + dx_ / 2.0 ) );
		distInfRight = sqrt( ( coord.x_ - position_.x_ - coordInfRight.x_ * dx_ + dx_ / 2.0 ) * ( coord.x_ - position_.x_ - coordInfRight.x_ * dx_ + dx_ / 2.0 ) + ( coord.y_ - position_.y_ - coordInfRight.y_ * dx_ + dx_ / 2.0 ) * ( coord.y_ - position_.y_ - coordInfRight.y_ * dx_ + dx_ / 2.0 ) );
		distSupLeft = sqrt( ( coord.x_ - position_.x_ - coordSupLeft.x_ * dx_ + dx_ / 2.0 ) * ( coord.x_ - position_.x_ - coordSupLeft.x_ * dx_ + dx_ / 2.0 ) + ( coord.y_ - position_.y_ - coordSupLeft.y_ * dx_ + dx_ / 2.0 ) * ( coord.y_ - position_.y_ - coordSupLeft.y_ * dx_ + dx_ / 2.0 ) );
		distSupRight = sqrt( ( coord.x_ - position_.x_ - coordSupRight.x_ * dx_ + dx_ / 2.0 ) * ( coord.x_ - position_.x_ - coordSupRight.x_ * dx_ + dx_ / 2.0 ) + ( coord.y_ - position_.y_ - coordSupRight.y_ * dx_ + dx_ / 2.0 ) * ( coord.y_ - position_.y_ - coordSupRight.y_ * dx_ + dx_ / 2.0 ) );


		/* std::cout<<"current coord :		"<<coord.x_<<","<<coord.y_<<std::endl;
		std::cout<<"coordInfLeft :		"<<coordInfLeft.x_<<","<<coordInfLeft.y_<<std::endl;
		std::cout<<"coordInfRight :		"<<coordInfRight.x_<<","<<coordInfRight.y_<<std::endl;
		std::cout<<"coordSupLeft:		"<<coordSupLeft.x_<<","<<coordSupLeft.y_<<std::endl;
		std::cout<<"coordSupRight :		"<<coordSupRight.x_<<","<<coordSupRight.y_<<std::endl;
		std::cout<<std::endl;

		assert( ( coordInfLeft.x_ < discreteSize_.x_ ) && ( coordInfLeft.x_ >= 0 ) ); assert( ( coordInfLeft.y_ < discreteSize_.y_ ) && ( coordInfLeft.y_ >= 0 ) );
		assert( ( coordInfRight.x_ < discreteSize_.x_ ) && ( coordInfRight.x_ >= 0 ) ); assert( ( coordInfRight.y_ < discreteSize_.y_ ) && ( coordInfRight.y_ >= 0 ) );
		assert( ( coordSupLeft.x_ < discreteSize_.x_ ) && ( coordSupLeft.x_ >= 0 ) ); assert( ( coordSupLeft.y_ < discreteSize_.y_ ) && ( coordSupLeft.y_ >= 0 ) );
		assert( ( coordSupRight.x_ < discreteSize_.x_ ) && ( coordSupRight.x_ >= 0 ) ); assert( ( coordSupRight.y_ < discreteSize_.y_ ) && ( coordSupRight.y_ >= 0 ) ); */

		//  take the closest point for the test
		coordClosest = coordInfLeft;
		if( distInfLeft <= distInfRight &&  distInfLeft <= distSupLeft && distInfLeft <= distSupRight ) coordClosest = coordInfLeft;
		else if( distInfRight <= distInfLeft &&  distInfRight <= distSupLeft && distInfRight <= distSupRight ) coordClosest = coordInfRight;
		else if( distSupLeft <= distInfRight &&  distSupLeft <= distInfLeft && distSupLeft <= distSupRight ) coordClosest = coordSupLeft;
		else if( distSupRight <= distInfRight &&  distSupRight <= distSupLeft && distSupRight <= distInfLeft ) coordClosest = coordSupRight;

		//  if coordClosest is out of bounds, take a point in bounds
		if( coordClosest.x_ < 0 ) coordClosest.x_ = 0;
		if( coordClosest.y_ < 0 ) coordClosest.y_ = 0;
		if( int( coordClosest.x_ ) >= discreteSize_.x_ ) coordClosest.x_ = discreteSize_.x_ - 1;
		if( int( coordClosest.y_ ) >= discreteSize_.y_ ) coordClosest.y_ = discreteSize_.y_ - 1;

		if( coord.z_ <= terrain_[ index2d( { int( coordClosest.x_ ), int( coordClosest.y_ ) } ) ] ) return true;
		else return false;

		//  check if coord.z is under the nearest point of the terrain
		/* if( distInfLeft <= distInfRight &&  distInfLeft <= distSupLeft && distInfLeft <= distSupRight ){
		if( coord.z_ <= terrain_[index2d( coordInfLeft )] ) return true;
		else return false;
	}
	else if( distInfRight <= distInfLeft &&  distInfRight <= distSupLeft && distInfRight <= distSupRight ){
	if( coord.z_ <= terrain_[index2d( coordInfRight )] ) return true;
	else return false;
}
else if( distSupLeft <= distInfRight &&  distSupLeft <= distInfLeft && distSupLeft <= distSupRight ){
if( coord.z_ <= terrain_[index2d( coordSupLeft )] ) return true;
else return false;
}
else if( distSupRight <= distInfRight &&  distSupRight <= distSupLeft && distSupRight <= distInfLeft ){
if( coord.z_ <= terrain_[index2d( coordSupRight )] ) return true;
else return false;
} */

return false;
}


std::vector<std::vector<int> > GridTerrain::getParticleDeposition( Int2 sampleSize ){
	std::vector< std::vector< int > > result;
	int maxId = 0;
	result.resize( sampleSize.x_ * sampleSize.y_ );
	//  get the highest family id in order to initialize each vector<int>
	for( auto pf : particleFamilies_ ) maxId = ( pf->familyId_ > maxId ) ? pf->familyId_ : maxId;
	//for( unsigned int i = 0; i < particleFamilies_.size(); i++ ) maxId = 0; //( (particleFamilies_[ i ]).get()->familyId_ > maxId ) ? (particleFamilies_[ i ]).get()->familyId_ : maxId;
	//for( std::shared_ptr< GenericParticleFamily > pf: particleFamilies_ ) maxId = ( pf->familyId_ > maxId ) ? pf->familyId_ : maxId;
	for( auto &site : result ) site.resize( maxId + 1 , 0 );
	//for( std::vector< int > &site: result ) site.resize( maxId + 1 , 0 );

	double dx = discreteSize_.x_  *  dx_  /  sampleSize.x_;
	int posX;
	int posY;

	for( auto&  part : storedParticles_ ){
		//for( Particle&  part: storedParticles_ ){
		posX = int( ( part.displacement_.x_ - position_.x_ ) / dx );
		posY = int( ( part.displacement_.y_ - position_.y_ ) / dx );
		if(  posX + posY * sampleSize.y_ > 0 && posX + posY * sampleSize.y_ < int( result.size() ) ){
			( result[ posX + posY * sampleSize.y_ ] )[ part.familyId_ ]++;
		}
	}
	return result;
}

std::vector<std::vector<int> > GridTerrain::getParticleDeposition(){
	return getParticleDeposition( discreteSize_ );
}

Int2 GridTerrain::getDiscreteSize(){
	return discreteSize_;
}

int GridTerrain::getXSize(){
	return discreteSize_.x_;
}

int GridTerrain::getYSize(){
	return discreteSize_.y_;
}

double GridTerrain::getDx(){
	return dx_;
}

Double2 GridTerrain::getPosition(){
	return position_;
}

int GridTerrain::getNumberOfParticles(){
	return storedParticles_.size() + garbage_.size();
}

int GridTerrain::countTerrainParticles(){
	return storedParticles_.size() + garbage_.size();
}

std::vector< int > GridTerrain::getStoredParticlesCounts(){
	std::vector< int > res( particleFamilies_.size(), 0 );
	for( auto part: storedParticles_ ) res[ part.familyId_ ]++;
	for( auto part: garbage_ ) res[ part.familyId_ ]++;
	return res;
}

} // namespace piaf
