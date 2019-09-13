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


#include "../../include/Simulator/ParticleRepository.hpp"

namespace piaf{

ParticleRepository::ParticleRepository( std::vector< std::shared_ptr< GenericParticleFamily> > particleFamilies ):
		isActive_			( true ),
		particleFamilies_	( particleFamilies )
{}

ParticleRepository::~ParticleRepository() {}

void ParticleRepository::depositParticle( Particle part ){
	if( isActive_ ){
		// just add the particle to the list of deposited particles...
		storedParticles_.push_back( part );
	}
}

void ParticleRepository::depositParticle( BoundaryParticle part ){
	if( isActive_ ){
		// just add the particle to the list of deposited particles...
		storedBoundaryParticles_.push_back( part );
	}
}

int ParticleRepository::getNumberOfParticles(){
	return storedParticles_.size();
}

void ParticleRepository::setActive( bool state ){
	isActive_ = state;
}

std::vector< std::shared_ptr< GenericParticleFamily > > ParticleRepository::getParticleFamilies(){
	return particleFamilies_;
}

std::vector< int > ParticleRepository::getStoredParticlesCounts(){
	std::vector< int > res( particleFamilies_.size(), 0 );
	for( auto part: storedParticles_ ) res[ part.familyId_ ]++;
	return res;
}

}// namespace piaf
