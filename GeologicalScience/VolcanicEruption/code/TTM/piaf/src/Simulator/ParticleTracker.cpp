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

#include "../../include/Simulator/ParticleTracker.hpp"

namespace piaf{

ParticleTracker::ParticleTracker( std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies ):
    particleFamilies_( particleFamilies ),
    interval_( -1.0 ),
    lastTrackTime(0.0)
{}

ParticleTracker::ParticleTracker( std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies, double interval ):
    particleFamilies_( particleFamilies ),
    interval_( interval ),
    lastTrackTime(0.0)
{}

ParticleTracker::~ParticleTracker(){}

void ParticleTracker::addParticle( Particle particle, double time ){
    if( times_.size()>0 && times_.back() == time ) particles_.back().push_back( particle );
    else if( time > times_.back() ){
        times_.push_back( time );
        lastTrackTime = time;
        particles_.push_back( std::vector< Particle >() );
        addParticle( particle, time );
    }
    else{
        throw std::runtime_error("Error in ParticleTracker, particles must be inserted in increasing time order");
    }
}

void ParticleTracker::clearParticles(){
    times_.clear();
    particles_.clear();
}

std::tuple< std::vector< double >, std::vector< std::vector< Particle > > > ParticleTracker::getParticles(){
    return make_tuple( times_, particles_ );
}

std::tuple< std::vector< double >, std::vector< std::vector< Particle > > > ParticleTracker::removeParticles(){
    auto res = getParticles();
    clearParticles();
    return res;
}

double ParticleTracker::getInterval(){
    return interval_;
}

void ParticleTracker::setInterval( double interval ){
    interval_ = interval;
}

double ParticleTracker::getLastTrackTime(){
    /*if( times_.size() > 0 ) return times_.back();
    return -1.0;*/
    return lastTrackTime;
}

std::vector< double > ParticleTracker::getTimes(){
    return times_;
}

} // namespace piaf
