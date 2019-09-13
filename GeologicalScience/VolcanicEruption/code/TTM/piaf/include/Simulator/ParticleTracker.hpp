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

#ifndef PARTICLETRACKER_H_
#define PARTICLETRACKER_H_

#include "GenericParticleFamily.hpp"
#include <vector>
#include <tuple>
#include <boost/shared_ptr.hpp>
#include <stdexcept>

namespace piaf{

class ParticleTracker {

protected :

  std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies_;

  std::vector< double > times_;

  std::vector< std::vector< Particle > > particles_;

  double interval_;

  double lastTrackTime;

public:

  ParticleTracker( std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies );

  ParticleTracker( std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies, double interval );

  virtual ~ParticleTracker();

  virtual void clearParticles();

  virtual void addParticle( Particle particle, double time );

  virtual std::tuple< std::vector< double >, std::vector< std::vector< Particle > > > getParticles();

  virtual std::tuple< std::vector< double >, std::vector< std::vector< Particle > > > removeParticles();

  virtual std::vector< double > getTimes();

  virtual double getInterval();

  virtual void setInterval( double interval );

  virtual double getLastTrackTime();



};

} // namespace piaf

#endif /* PARTICLEREPOSITORY_H_ */
