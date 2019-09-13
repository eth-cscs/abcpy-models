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

#ifndef MPIGRIDTERRAIN_H_
#define MPIGRIDTERRAIN_H_

#include "../GridTerrain.hpp"
#include "MPI3DGridTopology.hpp"
#include <boost/serialization/vector.hpp>
#include <sstream>

namespace piaf{

class MPIGridTerrain : public GridTerrain {
public:

	MPIGridTerrain( Double2 pos, Double2 size, double dx, std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies, MPI3DGridTopology* topology );
	virtual ~MPIGridTerrain();

	std::vector< std::vector< int > > getParticleDeposition( Int2 sampleSize );

	std::vector< std::vector< int > > getParticleDeposition();

	MPI3DGridTopology* getTopology();

	int getNumberOfParticles();

  int countTerrainParticles();

  /*!
  * return the number of particles stored for each class
  */
  virtual std::vector< int > getStoredParticlesCounts();

protected:

	void gather();

	/*! \brief The list of particles which fall outside the terrain ( result of the gather, relevant only for process 0 ) */
	std::vector< std::vector< Particle > > gatheredGarbage_;

	/*! \brief The list of stored particles ( result of the gather, relevant only for process 0 ) */
	std::vector< std::vector< Particle > > gatheredStoredParticles_;

	MPI3DGridTopology* topology_;



};

}// namespace piaf

#endif /* MPIGRIDTERRAIN_H_ */
