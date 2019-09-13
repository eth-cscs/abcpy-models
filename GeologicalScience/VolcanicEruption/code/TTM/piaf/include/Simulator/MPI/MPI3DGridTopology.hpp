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


#ifndef MPI3DGRIDTOPOLOGY_H_
#define MPI3DGRIDTOPOLOGY_H_



#include "../SimulatorTypes.hpp"
//#include "mpi.h"

#include "assert.h"

#include "MPITopology.hpp"
#include <algorithm>
#include "boost/foreach.hpp"

//using namespace piaf;

// includes for mpi
//#include <boost/serialization/vector.hpp>
//#include <boost/mpi.hpp>

namespace piaf{

class MPI3DGridTopology : public MPITopology {
public:
	//MPI3DGridTopology( int argc, char** argv );

    MPI3DGridTopology( int argc, char** argv, shared_ptr<MpiManager> mpiManager, std::string parallelDimensionPriority = "xyz" );

    //MPI3DGridTopology( );

	virtual ~MPI3DGridTopology();
    virtual void init();

	void to2DTopology();

	//int getSize();

	//int getRank();

	//boost::mpi::environment* getEnv();

	//boost::mpi::communicator* getCommWorld();

	Int3 getTopology();

	int getXSize();

	int getYSize();

	int getZSize();

	Int3 getProcessPos();

	int getXProcessPos();

	int getYProcessPos();

	int getZProcessPos();

	int gridPos2Rank( Int3 pos );

protected:

	// non copyable ( but not necessarily a singleton )
	MPI3DGridTopology( const MPI3DGridTopology& other );
	MPI3DGridTopology& operator=( const MPI3DGridTopology& );

	/*int rank_;

	int nProcess_;*/

	/*boost::mpi::environment* mpiEnv_;
	boost::mpi::communicator* mpiWorld_;*/

	Int3 topology_;

	Int3 processPos_;

   std::string parallelDimensionPriority;
};

} // namespace piaf

#endif /* MPI3DGRIDTOPOLOGY_H_ */
