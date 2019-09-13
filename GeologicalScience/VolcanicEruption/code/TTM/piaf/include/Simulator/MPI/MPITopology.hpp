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

/*
 * Updated By: Mohamed Ben Belgacem to use the  MPI class manager
 * date: 07.12.2015
 * */


#ifndef MPITOPOLOGY_H_
#define MPITOPOLOGY_H_

#include "mpi.h"
#include <memory>

#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
//#include <boost/mpi.hpp>
#include "musclehpc/parallelism/mpiManager.h"

using namespace std;

namespace piaf{

class MPITopology {
public:
    MPITopology( int argc, char** argv,  shared_ptr<MpiManager> mpiManager);
    //MPITopology( );
	virtual ~MPITopology();
    virtual void init();
	int getSize();
	int getRank();
    MpiManager & getMpiManager() const;

    //MpiManager * getMpiManagerPtr() const;
    void setMpiManager(shared_ptr<MpiManager> mpiManager);

protected:
	// non copyable ( but not necessarily a singleton )
	MPITopology( const MPITopology& other );
	MPITopology& operator=( const MPITopology& );
    shared_ptr<MpiManager> mpiManager;

};

} // namespace piaf

#endif /* MPITOPOLOGY_H_ */
