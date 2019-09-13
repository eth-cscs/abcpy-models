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


#include "../../../include/Simulator/MPI/MPI3DGridTopology.hpp"

namespace piaf{

MPI3DGridTopology::MPI3DGridTopology(int argc, char** argv, shared_ptr<MpiManager> mpiManager, std::string parallelDimensionPriority ):
                MPITopology( argc, argv, mpiManager )
{
    this->parallelDimensionPriority=parallelDimensionPriority;

}

void MPI3DGridTopology::init(){
	// convert priorities to indexes (reverse order)
	std::vector< int > ind;
	BOOST_FOREACH( char &c, parallelDimensionPriority ){
		if( c == 'x' ) ind.push_back( 2 );
		if( c == 'y' ) ind.push_back( 1 );
		if( c == 'z' ) ind.push_back( 0 );
	}
	// try to decompose the number of process in 3 factors
	int nTemp = getSize();//getSize();
	//int values[3] = { 1, 1, 1 };
	std::vector< int > values = { 1, 1, 1 };

	int val = ceil( pow( double( nTemp ), 1.0 / 3.0 ) );
	while( nTemp != ( nTemp / val ) * val ) val--;
	nTemp = nTemp / val;
	values[ 0 ] = val;
	val = pow( double( nTemp ), 1.0 / 2.0 );
	while( nTemp != ( nTemp / val ) * val )val--;
	nTemp = nTemp / val;
	values[ 1 ] = val;
	values[ 2 ] = getSize() / ( values[ 0 ] * values[ 1 ] );





	/*for( int ival = 0; ival < 3; ival++ ){
				if( ival != 2 ){
					for( int i = 2; i <= nTemp; i++ ){
						if( ( nTemp / i ) * i == nTemp ){
							values[ ival ] = i;
							break;
						}
					}
					if ( values[ ival ] == nTemp ){
						break;
					}
					nTemp = nTemp / values[ ival ];
				}
				else values[ ival ] = nTemp;
			}*/


	/* put the maximum of process on x axis... */
	//if(  ( values[ 0 ] >= values[ 1 ] ) && ( values[ 0 ] >= values[ 2 ] ) ) topology_ = { values[ 0 ], values[ 1 ], values[ 2 ]};
	//else if(  ( values[ 1 ] >= values[ 0 ] ) && ( values[ 1 ] >= values[ 2 ] ) ) topology_ = { values[ 1 ], values[ 2 ], values[ 0 ]};
	//else topology_ = { values[ 2 ], values[ 0 ], values[ 1 ]};
	//topology_ = { values[ 0 ], values[ 1 ], values[ 2 ]};

	std::sort( values.begin(), values.end() );

	topology_ = { values[ ind[ 0 ] ], values[ ind[ 1 ] ], values[ ind[ 2 ] ] };

	// test...
	if( getSize() % 2 == 0 ){
		topology_ = { getSize() / 2 , 2 , 1 };
		// make x and y dim as balanced as possible
		if( topology_.x_ > 4 ){
			while( ( abs( ( topology_.x_ / 2 ) - ( topology_.y_ * 2 ) ) < abs( ( topology_.x_ ) - ( topology_.y_  ) ) ) && ( ( topology_.x_ / 2 ) * 2 ) == topology_.x_ ){
				topology_.x_ /= 2;
				topology_.y_ *= 2;
			}
		}
	}

	//std::cout << "topology : " << topology_.x_ << ", " << topology_.y_ << ", " << topology_.z_ << std::endl;

	processPos_ = { getRank() % topology_.x_ , ( getRank() % ( topology_.x_ * topology_.y_ ) ) / topology_.x_ , getRank() / ( topology_.x_ * topology_.y_ )};

	assert( int( gridPos2Rank( processPos_ ) ) == getRank() );

}

//MPI3DGridTopology::MPI3DGridTopology( int argc, char** argv ):
//		MPI3DGridTopology( argc, argv,  )
//{}

/*MPI3DGridTopology::MPI3DGridTopology(  ):
						MPITopology()
//mpiEnv_		( NULL ),
//mpiWorld_	( NULL )
{}*/

MPI3DGridTopology::~MPI3DGridTopology() {
	//delete mpiEnv_;
	//delete mpiWorld_;
}

/*int MPI3DGridTopology::getSize(){
	//return mpiWorld.size();
	return getSize();
}

int MPI3DGridTopology::getRank(){
	//return getRank();
	return getRank();
}

boost::mpi::environment* MPI3DGridTopology::getEnv(){
	return mpiEnv_;
}

boost::mpi::communicator* MPI3DGridTopology::getCommWorld(){
	return mpiWorld_;
}*/

void MPI3DGridTopology::to2DTopology(){
	int nTemp = getSize();
	int values[3] = { 1, 1, 1};

	for( int ival = 0; ival < 2; ival++ ){
		if( ival != 1 ){
			for( int i = 2; i <= nTemp; i++ ){
				if( ( nTemp / i ) * i == nTemp ){
					values[ ival ] = i;
					break;
				}
			}
			if ( values[ ival ] == nTemp ){
				break;
			}
			nTemp = nTemp / values[ ival ];
		}
		else values[ ival ] = nTemp;
	}

	// for a volcanic simulation with volcanic column, its probably preferable that the z axis is the greatest axis
	if( values[ 1 ] > values[ 0 ] ) topology_ = { values[ 0 ], 1, values[ 1 ] };
	else topology_ = { values[ 1 ], 1, values[ 0 ] };

	processPos_ = { getRank() % topology_.x_ , ( getRank() % ( topology_.x_ * topology_.y_ ) ) / topology_.x_ , getRank() / ( topology_.x_ * topology_.y_ )};
}

int MPI3DGridTopology::gridPos2Rank( Int3 pos ){
	return pos.x_ + pos.y_ * topology_.x_ + pos.z_ * topology_.x_ * topology_.y_;
}

Int3 MPI3DGridTopology::getTopology(){ return topology_; }

int MPI3DGridTopology::getXSize(){ return topology_.x_; }

int MPI3DGridTopology::getYSize(){ return topology_.y_; }

int MPI3DGridTopology::getZSize(){ return topology_.z_; }

int MPI3DGridTopology::getXProcessPos(){ return processPos_.x_; }

int MPI3DGridTopology::getYProcessPos(){ return processPos_.y_; }

int MPI3DGridTopology::getZProcessPos(){ return processPos_.z_; }

} // namespace piaf
