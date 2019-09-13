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

#include "../../../include/Simulator/MPI/MPIGridTerrain.hpp"

namespace piaf{

	MPIGridTerrain::MPIGridTerrain( Double2 pos, Double2 size, double dx, std::vector< std::shared_ptr< GenericParticleFamily > > particleFamilies, MPI3DGridTopology* topology ):
	GridTerrain	( pos, size, dx, particleFamilies ),
	topology_	( topology )
	{

	}

		MPIGridTerrain::~MPIGridTerrain() {
			// TODO Auto-generated destructor stub
		}

		void MPIGridTerrain::gather(){


			gatheredStoredParticles_.clear();
			gatheredGarbage_.clear();


			int sendSize;
			int recvSize;
			//MPI_Status status;
			MpiManager & mpiManager= this->topology_->getMpiManager();
			// using a simple gather instead of mpi gather size is not the same for each process, see http://stackoverflow.com/questions/2546298/vector-usage-in-mpic for serialization
			/*! \todo replace by mpi_gatherv */
			if( topology_->getRank() == 0 ){
				gatheredStoredParticles_.resize( topology_->getSize() );
				gatheredGarbage_.resize( topology_->getSize() );
				gatheredStoredParticles_.push_back( storedParticles_ );
				gatheredGarbage_.push_back( garbage_ );
				for( int neighbor = 1; neighbor < topology_->getSize(); neighbor++ ){
					std::vector< Particle > receivedPart;
					std::vector< Particle > receivedGarb;
					int tag=0;
					// receive data size particles
					mpiManager.receive<int> (&recvSize, 1, neighbor, tag);
					receivedPart.resize(recvSize);
					// receive data particles
					//std::vector< char > incomingBuffer( recvSize );
					//mpiManager.receive<char> (&incomingBuffer[ 0 ], recvSize, neighbor, tag);
					mpiManager.receive<char> ((char*)receivedPart.data(), recvSize*sizeof(Particle), neighbor, tag);
					// unpack data particles
					//std::istringstream iss( std::string( &incomingBuffer[ 0 ], incomingBuffer.size() ) );
					//boost::archive::binary_iarchive ia(iss);
					//ia >> BOOST_SERIALIZATION_NVP( receivedPart );

					// receive data size garbage
					mpiManager.receive<int> ( &recvSize, 1, neighbor, tag );
					receivedGarb.resize(recvSize);
					// receive data garbage
					//incomingBuffer.resize( recvSize );
					//mpiManager.receive<char> (&incomingBuffer[ 0 ], recvSize, neighbor, tag);
					mpiManager.receive<char> ((char*)receivedGarb.data(), recvSize*sizeof(Particle), neighbor, tag);
					// unpack data garbage
					//std::istringstream iss2( std::string( &incomingBuffer[ 0 ], incomingBuffer.size() ) );
					//boost::archive::binary_iarchive ia2(iss2);
					//ia2 >> BOOST_SERIALIZATION_NVP( receivedGarb );

					// put received data in local structure
					gatheredStoredParticles_.push_back( receivedPart );
					gatheredGarbage_.push_back( receivedGarb );
				}

			}
			else{
				int tag=0;
				// prepare data particles to send

				/*std::ostringstream oss;
				{
					boost::archive::binary_oarchive oa( oss );
					oa << BOOST_SERIALIZATION_NVP( storedParticles_ );
				}
				sendSize = oss.str().size();*/

				sendSize = storedParticles_.size();

				// send particles size

				mpiManager.send<int>(&sendSize, 1,  mpiManager.bossId(), tag);
				// send particles data
				//mpiManager.send<char>((char *)oss.str().c_str(), sendSize, mpiManager.bossId(), tag );
				mpiManager.send<char>((char*)storedParticles_.data(), sendSize*sizeof(Particle), mpiManager.bossId(), tag );

				// prepare data garbage to send
			/*	std::ostringstream oss2;
				{
					boost::archive::binary_oarchive oa2( oss2 );
					oa2 << BOOST_SERIALIZATION_NVP( garbage_ );
				}
				sendSize = oss2.str().size();*/

				sendSize = garbage_.size();

				// send garbage size
				mpiManager.send<int>(&sendSize, 1,  mpiManager.bossId(), tag);
				// send garbage data
				//mpiManager.send<char>((char*)oss2.str().c_str(), sendSize, mpiManager.bossId(), tag);
				mpiManager.send<char>((char*)garbage_.data(), sendSize*sizeof(Particle), mpiManager.bossId(), tag);

			}


		}

		std::vector< std::vector< int > > MPIGridTerrain::getParticleDeposition( Int2 sampleSize ){
			/*! \todo: now that the gather has to be made by MPIGridTerrain, a gather must be done at each call... better to give this responsibility to the client ? (see in getNumberOfParticles, same thing) */
			gather();
			std::vector< std::vector< int > > result;
			int maxId = 0;
			result.resize( sampleSize.x_ * sampleSize.y_ );
			//  get the highest family id in order to initialize each vector<int>
			//for( std::shared_ptr< GenericParticleFamily > pf: particleFamilies_ ) maxId = ( pf->familyId_ > maxId ) ? pf->familyId_ : maxId;
			for( auto pf: particleFamilies_ ){ maxId = ( pf->familyId_ > maxId ) ? pf->familyId_ : maxId; }
			//for( std::vector< int > &site: result ) site.resize( maxId + 1 , 0 );
			for( auto &site: result ){ site.resize( maxId + 1 , 0 ); }

			double dx = discreteSize_.x_  *  dx_  /  sampleSize.x_;
			int posX;
			int posY;

			//for( std::vector< Particle >& vect: gatheredStoredParticles_ ){
			for( auto& vect: gatheredStoredParticles_ ){
				//for( Particle&  part: vect ){
				for( auto&  part: vect ){
					posX = int( ( part.displacement_.x_ - position_.x_ ) / dx );
					posY = int( ( part.displacement_.y_ - position_.y_ ) / dx );
					//if(  posX + posY * sampleSize.y_ > 0 && posX + posY * sampleSize.y_ < int( result.size() ) ){
					//	( result[ posX + posY * sampleSize.y_ ] )[ part.familyId_ ]++;
					//}
					if(  posX + posY * sampleSize.x_ > 0 && posX + posY * sampleSize.x_ < int( result.size() ) ){
						( result[ posX + posY * sampleSize.x_ ] )[ part.familyId_ ]++;
					}
				}
			}

			return result;
		}

		std::vector< std::vector< int > > MPIGridTerrain::getParticleDeposition(){
			return getParticleDeposition( discreteSize_ );
		}

		int MPIGridTerrain::getNumberOfParticles(){
			gather();

			int res = 0;
			//for( std::vector< Particle >& vect: gatheredStoredParticles_ ) res += vect.size();
			BOOST_FOREACH( std::vector< Particle >& vect, gatheredStoredParticles_ ) res += vect.size();
			//for( std::vector< Particle >& vect: gatheredGarbage_ ) res += vect.size();
			BOOST_FOREACH( std::vector< Particle >& vect, gatheredGarbage_ ) res += vect.size();

			return res;
		}

		int MPIGridTerrain::countTerrainParticles(){
			MpiManager & mpiManager= this->topology_->getMpiManager();
			int cpt = 0;
			int totalCpt;
			cpt = storedParticles_.size() + garbage_.size();
			mpiManager.reduce<int>(cpt,totalCpt, MPI_SUM,  mpiManager.bossId() );

			return totalCpt;
		}

		MPI3DGridTopology* MPIGridTerrain::getTopology(){
			return topology_;
		}

		std::vector< int > MPIGridTerrain::getStoredParticlesCounts(){
		        MpiManager &mpiManager=this->topology_->getMpiManager();

		        std::vector< int > resLoc( particleFamilies_.size(), 0 );
		        std::vector< int > res( particleFamilies_.size(), 0 );
		        for( auto part : storedParticles_ ) resLoc[ part.familyId_ ]++;
				for( auto part : garbage_ ) resLoc[ part.familyId_ ]++;

		        mpiManager.reduceVect<int>( resLoc, res, MPI_SUM );
		        return res;
		}

	} // namespace piaf
