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


#include "../../../include/Simulator/MPI/MPISimpleDomain.hpp"

namespace piaf{

void MPISimpleDomain::commonConstruct( Double3 pos, Double3 size, double dx , MPITopology* topology, int parallelDim ){
	parallelDim_ = parallelDim;

	// compute the number of elements per block
	if( parallelDim == 0 ) blockSize_ = { int( ceil( ( size.x_ / dx ) / (double)topology->getSize() ) ), int( size.y_ / dx), int( size.z_ / dx ) };
	else if( parallelDim == 1 ) blockSize_ = { int( size.x_ / dx), int( ceil( ( size.y_ / dx ) / (double)topology->getSize() ) ), int( size.z_ / dx ) };
	else if( parallelDim == 2 ) blockSize_ = { int( size.x_ / dx), int( size.y_ / dx), int( ceil( ( size.z_ / dx ) / (double)topology->getSize() ) ) };
	else{
		std::cout << "invalid parallel dimension, must be [0..2]" << std::endl;
		exit(0);
	}
	//blockSize_ = { (int)( ( ( size.x_ / dx ) / (double)topology->getSize() ) + 0.5 ), int( size.y_ / dx) , int( size.z_ / dx ) };
	//if( blockSize_.x_ ) blockSize_.x_ = 1;
	/*
			// the last process takes the remaining elemets
			if( ( ( blockSize_.x_ * topology->getSize() )  < ( size.x_ / dx ) ) && topology->getRank() == topology->getSize() - 1 ){
				blockSize_.x_ += ( size.x_ / dx ) - ( blockSize_.x_ * topology->getSize() );
			}*/
	// in the simple domain, there is only one block per process
	numBlockLocal_ = { 1, 1, 1 };

	if( parallelDim == 0 ) domainBlocks_.push_back( Domain( { position_.x_ + blockSize_.x_ * dx * topology->getRank(), position_.y_, position_.z_ }, { blockSize_.x_ * dx_ , blockSize_.y_ * dx_ , blockSize_.z_ * dx_ }, dx ) );
	else if( parallelDim == 1 ) domainBlocks_.push_back( Domain( { position_.x_ , position_.y_ + blockSize_.y_ * dx * topology->getRank(), position_.z_ }, { blockSize_.x_ * dx_ , blockSize_.y_ * dx_ , blockSize_.z_ * dx_ }, dx ) );
	else if( parallelDim == 2 ) domainBlocks_.push_back( Domain( { position_.x_ , position_.y_, position_.z_ + blockSize_.z_ * dx * topology->getRank() }, { blockSize_.x_ * dx_ , blockSize_.y_ * dx_ , blockSize_.z_ * dx_ }, dx ) );

	particlePackets_.resize( 2 );

	// we have to send and receive to / from 2 neighbors
	incomingBuffers_.resize( 2 );
	sendBuffers_.resize( 2 );
	//sendBuffers_ = std::vector< std::ostringstream >( 2 );
	recvSizes_.resize( 2 );
	sendSizes_.resize( 2 );
	requestsSend_.resize( 2 );
	statusesSend_.resize( 2 );
	requestsRecv_.resize( 2 );
	statusesRecv_.resize( 2 );
}

MPISimpleDomain::MPISimpleDomain( Double3 pos, Double3 size, double dx , MPITopology* topology ):
				MPIAbstractEulerianDomain( pos, size, dx, topology )
{
	commonConstruct( pos, size, dx , topology, 0 );
}

MPISimpleDomain::MPISimpleDomain( Double3 pos, Double3 size, double dx , MPITopology* topology, int parallelDim ):
				MPIAbstractEulerianDomain( pos, size, dx, topology )
{
	commonConstruct( pos, size, dx , topology, parallelDim );
}

MPISimpleDomain::MPISimpleDomain():
				MPIAbstractEulerianDomain()
{

}

MPISimpleDomain::~MPISimpleDomain() {
	// TODO Auto-generated destructor stub
}

Domain* MPISimpleDomain::getDomainBlockAtLocalIndex3d( Int3 pos ){
	/*! \todo handle with exception */
	if( ( pos.x_ != 0 ) || ( pos.y_ != 0 ) || ( pos.z_ != 0 ) ){
		std::cout << "wrong block request"<< std::endl;
		exit( 0 );
	}
	return &( domainBlocks_[ 0 ] );
}

void MPISimpleDomain::putParticleInPacket( Particle part, Domain* lastDomainBlock ){

	if( parallelDim_ == 0 ){
		if( part.displacement_.x_ < lastDomainBlock->getPosition().x_ ) particlePackets_[ 0 ].push_back( part );
		else particlePackets_[ 1 ].push_back( part );
	}
	else if( parallelDim_ == 1 ){
		if( part.displacement_.y_ < lastDomainBlock->getPosition().y_ ) particlePackets_[ 0 ].push_back( part );
		else particlePackets_[ 1 ].push_back( part );
	}
	else if( parallelDim_ == 2 ){
		if( part.displacement_.z_ < lastDomainBlock->getPosition().z_ ) particlePackets_[ 0 ].push_back( part );
		else particlePackets_[ 1 ].push_back( part );
	}


}

void MPISimpleDomain::exchangeParticles(){

	requestExchangeParticles();
	completeExchangeParticles();

}

/*! \todo
 * utiliser MPI_PROC_NULL
 */
void MPISimpleDomain::requestExchangeParticles(){


	int sourceNeighbor = topology_->getRank() - 1;
	int destNeighbor = topology_->getRank() + 1;
	std::vector < Particle > receivedPacket;

	//std::vector<char> incomingBuffer;
	//int recvSize, sendSize;
	MPI_Request requestLocal;
	MPI_Status statusLocal;
    MpiManager & mpiManager = this->topology_->getMpiManager();
    int tag=0;
    int tag_irec=1;
    int tag_isend=1;
	for( int i = 0; i < 2; i++ ){

		if( ( sourceNeighbor > -1 ) && ( sourceNeighbor < topology_->getSize() ) ){
			// begin asynchronous receive of size
            mpiManager.iRecv<int>(&recvSizes_[ i ], 1, sourceNeighbor, &requestLocal, tag);
		}

		std::ostringstream oss;
		if( ( destNeighbor > -1 ) && ( destNeighbor < topology_->getSize() ) ){
			// package (serialize) data to send
			//std::ostringstream oss;
			{
				boost::archive::binary_oarchive oa( oss );
				if( destNeighbor < topology_->getRank() ) oa << BOOST_SERIALIZATION_NVP( particlePackets_[ 0 ] );
				else oa << BOOST_SERIALIZATION_NVP( particlePackets_[ 1 ] );
			}
			// send size of data
			sendSizes_[ i ] = oss.str().size();
            mpiManager.send<int>(&sendSizes_[ i ], 1, destNeighbor, tag);
		}

		if( ( sourceNeighbor > -1 ) && ( sourceNeighbor < topology_->getSize() ) ){
			// wait for size to be received
            mpiManager.wait(&requestLocal, &statusLocal );
            //mohamed: MPI_Wait( &requestLocal, &statusLocal );
			// begin asynchronous receive of data
			incomingBuffers_[ i ].resize( recvSizes_[ i ] );
             mpiManager.iRecv<char>(&(incomingBuffers_[ i ])[ 0 ], recvSizes_[ i ], sourceNeighbor, &requestsRecv_[ i ], tag_irec);
		}

		if( ( destNeighbor > -1 ) && ( destNeighbor < topology_->getSize() ) ){
			// send data
			sendBuffers_[ i ] = oss.str();
            mpiManager.iSend<char>((char*)sendBuffers_[ i ].c_str(), sendSizes_[ i ], destNeighbor,  &requestsSend_[ i ] , tag_isend);
		}


		// second step, send to rank - 1 and receive from rank + 1
		sourceNeighbor = topology_->getRank() + 1;
		destNeighbor = topology_->getRank() - 1;
	}

}

void MPISimpleDomain::completeExchangeParticles(){

    MpiManager & mpiManager = this->topology_->getMpiManager();
	std::vector < Particle > receivedPacket;

	int sourceNeighbor = topology_->getRank() - 1;
	int destNeighbor = topology_->getRank() + 1;

	// wait for data to be sent and received
	for( int i = 0; i < 2; i++ ){
		// wait for data to be received
        if( ( sourceNeighbor > -1 ) && ( sourceNeighbor < topology_->getSize() ) ){
            mpiManager.wait( &requestsRecv_[ i ], &statusesRecv_[ i ] );
            //mohamed: MPI_Wait( &requestsRecv_[ i ], &statusesRecv_[ i ] );
        }
		// wait for data to be sent
        if( ( destNeighbor > -1 ) && ( destNeighbor < topology_->getSize() ) ) {
            mpiManager.wait( &requestsSend_[ i ], &statusesSend_[ i ] );
            //mohamed: MPI_Wait( &requestsSend_[ i ], &statusesSend_[ i ] );
        }

		sourceNeighbor = topology_->getRank() + 1;
		destNeighbor = topology_->getRank() - 1;
	}

	sourceNeighbor = topology_->getRank() - 1;
	destNeighbor = topology_->getRank() + 1;
	for( int i = 0; i < 2; i++ ){
		if( ( sourceNeighbor > -1 ) && ( sourceNeighbor < topology_->getSize() ) ){
			// unpack received data
			std::istringstream iss( std::string( &(incomingBuffers_[ i ])[ 0 ], incomingBuffers_[ i ].size() ) );
			boost::archive::binary_iarchive ia( iss );
			ia >> BOOST_SERIALIZATION_NVP( receivedPacket );

			for( auto& part: receivedPacket ){
				putParticle( part );
			}

			receivedPacket.clear();
		}

		if( ( destNeighbor > -1 ) && ( destNeighbor < topology_->getSize() ) ){
			if( destNeighbor < topology_->getRank() ) particlePackets_[ 0 ].clear();
			else particlePackets_[ 1 ].clear();
		}

		sourceNeighbor = topology_->getRank() + 1;
		destNeighbor = topology_->getRank() - 1;
	}

}

void MPISimpleDomain::insertParticles( Double3 coord, int incr, int familyId, double scaling ){
    containsParticles_ = 1;
    //if( ( domainBlocks_[ 0 ].getPosition().x_ <= coord.x_ ) && ( domainBlocks_[ 0 ].getPosition().x_ + domainBlocks_[ 0 ].getSize().x_ * dx_ > coord.x_) ){
    if( domainBlocks_[ 0 ].isInBounds( { int( round( ( coord.x_ - domainBlocks_[ 0 ].position_.x_ ) / dx_ ) ) , int( round( ( coord.y_ - domainBlocks_[ 0 ].position_.y_ ) / dx_  ) ) , int( round( ( coord.z_ - domainBlocks_[ 0 ].position_.z_ ) / dx_ ) ) } ) ){
        domainBlocks_[ 0 ].insertParticles( coord, incr, familyId, scaling );
    }
}

// void MPISimpleDomain::insertParticles( Double3 coord, int incr, int familyId ){
// 	containsParticles_ = 1;
// 	//if( ( domainBlocks_[ 0 ].getPosition().x_ <= coord.x_ ) && ( domainBlocks_[ 0 ].getPosition().x_ + domainBlocks_[ 0 ].getSize().x_ * dx_ > coord.x_) ){
// 	if( domainBlocks_[ 0 ].isInBounds( { int( round( ( coord.x_ - domainBlocks_[ 0 ].position_.x_ ) / dx_ ) ) , int( round( ( coord.y_ - domainBlocks_[ 0 ].position_.y_ ) / dx_  ) ) , int( round( ( coord.z_ - domainBlocks_[ 0 ].position_.z_ ) / dx_ ) ) } ) ){
// 		domainBlocks_[ 0 ].insertParticles( coord, incr, familyId );
// 	}
// }

Int3 MPISimpleDomain::global2Local( Int3 coord ){
	if( parallelDim_ == 0 ){
		return { coord.x_ % blockSize_.x_, coord.y_, coord.z_ };
	}
	else if( parallelDim_ == 1 ){
		return { coord.x_, coord.y_ % blockSize_.y_, coord.z_ };
	}
	else if( parallelDim_ == 2 ){
		return { coord.x_, coord.y_, coord.z_ % blockSize_.z_ };
	}
	return { -1, -1 , -1 };
}

void MPISimpleDomain::insertParticlesInBuffer( Int3 coord, int incr, int familyId ){
	containsParticles_ = 1;
	//if( ( topology_->getRank() * blockSize_.x_ <= coord.x_ ) && ( topology_->getRank() * blockSize_.x_ + blockSize_.x_ > coord.x_) ){
	if( domainBlocks_[ 0 ].isInBounds( global2Local( coord ) ) ){
		domainBlocks_[ 0 ].insertParticlesInBuffer( global2Local( coord ), incr, familyId );
	}
}

Double3 MPISimpleDomain::getSpeed( Int3 coord, int familyId ){
	//if( ( topology_->getRank() * blockSize_.x_ <= coord.x_ ) && ( topology_->getRank() * blockSize_.x_ + blockSize_.x_ > coord.x_) ){
	if( domainBlocks_[ 0 ].isInBounds( global2Local( coord ) ) ){
		return domainBlocks_[ 0 ].getSpeed( global2Local( coord ), familyId );
	}
	//log<LOG_DEBUG>(" warning - process %1% called getSpeed for coord %5% ") % topology_->getRank() % coord.x_ ;
	return { 0.0, 0.0, 0.0 };
}

bool MPISimpleDomain::isInGroundContactAt( Int3 coord ){
	//if( ( topology_->getRank() * blockSize_.x_ <= coord.x_ ) && ( topology_->getRank() * blockSize_.x_ + blockSize_.x_ > coord.x_) ){
	if( domainBlocks_[ 0 ].isInBounds( global2Local( coord ) ) ){
		return domainBlocks_[ 0 ].isInGroundContactAt( global2Local( coord ) );
	}
	//log<LOG_DEBUG>(" warning - process %1% called isInGroundContactAt for coord %5% ") % topology_->getRank() % coord.x_ ;
	return false;
}

void MPISimpleDomain::putParticle( Particle part, Int3 coord ){
	containsParticles_ = 1;
	//if( ( topology_->getRank() * blockSize_.x_ <= coord.x_ ) && ( topology_->getRank() * blockSize_.x_ + blockSize_.x_ > coord.x_) ){
	if( domainBlocks_[ 0 ].isInBounds( global2Local( coord ) ) ){
		domainBlocks_[ 0 ].putParticle( part, global2Local( coord ) );
	}
}

void MPISimpleDomain::putParticle( Particle part ){
	containsParticles_ = 1;

	//	if( ( domainBlocks_[ 0 ].getPosition().x_ <= part.displacement_.x_ ) && ( domainBlocks_[ 0 ].getPosition().x_ + domainBlocks_[ 0 ].getSize().x_ * dx_ > part.displacement_.x_) ){
	if( domainBlocks_[ 0 ].isInBounds( { int( round( ( part.displacement_.x_ - domainBlocks_[ 0 ].position_.x_ ) / dx_ ) ) , int( round( ( part.displacement_.y_ - domainBlocks_[ 0 ].position_.y_ ) / dx_  ) ) , int( round( ( part.displacement_.z_ - domainBlocks_[ 0 ].position_.z_ ) / dx_ ) ) } ) ){
		domainBlocks_[ 0 ].putParticle( part );
	}

}

void MPISimpleDomain::putParticleInBuffer( Particle part, Int3 coord ){
	containsParticles_ = 1;
	//if( ( topology_->getRank() * blockSize_.x_ <= coord.x_ ) && ( topology_->getRank() * blockSize_.x_ + blockSize_.x_ > coord.x_) ){
	if( domainBlocks_[ 0 ].isInBounds( global2Local( coord ) ) ){
		domainBlocks_[ 0 ].putParticleInBuffer( part, global2Local( coord ) );
	}
}

void MPISimpleDomain::putParticleInBuffer( Particle part ){
	containsParticles_ = 1;
	//if( ( domainBlocks_[ 0 ].getPosition().x_ <= part.displacement_.x_ ) && ( domainBlocks_[ 0 ].getPosition().x_ + domainBlocks_[ 0 ].getSize().x_ * dx_ > part.displacement_.x_) ){
	if( domainBlocks_[ 0 ].isInBounds( { int( round( ( part.displacement_.x_ - domainBlocks_[ 0 ].position_.x_ ) / dx_ ) ) , int( round( ( part.displacement_.y_ - domainBlocks_[ 0 ].position_.y_ ) / dx_  ) ) , int( round( ( part.displacement_.z_ - domainBlocks_[ 0 ].position_.z_ ) / dx_ ) ) } ) ){
		domainBlocks_[ 0 ].putParticleInBuffer( part );
	}
}

void MPISimpleDomain::init( Double3 pos, Double3 size, double dx ){
	dx_ = dx;
	size_ = { int( size.x_ / dx_ ) , int( size.y_ / dx_ ) , int( size.z_ / dx_ ) };
	position_ = pos;
}


void MPISimpleDomain::eraseAt( Int3 coord ){
	//if( ( topology_->getRank() * blockSize_.x_ <= coord.x_ ) && ( topology_->getRank() * blockSize_.x_ + blockSize_.x_ > coord.x_) ){
	if( domainBlocks_[ 0 ].isInBounds( global2Local( coord ) ) ){
		domainBlocks_[ 0 ].eraseAt( global2Local( coord ) );
	}
}

void MPISimpleDomain::eraseBufferAt( Int3 coord ){
	//if( ( topology_->getRank() * blockSize_.x_ <= coord.x_ ) && ( topology_->getRank() * blockSize_.x_ + blockSize_.x_ > coord.x_) ){
	if( domainBlocks_[ 0 ].isInBounds( global2Local( coord ) ) ){
		domainBlocks_[ 0 ].eraseBufferAt( global2Local( coord ) );
	}
}

void MPISimpleDomain::setTopology( MPITopology* topology ){
	topology_ = topology;

	blockSize_ = { int( ceil( ( size_.x_ / dx_ ) / (double)topology->getSize() ) ), int( size_.y_ / dx_), int( size_.z_ / dx_ ) };

	// in the simple domain, there is only one block per process
	numBlockLocal_ = { 1, 0, 0 };

	domainBlocks_.clear();
	domainBlocks_.push_back( Domain( { position_.x_ + blockSize_.x_ * dx_ * topology->getRank(), position_.y_, position_.z_ }, { blockSize_.x_ * dx_ , blockSize_.y_ * dx_ , blockSize_.z_ * dx_ }, dx_ ) );
}



} // namespace piaf
