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

#include "../../../include/Simulator/MPI/MPI3DBlockCyclicDomain.hpp"

namespace piaf
{

void MPI3DBlockCyclicDomain::commonConstruct(Double3 pos, Double3 size, double dx, MPI3DGridTopology *topology)
{

	gridTopology_ = (MPI3DGridTopology *)topology_;

	//        std::cout << "adjustmentFactor_ : " << adjustmentFactor_ << std::endl;
	//        std::cout << "size_.x_ : " << size_.x_  << ", gridTopology_->getXSize() : " << gridTopology_->getXSize() << std::endl;
	//        std::cout << "size_.y_ : " << size_.x_  << ", gridTopology_->getYSize() : " << gridTopology_->getYSize() << std::endl;
	//        std::cout << "size_.z_ : " << size_.x_  << ", gridTopology_->getZSize() : " << gridTopology_->getZSize() << std::endl;

	// if the simulation is 2D ( z size = 1 ), use a 2D processor topology
	/*! \todo faire une version "generique" de cette adaptation de topologie (pour toute taille et tout axe, pour le moment restreint a simulation 2D avec axe y = 1) */
	if (size_.y_ == 1)
		gridTopology_->to2DTopology();

	blockSize_ = {size_.x_ / (gridTopology_->getXSize() * adjustmentFactor_), 
	size_.y_ / (gridTopology_->getYSize() * adjustmentFactor_), 
	size_.z_ / (gridTopology_->getZSize() * adjustmentFactor_)};

	if (blockSize_.x_ == 0)
		blockSize_.x_ = 1;
	if (blockSize_.y_ == 0)
		blockSize_.y_ = 1;
	if (blockSize_.z_ == 0)
		blockSize_.z_ = 1;

	//blockSize_ = { 10, 10, 10 };

	numBlockLocal_.x_ = int(ceil(double(size_.x_) / double(gridTopology_->getXSize() * blockSize_.x_)));
	numBlockLocal_.y_ = int(ceil(double(size_.y_) / double(gridTopology_->getYSize() * blockSize_.y_)));
	numBlockLocal_.z_ = int(ceil(double(size_.z_) / double(gridTopology_->getZSize() * blockSize_.z_)));

	if (topology_->getRank() == 0)
	{
		std::cout << "blockSize_.x_ : " << blockSize_.x_ << std::endl;
		std::cout << "blockSize_.y_ : " << blockSize_.y_ << std::endl;
		std::cout << "blockSize_.z_ : " << blockSize_.z_ << std::endl;

		std::cout << "numBlockLocal_.x_ : " << numBlockLocal_.x_ << std::endl;
		std::cout << "numBlockLocal_.y_ : " << numBlockLocal_.y_ << std::endl;
		std::cout << "numBlockLocal_.z_ : " << numBlockLocal_.z_ << std::endl;
	}

	//domainBlocks_ = std::vector< Domain >;//( new Domain[ numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ ] );
	domainBlocks_.resize(numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_);

	for (int z = 0; z < numBlockLocal_.z_; z++)
	{
		for (int y = 0; y < numBlockLocal_.y_; y++)
		{
			for (int x = 0; x < numBlockLocal_.x_; x++)
			{
				domainBlocks_[blockLocalIndex3d({x, y, z})].init(blockLocalAddr2AbsolutePosition({x, y, z}), {blockSize_.x_ * dx_, blockSize_.y_ * dx_, blockSize_.z_ * dx_}, dx_);
			}
		}
	}

	//particlesPackets_.resize( numBlockLocal_.x_ * gridTopology_->getXSize() * numBlockLocal_.y_ * gridTopology_->getYSize() * numBlockLocal_.z_ * gridTopology_->getZSize() );

	particlePackets_.resize(gridTopology_->getXSize() * gridTopology_->getYSize() * gridTopology_->getZSize());
	//particlePackets_.resize( 27 );

	// we have to send and receive to / from 27 neighbors ( including ourselve... would be better to not communicate with ourselve... )
	incomingBuffers_.resize(27);
	//sendBuffers_ = std::vector< std::ostringstream >( 27 );
	sendBuffers_.resize(27);
	recvSizes_.resize(27);
	sendSizes_.resize(27);
	requestsSend_.resize(27);
	statusesSend_.resize(27);
	requestsRecv_.resize(27);
	statusesRecv_.resize(27);
}

MPI3DBlockCyclicDomain::MPI3DBlockCyclicDomain(Double3 pos, Double3 size, double dx, MPI3DGridTopology *topology) : MPIAbstractEulerianDomain(pos, size, dx, topology),
																													adjustmentFactor_(5)
{
	commonConstruct(pos, size, dx, topology);
}

MPI3DBlockCyclicDomain::MPI3DBlockCyclicDomain(Double3 pos, Double3 size, double dx, MPI3DGridTopology *topology, int nBlock) : MPIAbstractEulerianDomain(pos, size, dx, topology),
																																adjustmentFactor_(nBlock)
{
	commonConstruct(pos, size, dx, topology);
}

MPI3DBlockCyclicDomain::MPI3DBlockCyclicDomain() : MPIAbstractEulerianDomain(),
												   adjustmentFactor_(5) {}

MPI3DBlockCyclicDomain::~MPI3DBlockCyclicDomain()
{
	//delete[] containsParticlesArray_;
}

void MPI3DBlockCyclicDomain::test()
{
	//std::cout << "je suis un domaine mpi" << std::endl;
}

//std::vector< std::vector< Particle > > MPI3DBlockCyclicDomain::getDomainParticles(){
// gather and return the gathered domain ( relevant only if rank == 0 )
//gather();
//return gatheredDomain_.domain_;
//}

//Int3 MPI3DBlockCyclicDomain::getGlobalSize(){
//	return size_;
//}

//Int3 MPI3DBlockCyclicDomain::getBlockSize(){
//	return blockSize_;
//}

int MPI3DBlockCyclicDomain::blockLocalIndex3d(Int3 coord)
{
	return coord.x_ + coord.y_ * numBlockLocal_.x_ + coord.z_ * numBlockLocal_.x_ * numBlockLocal_.y_;
}

int MPI3DBlockCyclicDomain::blockGlobalIndex3d(Int3 coord)
{
	return coord.x_ + coord.y_ * numBlockLocal_.x_ * gridTopology_->getXSize() + coord.z_ * numBlockLocal_.x_ * gridTopology_->getXSize() * numBlockLocal_.y_ * gridTopology_->getYSize();
}

/*************************************************************************************************************************/

Double3 MPI3DBlockCyclicDomain::blockLocalAddr2AbsolutePosition(Int3 addr)
{
	Double3 res;
	res.x_ = (gridTopology_->getXProcessPos() * blockSize_.x_ + gridTopology_->getXSize() * blockSize_.x_ * addr.x_) * dx_ + position_.x_;
	res.y_ = (gridTopology_->getYProcessPos() * blockSize_.y_ + gridTopology_->getYSize() * blockSize_.y_ * addr.y_) * dx_ + position_.y_;
	res.z_ = (gridTopology_->getZProcessPos() * blockSize_.z_ + gridTopology_->getZSize() * blockSize_.z_ * addr.z_) * dx_ + position_.z_;
	return res;
}

Int3 MPI3DBlockCyclicDomain::absolutePosition2BlockGlobalAddr(Double3 coord)
{
	Int3 res;
	/* \todo : expliquer le + 0.05... 4.94 doit donner 4, mais 4.95 doit donner 5, d'ou le + 0.05 et le floor
		* quand taille bloc == 10 ca joue, mais soucis quand taille == 20 ? a voir
		*/
	/*res.x_ = floor( ( ( ( coord.x_ - position_.x_ ) / dx_ ) / blockSize_.x_ ) + 0.05 );
		res.y_ = floor( ( ( ( coord.y_ - position_.y_ ) / dx_ ) / blockSize_.y_ ) + 0.05 );
		res.z_ = floor( ( ( ( coord.z_ - position_.z_ ) / dx_ ) / blockSize_.z_ ) + 0.05 );*/
	/* \todo : la c'est louche... mais ca a l'air de marcher... a voir */
	res.x_ = floor((((coord.x_ - position_.x_) / dx_) / blockSize_.x_) + 5.0 / (blockSize_.x_ * 10.0));
	res.y_ = floor((((coord.y_ - position_.y_) / dx_) / blockSize_.y_) + 5.0 / (blockSize_.y_ * 10.0));
	res.z_ = floor((((coord.z_ - position_.z_) / dx_) / blockSize_.z_) + 5.0 / (blockSize_.z_ * 10.0));

	/*res.x_ = round( ( ( coord.x_ - position_.x_ ) / dx_ ) / blockSize_.x_ );
		res.y_ = round( ( ( coord.y_ - position_.y_ ) / dx_ ) / blockSize_.y_ );
		res.z_ = round( ( ( coord.z_ - position_.z_ ) / dx_ ) / blockSize_.z_ );*/
	/*res.x_ = roundhalfdown( ( ( coord.x_ - position_.x_ ) / dx_ ) / blockSize_.x_ );
		res.y_ = roundhalfdown( ( ( coord.y_ - position_.y_ ) / dx_ ) / blockSize_.y_ );
		res.z_ = roundhalfdown( ( ( coord.z_ - position_.z_ ) / dx_ ) / blockSize_.z_ );*/

	/*std::cout << "MPI3DBlockCyclicDomain::absolutePosition2BlockGlobalAddr( Double3 coord ) - coord : " << coord.x_ << ", " << coord.y_ << ", " << coord.z_ << std::endl;
		std::cout << "blockSize : " << blockSize_.x_ << ", " << blockSize_.y_ << ", " << blockSize_.z_ << std::endl;
		std::cout << "res : " << res.x_ << ", " << res.y_ << ", " << res.z_ << std::endl;*/

	return res;
}

// returns { numBlockLocal_.x_, numBlockLocal_.y_, numBlockLocal_.z_ } if the current process does not hold the global address
Int3 MPI3DBlockCyclicDomain::blockGlobalAddr2BlockLocalAddr(Int3 addr)
{
	Int3 res;
	res = {numBlockLocal_.x_, numBlockLocal_.y_, numBlockLocal_.z_};

	// 	if( ( ( int( floor( float( addr.x_ ) / float( numBlockLocal_.x_ ) ) ) ) % gridTopology_->getXSize()  == gridTopology_->getXProcessPos() ) &&
	// 	( int( ( floor( float( addr.y_ ) / float( numBlockLocal_.y_ ) ) ) ) % gridTopology_->getYSize()  == gridTopology_->getYProcessPos() ) &&
	// 	( int( ( floor( float( addr.z_ ) / float( numBlockLocal_.z_ ) ) ) ) % gridTopology_->getZSize()  == gridTopology_->getZProcessPos() ) ){
	// 	res.x_ = floor( float( addr.x_ ) / float( gridTopology_->getXSize() ) );
	// 	res.y_ = floor( float( addr.y_ ) / float( gridTopology_->getYSize() ) );
	// 	res.z_ = floor( float( addr.z_ ) / float( gridTopology_->getZSize() ) );
	// }

	if ((addr.x_ % gridTopology_->getXSize() == gridTopology_->getXProcessPos()) &&
		(addr.y_ % gridTopology_->getYSize() == gridTopology_->getYProcessPos()) &&
		(addr.z_ % gridTopology_->getZSize() == gridTopology_->getZProcessPos()))
	{
		res.x_ = floor(float(addr.x_) / float(gridTopology_->getXSize()));
		res.y_ = floor(float(addr.y_) / float(gridTopology_->getYSize()));
		res.z_ = floor(float(addr.z_) / float(gridTopology_->getZSize()));
	}

	return res;
}

Int3 MPI3DBlockCyclicDomain::domainGlobalAddr2BlockGlobalAddr(Int3 addr)
{
	Int3 res;
	res.x_ = addr.x_ / blockSize_.x_;
	res.y_ = addr.y_ / blockSize_.y_;
	res.z_ = addr.z_ / blockSize_.z_;
	return res;
}

Int3 MPI3DBlockCyclicDomain::domainGlobalAddr2DomainLocalAddr(Int3 addr)
{
	Int3 res;
	//res = { blockSize_.x_, blockSize_.y_, blockSize_.z_};
	//if( blockGlobalAddr2BlockLocalAddr( domainGlobalAddr2BlockGlobalAddr ( addr ) ).x_ != numBlockLocal_.x_ ){
	res.x_ = addr.x_ % blockSize_.x_;
	res.y_ = addr.y_ % blockSize_.y_;
	res.z_ = addr.z_ % blockSize_.z_;
	//}
	return res;
}

Int3 MPI3DBlockCyclicDomain::absolutePosition2BlockLocalAddr(Double3 coord)
{
	//auto blockGlobalCoord = absolutePosition2BlockGlobalAddr(coord);
	//std::cout << "Inserting in block global coord : " << blockGlobalCoord.x_ << ", " << blockGlobalCoord.y_ << ", " << blockGlobalCoord.z_ << std::endl;
	return blockGlobalAddr2BlockLocalAddr(absolutePosition2BlockGlobalAddr(coord));
}

/*************************************************************************************************************************/

Domain *MPI3DBlockCyclicDomain::getDomainBlockAtLocalIndex3d(Int3 addr)
{
	return &(domainBlocks_[blockLocalIndex3d(addr)]);
}

//void MPI3DBlockCyclicDomain::setTerrain( GridTerrain *terrain ){
//	for( int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++ ) domainBlocks_[ i ].setTerrain( terrain );
//}

void MPI3DBlockCyclicDomain::init(Double3 pos, Double3 size, double dx)
{
	dx_ = dx;
	size_ = {int(size.x_ / dx_), int(size.y_ / dx_), int(size.z_ / dx_)};
	position_ = pos;
}

void MPI3DBlockCyclicDomain::setTopology(MPITopology *topology)
{
	gridTopology_ = (MPI3DGridTopology *)topology;
	topology_ = topology;

	//containsParticlesArray_ = new int[ gridTopology_->getSize() ];

	blockSize_ = {size_.x_ / (gridTopology_->getXSize() * adjustmentFactor_), size_.y_ / (gridTopology_->getYSize() * adjustmentFactor_), size_.z_ / (gridTopology_->getZSize() * adjustmentFactor_)};
	if (blockSize_.x_ == 0)
		blockSize_.x_ = 1;
	if (blockSize_.y_ == 0)
		blockSize_.y_ = 1;
	if (blockSize_.z_ == 0)
		blockSize_.z_ = 1;

	numBlockLocal_.x_ = int(ceil(double(size_.x_) / double(gridTopology_->getXSize() * blockSize_.x_)));
	numBlockLocal_.y_ = int(ceil(double(size_.y_) / double(gridTopology_->getYSize() * blockSize_.y_)));
	numBlockLocal_.z_ = int(ceil(double(size_.z_) / double(gridTopology_->getZSize() * blockSize_.z_)));

	//domainBlocks_ = boost::shared_array< Domain >( new Domain[ numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ ] );
	domainBlocks_.resize(numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_);

	for (int z = 0; z < numBlockLocal_.z_; z++)
	{
		for (int y = 0; y < numBlockLocal_.y_; y++)
		{
			for (int x = 0; x < numBlockLocal_.x_; x++)
			{
				domainBlocks_[blockLocalIndex3d({x, y, z})].init(blockLocalAddr2AbsolutePosition({x, y, z}), {blockSize_.x_ * dx_, blockSize_.y_ * dx_, blockSize_.z_ * dx_}, dx_);
			}
		}
	}
}

//MPITopology* MPI3DBlockCyclicDomain::getTopology(){
//	return gridTopology_;
//}

//bool MPI3DBlockCyclicDomain::isInBounds( Int3 coord ){
//	if( coord.x_ < size_.x_ && coord.y_ < size_.y_ && coord.z_ < size_.z_ && coord.x_ >= 0 && coord.y_ >= 0 && coord.z_ >= 0) return true;
//	else return false;
//}

void MPI3DBlockCyclicDomain::insertParticles(Double3 coord, int incr, int familyId, double scaling)
{
	containsParticles_ = 1;
	int blockAddr = blockLocalIndex3d(absolutePosition2BlockLocalAddr(coord));
	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
		domainBlocks_[blockAddr].insertParticles(coord, incr, familyId, scaling);
}

// void MPI3DBlockCyclicDomain::insertParticles(Double3 coord, int incr, int familyId)
// {
// 	containsParticles_ = 1;
// 	int blockAddr = blockLocalIndex3d(absolutePosition2BlockLocalAddr(coord));
// 	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
// 		domainBlocks_[blockAddr].insertParticles(coord, incr, familyId);
// }

void MPI3DBlockCyclicDomain::insertParticlesInBuffer(Int3 coord, int incr, int familyId)
{
	containsParticles_ = 1;
	int blockAddr = blockLocalIndex3d(blockGlobalAddr2BlockLocalAddr(domainGlobalAddr2BlockGlobalAddr(coord)));
	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
		domainBlocks_[blockAddr].insertParticlesInBuffer(domainGlobalAddr2DomainLocalAddr(coord), incr, familyId);
}

Double3 MPI3DBlockCyclicDomain::getSpeed(Int3 coord, int familyId)
{
	int blockAddr = blockLocalIndex3d(blockGlobalAddr2BlockLocalAddr(domainGlobalAddr2BlockGlobalAddr(coord)));
	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
		return domainBlocks_[blockAddr].getSpeed(domainGlobalAddr2DomainLocalAddr(coord), familyId);
	//log<LOG_DEBUG>(" warning - process %1% located at %2% ; %3% ; %4% called getSpeed for coord %5% ; %6% ; %7") % gridTopology_->getRank() % gridTopology_->getXProcessPos() % gridTopology_->getYProcessPos() % gridTopology_->getZProcessPos() % coord.x_ % coord.y_ % coord.z_ ;
	return {0.0, 0.0, 0.0};
}

void MPI3DBlockCyclicDomain::putParticle(Particle part, Int3 coord)
{
	containsParticles_ = 1;
	int blockAddr = blockLocalIndex3d(blockGlobalAddr2BlockLocalAddr(domainGlobalAddr2BlockGlobalAddr(coord)));
	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
		domainBlocks_[blockAddr].putParticle(part, domainGlobalAddr2DomainLocalAddr(coord));
	//log<LOG_DEBUG>(" warning - process %1% located at %2% ; %3% ; %4% called putParticle for coord %5% ; %6% ; %7") % gridTopology_->getRank() % gridTopology_->getXProcessPos() % gridTopology_->getYProcessPos() % gridTopology_->getZProcessPos() % coord.x_ % coord.y_ % coord.z_ ;
}

void MPI3DBlockCyclicDomain::putParticle(Particle part)
{
	containsParticles_ = 1;
	int blockAddr = blockLocalIndex3d(absolutePosition2BlockLocalAddr(part.displacement_));
	//std::cout << "blockAddr = " << blockAddr << std::endl;
	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
		domainBlocks_[blockAddr].putParticle(part);
	//else{
	//	std::cout << "MPI3DBlockCyclicDomain::putParticle( Particle part ) -  la particule n'a pas pu etre placee" << std::endl;
	//	std::cout << "blockAddr : " << blockAddr << " - numBlockLocal_.x_ : " << numBlockLocal_.x_ << " - numBlockLocal_.y_ : " << numBlockLocal_.y_ << " - numBlockLocal_.z_ : " << numBlockLocal_.z_ << std::endl;
	//}
}

void MPI3DBlockCyclicDomain::putParticleInBuffer(Particle part, Int3 coord)
{
	containsParticles_ = 1;
	int blockAddr;
	blockAddr = blockLocalIndex3d(blockGlobalAddr2BlockLocalAddr(domainGlobalAddr2BlockGlobalAddr(coord)));
	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
		domainBlocks_[blockAddr].putParticleInBuffer(part, domainGlobalAddr2DomainLocalAddr(coord));
}

void MPI3DBlockCyclicDomain::putParticleInBuffer(Particle part)
{
	containsParticles_ = 1;
	int blockAddr = blockLocalIndex3d(absolutePosition2BlockLocalAddr(part.displacement_));
	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
		domainBlocks_[blockAddr].putParticleInBuffer(part);
}

void MPI3DBlockCyclicDomain::eraseAt(Int3 coord)
{
	int blockAddr;
	blockAddr = blockLocalIndex3d(blockGlobalAddr2BlockLocalAddr(domainGlobalAddr2BlockGlobalAddr(coord)));
	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
		domainBlocks_[blockAddr].eraseAt(domainGlobalAddr2DomainLocalAddr(coord));
}

void MPI3DBlockCyclicDomain::eraseBufferAt(Int3 coord)
{
	int blockAddr;
	blockAddr = blockLocalIndex3d(blockGlobalAddr2BlockLocalAddr(domainGlobalAddr2BlockGlobalAddr(coord)));
	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
		domainBlocks_[blockAddr].eraseBufferAt(domainGlobalAddr2DomainLocalAddr(coord));
}

//void MPI3DBlockCyclicDomain::swapBuffer(){
//	for(int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++) domainBlocks_[ i ].swapBuffer();
//}

//void MPI3DBlockCyclicDomain::computeSpeeds( double t, double dt ){
//	for(int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++) domainBlocks_[ i ].computeSpeeds( t, dt );
//}

//void MPI3DBlockCyclicDomain::computeStaticSpeeds(){
//	for(int i = 0; i < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_; i++) domainBlocks_[ i ].computeStaticSpeeds();
//}

bool MPI3DBlockCyclicDomain::isInGroundContactAt(Int3 coord)
{
	int blockAddr;
	blockAddr = blockLocalIndex3d(blockGlobalAddr2BlockLocalAddr(domainGlobalAddr2BlockGlobalAddr(coord)));
	if (blockAddr < numBlockLocal_.x_ * numBlockLocal_.y_ * numBlockLocal_.z_ && blockAddr >= 0)
		return domainBlocks_[blockAddr].isInGroundContactAt(domainGlobalAddr2DomainLocalAddr(coord));
	//log<LOG_DEBUG>(" warning - process %1% located at %2% ; %3% ; %4% called isInGroundContactAt for coord %5% ; %6% ; %7") % gridTopology_->getRank() % gridTopology_->getXProcessPos() % gridTopology_->getYProcessPos() % gridTopology_->getZProcessPos() % coord.x_ % coord.y_ % coord.z_ ;
	return false;
}

void MPI3DBlockCyclicDomain::putParticleInPacket(Particle part, Domain *lastDomainBlock)
{

	// with 3d grid topology and 3d bloc cyclic domain, we know that that the particle go in the direct neighbor processor
	Double3 partPosition;
	double boundXMax;
	double boundYMax;
	double boundZMax;
	double boundXMin;
	double boundYMin;
	double boundZMin;

	Int3 neighborPos;

	partPosition = part.displacement_;

	// compute max and min bound of the domain from which the particle come
	// bond max not included
	boundXMax = lastDomainBlock->position_.x_ + dx_ * lastDomainBlock->size_.x_ - dx_ / 2.0;
	boundYMax = lastDomainBlock->position_.y_ + dx_ * lastDomainBlock->size_.y_ - dx_ / 2.0;
	boundZMax = lastDomainBlock->position_.z_ + dx_ * lastDomainBlock->size_.z_ - dx_ / 2.0;

	// bound min included
	boundXMin = lastDomainBlock->position_.x_ - dx_ / 2.0;
	boundYMin = lastDomainBlock->position_.y_ - dx_ / 2.0;
	boundZMin = lastDomainBlock->position_.z_ - dx_ / 2.0;

	// if the coordinate of the particle is greater than the max bound in a given direction, then the particle must go to the neighbor in this direction
	// else, it means that the neighbor is in the opposit direction or in the same coordinate
	neighborPos.x_ = (partPosition.x_ >= boundXMax) ? (gridTopology_->getXProcessPos() + 1) % gridTopology_->getXSize() : (partPosition.x_ < boundXMin) ? (((int(gridTopology_->getXProcessPos()) - 1 >= 0) ? gridTopology_->getXProcessPos() - 1 : gridTopology_->getXSize() - 1)) : gridTopology_->getXProcessPos();
	neighborPos.y_ = (partPosition.y_ >= boundYMax) ? (gridTopology_->getYProcessPos() + 1) % gridTopology_->getYSize() : (partPosition.y_ < boundYMin) ? (((int(gridTopology_->getYProcessPos()) - 1 >= 0) ? gridTopology_->getYProcessPos() - 1 : gridTopology_->getYSize() - 1)) : gridTopology_->getYProcessPos();
	neighborPos.z_ = (partPosition.z_ >= boundZMax) ? (gridTopology_->getZProcessPos() + 1) % gridTopology_->getZSize() : (partPosition.z_ < boundZMin) ? (((int(gridTopology_->getZProcessPos()) - 1 >= 0) ? gridTopology_->getZProcessPos() - 1 : gridTopology_->getZSize() - 1)) : gridTopology_->getZProcessPos();

	particlePackets_[gridTopology_->gridPos2Rank(neighborPos)].push_back(part);
}

void MPI3DBlockCyclicDomain::exchangeParticles()
{
	requestExchangeParticles();
	completeExchangeParticles();
}


void MPI3DBlockCyclicDomain::requestExchangeParticles()
{

	const int xSize = gridTopology_->getXSize();
	const int ySize = gridTopology_->getYSize();
	const int zSize = gridTopology_->getZSize();
	const int xProcPos = gridTopology_->getXProcessPos();
	const int yProcPos = gridTopology_->getYProcessPos();
	const int zProcPos = gridTopology_->getZProcessPos();
	Int3 destNeighbor, sourceNeighbor;
	Int3 move;

	//std::vector < Particle > receivedPacket;

	//std::vector<char> incomingBuffer;

	//int recvSize, sendSize;
	MPI_Request requestLocal;
	MPI_Status statusLocal;

	int rankSource, rankDest;
	int i = 0;
	int tag_send = 0;
	int tag_rec = 1;
	MpiManager &mpiManager = gridTopology_->getMpiManager();

	for (move.z_ = -1; move.z_ <= 1; move.z_++)
	{
		for (move.y_ = -1; move.y_ <= 1; move.y_++)
		{
			for (move.x_ = -1; move.x_ <= 1; move.x_++)
			{

				destNeighbor.x_ = (xProcPos + move.x_ >= 0) ? (xProcPos + move.x_) % xSize : xSize + xProcPos + move.x_;
				destNeighbor.y_ = (yProcPos + move.y_ >= 0) ? (yProcPos + move.y_) % ySize : ySize + yProcPos + move.y_;
				destNeighbor.z_ = (zProcPos + move.z_ >= 0) ? (zProcPos + move.z_) % zSize : zSize + zProcPos + move.z_;

				sourceNeighbor.x_ = (xProcPos - move.x_ >= 0) ? (xProcPos - move.x_) % xSize : xSize + xProcPos - move.x_;
				sourceNeighbor.y_ = (yProcPos - move.y_ >= 0) ? (yProcPos - move.y_) % ySize : ySize + yProcPos - move.y_;
				sourceNeighbor.z_ = (zProcPos - move.z_ >= 0) ? (zProcPos - move.z_) % zSize : zSize + zProcPos - move.z_;

				rankSource = gridTopology_->gridPos2Rank(sourceNeighbor);
				rankDest = gridTopology_->gridPos2Rank(destNeighbor);

				// begin asynchronous receive size of the data
				mpiManager.iRecv<int>(&recvSizes_[i], 1, rankSource, &requestLocal, tag_send);

				// prepare data to send ( taken from http://stackoverflow.com/questions/2546298/vector-usage-in-mpic )
				std::ostringstream oss;
				{
					boost::archive::binary_oarchive oa(oss);
					oa << BOOST_SERIALIZATION_NVP(particlePackets_[rankDest]);
				}

				particlePackets_[rankDest].clear();

				sendBuffers_[i] = oss.str();
				// send size of data
				sendSizes_[i] = sendBuffers_[i].size();
				mpiManager.send<int>(&sendSizes_[i], 1, rankDest, tag_send);
				// wait for size to be received
				mpiManager.wait(&requestLocal, &statusLocal);
				//mohamed: MPI_Wait( &requestLocal, &statusLocal );

				incomingBuffers_[i].resize(recvSizes_[i]);
				// begin asynchronous receive of data
				mpiManager.iRecv<char>(&(incomingBuffers_[i])[0], recvSizes_[i], rankSource, &requestsRecv_[i], tag_rec);
				// begin asynchronous send of data
				mpiManager.iSend<char>((char *)sendBuffers_[i].c_str(), sendSizes_[i], rankDest, &requestsSend_[i], tag_rec);

				i++;
			}
		}
	}
}

void MPI3DBlockCyclicDomain::completeExchangeParticles()
{

	//const int xSize = gridTopology_->getXSize();
	//const int ySize = gridTopology_->getYSize();
	//const int zSize = gridTopology_->getZSize();
	//const int xProcPos = gridTopology_->getXProcessPos();
	//const int yProcPos = gridTopology_->getYProcessPos();
	//const int zProcPos = gridTopology_->getZProcessPos();
	//Int3 destNeighbor; //sourceNeighbor;
	Int3 move;

	std::vector<Particle> receivedPacket;

	//int rankDest;
	int i = 0;
	MpiManager &mpiManager = gridTopology_->getMpiManager();
	// wait for data to be sent and received
	for (int i = 0; i < 27; i++)
	{
		// wait for data to be received
		mpiManager.wait(&requestsRecv_[i], &statusesRecv_[i]);
		//mohamed: MPI_Wait( &requestsRecv_[ i ], &statusesRecv_[ i ] );
		// wait for data to be sent
		mpiManager.wait(&requestsSend_[i], &statusesSend_[i]);
		//mohamed: MPI_Wait( &requestsSend_[ i ], &statusesSend_[ i ] );
	}

	for (move.z_ = -1; move.z_ <= 1; move.z_++)
	{
		for (move.y_ = -1; move.y_ <= 1; move.y_++)
		{
			for (move.x_ = -1; move.x_ <= 1; move.x_++)
			{
				//destNeighbor.x_ =  ( xProcPos + move.x_ >= 0 ) ? ( xProcPos + move.x_ ) % xSize : xSize + xProcPos + move.x_;
				//destNeighbor.y_ =  ( yProcPos + move.y_ >= 0 ) ? ( yProcPos + move.y_ ) % ySize : ySize + yProcPos + move.y_;
				//destNeighbor.z_ =  ( zProcPos + move.z_ >= 0 ) ? ( zProcPos + move.z_ ) % zSize : zSize + zProcPos + move.z_;

				//rankSource = gridTopology_->gridPos2Rank( sourceNeighbor );
				//rankDest = gridTopology_->gridPos2Rank( destNeighbor );

				// wait for data to be received
				//MPI_Wait( &requestsRecv_[ i ], &statusesRecv_[ i ] );

				// wait for data to be sent
				//MPI_Wait( &requestsSend_[ i ], &statusesSend_[ i ] );

				// unpack received data
				std::istringstream iss(std::string(&(incomingBuffers_[i])[0], incomingBuffers_[i].size()));
				boost::archive::binary_iarchive ia(iss);
				ia >> BOOST_SERIALIZATION_NVP(receivedPacket);

				//particlePackets_[ rankDest ].clear();
				//sendBuffers_[ i ].clear();
				//incomingBuffers_[ i ].clear();

				//if( receivedPacket.size() > 0 ) std::cout << "le processus " << gridTopology_->getRank() << " a recu " << receivedPacket.size() << " particules dans un paquet" << std::endl;

				//for( Particle& part: receivedPacket ) {
				BOOST_FOREACH (Particle &part, receivedPacket)
				{
					//std::cout << "une particule va etre placee dans le domaine" << std::endl;
					putParticle(part);
					//std::cout << "une particule placee dans le domaine" << std::endl;
				}

				receivedPacket.clear();

				i++;
			}
		}
	}

	/*incomingBuffers_ = std::vector< std::vector<char> >();
	sendBuffers_ = std::vector< std::string >();
	recvSizes_ = std::vector< int >();
	sendSizes_ = std::vector< int >();
	requestsSend_ = std::vector < MPI_Request >();
	statusesSend_ = std::vector < MPI_Status >();
	requestsRecv_ = std::vector < MPI_Request >();
	statusesRecv_ = std::vector < MPI_Status >();
	particlePackets_ = std::vector< std::vector < Particle > > ();


	incomingBuffers_.resize( 27 );
	sendBuffers_.resize( 27 );
	recvSizes_.resize( 27 );
	sendSizes_.resize( 27 );
	requestsSend_.resize( 27 );
	statusesSend_.resize( 27 );
	requestsRecv_.resize( 27 );
	statusesRecv_.resize( 27 );
	particlePackets_.resize( 27 );*/
}

} // namespace piaf
