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


#include "../../include/Simulator/AbstractSimulator.hpp"

namespace piaf{

    AbstractSimulator::AbstractSimulator( AbstractEulerianDomain *d, double initialTime, double dt, GridTerrain *terrain, ParticleRepository *repository, bool ddt ):
    domain_			( d ),
    currentTime_	( initialTime ),
    dt_				( dt ),
    ddt_			( ddt ),
    terrain_		( terrain ),
    repository_		( repository ),
    partTracker_ (NULL)
    {

        boost::posix_time::ptime time = boost::posix_time::microsec_clock::local_time();
        boost::posix_time::time_duration duration( time.time_of_day() );
        seed_ = duration.total_microseconds();

        initRand(seed_);
        /*!
        * \todo do a better exception handling...
        */
        if( domain_ == NULL ){
            //std::cout << "domain == null, can't simulate" << std::endl;
            exit( 0 );
        }
        if( terrain_ == NULL ){
            //std::cout << " warning : terrain == null " << std::endl;
        }
        if( repository_ == NULL ){
            //std::cout << "warning : repository == null" << std::endl;
        }

        if( terrain_ != NULL ){
            d->setTerrain( terrain_ );
        }
    }




    // Simulator::Simulator(Domain &d, double dt, GridTerrain &terrain):domain_(d),
    // currentTime_(0),
    // dt_(dt),
    // terrain_(terrain),
    // isInGroundContact_(boost::shared_array<bool>(new bool[d.getXSize()*d.getYSize()*d.getZSize()]))
    // {
    // 	for(int z=0;z<d.getZSize();z++){
    // 		for(int y=0;y<d.getYSize();y++){
    // 			for(int x=0;x<d.getXSize();x++){
    // 				//if(y==0)isInGroundContact_[d.index3d({x,y,z})]=false;
    // 				//else isInGroundContact_[d.index3d({x,y,z})]=false;
    // 				isInGroundContact_[d.index3d({x,y,z})]=false;
    // 			}
    // 		}
    // 	}
    // }

    void AbstractSimulator::setParticleTracker( ParticleTracker *partTracker ){
        partTracker_ = partTracker;
    }

    double AbstractSimulator::getDt(){
        return dt_;
    }

    void AbstractSimulator::enableDdt(){
        ddt_ = true;
    }

    void AbstractSimulator::disableDdt(){

        ddt_ = false;
    }

    void AbstractSimulator::setDdtValue( double value ){
        ddtValue_ = value;
    }

    double AbstractSimulator::getDdtValue(){
        return ddtValue_;
    }

    void AbstractSimulator::initRand( int seed ){
        seed_ = seed;
        srand( seed_ );
    }

    int AbstractSimulator::getSeed(){
        return seed_;
    }

    double AbstractSimulator::getTime(){
        return currentTime_;
    }

    AbstractSimulator::~AbstractSimulator() {}

} // namespace piaf
