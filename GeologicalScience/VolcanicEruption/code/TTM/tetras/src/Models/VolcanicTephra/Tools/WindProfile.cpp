/*
TEphra TRAsport Simulator (TETRAS)
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

#include "WindProfile.hpp"

WindProfile::WindProfile(): linearSamplingStep_(0.0),
samplingConsistent_ (false)
{}

WindProfile::~WindProfile(){}

std::vector< piaf::Double3 >* WindProfile::getLinearSampling( double time, double step ){
  if( (step != linearSamplingStep_) || (!samplingConsistent_) ) computeLinearSampling( step );
  return &(linearSampling_[ indexGet( time ) ]);
}

void WindProfile::computeLinearSampling( double step ){

  const double pi = boost::math::constants::pi<double>();

  linearSamplingStep_ = step;


  for( uint i=0; i<time_.size(); i++ ){

    linearSampling_[i].clear();

    uint indexClosestHeight = 0;
    double currentHeight = 0.0;

    //int c = 0;

  //  std::cout << "height[i] size : "<< height_[i].size() << std::endl;

    while( indexClosestHeight < height_[i].size() ){
      //if( indexClosestHeight < height_[i].size() ){

      //std::cout << "height_[i][indexClosestHeight] = " << height_[i][indexClosestHeight] << std::endl;
      //std::cout << "currentHeight = " << currentHeight << std::endl;

        if( abs(height_[i][indexClosestHeight]-currentHeight) > abs(height_[i][indexClosestHeight+1]-currentHeight) ){
          indexClosestHeight++;


        }
      //}



      double speed = speed_[i][indexClosestHeight];
      double directionTo = directionTo_[i][indexClosestHeight];
      double beta = (90.0 - directionTo) * (pi/180.0);

      linearSampling_[i].push_back ({ speed * cos(beta) , speed * sin(beta) , 0.0 } );

      currentHeight += linearSamplingStep_;
      //std::cout << currentHeight << std::endl;
      //c++;

      if( currentHeight >= 100000.0 ) break;

    }

    //std::cout << "c : " << c << std::endl;

  }

  samplingConsistent_ = true;

}

std::vector< double >* WindProfile::getTime(){
  return &time_;
}
std::vector< double >* WindProfile::getDirectionTo( double time ){
  return &(directionTo_[ indexGet( time ) ]);
}

std::vector< double >* WindProfile::getHeight( double time ){
  return &(height_[ indexGet( time ) ]);
}

void WindProfile::addDirectionTo( std::vector< double > directionTo, double time ){
  samplingConsistent_ = false;
  directionTo_[ indexAdd( time ) ] = directionTo;
}

void WindProfile::addHeight( std::vector< double > height, double time ){
  samplingConsistent_ = false;
  height_[ indexAdd( time ) ] = height;
}

void WindProfile::addSpeed( std::vector< double > speed, double time ){
  samplingConsistent_ = false;
  speed_[ indexAdd( time ) ] = speed;
}

int WindProfile::indexGet( double time ){
  if( time_.size() == 0 ) throw std::runtime_error("trying to access data from uninitialised wind profile");

  if( time_[ 0 ] > time ) return 0;
  if( time_.back() < time ) return time_.size()-1;

  for( uint i = 0; i < time_.size(); i++ ){
    if( time_[ i ] > time ) return i-1;
    if( time_[ i ] == time ) return i;
  }
  throw("index not found");
}

int WindProfile::indexAdd( double time ){

  samplingConsistent_ = false;
  uint i = 0;
  for( i=0; i<time_.size(); i++ ){
    if( time_[i] == time ) return i;
    if( time_[i] > time ) break;
  }

  // create elements at the end of vectors
  if( time_.size() == 0 || time_[i] < time ){
    time_.push_back( time );
    directionTo_.push_back( std::vector< double >() );
    height_.push_back( std::vector< double >() );
    speed_.push_back( std::vector< double >() );
    linearSampling_.push_back( std::vector< piaf::Double3 >() );
    return time_.size()-1;
  }
  // create elements before the indicated one
  else{
    time_.insert( time_.begin()+i, time );
    directionTo_.insert( directionTo_.begin()+i, std::vector< double >() );
    height_.insert( height_.begin()+i, std::vector< double >() );
    speed_.insert( speed_.begin()+i, std::vector< double >() );
    linearSampling_.insert( linearSampling_.begin()+i, std::vector< piaf::Double3 >() );
    return i;
  }

}
