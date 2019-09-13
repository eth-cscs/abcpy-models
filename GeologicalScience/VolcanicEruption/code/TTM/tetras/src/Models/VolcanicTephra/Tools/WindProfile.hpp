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

#ifndef WINDPROFILE_HPP_
#define WINDPROFILE_HPP_

#include <vector>
#include <stdexcept>
#include "piaf.hpp"

class WindProfile {
public:

  WindProfile();

  ~WindProfile();

  std::vector< double >* getTime();

  std::vector< double >* getDirectionTo( double time );

  std::vector< double >* getHeight( double time );

  std::vector< double >* getSpeed( double time );

  std::vector< piaf::Double3 >* getLinearSampling( double time, double step);

  void addDirectionTo( std::vector< double > directionTo, double time );

  void addHeight( std::vector< double > height, double time );

  void addSpeed( std::vector< double > speed, double time );




private:

  // return index for getting data at specified time
  int indexGet( double time );

  // return index for inserting data for specidied time
  // data structures are extended if needed
  int indexAdd( double time );

  std::vector< double > time_;

  std::vector< std::vector< double > > directionTo_;

  std::vector< std::vector< double > > height_;

  std::vector< std::vector< double > > speed_;

  std::vector< std::vector< piaf::Double3 > > linearSampling_;

  void computeLinearSampling( double step );

  double linearSamplingStep_;

  bool samplingConsistent_;

};

#endif /* WINDPROFILE_HPP_ */
