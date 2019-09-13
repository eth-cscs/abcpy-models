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

#ifndef WIMPLUME_HPP_
#define WIMPLUME_HPP_

#include <vector>
#include <stdexcept>
#include <cmath>
#include "piaf.hpp"

class WimPlume {
public:

  WimPlume( piaf::Double3 ventPosition );

  ~WimPlume();

  std::vector< double >* getTime();

  double getDirectionTo( double time );

  std::vector< double >* getAngle( double time );

  std::vector< double >* getZ( double time );

  std::vector< double >* getX( double time );

  double getMs( double time );

  std::vector< double >* getR( double time );

  std::vector< double >* getU( double time );

  std::vector< piaf::Double3 >* getUVect( double time );

  std::vector< double >* getM( double time );

  std::vector< piaf::Double3 >* getCenterLinePoints( double time );

  // double getHt(double time);

  // piaf::Double3 getVentPosition();

  // return the index of the closest point of the center line of the plume
  // from the given point, first element is the index, second is a boolean telling
  // if the position of otherPoint is in the plume (true if in the plume, false otherwise)
  std::tuple<int, bool> indexClosestPoint( piaf::Double3 otherPoint, double time, bool trace );

  void addDirectionTo( double directionTo ,double time );

  void addAngle( std::vector< double > angle ,double time );

  void addZ( std::vector< double > z ,double time );

  void addX( std::vector< double > x ,double time );

  void addMs( double ms ,double time );

  void addR( std::vector< double > r ,double time );

  void addU( std::vector< double > u ,double time );

  void addM( std::vector< double > m ,double time );

  void setRadiusLimit( double rl );

private:

  // compute coordinates of centerline points for all times from x, z and direction to
  void computeCenterLinePoints();

  // compute u in vector form from direction and speed
  void computeUVect();

  void limitRadius();

  void refreshData();

  // return index for getting data at specified time
  int indexGet( double time );

  // return index for inserting data for specidied time
  // data structures are extended if needed
  int indexAdd( double time );

  std::vector< double > time_;

  bool consistentData_;

  // directions of the plume (from north, in degrees, with anti-trigo direction)
  std::vector< double > directionTo_;

  std::vector< std::vector< double > > angle_;

  std::vector< std::vector< double > > z_;

  std::vector< std::vector< double > > x_;

  std::vector< double > ms_;

  std::vector< std::vector< double > > r_;

  std::vector< std::vector< double > > u_;

  std::vector< std::vector< piaf::Double3 > > uVect_;

  std::vector< std::vector< double > > m_;

  std::vector< std::vector< piaf::Double3 > > centerLinePoints_;

  piaf::Double3 ventPosition_;

  double radiusLimit_;


};

#endif /* WimPlume_HPP_ */
