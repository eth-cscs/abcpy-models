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

#include "WimPlume.hpp"

WimPlume::WimPlume(piaf::Double3 ventPosition):
    consistentData_(false),
    ventPosition_ (ventPosition),
    radiusLimit_ (-1.0)
{}

WimPlume::~WimPlume(){}

void WimPlume::computeCenterLinePoints(){
    const double pi = boost::math::constants::pi<double>();
    double beta, r, x, y;
    centerLinePoints_.clear();
    for( uint i=0; i<directionTo_.size(); i++ ){
        beta = directionTo_[i];
        centerLinePoints_.push_back( std::vector< piaf::Double3 >() );

        //std::cout << "center line points at t = " << time_[i] << std::endl;

        for( uint j=0; j<x_[i].size(); j++ ){
            r = x_[i][j];
            x = r * cos( (90.0 - beta) * (pi/180.0) );
            y = r * sin( (90.0 - beta) * (pi/180.0) );
            centerLinePoints_.back().push_back( { x + ventPosition_.x_, y+ventPosition_.y_, z_[i][j]+ventPosition_.z_ } );
            //centerLinePoints_.back().push_back( { x + ventPosition_.x_, y+ventPosition_.y_, z_[i][j] } );

            //std::cout << "("<<x + ventPosition_.x_<<","<< y+ventPosition_.y_<<","<< z_[i][j]+ventPosition_.z_<<")"<<","<<std::endl;

        }
    }
}

void WimPlume::computeUVect(){
    const double pi = boost::math::constants::pi<double>();
    double alpha, beta, norm, x, y, z;
    // beta : rotation around z
    // alpha : rotation around y
    uVect_.clear();
    for( uint i=0; i<directionTo_.size(); i++ ){
        beta = (90.0 - directionTo_[i]) * (pi/180.0);
        uVect_.push_back( std::vector< piaf::Double3 >() );
        for( uint j=0; j<u_[i].size(); j++ ){
            norm = u_[i][j];
            alpha = angle_[i][j];
            x = norm * cos(alpha);
            z = norm * sin(alpha);
            y = x * sin(beta);
            x = x * cos(beta);
            uVect_.back().push_back( { x, y, z } );
        }
    }
}

void WimPlume::limitRadius(){
    if(radiusLimit_ != -1.0){
        for(auto& radiusVector : r_){
            for(auto& radius : radiusVector){
                if(radius > radiusLimit_) radius = radiusLimit_;
            }
        }
    }
}

void WimPlume::refreshData(){
    computeCenterLinePoints();
    computeUVect();
    limitRadius();
    consistentData_ = true;
}

std::vector< double >* WimPlume::getTime(){
    return &time_;
}

double WimPlume::getDirectionTo( double time ){
    return directionTo_[ indexGet( time ) ];
}

std::vector< double >* WimPlume::getAngle( double time ){
    return &(angle_[ indexGet( time ) ]);
}

std::vector< double >* WimPlume::getZ( double time ){
    return &(z_[ indexGet( time ) ]);
}

std::vector< double >* WimPlume::getX( double time ){
    return &(x_[ indexGet( time ) ]);
}

double WimPlume::getMs( double time ){
    //return ms_[ indexGet( time ) ];
    return ms_[ indexGet( time ) ]*3.14;
}

std::vector< double >* WimPlume::getR( double time ){
    return &(r_[ indexGet( time ) ]);
}

std::vector< double >* WimPlume::getU( double time ){
    return &(u_[ indexGet( time ) ]);
}

std::vector< piaf::Double3 >* WimPlume::getUVect( double time ){
    if( !consistentData_ ) refreshData();
    return &(uVect_[ indexGet( time ) ]);
}

std::vector< double >* WimPlume::getM( double time ){
    return &(m_[ indexGet( time ) ]);
}

std::vector< piaf::Double3 >* WimPlume::getCenterLinePoints( double time ){
    if( !consistentData_ ) refreshData();
    return &(centerLinePoints_[ indexGet( time ) ]);
}

// double WimPlume::getHt( double time ){
//     return (z_[ indexGet( time ) ]).back();
// }

// piaf::Double3 WimPlume::getVentPosition(){
//     return ventPosition_;
// }

std::tuple<int, bool> WimPlume::indexClosestPoint( piaf::Double3 otherPoint, double time, bool trace ){

    if(trace) std::cout << "searching the closest point of " << otherPoint.x_ << "," << otherPoint.y_ << "," << otherPoint.z_ << " in centerline points" << std::endl;

    if( !consistentData_ ) refreshData();

    //http://comproguide.blogspot.ch/2015/02/finding-maximum-element-in-bitonic-array.html
    int it = indexGet( time );

    std::vector< piaf::Double3 >& points = centerLinePoints_[ it ];

    if(trace){
        std::cout << "centerline points : " << std::endl;
        for( piaf::Double3 p : points ){
            std::cout << "(" << p.x_ << "," << p.y_ << "," << p.z_ << ") , " ;
        }
        std::cout << std::endl;
    }

    int size = points.size();
    int left = 0;
    int right = size-1;
    int middle;

    while( left < right ){ // boundary not reached

        middle = left + (right-left) / 2;

        if(trace) std::cout << "middle : " << middle << ", left : " << left << ", right : " << right << std::endl;

        if( middle-1 >= 0 && middle+1 < size ){

            // if dist with middle point is smaller than with its neighbors, we found the closest point
            if( dist3D( points[ middle ], otherPoint ) < dist3D( points[ middle-1 ], otherPoint ) &&
                    dist3D( points[ middle ], otherPoint ) < dist3D( points[ middle+1 ], otherPoint ) ){
                if( dist3D( points[ middle ], otherPoint ) > r_[it][middle]  ){
                    if(trace) std::cout << "closest point found " << points[ middle ].x_ << "," << points[ middle ].y_ << "," << points[ middle ].z_ << std::endl;
                    return std::make_tuple(middle, false);
                }
                if(trace) std::cout << "closest point found " << points[ middle ].x_ << "," << points[ middle ].y_ << "," << points[ middle ].z_ << std::endl;
                return std::make_tuple(middle, true);
            }

            // search in the left part
            else if( dist3D( points[ middle ], otherPoint ) > dist3D( points[ middle-1 ], otherPoint ) ){
                right = middle-1;
            }

            // search in the right part
            else{
                left = middle+1;
            }

        }

        else if( middle > 0 ){ // reached right end
            if( dist3D( points[ middle ], otherPoint ) <= dist3D( points[ middle-1 ], otherPoint ) ){
                left = middle;
            }
            else{
                right = middle-1;
            }
        }

        else{ // reached left end
            if( dist3D( points[ middle ], otherPoint ) < dist3D( points[ middle+1 ], otherPoint ) ){
                right = middle;
            }
            else{
                left = middle+1;
            }
        }
    }

    // if the closest point is the last point of the plume, we consider the point to be outside
    if( dist3D( points[ left ], otherPoint ) > r_[it][left] || left==int(points.size())-1  ){
        if(trace) std::cout << "closest point found " << points[ left ].x_ << "," << points[ left ].y_ << "," << points[ left ].z_ << std::endl;
        return std::make_tuple(left, false);
    }

    if(trace) std::cout << "closest point found " << points[ left ].x_ << "," << points[ left ].y_ << "," << points[ left ].z_ << std::endl;
    return std::make_tuple(left, true);

}

void WimPlume::addDirectionTo( double directionTo ,double time ){
    consistentData_ = false;
    directionTo_[ indexAdd( time ) ] = directionTo;
}

void WimPlume::addAngle( std::vector< double > angle ,double time ){
    consistentData_ = false;
    angle_[ indexAdd( time ) ] = angle;
}

void WimPlume::addZ( std::vector< double > z ,double time ){
    consistentData_ = false;
    z_[ indexAdd( time ) ] = z;
}

void WimPlume::addX( std::vector< double > x ,double time ){
    consistentData_ = false;
    x_[ indexAdd( time ) ] = x;
}

void WimPlume::addMs( double ms ,double time ){
    consistentData_ = false;
    ms_[ indexAdd( time ) ] = ms;
}

void WimPlume::addR( std::vector< double > r ,double time ){
    consistentData_ = false;
    r_[ indexAdd( time ) ] = r;
}

void WimPlume::addU( std::vector< double > u ,double time ){
    consistentData_ = false;
    u_[ indexAdd( time ) ] = u;
}

void WimPlume::addM( std::vector< double > m ,double time ){
    consistentData_ = false;
    m_[ indexAdd( time ) ] = m;
}

int WimPlume::indexGet( double time ){
    if( time_.size() == 0 ) throw std::runtime_error("trying to access data from uninitialised plume model");

    if( time_[ 0 ] > time ) return 0;
    if( time_.back() < time ) return time_.size()-1;

    for( uint i = 0; i < time_.size(); i++ ){
        if( time_[ i ] > time ) return i-1;
        if( time_[ i ] == time ) return i;
    }
    throw std::runtime_error("index not found");
}

int WimPlume::indexAdd( double time ){

    uint i = 0;
    for( i=0; i<time_.size(); i++ ){
        if( time_[i] == time ) return i;
        if( time_[i] > time ) break;
    }

    // create elements at the end of vectors
    if( time_.size() == 0 || time_[i] < time ){
        time_.push_back( time );
        directionTo_.push_back( 0.0 );
        ms_.push_back( 0.0 );
        angle_.push_back( std::vector< double >() );
        z_.push_back( std::vector< double >() );
        x_.push_back( std::vector< double >() );
        r_.push_back( std::vector< double >() );
        u_.push_back( std::vector< double >() );
        m_.push_back( std::vector< double >() );
        return time_.size()-1;
    }
    // create elements before the indicated one (i)
    else{
        time_.insert( time_.begin()+i, time );
        directionTo_.insert( directionTo_.begin()+i, 0.0 );
        ms_.insert( ms_.begin()+i, 0.0 );
        angle_.insert( angle_.begin()+i, std::vector< double >() );
        z_.insert( z_.begin()+i, std::vector< double >() );
        x_.insert( x_.begin()+i, std::vector< double >() );
        r_.insert( r_.begin()+i, std::vector< double >() );
        u_.insert( u_.begin()+i, std::vector< double >() );
        m_.insert( m_.begin()+i, std::vector< double >() );
        return i;
    }

}

void WimPlume::setRadiusLimit( double rl ){
    consistentData_ = false;
    radiusLimit_ = rl;
}
