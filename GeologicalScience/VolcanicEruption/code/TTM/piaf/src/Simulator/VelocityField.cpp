/*
Particles in Advection Field (PIAF)
Copyright (C) 2018  University of Geneva, Switzerland

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

#include "../../include/Simulator/VelocityField.hpp"

namespace piaf{

    VelocityField::VelocityField( std::vector< std::shared_ptr< SpeedFunctor > > speedFunctors, Double3 origin, Double3 size, double dx, std::shared_ptr<GenericParticleFamily> particleFamily ):
    speedFunctors_ (speedFunctors),
    discreteSize_ ({uint(ceil(size.x_/dx)), uint(ceil(size.y_/dx)), uint(ceil(size.z_/dx))}),
    origin_ (origin),
    size_ (size),
    dx_ (dx),
    particleFamily_ (particleFamily)
    {}

    Double3 VelocityField::sumSpeedFunctors(Double3 position, double t, double dt){
        Double3 res = {0.0, 0.0, 0.0};
        for(auto f : speedFunctors_){
            if(!f->isDiffusion()){
                res = res + (*f)(position, t, dt, particleFamily_);
            }
        }
        return res;
    }

    std::shared_ptr<Double3> VelocityField::getHdf5Ordering( double t, double dt ){
        if(!velocityField_){
            velocityField_ = std::shared_ptr<Double3>(new Double3[discreteSize_.x_ * discreteSize_.y_ * discreteSize_.z_], [](Double3* p){delete[] p;});
        }
        for(Uint x=0; x<discreteSize_.x_; x++){
            for(Uint y=0; y<discreteSize_.y_; y++){
                for(Uint z=0; z<discreteSize_.z_; z++){
                    velocityField_.get()[x*discreteSize_.z_*discreteSize_.y_ + y*discreteSize_.z_ + z] = sumSpeedFunctors({x*dx_+origin_.x_, y*dx_+origin_.y_, z*dx_+origin_.z_}, t, dt);
                }
            }
        }
        return velocityField_;
    }

    std::shared_ptr<double> VelocityField::getXHdf5Ordering( double t, double dt ){
        if(!velocityFieldX_){
            velocityFieldX_ = std::shared_ptr<double>(new double[discreteSize_.x_ * discreteSize_.y_ * discreteSize_.z_], [](double* p){delete[] p;});
        }
        for(Uint x=0; x<discreteSize_.x_; x++){
            for(Uint y=0; y<discreteSize_.y_; y++){
                for(Uint z=0; z<discreteSize_.z_; z++){
                    velocityFieldX_.get()[x*discreteSize_.z_*discreteSize_.y_ + y*discreteSize_.z_ + z] = sumSpeedFunctors({x*dx_+origin_.x_, y*dx_+origin_.y_, z*dx_+origin_.z_}, t, dt).x_;
                }
            }
        }
        return velocityFieldX_;
    }


    std::shared_ptr<double> VelocityField::getYHdf5Ordering( double t, double dt ){
        if(!velocityFieldY_){
            velocityFieldY_ = std::shared_ptr<double>(new double[discreteSize_.x_ * discreteSize_.y_ * discreteSize_.z_], [](double* p){delete[] p;});
        }
        for(Uint x=0; x<discreteSize_.x_; x++){
            for(Uint y=0; y<discreteSize_.y_; y++){
                for(Uint z=0; z<discreteSize_.z_; z++){
                    velocityFieldY_.get()[x*discreteSize_.z_*discreteSize_.y_ + y*discreteSize_.z_ + z] = sumSpeedFunctors({x*dx_+origin_.x_, y*dx_+origin_.y_, z*dx_+origin_.z_}, t, dt).y_;
                }
            }
        }
        return velocityFieldY_;
    }


    std::shared_ptr<double> VelocityField::getZHdf5Ordering( double t, double dt ){
        if(!velocityFieldZ_){
            velocityFieldZ_ = std::shared_ptr<double>(new double[discreteSize_.x_ * discreteSize_.y_ * discreteSize_.z_], [](double* p){delete[] p;});
        }
        for(Uint x=0; x<discreteSize_.x_; x++){
            for(Uint y=0; y<discreteSize_.y_; y++){
                for(Uint z=0; z<discreteSize_.z_; z++){
                    velocityFieldZ_.get()[x*discreteSize_.z_*discreteSize_.y_ + y*discreteSize_.z_ + z] = sumSpeedFunctors({x*dx_+origin_.x_, y*dx_+origin_.y_, z*dx_+origin_.z_}, t, dt).z_;
                }
            }
        }
        return velocityFieldZ_;
    }

    Double3 VelocityField::getSize(){
        return size_;
    }

    Uint3 VelocityField::getDiscreteSize(){
        return discreteSize_;
    }

    double VelocityField::getDx(){
        return dx_;
    }

}