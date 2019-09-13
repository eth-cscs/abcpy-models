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


#ifndef VELOCITYFIELD_H_
#define VELOCITYFIELD_H_

#include "SpeedFunctor.hpp"
#include <vector>

namespace piaf{

    class VelocityField{
    public:
        
        VelocityField( std::vector< std::shared_ptr< SpeedFunctor > > speedFunctors, Double3 origin, Double3 size, double dx, std::shared_ptr<GenericParticleFamily> particleFamily );

        std::shared_ptr<Double3> getHdf5Ordering( double t, double dt );
        std::shared_ptr<double> getXHdf5Ordering( double t, double dt );
        std::shared_ptr<double> getYHdf5Ordering( double t, double dt );
        std::shared_ptr<double> getZHdf5Ordering( double t, double dt );
        Double3 getSize();
        Uint3 getDiscreteSize();
        double getDx();
    
    protected:

        Double3 sumSpeedFunctors(Double3 position, double t, double dt);

        std::vector< std::shared_ptr< SpeedFunctor > > speedFunctors_;
        Uint3 discreteSize_;
        Double3 origin_;
        Double3 size_;
        double dx_;
        std::shared_ptr<Double3> velocityField_;
        std::shared_ptr<double> velocityFieldX_;
        std::shared_ptr<double> velocityFieldY_;
        std::shared_ptr<double> velocityFieldZ_;
        std::shared_ptr<GenericParticleFamily> particleFamily_;

    };

}

#endif