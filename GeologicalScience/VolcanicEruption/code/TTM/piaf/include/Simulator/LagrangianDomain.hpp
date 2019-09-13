#ifndef LAGRANGIANDOMAIN_H_
#define LAGRANGIANDOMAIN_H_

#include <vector>
#include "SpeedFunctor.hpp"
#include "GenericParticleFamily.hpp"
#include "SimulatorTypes.hpp"

namespace piaf{
class LagrangianDomain {
	public:

        LagrangianDomain();

        virtual ~LagrangianDomain();

    protected:

        std::vector<Particle> domain_;

};
}

#endif