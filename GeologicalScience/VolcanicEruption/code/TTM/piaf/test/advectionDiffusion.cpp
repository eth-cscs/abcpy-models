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


#include "piaf.hpp"

int main( int argc, char **argv ) {

	const double D = 0.5;
	const double dx = 1.0;
	const double dt = 0.1;
	const int nPart = 1000000;
	const piaf::Double3 speed = { 1.0, 0.0, 0.0 };
	//const piaf::Double3 speed = { 0.0, 0.0, 0.0 };
	const double simTime = 10.0;

	//std::cout << "ur : " << sqrt( ( 4.0 * D ) / dt ) << std::endl;

	std::shared_ptr< piaf::Diffusion > diffusion(new piaf::Diffusion( D ) );
	std::shared_ptr< piaf::ConstantSpeed > constantSpeed( new piaf::ConstantSpeed( speed ) );

	piaf::Domain domain = piaf::Domain( { -50.0, -50.0, -50.0 }, { 101.0, 101.0, 101.0 }, dx );
	piaf::ExactSimulator simulator = piaf::ExactSimulator( &domain, 0.0, dt, NULL, NULL, false );

	domain.addSpeed( diffusion );
	domain.addSpeed( constantSpeed );

	std::shared_ptr< piaf::GenericParticleFamily > pf( new piaf::GenericParticleFamily( 0 ) );
	domain.addParticleFamily( pf );
	domain.insertParticles( { 0.0, 0.0, 0.0 }, nPart, 0 );

	double cpt = 0.0;

	for( int i = 0; i < int( simTime / dt ); i++ ){
		simulator.step();
		cpt += 1.0;
	}

	piaf::FileManager fm;
	fm.write( "test-3d.h5", &domain );

	if( domain.getNParticles() == nPart ) return 0;
	return 1;

}
