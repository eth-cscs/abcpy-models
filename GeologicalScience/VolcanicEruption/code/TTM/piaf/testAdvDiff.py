"""
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
"""

#!/usr/bin/python

import os, sys

print "Executing simple advection-diffusion simulation... this takes about 30 sec in a core i5 CPU"
success = os.system( "bin/test/advectionDiffusion" )
if( success != 0 ):
    print "Test failed, number of particles varies along the simulation"
    sys.exit( 1 )

print "Checking particles distribution in the domain"
success = os.system( "bin/test/advectionDiffusionValid3d.py test-3d.h5" )
if( success != 0 ):
    print "Test failed, invalid particles distribution"
    sys.exit( 1 )

os.system( "rm test-3d.h5" )
print "Test succeeded"
sys.exit( 0 )
