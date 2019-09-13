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

#import mayavi.mlab as mlab
from numpy import *
#import matplotlib.pyplot as plt
import h5py
import sys

#def gauss( x , scale, loc ):
#    return ( 1 / ( scale * sqrt( 2 * pi ) ) ) * exp( - ( ( x - loc )**2 ) / ( 2 * scale**2 ) )

def gauss2d( x, d, loc, t ):
    return ( 1000000.0 / ( 4.0 * pi * t * sqrt( d * d ) ) ) * exp( -( ( x - loc )**2.0 / (  4.0 * d * t ) ) )

def gauss3d( x, d, loc, t ):
    return ( 1000000.0 / ( pow( 4.0 * pi * t, 3.0/2.0 ) * sqrt( d * d * d ) ) ) * exp( -( ( x - loc )**2.0 / (  4.0 * d * t ) ) )

def main():
    d = 0.5 
    v = 1.0
    t = 10.0
    
    #scale = sqrt( D * 2 * t )
    #scale = sqrt( D * 3/2 * t )
    #scale = sqrt( D * 6 * t )
    #scale = sqrt(8)
    #scale = sqrt( ( 5.65685**2 * 0.05 * 10 ) / 2 ) 
    #scale = sqrt( ( 7.74597**2 * 0.05 * 10 ) / 4) 
    
    loc = v * t
    #loc = 10

    infile = h5py.File( sys.argv[1], 'r' )
    dataset = array( infile[ '/f0' ] )
    
    #yCut = dataset[ 50 , :: , 50 ]
    yCut = dataset[ 60 , :: , 60 ]
    xCut = dataset[ :: , 60 , 60 ]
    zCut = dataset[ 60 , 60 , :: ]
   
    theoricalCut = gauss3d( arange( -50, 51, 1 ), d, loc, t )

    diff = abs( xCut - theoricalCut)
    if( sum( diff ) / sum( theoricalCut ) > 0.05 ):
        return 1
    
    diff = abs( yCut - theoricalCut)
    if( sum( diff ) / sum( theoricalCut ) > 0.05 ):
        return 1

    diff = abs( zCut - theoricalCut)
    if( sum( diff ) / sum( theoricalCut ) > 0.05 ):
        return 1
    
    return 0

sys.exit( main() )
