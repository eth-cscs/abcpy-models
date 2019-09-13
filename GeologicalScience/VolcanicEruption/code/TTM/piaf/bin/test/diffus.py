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

from numpy import *
import matplotlib.pyplot as plt
import h5py
import sys


def main():
  infile = h5py.File( sys.argv[1], 'r' )
  dataset = array( infile[ '/f0' ] )
  cut = dataset[ :: , :: , 60 ]
  plt.rc('xtick', labelsize=25) 
  plt.rc('ytick', labelsize=25) 
  plt.imshow( cut, origin = 'lower', extent = [-50,50,-50,50] )
  #plt.imshow( cut, origin = 'lower')
  plt.xlabel( '$x$', fontsize=30 )
  plt.ylabel( '$y$', fontsize=30 )
  plt.xlim( [ -25, 25 ] )
  plt.ylim( [ -25, 25 ] )
  plt.grid(ls='solid')

  fig = plt.gcf()
  fig.set_size_inches(12,12)
  fig.savefig("fig1.pdf") 
  
  plt.show()



main()
