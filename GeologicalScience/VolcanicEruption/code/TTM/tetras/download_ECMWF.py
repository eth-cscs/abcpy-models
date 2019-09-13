#!/usr/bin/env python
import sys
from ecmwfapi import ECMWFDataServer
import calendar
import os
import socket

# date format : yyyymmdd
# bd : begin date, ed : end date, n : north, s : south, w : west, e : east
def parseargs( argv ):
  bdate, edate = "", ""
  north, south, west, east = -1000.0, -1000.0, -1000.0, -1000.0
  uid = ""
  for e in range(len(argv)):
    if argv[e] == "-id":
        uid = sys.argv[e+1]
    elif argv[e] == "-bd":
      bdate = sys.argv[e+1]
    elif argv[e] == "-ed":
      edate = sys.argv[e+1]
    elif argv[e] == "-n":
      north = float(sys.argv[e+1])
    elif argv[e] == "-s":
      south = float(sys.argv[e+1])
    elif argv[e] == "-w":
      west = float(sys.argv[e+1])
    elif argv[e] == "-e":
      east = float(sys.argv[e+1])

  if bdate == "" or edate == "" or north == -1000.0 or south == -1000.0 or west == -1000.0 or east == -1000.0:
    print "Missing argument, usage : python download_ECMWF.py -bd yyyymmdd -ed yyyymmdd -n 0.0 -s 0.0 -w 0.0 -e 0.0"
    sys.exit(0)

  return (uid, bdate, edate, north, south, west, east)

# use of geopotential height or geometrical height ? geometrical height is only available through Grib1 and Grib2 file formats

out_path    = './'

(uid, bdate, edate, north, south, west, east) = parseargs(sys.argv)

server = ECMWFDataServer()

print "######### ERA-interim  #########"
print 'Accessing wind data from ', bdate,' to ',edate,' (YYYYMMDD)'
print "################################"

server.retrieve({
'dataset'   : "interim",
'date'      : "%s/to/%s"%(bdate,edate),
'time'      : "00/03/06/09/12/15/18/21",
'stream'    : "oper",
'levtype'   : "pl",
'levelist'  : "all",
'type'      : "an",
'class'     : "ei",
'grid'      : "0.25/0.25",
# http://apps.ecmwf.int/codes/grib/param-db/
# 129 geopotential m^2*s^-2 / 130 Temperature K / 131 U component of wind m*s^-1 / 132 V component of wind m*s^-1 / 135 vertical velocity Pa*s^-1 / 156 geopotential height gpm ??? ( geopotential/g) / 157 Relative humidity %
'param'     : "129/130/131/132/135/156/157",
# 156 geopotential height doesn't work -> using geopotential divided by 9.80665 (http://apps.ecmwf.int/codes/grib/param-db/?id=156)
'area'      : "%d/%d/%d/%d"%(north, west, south, east),
'format'	  : 'netcdf',
'target'    : out_path+uid+"_"+bdate+"_"+edate+"_"+str(int(north))+"_"+str(int(south))+"_"+str(int(west))+"_"+str(int(east))+".nc"
})
