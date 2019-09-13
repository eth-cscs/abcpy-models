import h5py
from numpy import *
import scipy.ndimage as ndimage
from smoothing import *

# read physical parameters file
# UTM x , UTM y, distance vent, altitude, thickness, mass / area, ( weight pct size -7 ( id 0 ) to size 11 ( id 18 ) )
def readPhysicalParam( filename ):
	measures = { }
	physicalParamFile = open( filename )
	for line in physicalParamFile:
		tokens = line.split( ";" )
		if( tokens[ 0 ][ 0 : 2 ] == "PL" ):
			data = [ float( tokens[ 3 ] ), float( tokens[ 4 ] ), float( tokens[ 5 ] ), float( tokens[ 6 ] ) ]
			if( tokens[ 7 ] != "" ):
				data = data + [ float( tokens[ 7 ] ), float( tokens[ 8 ] ) ]
				if( tokens[ 9 ] != "" ):
					for i in range( 9, 28 ):
						data = data + [ float( tokens[ i ] ) ]
			measures[ tokens[ 0 ] ] = data
	return measures

def getTotalDepositPoints( measures ):
	depositPoints = []
	for elem in measures:
		if( len( measures[ elem ] ) >= 6 ):
			position = ( measures[ elem ][ 0 ], measures[ elem ][ 1 ] )
			measuredMassarea = measures[ elem ][ 5 ]
			depositPoints = depositPoints + [ ( position, measuredMassarea ) ]
	return depositPoints



def massareaPoints( filename, measures, particleId, sigma_ ):
	depositfile = h5py.File(filename,'r')
	if( particleId == "all" ):

		computedDeposit = zeros( array(depositfile['/f0']).shape )
		for elem in depositfile:
			computedDeposit = add( computedDeposit, array( depositfile[ elem ] ) )
		#computedDeposit = laplacian_smoothing( computedDeposit, 100 )
		computedDeposit = blur_image( computedDeposit, sigma_ )
		#computedDeposit = ndimage.gaussian_filter(computedDeposit, sigma=sigma_, order=0)
		terrainPosition = depositfile['/'].attrs.get("terrain_position")
    
		terrainPosition[0] = 780000.0-30000.0
		terrainPosition[1] = 10005500.0-30000.0

		dx = depositfile['/'].attrs.get("terrain_dx")[ 0 ]   
		#print "position terrain : " + str( terrainPosition )
		points = [] 
		# construct data measured mass area / computed mass area
		for elem in measures:
			if( len( measures[ elem ] ) >= 6 ):
				position = [ measures[ elem ][ 0 ], measures[ elem ][ 1 ] ]
				measuredMassarea = measures[ elem ][ 5 ]
				computedMassarea = computedDeposit[ int(( position[ 0 ] - terrainPosition[ 0 ] ) / dx), int(( position[ 1 ] - terrainPosition[ 1 ] ) / dx) ]
				#print str( measuredMassarea ) + " - " + str( computedMassarea )
				points = points + [ [ measuredMassarea, computedMassarea ] ]

		#print points
          

	else:
		computedDeposit = array(depositfile['/f' + particleId])
		#computedDeposit = laplacian_smoothing( computedDeposit, 10 )
		computedDeposit = blur_image( computedDeposit, sigma_ )
		#computedDeposit = ndimage.gaussian_filter(computedDeposit, sigma=sigma_, order=0)		
		phi = depositfile['/f' + particleId].attrs.get("phi")[ 0 ]
		terrainPosition = depositfile['/'].attrs.get("terrain_position")

		terrainPosition[0] = 780000.0-30000.0
		terrainPosition[1] = 10005500.0-30000.0

		dx = depositfile['/'].attrs.get("terrain_dx")[ 0 ]

		#print "position terrain : " + str( terrainPosition )
		points = []

		# construct data measured mass area / computed mass area
		for elem in measures:
			if( len( measures[ elem ] ) == 25 ):
				weightpct = measures[ elem ][ 13 + int( phi ) ]
				position = [ measures[ elem ][ 0 ], measures[ elem ][ 1 ] ]
				measuredMassarea = ( measures[ elem ][ 5 ] / 100.0 ) * weightpct
				#print "position : " + str( position )
				computedMassarea = computedDeposit[ ( position[ 0 ] - terrainPosition[ 0 ] ) / dx, ( position[ 1 ] - terrainPosition[ 1 ] ) / dx ]
				#print str( measuredMassarea ) + " - " + str( computedMassarea )
				points = points + [ [ measuredMassarea, computedMassarea ] ]

	return points


def readDeposit( filename, sigma_ ):

	depositfile = h5py.File(filename,'r')
	computedDeposit = zeros( array(depositfile['/f0']).shape )
	for elem in depositfile:
		computedDeposit = add( computedDeposit, array( depositfile[ elem ] ) )
	dx = depositfile['/'].attrs.get("terrain_dx")[ 0 ]   
	terrainPosition = depositfile['/'].attrs.get("terrain_position")
	#computedDeposit = ndimage.gaussian_filter(computedDeposit, sigma=sigma_, order=0)
	computedDeposit = blur_image( computedDeposit, sigma_ )
	return computedDeposit, dx, terrainPosition




