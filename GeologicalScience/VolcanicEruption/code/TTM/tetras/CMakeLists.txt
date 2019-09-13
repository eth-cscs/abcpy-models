cmake_minimum_required(VERSION 2.6)
project( tetras )

set( CMAKE_INCLUDE_DIRECTORIES_BEFORE ON )


# the program needs at least gcc/g++ 4.4
set( CMAKE_CXX_COMPILER "mpicxx" )
set( CMAKE_C_COMPILER "mpiCC" )

#choose compiling options
set( TRACE "no" )
set( DISPLAY "no" )
set( CLI "no" )

# set here include directories of mpi, boost, hdf5 ( optionally opencv and ncurses ) if the system is not able to find them itself
set( PIAF_INCLUDE_DIR "../piaf/include/" CACHE FILEPATH "include dir for piaf" )
set( BOOST_INCLUDE_DIR " " CACHE FILEPATH "path to the boost include dir" )
set( HDF5_INCLUDE_DIR " " CACHE FILEPATH "path to the hdf5 include dir" )
set( NCURSES_INCLUDE_DIR " " CACHE FILEPATH "path to the ncurses include dir" )
set( OPENCV_INCLUDE_DIR " " CACHE FILEPATH "path to the opencv include dir" )

# set here lib directories of mpi, boost, hdf5 ( optionally opencv and ncurses ) if the system is not able to find them itself
set( PIAF_LIB_DIR "../piaf/lib/" CACHE FILEPATH "library dir for piaf"  )
set( BOOST_LIB_DIR " " CACHE FILEPATH "path to the boost lib dir" )
set( HDF5_LIB_DIR " " CACHE FILEPATH "path to the hdf5 lib dir" )
set( NCURSES_LIB_DIR " " CACHE FILEPATH "path to the ncurses lib dir" )
set( OPENCV_LIB_DIR " " CACHE FILEPATH "path to the opencv lib dir" )

#=====================================================================================================

#for version number
execute_process(COMMAND git describe --long --abbrev=12 OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)
file(WRITE src/version.h "#define VERSION \"${_output}\"")


set( TRACE_FLAG "" )
set( DISPLAY_FLAG "" )
set( CLI_FLAG "" )

set( MUSCLE_FLAG "" )
set( MUSCLE_FILES "" )
set( MUSCLE_LIB_DIR "" )
set( MUSCLE_LIBRARY_NAME "" )
set( MUSCLE_INCLUDE_DIR "" )

# -DUSE_PARALLEL_MPI is required to compile util.h with mpi
add_definitions(-DUSE_PARALLEL_MPI)
#===================================== MUSCLE FLAG COMPILATION =======================================
option(BUILD_MUSCLE "Build C++ part related to MUSCLE" ON) # or OFF
#=====================================================================================================
IF( TRACE MATCHES "yes" )
	set( TRACE_FLAG "-DTRACE " )
ENDIF( TRACE MATCHES "yes" )
IF( DISPLAY MATCHES "yes" )
	set( DISPLAY_FLAG "-DDISPLAY " )
ENDIF( DISPLAY MATCHES "yes" )
IF( CLI MATCHES "yes" )
	set( CLI_FLAG "-DCLI " )
ENDIF( CLI MATCHES "yes" )




# for gcc versions from 4.4 to 4.6, the param -std to use is -std=c++0x, otherwise -std=c++11
# from
# http://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake
if(CMAKE_COMPILER_IS_GNUCXX)
   execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE
GCC_VERSION)
   if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
        message(STATUS "C++11 activated.")
        #add_definitions("-std=c++11")
        set( CPP11FLAG "-std=c++11" )
   elseif(GCC_VERSION VERSION_GREATER 4.3)
        message(WARNING "C++0x activated. If you get any errors update to a
compiler which fully supports C++11")
        #add_definitions("-std=c++0x")
        set( CPP11FLAG "-std=c++0x" )
   else ()
        message(FATAL_ERROR "C++11 needed. Therefore a gcc compiler with a
version higher or equal to 4.4 is needed.")
   endif()
else(CMAKE_COMPILER_IS_GNUCXX)
   #add_definitions("-std=c++0x")
   set( CPP11FLAG "-std=c++0x" )
endif(CMAKE_COMPILER_IS_GNUCXX)

#=====================================================================================================

if(BUILD_MUSCLE) 
	if (NOT DEFINED $ENV{MUSCLE_HOME})
		message ( STATUS "May be you forgot to export 'muscle.profile' to define MUSCLE_HOME" )
	endif() 
	set (MUSCLE_HOME_DIR $ENV{MUSCLE_HOME}  CACHE FILEPATH "Path to muscle home" FORCE)
	set (MUSCLE_INCLUDE_DIR "${MUSCLE_HOME_DIR}/include/muscle2/"CACHE FILEPATH "path to headers file of muscle" FORCE)
	set (MUSCLE_LIB_DIR "${MUSCLE_HOME_DIR}/lib/" CACHE FILEPATH "path to .so libs" FORCE)
	set (MUSCLE_LIBRARY_NAME "muscle2" CACHE STRING "name of the muscle library" FORCE)
	#set( MUSCLE_FLAG "-DMUSCLE" CACHE STRING "Flag to compile or not the code with muscle extention" FORCE )
	add_definitions(-DMUSCLE)
	set( MUSCLE_FILES "../communication/mapper/*;src/Tools/Simulation/distributedsimulation.*"	)
	
endif (BUILD_MUSCLE)
#=====================================================================================================

set( CMAKE_CXX_FLAGS "${CPP11FLAG} -O3 -ffast-math -funroll-loops ${TRACE_FLAG} ${DISPLAY_FLAG} ${CLI_FLAG} ${MUSCLE_FLAG}" )
set( CMAKE_C_FLAGS "-O3 -ffast-math -funroll-loops ${TRACE_FLAG} ${DISPLAY_FLAG} ${CLI_FLAG} ${MUSCLE_FLAG}" )

include_directories( ${PIAF_INCLUDE_DIR} ${BOOST_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${MUSCLE_INCLUDE_DIR})
link_directories( ${PIAF_LIB_DIR} ${BOOST_LIB_DIR} ${HDF5_LIB_DIR} ${MUSCLE_LIB_DIR})


# include display and cli files only if needed
set( LIB_DISPLAY "" )
set( LIB_CLI "" )
set( DISPLAY_FILES "" )
set( CLI_FILES "" )

set( BASE_FILES "src/Models/*;src/Simulator/*;src/Tools/Eruption/*;src/Tools/FileManager/*;src/Tools/FormattedLog/*;src/Tools/Simulation/Simulation.*;../communication/mapper/util.*;src/Tools/Test/*;src/Tools/ToolsTypes.h" )
IF( DISPLAY MATCHES "yes" )
	set( DISPLAY_FILES "src/Tools/Display2D/*" )
	set( LIB_DISPLAY "opencv_core" "opencv_highgui" )
ENDIF( DISPLAY MATCHES "yes" )
IF( CLI MATCHES "yes" )
	set( CLI_FILES "src/Tools/ConsoleInterface/*" )
	set( LIB_CLI "boost_thread" "ncurses" )
ENDIF( CLI MATCHES "yes" )

# include selected files
file( GLOB_RECURSE source_files ${BASE_FILES} ${DISPLAY_FILES} ${CLI_FILES} ${MUSCLE_FILES})


add_executable( tetras ${source_files} "src/tetras.cpp")
add_executable( tetrascontinent ${source_files} "src/tetrascontinent.cpp" )
target_link_libraries( tetrascontinent piaf boost_serialization boost_system hdf5 hdf5_hl ${LIB_DISPLAY} ${LIB_CLI} ${MUSCLE_LIBRARY_NAME})
target_link_libraries( tetras piaf boost_serialization boost_system hdf5 hdf5_hl ${LIB_DISPLAY} ${LIB_CLI} ${MUSCLE_LIBRARY_NAME})
include_directories( ${PIAF_INCLUDE_DIR} ${BOOST_INCLUDE_DIR} ${HDF5_INCLUDE_DIR})



