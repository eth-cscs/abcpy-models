cmake_minimum_required(VERSION 2.6)

project(piaf)

set( CMAKE_INCLUDE_DIRECTORIES_BEFORE ON )

set( CMAKE_CXX_COMPILER "mpicxx" )
set( CMAKE_C_COMPILER "mpicc" )


#exec_program(
#	${CMAKE_CXX_COMPILER}
#	ARGS                    --version
#	OUTPUT_VARIABLE _compiler_output)
#string(REGEX REPLACE ".* ([0-9]\\.[0-9]\\.[0-9]) .*" "\\1" gcc_compiler_version ${_compiler_output})
#message(STATUS "C++ compiler version: ${gcc_compiler_version} [${CMAKE_CXX_COMPILER}]")
#
#if(gcc_compiler_version MATCHES "4\\.[4-6]\\.[0-9]")
#	message("version ok c++0x")
#	set( CPP11FLAG "-std=c++0x" )
#elseif(gcc_compiler_version VERSION_GREATER "4.6.99")
#	message("version ok c++11")
#	set( CPP11FLAG "-std=c++11" )
#else()
#	message(FATAL_ERROR "compilation need at least gcc 4.4")
#endif()

set( CPP11FLAG "-std=c++11" )

#set( CPP11FLAG "-std=c++11" )

# -std=c++11 : use c++11 standard
# -O3 : enable optimizations (and vectorization)
# -ftree-vectorizer-verbose=2 : show informations about vectorization
# -msse2 : allow vectorization on x86_64 platform with sse ( works for AMD or only intel ? possible to use AVX ? and sse or sse2 ? and regarding other CPU types ? )
# -ffast-math : allow vectorization of float reductions
# for gcc versions from 4.4 to 4.6, the param -std to use is -std=c++0x
#set( CMAKE_CXX_FLAGS "-O3 -ftree-vectorizer-verbose=2 -msse2 -ffast-math ${CPP11FLAG}"  )
#set( CMAKE_C_FLAGS "-O3 -ftree-vectorizer-verbose=2 -msse2 -ffast-math" )

set( CMAKE_CXX_FLAGS "-O3 -funroll-loops ${CPP11FLAG} -DPROFILE"  )
set( CMAKE_C_FLAGS "-O3 -funroll-loops -DPROFILE" )

set( BOOST_INCLUDE_DIR " " CACHE FILEPATH "path to the boost include dir" )
set( BOOST_LIB_DIR " " CACHE FILEPATH "path to the boost lib dir" )

#set( MPI_INCLUDE_DIR " " CACHE FILEPATH "path to the boost include dir" )
#set( MPI_LIB_DIR " " CACHE FILEPATH "path to the boost lib dir" )

set( HDF5_INCLUDE_DIR " " CACHE FILEPATH "path to the boost include dir" )
set( HDF5_LIB_DIR " " CACHE FILEPATH "path to the boost lib dir" )

set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY bin )

include_directories(/usr/include/hdf5/serial/ /home/pierre/GIT/mm/muscleHPC/include/ include include/Simulator include/Simulator/MPI include/Tools/FormattedLog ${BOOST_INCLUDE_DIR} ${HDF5_INCLUDE_DIR})
link_directories( ${BOOST_LIB_DIR} ${HDF5_LIB_DIR} /usr/lib/x86_64-linux-gnu/hdf5/serial/ )

file(

        GLOB_RECURSE

        source_files

        src/*

        include/*

)

file(

	GLOB_RECURSE

	source_files_test

	test/advectionDiffusion.cpp

)

file(

        GLOB_RECURSE

        source_files_test2

        test/advectionDirectionalDiffusion.cpp

)


file(

        GLOB_RECURSE

        source_files_benchmark2

        test/benchmark2.cpp

)


file(

  GLOB_RECURSE

  source_files_benchmark

  test/benchmark.cpp

)

file(

  GLOB_RECURSE

  source_files_testeventlogger

  test/testEventLogger.cpp

)

add_library(

        piaf

        STATIC

        ${source_files}

)

target_link_libraries( piaf boost_serialization boost_system hdf5 hdf5_hl )

#add_executable( benchmark ${source_files_benchmark} )
#target_link_libraries( benchmark boost_serialization boost_system piaf hdf5 hdf5_hl )
#set_target_properties( benchmark PROPERTIES RUNTIME_OUTPUT_DIRECTORY bin/test )

#add_executable( benchmark2 ${source_files_benchmark2} )
#target_link_libraries( benchmark2 boost_serialization boost_system piaf hdf5 hdf5_hl )
#set_target_properties( benchmark2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY bin/test )

#add_executable( advectionDiffusion ${source_files_test} )
#target_link_libraries( advectionDiffusion boost_serialization boost_system piaf hdf5 hdf5_hl )
#set_target_properties( advectionDiffusion PROPERTIES RUNTIME_OUTPUT_DIRECTORY bin/test )

#add_executable( advectionDirectionalDiffusion ${source_files_test2} )
#target_link_libraries( advectionDirectionalDiffusion boost_serialization boost_system piaf hdf5 hdf5_hl )
#set_target_properties( advectionDirectionalDiffusion PROPERTIES RUNTIME_OUTPUT_DIRECTORY bin/test )

add_executable( testEventLogger ${source_files_testeventlogger} )
target_link_libraries( testEventLogger boost_serialization boost_system piaf hdf5 hdf5_hl )
set_target_properties( testEventLogger PROPERTIES RUNTIME_OUTPUT_DIRECTORY bin/test )

enable_testing()
add_test( NAME AdvectionDiffusion COMMAND "./testAdvDiff.py"  )
#add_test( NAME testtest COMMAND "/home/pierre/workspace/piaf/trunk/bin/test/testAdvDiff.py" )
