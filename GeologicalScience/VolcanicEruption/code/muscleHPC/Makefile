##########################################################################
## Makefile.
##
## The present Makefile is a pure configuration file, in which 
## you can select compilation options. Compilation dependencies
## are managed automatically through the Python library SConstruct.
##
## If you don't have Python, or if compilation doesn't work for other
## reasons, consult the tetras user's guide for instructions on manual
## compilation.
# (Taken fron Palabos)
#
# Author: Mohamed Ben Belgacem
# date: 15.11.2015
##########################################################################

# Leading directory of the tetras source code
muscleHPCRoot  = ./
# Set optimization flags on/off
optimize     = true
# Set debug mode and debug flags on/off
debug        = true
# Set profiling flags on/off
profile      = false
# Path to external libraries (other than tetras)
libraryPaths =
# Path to inlude directories (other than tetras)
includePaths =
# Dynamic and static libraries (other than tetras)
libraries    = boost_regex boost_program_options 
 # Compiler to use with MPI parallelism
parallelCXX  = CC
# General compiler flags (e.g. -Wall to turn on all warnings on g++)
compileFlags = -Wall -Wnon-virtual-dtor  -O3 -ffast-math -funroll-loops -std=c++11 -DUSE_PARALLEL_MPI
# General linker flags (don't put library includes into this flag)
linkFlags    = -ffast-math -funroll-loops -std=c++11 
# Compiler flags to use when optimization mode is on
optimFlags   = -O3
# Compiler flags to use when debug mode is on
debugFlags   = -g
# Compiler flags to use when profile mode is on
profileFlags = -pg
# names of directories to execlude fom compilation
excludedDirs= Display2D ConsoleInterface SConstruct
# Whether use static or dynamic generated libraries: true means generte static library and link them statically 
useStaticLink= true

##########################################################################
# All code below this line is just about forwarding the options
# to SConstruct. It is recommended not to modify anything there.
##########################################################################

SCons     = $(muscleHPCRoot)/scons/scons.py -j 8 -f $(muscleHPCRoot)/SConstruct

SConsArgs = muscleHPCRoot=$(muscleHPCRoot) \
            optimize=$(optimize) \
            debug=$(debug) \
            profile=$(profile) \
            parallelCXX=$(parallelCXX) \
            compileFlags="$(compileFlags)" \
            linkFlags="$(linkFlags)" \
            optimFlags="$(optimFlags)" \
            debugFlags="$(debugFlags)" \
            profileFlags="$(profileFlags)" \
            libraryPaths="$(libraryPaths)" \
            includePaths="$(includePaths)" \
            libraries="$(libraries)"\
	    excludedDirs="$(excludedDirs)"\
	    useStaticLink="$(useStaticLink)"

compile:
	python $(SCons) $(SConsArgs)

clean:
	python $(SCons) -c $(SConsArgs)
	/bin/rm -vf `find $(muscleHPCRoot) -name '*~'`
