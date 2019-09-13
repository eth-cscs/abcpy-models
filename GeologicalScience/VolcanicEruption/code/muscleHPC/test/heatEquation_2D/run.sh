#!/bin/sh
cd "`dirname $0`"
root_folder="."
executable="${root_folder}/2d_heat_equation"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../lib

# add -p option to generate profiling files
#mpiexec -np 9  valgrind --leak-check=yes --track-origins=yes  ${executable} -a -m -c coupling.txt 
mpiexec -np 9 ${executable} -a -m -c coupling.txt 

cd -
