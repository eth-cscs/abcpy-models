#!/bin/sh
cd "`dirname $0`"
root_folder="./"
executable="${root_folder}/sample3submodels"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../lib

# add -p option to generate profiling files
mpiexec -np 11 ${executable} -a -m -c coupling.txt 

cd -
