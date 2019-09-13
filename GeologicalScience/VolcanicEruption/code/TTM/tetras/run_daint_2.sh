#!/bin/bash -l

#SBATCH --account=cad01
#SBATCH --job-name=volcano
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=35
#SBATCH --partition=debug
#SBATCH --constraint=mc

module load daint-mc
module load cray-python
module load PyExtensions/3.6.5.1-CrayGNU-18.08
module load Boost/1.67.0-CrayGNU-18.08
module load GSL/2.5-CrayGNU-18.08
module load h5py/2.8.0-CrayGNU-18.08-python3-serial
module load netcdf-python/1.4.1-CrayGNU-18.08-python3

srun tetras 158.41431907 34.23246676 738 output.h5
