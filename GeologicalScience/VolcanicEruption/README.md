



**To compile and run the abcpy version :**

Be sure that your working directory is in scratch (cd $SCRATCH)

* Load following modules

module load daint-mc

module load cray-python

module load PyExtensions/3.6.5.1-CrayGNU-18.08 

module load Boost/1.67.0-CrayGNU-18.08

module load GSL/2.5-CrayGNU-18.08

module load h5py/2.8.0-CrayGNU-18.08-python3-serial

module load netcdf-python/1.4.1-CrayGNU-18.08-python3

* inside your python venv install libraries :
`source venv-3.6/bin/activate`
`pip3 install --upgrade --force-reinstall numpy`
`pip3 install tqdm`
`pip3 install torch`
`pip3 install metric_learn`

* Inside muscleHPC directory run :
`make`

* Inside TTM directory, correct TTM/tetras/src/tetras.cpp by providing the correct location of 'eruption-data' folder 

* Inside TTM directory, run :
`make binaries`

* Add following lines to .bashrc file

export LD_LIBRARY_PATH=/users/duttar/abcvolcano/code/TTM/lib/:/users/duttar/abcvolcano/code/muscleHPC/lib/:$LD_LIBRARY_PATH

export CRAYPE_LINK_TYPE=dynamic

* Inside volcanoinference directory

* Correct compile_daint.sh by mentioning path of swig 

* run :
`sh compile_daint.sh`

* Then you should be able to run :
`sbatch run_daint_test.sh`
to just test the model execution and :
`sbatch run_daint_*.sh`
to run the inference scheme with the chosen statistics calculator 
