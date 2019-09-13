#wget "https://raw.githubusercontent.com/numpy/numpy/master/tools/swig/numpy.i"
~/software/swig/bin/swig -python -c++ -o model_wrap.cpp model.i
#mpiCC -fPIC -c ../src/tetras_model.cpp -o ../src/tetras_model.o
CC -fPIC -O3 -I /opt/python/3.6.5.1/include/python3.6m/ -I /opt/python/3.6.5.1/lib/python3.6/site-packages/mpi4py/include/ -I /opt/python/3.6.5.1/lib/python3.6/site-packages/numpy/core/include/ -c model_wrap.cpp -o model_wrap.o
# mpiCC -fPIC -shared ../src/tetras_model.o model_wrap.o -o _model.so
CC -fPIC -shared ../TTM/tetras/src/tetras.o model_wrap.o -o _model.so -L../TTM/lib/ -lpiaf_mpi -L../TTM/lib/ -ltetras_mpi -L../muscleHPC/lib/ -lmusclehpc
