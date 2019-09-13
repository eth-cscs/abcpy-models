#wget "https://raw.githubusercontent.com/numpy/numpy/master/tools/swig/numpy.i"
swig -python -c++ -o model_wrap.cpp model.i
#mpiCC -fPIC -c ../src/tetras_model.cpp -o ../src/tetras_model.o
mpiCC -fPIC -O3 -I /usr/include/python3.6/ -I /usr/local/lib/python3.6/dist-packages/mpi4py/include/ -c model_wrap.cpp -o model_wrap.o
# mpiCC -fPIC -shared ../src/tetras_model.o model_wrap.o -o _model.so
mpiCC -fPIC -shared ../src/tetras.o model_wrap.o -o _model.so -L/home/pierre/GIT/mm/TTM/lib/ -lpiaf_mpi -L/home/pierre/GIT/mm/TTM/lib/ -ltetras_mpi -L/home/pierre/GIT/mm/muscleHPC/lib/ -lmusclehpc
