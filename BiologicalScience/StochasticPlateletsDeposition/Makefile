SWIG=/users/pkunzli/software/swig/bin/swig
SWIGFLAGS=-python -c++
WGET=wget -q

CC=g++
# -D:
CPPFLAGS=-fPIC -O3 -std=c++11 -DGAUSSIAN

INCLUDEPATH=/opt/python/3.6.5.7/include/python3.6m/
# The line below is system dependent, and so you have to search it based on the machine you are working on
INCLUDEPATH+=-I /opt/python/3.6.5.7/lib/python3.6/site-packages/numpy/core/include/

cpp: _model.so model.py

%.py: %.i
	$(SWIG) $(SWIGFLAGS) -o $@ $<

%_wrap.cpp: %.i numpy.i
	$(SWIG) $(SWIGFLAGS) -o $@ $<

%.o: %.cpp
	$(CC) $(CPPFLAGS) -I $(INCLUDEPATH) -c $< -o $@

_%.so: %.o %_wrap.o
	$(CC) -shared $^ -o $@

%.i:
	$(WGET) "https://raw.githubusercontent.com/numpy/numpy/master/tools/swig/numpy.i"
