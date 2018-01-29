%module model
%{
  #define SWIG_FILE_WITH_INIT
  
  #include <iostream>
  #include <fstream>
  #include <cstdlib>
  #include <string>
  #include <cmath>
  #include <sstream>

  #include <boost/random.hpp>
  
  extern void model(double* results, unsigned int rsize, unsigned int k, double pAd, double pAg, double pT, double pF, double aT, int seed);
%}

%include "numpy.i"

%init %{
  import_array();
%}

%apply (double* ARGOUT_ARRAY1, int DIM1 ) {(double* results, unsigned int rsize)};

extern void model(double* results, unsigned int rsize, unsigned int k, double pAd, double pAg, double pT, double pF, double aT, int seed);
