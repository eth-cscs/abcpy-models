%module model
%{
  #define SWIG_FILE_WITH_INIT
  
  #include <iostream>
  #include <fstream>
  #include <cstdlib>
  #include <string>
  #include <cmath>
  #include <sstream>
  #include <random>
  #include <list>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <unistd.h>
  
  extern void model( double *results, unsigned int rsize, unsigned int k,
                     int noAP, int noNAP,
                     double SR_x,
                     double pAd, double pAg, double pT, double pF, double aT, 
                     double v_z_AP, double v_z_NAP,
                     int seed );
%}

%include "numpy.i"

%init %{
  import_array();
%}

%apply (double* ARGOUT_ARRAY1, int DIM1 ) {(double* results, unsigned int rsize)};

extern void model( double *results, unsigned int rsize, unsigned int k,
                   int noAP, int noNAP,
                   double SR_x,
                   double pAd, double pAg, double pT, double pF, double aT, 
                   double v_z_AP, double v_z_NAP,
                   int seed );
