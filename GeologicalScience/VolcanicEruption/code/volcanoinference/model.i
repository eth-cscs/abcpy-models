%apply double *OUTPUT { int* success };

%module model
%{
  #define SWIG_FILE_WITH_INIT

  #include "mpi.h"

  extern void model( MPI_Comm communicator, double U0, double L0, int seed, char* outputfile, int len, int* success );
%}

%include mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%include "numpy.i"

%init %{
  import_array();
%}


%apply (char *STRING, int LENGTH) { (char* outputfile, int len) };

extern void model( MPI_Comm communicator, double U0, double L0, int seed, char* outputfile, int len, int* success );
