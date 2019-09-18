/* File: TRKwebpageutils.i */
%module TRKwebpageutils

%{
#define SWIG_FILE_WITH_INIT
#include "TRKwebpageutils.h"
%}

%include "std_vector.i"

%include "TRKwebpageutils.h"
%include "exampleModels.h"
%include "TRK.h"

/* creates DoubleVector object */
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}


std::vector <double> requestHandler(int fType, std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> allparamsguess, int dataSize, int pivotCheck, int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec, int opScale, int findUncertainties);
