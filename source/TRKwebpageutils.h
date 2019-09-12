#include "exampleModels.h"

Priors getPriors(int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec, int paramCount);

std::vector <double> requestHandler(int fType, std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> allparamsguess, int dataSize, int pivotCheck, int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec);
