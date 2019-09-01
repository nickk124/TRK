#pragma once
#include <vector>
#include "TRK.h"


//double pivot; //this is the pivot point(s) that your model uses. It must be defined globally so that TRK can access it.

// LINEAR
double linear(double x, std::vector <double> params);
double dLinear(double x, std::vector <double> params);
double ddLinear(double x, std::vector <double> params);
double pivotLinear(std::vector <double> params1, std::vector <double> params2);

// QUADRATIC
double quadratic(double x, std::vector <double> params);
double dQuadratic(double x, std::vector <double> params);
double ddQuadratic(double x, std::vector <double> params);
double pivotQuadratic(std::vector <double> params1, std::vector <double> params2);

// CUBIC
double cubic(double x, std::vector <double> params);
double dCubic(double x, std::vector <double> params);
double ddCubic(double x, std::vector <double> params);
double pivotCubic(std::vector <double> params1, std::vector <double> params2);

// POWER LAW
double powerlaw(double x, std::vector <double> params);
double dPowerlaw(double x, std::vector <double> params);
double ddPowerlaw(double x, std::vector <double> params);
double pivotPowerLaw(std::vector <double> params1, std::vector <double> params2);

// EXPONENTIAL
double exponential(double x, std::vector <double> params);
double dExponential(double x, std::vector <double> params);
double ddExponential(double x, std::vector <double> params);
double pivotExponential(std::vector <double> params1, std::vector <double> params2);


// LOGARITHMIC
double logarithmic(double x, std::vector <double> params);
double dLogarithmic(double x, std::vector <double> params);
double ddLogarithmic(double x, std::vector <double> params);
double pivotLogarithmic(std::vector <double> params1, std::vector <double> params2);