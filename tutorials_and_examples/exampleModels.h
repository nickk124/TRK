/*
Trotter Reichart Konz (TRK) Official Codebase
Author: Nick C. Konz
See license at https://github.com/nickk124/TRK
*/

//#pragma once
#include <vector>
#include "TRK.h"


//double pivot; //this is the pivot point(s) that your model uses. It must be defined globally so that TRK can access it.

// LINEAR
double linear(double x, std::vector <double> params);
double dLinear(double x, std::vector <double> params);
double ddLinear(double x, std::vector <double> params);
//double pivotLinear(std::vector <double> params1, std::vector <double> params2);
double linearIntercept(std::vector <double> params);
double linearSlope(std::vector <double> params);

// QUADRATIC
double quadratic(double x, std::vector <double> params);
double dQuadratic(double x, std::vector <double> params);
double ddQuadratic(double x, std::vector <double> params);
//double pivotQuadratic(std::vector <double> params1, std::vector <double> params2);

// CUBIC
double cubic(double x, std::vector <double> params);
double dCubic(double x, std::vector <double> params);
double ddCubic(double x, std::vector <double> params);
//double pivotCubic(std::vector <double> params1, std::vector <double> params2);

// POWER LAW
double powerlaw(double x, std::vector <double> params);
double dPowerlaw(double x, std::vector <double> params);
double ddPowerlaw(double x, std::vector <double> params);
//double pivotPowerLaw(std::vector <double> params1, std::vector <double> params2);
double powerlawIntercept(std::vector <double> params);
double powerlawSlope(std::vector <double> params);

// EXPONENTIAL
double exponential(double x, std::vector <double> params);
double dExponential(double x, std::vector <double> params);
double ddExponential(double x, std::vector <double> params);
//double pivotExponential(std::vector <double> params1, std::vector <double> params2);
double exponentialIntercept(std::vector <double> params);
double exponentialSlope(std::vector <double> params);


// LOGARITHMIC
double logarithmic(double x, std::vector <double> params);
double dLogarithmic(double x, std::vector <double> params);
double ddLogarithmic(double x, std::vector <double> params);
//double pivotLogarithmic(std::vector <double> params1, std::vector <double> params2);
double logarithmicIntercept(std::vector <double> params);
double logarithmicSlope(std::vector <double> params);


// DUST EXTINCTION FITS

//c1 vs c2
double c1c2(double c2, std::vector <double> params);
double dc1c2(double c2, std::vector <double> params);
double ddc1c2(double c2, std::vector <double> params);


//BH vs c2
double bhc2(double c2, std::vector <double> params);
double dbhc2(double c2, std::vector <double> params);
double ddbhc2(double c2, std::vector <double> params);
double bhc2Intercept1(std::vector <double> params);
double bhc2Slope1(std::vector <double> params);
double bhc2Intercept2(std::vector <double> params);
double bhc2Slope2(std::vector <double> params);


//RV vs c2
double rvc2(double c2, std::vector <double> params);
double drvc2(double c2, std::vector <double> params);
double ddrvc2(double c2, std::vector <double> params);
double rvc2Intercept1(std::vector <double> params);
double rvc2Slope1(std::vector <double> params);
double rvc2Intercept2(std::vector <double> params);
double rvc2Slope2(std::vector <double> params);


// x0/gamma
double x0gam(double c2, std::vector <double> params);

// COVID-19 MODELS
static bool covid_fitInLogSpace = false;

// 1-segment
double covid19_1BLPO(double t, std::vector <double> params);

// broken-linear (power law form)
double covid19_BL(double t, std::vector <double> params);
double covid19_BL_Intercept1(std::vector <double> params);
double covid19_BL_Slope1(std::vector <double> params);
double covid19_BL_Intercept2(std::vector <double> params);
double covid19_BL_Slope2(std::vector <double> params);
double covid19_BL_Intercept3(std::vector <double> params);
double covid19_BL_Slope3(std::vector <double> params);
double covid19_BL_Intercept4(std::vector <double> params);
double covid19_BL_Slope4(std::vector <double> params);
double covid19_BL_Intercept5(std::vector <double> params);
double covid19_BL_Slope5(std::vector <double> params);
double covid19_BL_Intercept6(std::vector <double> params);
double covid19_BL_Slope6(std::vector <double> params);

// with fixed s
double covid19_BL_fixed(double t, std::vector <double> params);

// with fixed s and 2nd line
double covid19_BL_line2fixed(double t, std::vector <double> params);

// piecewise
// piecewise
double covid19_PW(double t, std::vector <double> params);
double covid19_PW_Intercept1(std::vector <double> params);
double covid19_PW_Slope1(std::vector <double> params);
double covid19_PW_Intercept2(std::vector <double> params);
double covid19_PW_Slope2(std::vector <double> params);


// oscillating model to account for weekend variability
double covid19_BL_oscil(double t, std::vector <double> params);
double covid19_BL_oscil_fixed(double t, std::vector <double> params);
double covid19_BL_oscil_fixed_s(double t, std::vector <double> params);
double covid19_BL_oscil_split(double t, std::vector <double> params);
double covid19_BL_oscil_split_fixed(double t, std::vector <double> params);
double covid19_BL_oscil_split_all_fixed(double t, std::vector <double> params);

double covid19_polynomial(double t, double t0, std::vector <double> p);
double covid19_BL_oscil_poly(double t, std::vector <double> params);
double covid19_BL_oscil_poly_fixed(double t, std::vector <double> params);
double covid19_BL_oscil_poly_fixedoscil(double t, std::vector <double> params);

// 2-broken line model (most recent)
double covid19_2BLPO(double t, std::vector <double> params);
// remove weekly cycle
double covid19_2BLPO_fixed(double t, std::vector <double> params);

// 3 Broken Line model with Polynomial Oscillation frequency
double covid19_3BLPO(double t, std::vector <double> params);
double covid19_3BLPO_fixed(double t, std::vector <double> params);

// Same, but with 4 lines
// 4-broken line model
double covid19_4BLPO(double t, std::vector <double> params);
double covid19_4BLPO_fixed(double t, std::vector <double> params);

// 5+ lines...we need to go deeper
double covid19_5BLPO(double t, std::vector <double> params);
double covid19_6BLPO(double t, std::vector <double> params);
