/*
Trotter Reichart Konz (TRK) Official Codebase
Author: Nick C. Konz
See license at https://github.com/nickk124/TRK
*/

#include "exampleModels.h"

// model functions given, then their first two derivatives, then the functions needed to compute their linearized intercept(s) and slope(s) for pivot point finding.
// if you have a pivot point, it must be called as TRK::pivot.

// LINEAR
double linear(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

    return a0 + a1 * (x - TRK::CorrelationRemoval::pivots[0]);
}

double dLinear(double x, std::vector <double> params) {
	double a1 = params[1];

	return a1;
}

double ddLinear(double x, std::vector <double> params) {

	return 0;
}

double linearIntercept(std::vector <double> params){
    return params[0];
}

double linearSlope(std::vector <double> params){
    return params[1];
}


// QUADRATIC
double quadratic(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];
	double a2 = params[2];

	return a0 + a1 * (x) + a2 * std::pow((x), 2.0);
}

double dQuadratic(double x, std::vector <double> params) {
	double a1 = params[1];
	double a2 = params[2];

	return a1 + 2.0 * a2 * (x);
}

double ddQuadratic(double x, std::vector <double> params) {
	double a2 = params[2];

	return 2.0*a2;
}



// CUBIC
double cubic(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];
	double a2 = params[2];
	double a3 = params[3];

	return a0 + a1 * (x) + a2 * std::pow((x), 2.0) + a3 * std::pow((x), 3.0);
}

double dCubic(double x, std::vector <double> params) {
	double a1 = params[1];
	double a2 = params[2];
	double a3 = params[3];

	return a1 + 2.0 * a2 * (x) + 3.0 * a3 * std::pow((x), 2.0);
}

double ddCubic(double x, std::vector <double> params) {
	double a2 = params[2];
	double a3 = params[3];

	return 2.0 * a2 + 6.0 * a3 * (x);
}


// POWER LAW
double powerlaw(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * std::pow((x / std::pow(10.0, TRK::CorrelationRemoval::pivots[0])), a1);
}

double dPowerlaw(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return std::pow(10.0, -1.0*TRK::CorrelationRemoval::pivots[0]) * a0 * a1 * std::pow(std::pow(10.0, -TRK::CorrelationRemoval::pivots[0]) * x, a1 - 1.0);
}

double ddPowerlaw(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return std::pow(10.0, -2.0*TRK::CorrelationRemoval::pivots[0]) * a0 * (a1 - 1.0) * a1 * std::pow(std::pow(10.0, -TRK::CorrelationRemoval::pivots[0]) * x, a1 - 2.0);
}


double powerlawIntercept(std::vector <double> params){
    return std::log10(params[0]);
}

double powerlawSlope(std::vector <double> params){
    return params[1];
}


// EXPONENTIAL
double exponential(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * std::pow(10, a1*(x - TRK::CorrelationRemoval::pivots[0]));
}

double dExponential(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return std::pow(10, a1*(x - TRK::CorrelationRemoval::pivots[0])) * a0 * a1 * std::log(10.0);
}

double ddExponential(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return std::pow(10, a1*(x - TRK::CorrelationRemoval::pivots[0])) * a0 * std::pow(a1, 2.0) * std::pow(std::log(10.0) , 2.0);
}

double exponentialIntercept(std::vector <double> params){
    return std::log10(params[0]);
}

double exponentialSlope(std::vector <double> params){
    return params[1];
}


// LOGARITHMIC
double logarithmic(double x, std::vector <double> params) {
	double a0 = params[0];
    double a1 = params[0];

	return a0 + a1 * std::log10(x/TRK::CorrelationRemoval::pivots[0]);
}

double dLogarithmic(double x, std::vector <double> params) {
	double a1 = params[1];

    return a1 / (x*std::log(10.0));
}

double ddLogarithmic(double x, std::vector <double> params) {
	double a1 = params[1];

	return -1.0 * a1 / (std::pow(x, 2.0) * std::log(10.0));
}


double logarithmicIntercept(std::vector <double> params){
    return params[0];
}

double logarithmicSlope(std::vector <double> params){
    return params[1];
}


// DUST EXTINCTION FITS

//c1 vs c2
double c1c2(double c2, std::vector <double> params) {
    double bc1 = params[0];
    double mc1 = params[1];

    return bc1 + mc1*(c2 - TRK::CorrelationRemoval::pivots[0]);
}

double dc1c2(double c2, std::vector <double> params) {
    double mc1 = params[1];

    return mc1;
}

double ddc1c2(double c2, std::vector <double> params) {
    return 0.0;
}


//BH vs c2
double bhc2(double c2, std::vector <double> params) {
    double b1BH = params[0];
    double theta1BH = params[1];

    double b2BH = params[2];
    double theta2BH = params[3];

    return -std::log(std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::CorrelationRemoval::pivots[0])) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::CorrelationRemoval::pivots[1])));
}

double dbhc2(double c2, std::vector <double> params) {
    double b1BH = params[0];
    double theta1BH = params[1];

    double b2BH = params[2];
    double theta2BH = params[3];

    double top = -std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::CorrelationRemoval::pivots[0]))*std::tan(theta1BH) - std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::CorrelationRemoval::pivots[1]))*std::tan(theta2BH);
    double bottom = std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::CorrelationRemoval::pivots[0])) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::CorrelationRemoval::pivots[1]));

    return -top/bottom;
}

double ddbhc2(double c2, std::vector <double> params) {
    double b1BH = params[0];
    double theta1BH = params[1];

    double b2BH = params[2];
    double theta2BH = params[3];

    double top1 = std::pow(-std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::CorrelationRemoval::pivots[0]))*std::tan(theta1BH) - std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::CorrelationRemoval::pivots[1]))*std::tan(theta2BH), 2.0);
    double bottom1 = std::pow(std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::CorrelationRemoval::pivots[0])) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::CorrelationRemoval::pivots[1])), 2.0);

    double top2 = std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::CorrelationRemoval::pivots[0]))*std::pow(std::tan(theta1BH), 2.0) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::CorrelationRemoval::pivots[1]))*std::pow(std::tan(theta2BH), 2.0);
    double bottom2 = std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::CorrelationRemoval::pivots[0])) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::CorrelationRemoval::pivots[1]));

    return top1 / bottom1 - top2 / bottom2;
}

double bhc2Intercept1(std::vector <double> params){
    return params[0];
}

double bhc2Slope1(std::vector <double> params){
    return std::tan(params[1]);
}

double bhc2Intercept2(std::vector <double> params){
    return params[2];
}

double bhc2Slope2(std::vector <double> params){
    return std::tan(params[3]);
}


//RV vs c2
double rvc2(double c2, std::vector <double> params) {
    double b1RV = params[0];
    double theta1RV = params[1];

    double b2RV = params[2];
    double theta2RV = params[3];

    return std::log(std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::CorrelationRemoval::pivots[0])) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::CorrelationRemoval::pivots[1])));
}

double drvc2(double c2, std::vector <double> params) {
    double b1RV = params[0];
    double theta1RV = params[1];

    double b2RV = params[2];
    double theta2RV = params[3];

    double top = std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::CorrelationRemoval::pivots[0]))*std::tan(theta1RV) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::CorrelationRemoval::pivots[1]))*std::tan(theta2RV);
    double bottom = std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::CorrelationRemoval::pivots[0])) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::CorrelationRemoval::pivots[1]));

    return top / bottom;
}

double ddrvc2(double c2, std::vector <double> params) {
    double b1RV = params[0];
    double theta1RV = params[1];

    double b2RV = params[2];
    double theta2RV = params[3];

    double top1 = std::pow(std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::CorrelationRemoval::pivots[0]))*std::tan(theta1RV) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::CorrelationRemoval::pivots[1]))*std::tan(theta2RV), 2.0);
    double bottom1 = std::pow(std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::CorrelationRemoval::pivots[0])) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::CorrelationRemoval::pivots[1])), 2.0);

    double top2 = std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::CorrelationRemoval::pivots[0]))*std::pow(std::tan(theta1RV), 2.0) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::CorrelationRemoval::pivots[1]))*std::pow(std::tan(theta2RV), 2.0);
    double bottom2 = std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::CorrelationRemoval::pivots[0])) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::CorrelationRemoval::pivots[1]));

    return -(top1 / bottom1) + top2 / bottom2;
}

double rvc2Intercept1(std::vector <double> params){
    return params[0];
}

double rvc2Slope1(std::vector <double> params){
    return std::tan(params[1]);
}

double rvc2Intercept2(std::vector <double> params){
    return params[2];
}

double rvc2Slope2(std::vector <double> params){
    return std::tan(params[3]);
}


// COVID-19 MODELS

// broken-linear (power law form)
double covid19_BL(double t, std::vector <double> params) {
    double a1 = params[0];
    double b1 = params[1];

    double a2 = params[2];
    double b2 = params[3];
    
    double s = params[4];
    
    double y1 = a1 + b1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = a2 + b2*(t - TRK::CorrelationRemoval::pivots[1]);

//    return std::pow(std::pow(y1,s) + std::pow(y2, s), 1/s);
    
    if (TRK::COVID19::logModel){
        return std::log10((1/s) * std::log(std::exp(s*y1) + std::exp(s*y2)));
    }
    else {
        return (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    }
}

double covid19_BL_Intercept1(std::vector <double> params){
    return params[0];
}

double covid19_BL_Slope1(std::vector <double> params){
    return params[1];
}

double covid19_BL_Intercept2(std::vector <double> params){
    return params[2];
}

double covid19_BL_Slope2(std::vector <double> params){
    return params[3];
}

double covid19_BL_Intercept3(std::vector <double> params){
    return params[4];
}

double covid19_BL_Slope3(std::vector <double> params){
    return params[5];
}

// fixed at s=infty
double covid19_BL_fixed(double t, std::vector <double> params) {
    double a1 = params[0];
    double b1 = params[1];

    double a2 = params[2];
    double b2 = params[3];
    
    double y1 = a1 + b1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = a2 + b2*(t - TRK::CorrelationRemoval::pivots[1]);

//    return std::pow(std::pow(y1,s) + std::pow(y2, s), 1/s);
    
    double s = TRK::COVID19::s;
    
    if (TRK::COVID19::logModel){
        return std::log10((1/s) * std::log(std::exp(s*y1) + std::exp(s*y2)));
    }
    else {
        return (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    }
}


// 2nd line fixed
double covid19_BL_line2fixed(double t, std::vector <double> params) {
    double a1 = params[0];
    double b1 = params[1];

    double a2 = 4.0;
    double b2 = 0.6;
    
    double s = 1.0;
    
    double y1 = a1 + b1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = a2 + b2*(t - TRK::CorrelationRemoval::pivots[1]);

//    return std::pow(std::pow(y1,s) + std::pow(y2, s), 1/s);
    
    if (TRK::COVID19::logModel){
        return std::log10((1/s) * std::log(std::exp(s*y1) + std::exp(s*y2)));
    }
    else {
        return (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    }
}


// piecewise
double covid19_PW(double t, std::vector <double> params) {
    double b = params[0];
    double m1 = params[1];
    double m2 = params[2];
    if ( t <=12 ){
        return b + m1 * (t - TRK::CorrelationRemoval::pivots[0]);
    }
    else {
        return b + m1 * (12 - TRK::CorrelationRemoval::pivots[0]) + m2 * (t-12);
    }
}

double covid19_PW_Intercept1(std::vector <double> params){
    return params[0];
}

double covid19_PW_Slope1(std::vector <double> params){
    return params[1];
}

double covid19_PW_Intercept2(std::vector <double> params){
    return TRK::COVID19::y12;
}

double covid19_PW_Slope2(std::vector <double> params){
    return params[2];
}

// oscillating model to account for weekend variability
double covid19_BL_oscil(double t, std::vector <double> params) {
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double s = params[4];
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);

    double A = params[5];
    double t0 = params[6];
    double c = params[7];
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
//    double c = 1.0;
//    double a2 = 0.0;
//    double a3 = 0.0;
    
    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    
    double h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    double p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    double g = std::cos(2.0 * PI * p / 7.0);
    double f = 1.0 + A * g * h;
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
    return ret;
}

double covid19_BL_oscil_fixed(double t, std::vector <double> params) {
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);
    
    double s = TRK::COVID19::fixed_params[0];
    double A = TRK::COVID19::fixed_params[1];
    double t0 = TRK::COVID19::fixed_params[2];
    double c = TRK::COVID19::fixed_params[3];
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    
    double h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    double p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    double g = std::cos(2.0 * PI * p / 7.0);
    double f = 1.0 + A * g * h;
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
//    if (isnan(ret)){
//        std::cout << std::endl;
//    }
    
    return ret;
}

double covid19_BL_oscil_fixed_s(double t, std::vector <double> params) {
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);
    
    double s = TRK::COVID19::s;
    
    double A = params[4];
    double t0 = params[5];
    double c = params[6];
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    
    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    
    double h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    double p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    double g = std::cos(2.0 * PI * p / 7.0);
    double f = 1.0 + A * g * h;
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
//    if (isnan(ret)){
//        std::cout << std::endl;
//    }
    
    return ret;
}

double covid19_BL_oscil_split(double t, std::vector <double> params) {
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double s = params[4];
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);

    double A = params[5];
    double t0 = params[6];
    double c = params[7];
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    double A_2 = params[8];
    double t0_2 = params[9];
    double c_2 = params[10];
    double a2_2 = 3.0/7.0 * (1 - c_2);
    double a3_2 = a2_2 * -2.0 / 21.0;
    
    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    
    double h, p, g, f;
    if (t <= TRK::COVID19::t_split){
        h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
        p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
        g = std::cos(2.0 * PI * p / 7.0);
        f = 1.0 + A * g * h;
    } else {
        h = c_2 + 2.0 * a2_2 * std::fmod(t - t0_2, 7) + 3.0 * a3_2 * std::pow(std::fmod(t - t0_2, 7), 2.0);
        p = c_2 * std::fmod(t - t0_2, 7) + a2_2 * std::pow(std::fmod(t - t0_2, 7), 2.0) + a3_2 * std::pow(std::fmod(t - t0_2, 7), 3.0);
        g = std::cos(2.0 * PI * p / 7.0);
        f = 1.0 + A_2 * g * h;
    }
    
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
    return ret;
}

double covid19_BL_oscil_split_all_fixed(double t, std::vector <double> params) {
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double s = TRK::COVID19::s;
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);

    double A = 0.1;
    double t0 = 3.25;
    double c = 1.7;
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    double A_2 = 0.15;
    double t0_2 = 4.25;
    double c_2 = 1.0;
    double a2_2 = 3.0/7.0 * (1 - c_2);
    double a3_2 = a2_2 * -2.0 / 21.0;
    
    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    
    double h, p, g, f;
    if (t <= TRK::COVID19::t_split){
        h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
        p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
        g = std::cos(2.0 * PI * p / 7.0);
        f = 1.0 + A * g * h;
    } else {
        h = c_2 + 2.0 * a2_2 * std::fmod(t - t0_2, 7) + 3.0 * a3_2 * std::pow(std::fmod(t - t0_2, 7), 2.0);
        p = c_2 * std::fmod(t - t0_2, 7) + a2_2 * std::pow(std::fmod(t - t0_2, 7), 2.0) + a3_2 * std::pow(std::fmod(t - t0_2, 7), 3.0);
        g = std::cos(2.0 * PI * p / 7.0);
        f = 1.0 + A_2 * g * h;
    }
    
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
    return ret;
}

double covid19_BL_oscil_split_fixed(double t, std::vector <double> params) {
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];

    double s = TRK::COVID19::s;

    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);

    double A = params[4];
    double t0 = params[5];
    double c = params[6];
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;

    double A_2 = params[7];
    double t0_2 = params[8];
    double c_2 = params[9];
    double a2_2 = 3.0/7.0 * (1 - c_2);
    double a3_2 = a2_2 * -2.0 / 21.0;

    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));

    double h, p, g, f;
    if (t <= TRK::COVID19::t_split){
        h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
        p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
        g = std::cos(2.0 * PI * p / 7.0);
        f = 1.0 + A * g * h;
    } else {
        h = c_2 + 2.0 * a2_2 * std::fmod(t - t0_2, 7) + 3.0 * a3_2 * std::pow(std::fmod(t - t0_2, 7), 2.0);
        p = c_2 * std::fmod(t - t0_2, 7) + a2_2 * std::pow(std::fmod(t - t0_2, 7), 2.0) + a3_2 * std::pow(std::fmod(t - t0_2, 7), 3.0);
        g = std::cos(2.0 * PI * p / 7.0);
        f = 1.0 + A_2 * g * h;
    }


    double ret = line*f;

    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }

    return ret;
}

double covid19_polynomial(double t, double t0, std::vector <double> p){
    
    double P = p[0];
    
    for (int o = 1; o < (int) p.size(); o++){
        P += p[o] * std::pow(t - t0, (double) o);
    }
    
    return P;
}

double covid19_BL_oscil_poly(double t, std::vector <double> params) {
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double s = params[4];
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);

    std::vector <double> A_coefs, t0_coefs, c_coefs;
    int D = (int) (params.size() - 5) / 3; // polynomial order
    for (int d = 0; d < D; d++){
        A_coefs.push_back(params[5 + d]);
        t0_coefs.push_back(params[5 + D + d]);
        c_coefs.push_back(params[5 + 2*D + d]);
    }
    
    double A = covid19_polynomial(t, TRK::COVID19::tmed, A_coefs);
    double t0 = covid19_polynomial(t, TRK::COVID19::tmed, t0_coefs);
    double c = covid19_polynomial(t, TRK::COVID19::tmed, c_coefs);
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    
    double h, p, g, f;
    h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    g = std::cos(2.0 * PI * p / 7.0);
    f = 1.0 + A * g * h;
    
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
    if (std::isnan(ret)){
        double x = 0.0;
        x += 1.0;
    }
    
    return ret;
}

double covid19_BL_oscil_poly_fixed(double t, std::vector <double> params) {
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double s = TRK::COVID19::s;
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);

    std::vector <double> A_coefs, t0_coefs, c_coefs;
    int D = (int) (params.size() - 4) / 3; // polynomial order
    for (int d = 0; d < D; d++){
        A_coefs.push_back(params[4 + d]);
        t0_coefs.push_back(params[4 + D + d]);
        c_coefs.push_back(params[4 + 2*D + d]);
    }
    
    double A = covid19_polynomial(t, TRK::COVID19::tmed, A_coefs);
    double t0 = covid19_polynomial(t, TRK::COVID19::tmed, t0_coefs);
    double c = covid19_polynomial(t, TRK::COVID19::tmed, c_coefs);
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    
    double h, p, g, f;
    h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    g = std::cos(2.0 * PI * p / 7.0);
    f = 1.0 + A * g * h;
    
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
    return ret;
}

double covid19_BL_oscil_poly_fixedoscil(double t, std::vector <double> params) {
    params = concat(params, TRK::COVID19::fixed_params);
    
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double s = params[4];
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);

    std::vector <double> A_coefs, t0_coefs, c_coefs;
    int D = (int) (params.size() - 5) / 3; // polynomial order
    for (int d = 0; d < D; d++){
        A_coefs.push_back(params[5 + d]);
        t0_coefs.push_back(params[5 + D + d]);
        c_coefs.push_back(params[5 + 2*D + d]);
    }
    
    double A = covid19_polynomial(t, TRK::COVID19::tmed, A_coefs);
    double t0 = covid19_polynomial(t, TRK::COVID19::tmed, t0_coefs);
    double c = covid19_polynomial(t, TRK::COVID19::tmed, c_coefs);
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    
    double h, p, g, f;
    h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    g = std::cos(2.0 * PI * p / 7.0);
    f = 1.0 + A * g * h;
    
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
//    if (std::isnan(ret)){
//        printf("ouch\n");
//    }
    
    return ret;
}

// 3 broken line model with polynomial oscillation frequency
double covid19_3BLPO(double t, std::vector <double> params){
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double b3 = params[4];
    double m3 = params[5];
    
    double s12 = params[6];
    double s23 = params[7];
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);
    double y3 = b3 + m3*(t - TRK::CorrelationRemoval::pivots[2]);

    std::vector <double> A_coefs, t0_coefs, c_coefs;
    int D = (int) (params.size() - 8) / 3; // polynomial order
    for (int d = 0; d < D; d++){
        A_coefs.push_back(params[8 + d]);
        t0_coefs.push_back(params[8 + D + d]);
        c_coefs.push_back(params[8 + 2*D + d]);
    }
    
    double A = covid19_polynomial(t, TRK::COVID19::tmed, A_coefs);
    double t0 = covid19_polynomial(t, TRK::COVID19::tmed, t0_coefs);
    double c = covid19_polynomial(t, TRK::COVID19::tmed, c_coefs);
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    double line = (1.0/s23) * std::log(std::exp((s23/s12) * std::log(std::exp(s12 * y1) + std::exp(s12 * y2))) + std::exp(s23 * y3));
    
    double h, p, g, f;
    h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    g = std::cos(2.0 * PI * p / 7.0);
    f = 1.0 + A * g * h;
    
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
    if (std::isnan(ret)){
        double x = 0.0;
        x += 1.0;
    }
    
    return ret;
}

double covid19_3BLPO_fixed(double t, std::vector <double> params){
//    std::vector<double> _2BLparams = slice(params, 0, 6);
//    std::vector <double> smoothingparams = TRK::COVID19::fixed_params;
//    std::vector <double> polyparams = slice(params, 6, 18);
//    params = concat(_2BLparams, concat(smoothingparams, polyparams));
    
    params = concat(TRK::COVID19::fixed_params, params);
    
//    params = concat(params, TRK::COVID19::fixed_params);
    
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double b3 = params[4];
    double m3 = params[5];
    
    double s12 = params[6];
    double s23 = params[7];
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);
    double y3 = b3 + m3*(t - TRK::CorrelationRemoval::pivots[2]);

    std::vector <double> A_coefs, t0_coefs, c_coefs;
    int D = (int) (params.size() - 8) / 3; // polynomial order
    for (int d = 0; d < D; d++){
        A_coefs.push_back(params[8 + d]);
        t0_coefs.push_back(params[8 + D + d]);
        c_coefs.push_back(params[8 + 2*D + d]);
    }
    
    double A = covid19_polynomial(t, TRK::COVID19::tmed, A_coefs);
    double t0 = covid19_polynomial(t, TRK::COVID19::tmed, t0_coefs);
    double c = covid19_polynomial(t, TRK::COVID19::tmed, c_coefs);
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    double line = (1.0/s23) * std::log(std::exp((s23/s12) * std::log(std::exp(s12 * y1) + std::exp(s12 * y2))) + std::exp(s23 * y3));
    
    double h, p, g, f;
    h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    g = std::cos(2.0 * PI * p / 7.0);
    f = 1.0 + A * g * h;
    
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
    if (std::isnan(ret)){
        double x = 0.0;
        x += 1.0;
    }
    
    return ret;
}

// 4-broken line model
double covid19_4BLPO(double t, std::vector <double> params){
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double b3 = params[4];
    double m3 = params[5];
    
    double b4 = params[6];
    double m4 = params[7];
    
    double s12 = params[8];
    double s23 = params[9];
    double s34 = params[10];
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);
    double y3 = b3 + m3*(t - TRK::CorrelationRemoval::pivots[2]);
    double y4 = b4 + m4*(t - TRK::CorrelationRemoval::pivots[3]);

    std::vector <double> A_coefs, t0_coefs, c_coefs;
    int D = (int) (params.size() - 11) / 3; // polynomial order
    for (int d = 0; d < D; d++){
        A_coefs.push_back(params[11 + d]);
        t0_coefs.push_back(params[11 + D + d]);
        c_coefs.push_back(params[11 + 2*D + d]);
    }
    
    double A = covid19_polynomial(t, TRK::COVID19::tmed, A_coefs);
    double t0 = covid19_polynomial(t, TRK::COVID19::tmed, t0_coefs);
    double c = covid19_polynomial(t, TRK::COVID19::tmed, c_coefs);
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    double line = (1.0/s23) * std::log(std::exp((s23/s12)*std::log(exp(s12*y1) + std::exp(s12*y2))) +  std::exp((s23/s34)*std::log(exp(s34*y3) + std::exp(s34*y4))));
    
    double h, p, g, f;
    h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    g = std::cos(2.0 * PI * p / 7.0);
    f = 1.0 + A * g * h;
    
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
    if (std::isnan(ret)){
        double x = 0.0;
        x += 1.0;
    }
    
    return ret;
}

double covid19_4BLPO_fixed(double t, std::vector <double> params){
    
    params = concat(TRK::COVID19::fixed_params, params);
    
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double b3 = params[4];
    double m3 = params[5];
    
    double b4 = params[6];
    double m4 = params[7];
    
    double s12 = params[8];
    double s23 = params[9];
    double s34 = params[10];
    
    double y1 = b1 + m1*(t - TRK::CorrelationRemoval::pivots[0]);
    double y2 = b2 + m2*(t - TRK::CorrelationRemoval::pivots[1]);
    double y3 = b3 + m3*(t - TRK::CorrelationRemoval::pivots[2]);
    double y4 = b4 + m4*(t - TRK::CorrelationRemoval::pivots[3]);

    std::vector <double> A_coefs, t0_coefs, c_coefs;
    int D = (int) (params.size() - 11) / 3; // polynomial order
    for (int d = 0; d < D; d++){
        A_coefs.push_back(params[11 + d]);
        t0_coefs.push_back(params[11 + D + d]);
        c_coefs.push_back(params[11 + 2*D + d]);
    }
    
    double A = covid19_polynomial(t, TRK::COVID19::tmed, A_coefs);
    double t0 = covid19_polynomial(t, TRK::COVID19::tmed, t0_coefs);
    double c = covid19_polynomial(t, TRK::COVID19::tmed, c_coefs);
    double a2 = 3.0/7.0 * (1 - c);
    double a3 = a2 * -2.0 / 21.0;
    
    double line = (1.0/s23) * std::log(std::exp((s23/s12)*std::log(exp(s12*y1) + std::exp(s12*y2))) +  std::exp((s23/s34)*std::log(exp(s34*y3) + std::exp(s34*y4))));
    
    double h, p, g, f;
    h = c + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    p = c * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    g = std::cos(2.0 * PI * p / 7.0);
    f = 1.0 + A * g * h;
    
    
    double ret = line*f;
    
    if (TRK::COVID19::logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
    if (std::isnan(ret)){
        double x = 0.0;
        x += 1.0;
    }
    
    return ret;
}

