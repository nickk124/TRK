#include "exampleModels.h"

// model functions given, then their first two derivatives, then the functions needed to compute their linearized intercept(s) and slope(s) for pivot point finding.
// if you have a pivot point, it must be called as TRK::pivot.

// LINEAR
double linear(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 + a1 * (x - TRK::pivot);
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

	return a0 * std::pow((x / std::pow(10.0, TRK::pivot)), a1);
}

double dPowerlaw(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return std::pow(10.0, -1.0*TRK::pivot) * a0 * a1 * std::pow(std::pow(10.0, -TRK::pivot) * x, a1 - 1.0);
}

double ddPowerlaw(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return std::pow(10.0, -2.0*TRK::pivot) * a0 * (a1 - 1.0) * a1 * std::pow(std::pow(10.0, -TRK::pivot) * x, a1 - 2.0);
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

	return a0 * std::pow(10, a1*(x - TRK::pivot));
}

double dExponential(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return std::pow(10, a1*(x - TRK::pivot)) * a0 * a1 * std::log(10.0);
}

double ddExponential(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return std::pow(10, a1*(x - TRK::pivot)) * a0 * std::pow(a1, 2.0) * std::pow(std::log(10.0) , 2.0);
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

	return a0 + a1 * std::log10(x/TRK::pivot);
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

    return bc1 + mc1*(c2 - TRK::pivot);
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

    return -std::log(std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::pivot)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::pivot2)));
}

double dbhc2(double c2, std::vector <double> params) {
    double b1BH = params[0];
    double theta1BH = params[1];

    double b2BH = params[2];
    double theta2BH = params[3];

    double top = -std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::pivot))*std::tan(theta1BH) - std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::pivot2))*std::tan(theta2BH);
    double bottom = std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::pivot)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::pivot2));

    return -top/bottom;
}

double ddbhc2(double c2, std::vector <double> params) {
    double b1BH = params[0];
    double theta1BH = params[1];

    double b2BH = params[2];
    double theta2BH = params[3];

    double top1 = std::pow(-std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::pivot))*std::tan(theta1BH) - std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::pivot2))*std::tan(theta2BH), 2.0);
    double bottom1 = std::pow(std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::pivot)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::pivot2)), 2.0);

    double top2 = std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::pivot))*std::pow(std::tan(theta1BH), 2.0) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::pivot2))*std::pow(std::tan(theta2BH), 2.0);
    double bottom2 = std::exp(-b1BH - std::tan(theta1BH)*(c2 - TRK::pivot)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - TRK::pivot2));

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

    return std::log(std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::pivot)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::pivot2)));
}

double drvc2(double c2, std::vector <double> params) {
    double b1RV = params[0];
    double theta1RV = params[1];

    double b2RV = params[2];
    double theta2RV = params[3];

    double top = std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::pivot))*std::tan(theta1RV) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::pivot2))*std::tan(theta2RV);
    double bottom = std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::pivot)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::pivot2));

    return top / bottom;
}

double ddrvc2(double c2, std::vector <double> params) {
    double b1RV = params[0];
    double theta1RV = params[1];

    double b2RV = params[2];
    double theta2RV = params[3];

    double top1 = std::pow(std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::pivot))*std::tan(theta1RV) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::pivot2))*std::tan(theta2RV), 2.0);
    double bottom1 = std::pow(std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::pivot)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::pivot2)), 2.0);

    double top2 = std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::pivot))*std::pow(std::tan(theta1RV), 2.0) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::pivot2))*std::pow(std::tan(theta2RV), 2.0);
    double bottom2 = std::exp(b1RV + std::tan(theta1RV)*(c2 - TRK::pivot)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - TRK::pivot2));

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
    
    double y1 = a1 + b1*(t - TRK::pivot);
    double y2 = a2 + b2*(t - TRK::pivot2);

//    return std::pow(std::pow(y1,s) + std::pow(y2, s), 1/s);
    
    if (TRK::covid_logModel){
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

// fixed at s=infty
double covid19_BL_fixed(double t, std::vector <double> params) {
    double a1 = params[0];
    double b1 = params[1];

    double a2 = params[2];
    double b2 = params[3];
    
    double y1 = a1 + b1*(t - TRK::pivot);
    double y2 = a2 + b2*(t - TRK::pivot2);

//    return std::pow(std::pow(y1,s) + std::pow(y2, s), 1/s);
    
    double s = TRK::covid_s;
    
    if (TRK::covid_logModel){
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
    
    double y1 = a1 + b1*(t - TRK::pivot);
    double y2 = a2 + b2*(t - TRK::pivot2);

//    return std::pow(std::pow(y1,s) + std::pow(y2, s), 1/s);
    
    if (TRK::covid_logModel){
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
        return b + m1 * (t - TRK::pivot);
    }
    else {
        return b + m1 * (12 - TRK::pivot) + m2 * (t-12);
    }
}

double covid19_PW_Intercept1(std::vector <double> params){
    return params[0];
}

double covid19_PW_Slope1(std::vector <double> params){
    return params[1];
}

double covid19_PW_Intercept2(std::vector <double> params){
    return TRK::covid_y12;
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
    
    double y1 = b1 + m1*(t - TRK::pivot);
    double y2 = b2 + m2*(t - TRK::pivot2);

    double A = params[5];
    double t0 = params[6];
    double a1 = params[7];
    double a2 = params[8];
    double a3 = params[9];
    
    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    
    double h = a1 + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    double p = a1 * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    double g = std::cos(2.0 * PI * p / 7.0);
    double f = 1.0 + A * g * h;
    
    double ret = line + f;
    
    if (TRK::covid_logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
    return ret;
}

double covid19_BL_oscil_fixed(double t, std::vector <double> params) {
    double b1 = params[0];
    double m1 = params[1];

    double b2 = params[2];
    double m2 = params[3];
    
    double y1 = b1 + m1*(t - TRK::pivot);
    double y2 = b2 + m2*(t - TRK::pivot2);
    
    double s = TRK::covid_fixed_params[0];
    double A = TRK::covid_fixed_params[1];
    double t0 = TRK::covid_fixed_params[2];
    double a1 = TRK::covid_fixed_params[3];
    double a2 = TRK::covid_fixed_params[4];
    double a3 = TRK::covid_fixed_params[5];
    
    double line = (1/s) * std::log(std::exp(s*y1) + std::exp(s*y2));
    
    double h = a1 + 2.0 * a2 * std::fmod(t - t0, 7) + 3.0 * a3 * std::pow(std::fmod(t - t0, 7), 2.0);
    double p = a1 * std::fmod(t - t0, 7) + a2 * std::pow(std::fmod(t - t0, 7), 2.0) + a3 * std::pow(std::fmod(t - t0, 7), 3.0);
    double g = std::cos(2.0 * PI * p / 7.0);
    double f = 1.0 + A * g * h;
    
    double ret = line + f;
    
    if (TRK::covid_logModel){
        ret = std::log10(line) + std::log10(f);
    }
    
//    if (isnan(ret)){
//        std::cout << std::endl;
//    }
    
    return ret;
}
