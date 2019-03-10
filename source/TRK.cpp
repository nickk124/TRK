#include "pch.h"
#include "TRK.h"
#include <iostream>
#include <cmath>

double yC(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * std::sin(a1 * x);
}

double dyC(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * a1 * std::cos(a1 * x);
}

double ddyC(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return -1.0 * a0 * a1 * a1 * std::sin(a1 * x);
}

double newtonRaphson(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess) {
	double x0 = xguess;
	int itercount = 0;

	//initial iteration
	double f = (yC(x0, params) - y_n) * dyC(x0, params) * Sig_xn2 + (x0 - x_n) * Sig_yn2;  //function we're finding zero of
	double df = (std::pow(dyC(x0, params), 2.0) + (yC(x0, params) - y_n)*ddyC(x0, params)) * Sig_xn2 + Sig_yn2; //derivative of above

	double x1 = x0 - f / df;

	double tol = 1e-9;

	while (std::abs(x1 - x0) > tol) {
		x0 = x1;

		f = (yC(x0, params) - y_n) * dyC(x0, params) * Sig_xn2 + (x0 - x_n) * Sig_yn2;  //function we're finding zero of
		df = (std::pow(dyC(x0, params), 2.0) + (yC(x0, params) - y_n)*ddyC(x0, params)) * Sig_xn2 + Sig_yn2; //derivative of above

		x1 = x0 - f / df;
		std::cout << x1 << std::endl;

		itercount += 1;
	}

	std::cout << itercount << " iterations." << std::endl;
	return x1;
}

double twoPointNR(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess)
{
	double tol = 1e-9;
	double xkm1 = xguess;
	//double xk = xguess + xguess / 100.0;
	double xk = xguess + Sig_xn2;

	double xkp1;
	double r;
	double ykm1;
	double yk;
	double dyk;

	int itercount = 0;

	while (std::abs(xk - xkm1) > tol) {
		ykm1 = (yC(xkm1, params) - y_n) * dyC(xkm1, params) * Sig_xn2 + (xkm1 - x_n) * Sig_yn2;
		yk = (yC(xk, params) - y_n) * dyC(xk, params) * Sig_xn2 + (xk - x_n) * Sig_yn2;  //function we're finding zero of
		dyk = (std::pow(dyC(xk, params), 2.0) + (yC(xk, params) - y_n)*ddyC(xk, params)) * Sig_xn2 + Sig_yn2; //derivative of above

		r = 1.0 - ((yk / ykm1) * (((yk - ykm1) / (xk - xkm1)) / dyk));

		xkp1 = (1.0 - 1.0 / r)*xkm1 + xk / r;
		std::cout << xkp1 << std::endl;

		xkm1 = xk;
		xk = xkp1;

		itercount += 1;
	}

	std::cout << itercount << " iterations." << std::endl;
	return xkp1;
}
int main()
{
    std::cout << "Root Finder Testing\n"; 

	double Sig_xn2 = std::pow(0.1, 2.0) + std::pow(0.1, 2.0);
	double Sig_yn2 = std::pow(0.2, 2.0) + std::pow(0.1, 2.0);

	double x_n = 0.97;
	double y_n = 0.9;

	std::vector <double> params = {1.0, 8.3};

	//double xguess = 0.97;
	double xguess = 0.97;

	double xres = twoPointNR(yC, dyC, ddyC, params, x_n, y_n, Sig_xn2, Sig_yn2, xguess);
	
	std::cout << xres << std::endl;
}
