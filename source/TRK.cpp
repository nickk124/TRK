#include "pch.h"
#include "TRK.h"
#include <iostream>
#include <cmath>
#include <algorithm>

// CONSTRUCTORS

TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>)) {
	this->yc = (*yc);
	this->dyc = (*dyc);
	this->ddyc = (*ddyc);
}

//default
TRK::TRK() {

}

// TOOLS
double TRK::newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess) {
	double x0 = xguess;
	int itercount = 0;

	//initial iteration
	double f = (yc(x0, params) - y_n) * dyc(x0, params) * Sig_xn2 + (x0 - x_n) * Sig_yn2;  //function we're finding zero of
	double df = (std::pow(dyc(x0, params), 2.0) + (yc(x0, params) - y_n)*ddyc(x0, params)) * Sig_xn2 + Sig_yn2; //derivative of above

	double x1 = x0 - f / df;

	double tol = 1e-9;

	while (std::abs(x1 - x0) > tol) {
		x0 = x1;

		f = (yc(x0, params) - y_n) * dyc(x0, params) * Sig_xn2 + (x0 - x_n) * Sig_yn2;  //function we're finding zero of
		df = (std::pow(dyc(x0, params), 2.0) + (yc(x0, params) - y_n)*ddyc(x0, params)) * Sig_xn2 + Sig_yn2; //derivative of above

		x1 = x0 - f / df;
		std::cout << x1 << std::endl;

		itercount += 1;
	}

	std::cout << itercount << " iterations." << std::endl;
	return x1;
}

double TRK::twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1)
{
	double tol = 1e-9;
	double xkm1 = xguess;
	//double xk = xguess + xguess / 100.0;
	double xk = xguessp1;

	double xkp1;
	double r;
	double ykm1;
	double yk;
	double dyk;

	int itercount = 0;

	while (std::abs(xk - xkm1) > tol) {
		ykm1 = (yc(xkm1, params) - y_n) * dyc(xkm1, params) * Sig_xn2 + (xkm1 - x_n) * Sig_yn2;
		yk = (yc(xk, params) - y_n) * dyc(xk, params) * Sig_xn2 + (xk - x_n) * Sig_yn2;  //function we're finding zero of
		dyk = (std::pow(dyc(xk, params), 2.0) + (yc(xk, params) - y_n)*ddyc(xk, params)) * Sig_xn2 + Sig_yn2; //derivative of above

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

std::vector <double> TRK::cubicSolver(double A, double B, double C, double D) {
	//cubic solver for three real and distinct roots
	std::vector <double> roots;

	double a1 = B / A;
	double a2 = C / A;
	double a3 = D / A;

	double Q = (a1 * a1 - 3.0 * a2) / 9.0;
	double R = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
	double Qc = std::pow(Q, 3.0);
	double d = Qc - std::pow(R, 2.0);

	double theta = std::acos(R / sqrt(Qc));

	roots.push_back( -2 * std::sqrt(Q) * std::cos(theta / 3) - a1 / 3 );
	roots.push_back(-2 * std::sqrt(Q) * std::cos((theta + 2 * PI) / 3) - a1 / 3 );
	roots.push_back(-2 * std::sqrt(Q) * std::cos((theta + 4 * PI) / 3) - a1 / 3 );

	return roots;
}

std::vector <double> TRK::approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg1) {
	//using board derivation notation:
	double b = yc(xg1, params) - dyc(xg1, params) * xg1 + (ddyc(xg1, params) / 2.0) * std::pow(xg1, 2.0); //coefficients of quadratic approximation from taylor expansion
	double m = dyc(xg1, params) - ddyc(xg1, params) * xg1;
	double a = ddyc(xg1, params) / 2.0;


	//DIFFERENT FROM NOTATION OF BOARD DERIVATION!
	double A = 2.0 * std::pow(a, 2.0);											// coef of x^3
	double B = 3.0 * a * m;														//coef of x^2
	double C = (2.0 * a * (b - y_n) + std::pow(m, 2.0)) * Sig_xn2 + Sig_yn2;	//coef of x
	double D = m * (b - y_n) * Sig_xn2 - x_n * Sig_yn2;							//coef of 1

	double discriminant = 18.0*A*B*C*D - 4.0*std::pow(B, 3.0)*D + std::pow(B, 2.0)*std::pow(C, 2.0) - 4.0*A*std::pow(C, 3.0) - 27.0*std::pow(A, 2.0)*std::pow(D, 2.0);

	std::vector <double> roots;

	if (discriminant > 0) {
		roots = cubicSolver(A, B, C, D);
		return roots;
	}
	//returns no extra roots (empty vector) if the other two roots are 
	return roots;
}

std::vector <double> TRK::tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg) {
	
	std::vector <double> result;

	double xg1 = xg;

	while (true) {
		result.clear();

		double xr1 = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg1, xg1 + std::sqrt(Sig_xn2));
		double xg2, xg3;

		//Quadratic Approximation

		std::vector <double> allRoots = approxQuadraticRoots(params, x_n, y_n, Sig_xn2, Sig_yn2, xg1); //get approx. roots from quadratic taylor approximation
		std::vector <double> extraRoots;

		if (allRoots.size() == 3) { //either add the two other roots, or no more roots (depending on discriminant of cubic equation)
			//grab to new roots 
			for (int i = 0; i < 3; i++) {
				double root = allRoots[i];
				if (std::abs(root - xr1) >= 1e-9) { //checks if roots are (numerically) identical or not
					extraRoots.push_back(root);
				}
			}
		}

		if (extraRoots.size() == 2) { //if have 2 additional, real roots
			double xr2, xr3;
			xg2 = extraRoots[1];
			xg3 = extraRoots[2];

			if (xg2 < xr1 && xr1 < xg3) {
				xr2 = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg2, xg2 + std::sqrt(Sig_xn2));
				xr3 = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg3, xg3 + std::sqrt(Sig_xn2));

				result.push_back(xr1);
				result.push_back(xr2);
				result.push_back(xr3);

				break;
			}
			else {
				std::vector <double> rootVec = { xr1, xg2, xg3 };
				std::sort(rootVec.begin(), rootVec.end());

				xg1 = rootVec[1];
			}
		} else if (extraRoots.size() == 0) { //if only have one root (the one initially found)
			result.push_back(xr1);
			break;
		}
	}

	return result;
}