#include <vector>
#include <iostream>
#include <cmath>

const double PI = 3.1415926535897932384626434;

class TRK
{
	public:
		std::vector <double> downhillSimplex(double(*f)(std::vector <double>), std::vector <double> allparams_guess);

		//function pointers
		double (*yc)(double, std::vector <double>);
		double (*dyc)(double, std::vector <double>);
		double (*ddyc)(double, std::vector <double>);

		//dataset
		std::vector <double> x, y, sx, sy, w; //datapoints; errorbars
		double N, M;

		//parameter guesses
		std::vector <double> params_guess;
		double slop_x_guess, slop_y_guess;

		//constructors
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess);

		//default constructor:
		TRK();

	private:
		// OTHER ALGORITHMS

		// FITTING TOOLS
		double newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess);
		double twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1);
		std::vector <double> cubicSolver(double A, double B, double C, double D);
		// STATISTICS
		double singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn);
		double modifiedChiSquared(std::vector <double> allparams);
		// TANGENT_FINDING ALGORITHMS
		std::vector <double> approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xr1);
		std::vector <double> tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg);
		double findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec);
		// SCALE OPTIMIZATION ALGORITHMS
		std::vector <double> findCentroid(std::vector <std::vector <double> > vertices);
};