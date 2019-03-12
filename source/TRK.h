#include <vector>
#include <iostream>
#include <cmath>

const double PI = 3.1415926535897932384626434;

enum whichScaleExtrema{ slopx, slopy, none };

class TRK
{
	public:
		//function pointers
		double (*yc)(double, std::vector <double>);
		double (*dyc)(double, std::vector <double>);
		double (*ddyc)(double, std::vector <double>);

		//dataset
		std::vector <double> x, y, sx, sy, w; //datapoints; errorbars
		double N, M;

		//scaling
		double s, a, b;
		std::vector <double> x_t_slopx, x_t_slopy, x_t_a, x_t_b;

		whichScaleExtrema whichExtrema;

		//parameter guesses
		std::vector <double> params_guess;
		double slop_x_guess, slop_y_guess;

		std::vector <double> allparams_guess;
		std::vector <double> iterative_allparams_guess;

		//constructors
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess);
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess);
		
		//default constructor:
		TRK();

	private:

		// FITTING TOOLS
		double newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess);
		double twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1);
		std::vector <double> cubicSolver(double A, double B, double C, double D);
		std::vector <double> downhillSimplex(double(TRK::*f)(std::vector <double>), std::vector <double> allparams_guess);
		std::vector <double> downhillSimplex(double(*f)(std::vector <double>), std::vector <double> allparams_guess);
		
		// STATISTICS
		double singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn);
		double likelihood(std::vector <double> allparams);
		double stDevUnweighted(std::vector <double> x);
		double modifiedChiSquared(std::vector <double> allparams);
		
		// TANGENT_FINDING ALGORITHMS
		std::vector <double> approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xr1);
		std::vector <double> tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg);
		double findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec);
		
		// SCALE OPTIMIZATION ALGORITHMS
		double innerSlopX_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
		double innerSlopY_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
		std::vector <double> findCentroid(std::vector <std::vector <double> > vertices);
		double optimize_s_SlopX(); //returns the scale s that minimizes slopx
		double optimize_s_SlopY(); //"" slopy
		double R2TRK_prime();
		void optimizeScale();
};