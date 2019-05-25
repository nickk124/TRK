#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <string>

const double PI = 3.1415926535897932384626434;

const double PHI = (std::sqrt(5.0) + 1.0) / 2.0;

enum whichScaleExtrema{ S, slopx, slopy, none };

class TRK
{
	public:

		//public only for testing
		std::vector <double> pegToZeroSlop(std::vector <double> vertex);
		double s;
		double a, b;
		whichScaleExtrema whichExtrema;
		std::vector <double> iterative_allparams_guess;

		std::vector <double> x_t_slopx, x_t_slopy, x_t_a, x_t_b, x_t_s;
		std::vector <double> params_slopx, params_slopy, params_a, params_b, params_s;

		double innerSlopX_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
		double innerSlopY_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
		double innerR2_iter_Simplex(std::vector <double> ss, std::vector <double> allparams_guess, double s0);
		double iterateR2_OptimumScale(double s0);

		std::vector <double> tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg);
		double twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1);
		double modifiedChiSquared(std::vector <double> allparams);
		std::vector <double> downhillSimplex(double(TRK::*f)(std::vector <double>), std::vector <double> allparams_guess);
		void optimizeScale();
		double optimize_s_SlopX();
		double optimize_s_SlopY();
		double optimize_s0_R2();
		double optimize_s_prime_R2(double s0);
		double regularChiSquared(std::vector <double> params);

		//function pointers
		double (*yc)(double, std::vector <double>);
		double (*dyc)(double, std::vector <double>);
		double (*ddyc)(double, std::vector <double>);

		//dataset
		std::vector <double> x, y, sx, sy, w; //datapoints; errorbars
		double datawidth;
		void getDataWidth();
		double root_bound = 1000;

		//parameter guesses
		std::vector <double> params_guess;
		double slop_x_guess, slop_y_guess;

		//constructors
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess);
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess);
		
		//default constructor:
		TRK();

	private:

		// dataset
		double N, M;

		//scaling
		
		
		

		

		//parameter guesses
		std::vector <double> allparams_guess;
		

		// OTHER TOOLS
		std::vector <double> minMax(std::vector <double> vec);

		// FITTING TOOLS
		double newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess);
		std::vector <double> cubicSolver(double A, double B, double C, double D);
		
		// STATISTICS
		double singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn);
		double likelihood(std::vector <double> allparams);
		double stDevUnweighted(std::vector <double> x);
		
		// TANGENT_FINDING ALGORITHMS
		std::vector <double> approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xr1);
		double findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec);
		
		// SCALE OPTIMIZATION ALGORITHMS
		
		double innerR2_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
		std::vector <double> findCentroid(std::vector <std::vector <double> > vertices);

		double R2TRK_prime_as();
		double R2TRK_prime_sb();
		double R2TRK_prime_as0(double s0, std::vector <double> x_t_s0, std::vector <double> params_s0);
		double R2TRK_prime_s0b(double s0, std::vector <double> x_t_s0, std::vector <double> params_s0);
};


//for testing only

clock_t startTimer();

double secElapsed(clock_t t_i);

void writeResults(TRK TRKobj, double t_sec, std::string filename);