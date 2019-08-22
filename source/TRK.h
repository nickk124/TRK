#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <string>
#include <random>
#include <cfloat>
#include <execution>
#include <future>
#include <functional>
#include <thread>

const double PI = 3.1415926535897932384626434;
const std::vector <double> SIGMAS = { 0.682639, 0.954500, 0.997300 };

enum whichScaleExtrema{ S, slopx, slopy, none };

enum priorTypes { CUSTOM, GAUSSIAN, CONSTRAINED, MIXED};

class Priors
{
public:
	//constructors:
	Priors(priorTypes priorType, std::vector < std::vector <double> > params); //Only Gaussian or only bounded/constrained
	Priors(priorTypes priorType, std::vector < std::vector <double> > gaussianParams, std::vector < std::vector <double> > paramBounds); //mixed
	Priors(priorTypes priorType, std::vector <double(*)(double)> priorsPDFs); //custom
	Priors();

	priorTypes priorType;
	std::vector < std::vector <double> > gaussianParams; // a vector that contains a vector of mu and sigma for the guassian prior of each param. If no prior, then just use NANs for one or both mu and sigma
	std::vector < std::vector <double> > paramBounds; // a vector that contains vectors of the bounds of each param. If not bounded, use NANs, and if there's only one bound, use NAN for the other "bound".
	std::vector <double(*)(double)> priorsPDFs; //a vector that contains (pointers to) the custom prior probability distribution functions for each parameter. If no prior (uninformative) for a parameter, use the function noPrior()

};

struct Results
{
	public:
		std::vector <double> bestFitParams;
		double slop_x;
		double slop_y;
		double optimumScale, minimumScale, maximumScale;
		std::vector < std::vector < std::vector <double> > > bestFit_123Sigmas;
		std::vector < std::vector <double> > slopX_123Sigmas;
		std::vector < std::vector <double> > slopY_123Sigmas;
		std::vector < std::vector < std::vector <double> > > paramDistributionHistograms; // vector: {bins, edges}

};

class TRK
{
	public:
		//constructors
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, std::vector <double> params_sigmas_guess, double slop_x_sigma_guess, double slop_y_sigma_guess);
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, std::vector <double> params_sigmas_guess, double slop_x_sigma_guess, double slop_y_sigma_guess);

		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, std::vector <double> params_sigmas_guess, double slop_x_sigma_guess, double slop_y_sigma_guess, Priors priorsObject);
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, std::vector <double> params_sigmas_guess, double slop_x_sigma_guess, double slop_y_sigma_guess, Priors priorsObject);

		//default constructor:
		TRK();

		//dataset 
		int N, M;
		std::vector <double> SigXVec, SigYVec;
		std::vector <double> x, y, sx, sy, w; //datapoints; errorbars
		double datawidth, x_min, x_max;
		void getDataWidth();

		//core algorithms
		void performTRKFit(); //finds optimum scale AND calculates uncertainties
		void performTRKFit(double scale); //perform fit on some provided scale (for example, if they already know optimum scale, they can just start with this)
		void performSimpleTRKFit(); //finds optimum scale and and performs TRK fit but without finding uncertainties

		//pivot points
		void findPivots();

		//results
		Results results;

		//simplex tools

		std::vector <double> avoidNegativeSlop(std::vector <double> vertex, int n);
		std::vector <double> pegToZeroSlop(std::vector <double> vertex);
		double evalWPriors(double(TRK::*f)(std::vector <double>, double), std::vector <double> vertex, double s);
		double simplex_size = 0.1;

		//scales

		double s, s_sx, s_sy;
		double a, b;
		whichScaleExtrema whichExtrema;

		//tolerances

		double pegToZeroTol = 0.004;
		double slopTol = 1e-4;
		double simplexTol = 1e-5;

		//tangent points, parameters, and guesses
		std::vector <double> x_t_slopx, x_t_slopy, x_t_a, x_t_b, x_t_s;
		std::vector <double> params_slopx, params_slopy, params_a, params_b, params_s, allparams_s, iterative_allparams_guess;
		std::vector <double> params_guess, params_sigmas_guess;
		double slop_x_guess, slop_y_guess, slop_x_sigma_guess, slop_y_sigma_guess;
		std::vector <double> allparams_guess, allparams_sigmas_guess;

		//scale optimization
		std::vector <double> downhillSimplex(double(TRK::*f)(std::vector <double>, double), std::vector <double> allparams_guess, double s);
		double innerSlopX_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
		double innerSlopY_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
		double innerR2_iter_Simplex(std::vector <double> ss, std::vector <double> allparams_guess, double s0);
		double iterateR2_OptimumScale(double s0);

		bool firstGuess = true;
		double slopYGuess;
		double slopYScaleGuess = 1.0;
		void getBetterSlopYGuess(double slop_y, double s);
		void optimizeScale();
		double optimize_s_SlopX();
		double optimize_s_SlopY();
		std::vector <double (TRK::*)()> optimizeList = {&TRK::optimize_s_SlopX, &TRK::optimize_s_SlopY};
		double optimize_s0_R2();
		double optimize_s_prime_R2(double s0);

		double innerR2_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
		std::vector <double> findCentroid(std::vector <std::vector <double> > vertices);

		double R2TRK_prime_as();
		double R2TRK_prime_sb();
		double R2TRK_prime_as0(double s0, std::vector <double> x_t_s0, std::vector <double> params_s0);
		double R2TRK_prime_s0b(double s0, std::vector <double> x_t_s0, std::vector <double> params_s0);


		//tangent finding
		double findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec);
		std::vector <double> approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xr1);
		std::vector <double> tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg);
		double twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1);
		double newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess);
		std::vector <double> cubicSolver(double A, double B, double C, double D);
		double root_bound = 10;
		std::vector <double> tangentParallel(std::vector<double> params, double slop_x, double slop_y, int n, double s);

		//statistics
		bool hasPriors;
		Priors priorsObject;
		double regularChiSquared(std::vector <double> params);
		double modifiedChiSquared(std::vector <double> allparams, double s);
		double normal(double x, double mu, double sig);
		double singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn);
		double likelihood(std::vector <double> allparams);
		double priors(std::vector <double> allparams);
		double posterior(std::vector <double> allparams);
		double stDevUnweighted(std::vector <double> x);
		std::vector <double> tangentParallelLikelihood(std::vector<double> params, double slop_x, double slop_y, int n);

		//function pointers
		double (*yc)(double, std::vector <double>);
		double (*dyc)(double, std::vector <double>);
		double (*ddyc)(double, std::vector <double>);

		// MCMC/uncertainty calculation
		std::vector <std::vector <double >> methastPosterior(int R, int burncount, std::vector <double> sigmas_guess);
		std::vector <std::vector <double >> checkSlopSignMCMC(std::vector <std::vector <double >> result_final);
		std::vector <double> pegToNonZeroDelta(std::vector <double> vertex, std::vector <double> lastvertex);
		std::vector <double> optimizeMetHastDeltas(int burncount, std::vector <double> delta_guess);
		double innerMetHastSimplex(int burncount, std::vector <double> delta, double best_ratio);
		double rnorm(double mu, double sig);
		double runiform(double a, double b);
		std::vector <std::vector <std::vector <double> > >  lowerBar(std::vector <std::vector <double> > allparam_samples);
		std::vector <std::vector <double> > getHistogram(std::vector <double> data);
		void calculateUncertainties();

		int R = 100000; //adjustable; could make accessible by users later
		int burncount = 10000;

		// OTHER TOOLS
		std::vector <double> minMax(std::vector <double> vec);
		std::vector <double> slice(std::vector <double> vec, int l, int r);
		double getAverage(std::vector <double> x);

		// SETTINGS
		bool outputDistributionToFile = false;
		bool cpp17MultiThread = false;
		bool openMPMultiThread = false;
		bool findPivotPoints = false;
		int maxThreads = 8;
};

double noPrior(double param);


//

//for testing only

clock_t startTimer();

double secElapsed(clock_t t_i);

void writeResults(TRK TRKobj, double t_sec, std::string filename);

double toRad(double deg);

double toDeg(double rad);