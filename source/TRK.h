#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <string>
#include <random>
#include <cfloat>
//#include <execution>
#include <future>
#include <functional>
#include <thread>
#include <sstream>

const double PI = 3.1415926535897932384626434;
const std::vector <double> SIGMAS = { 0.682639, 0.954500, 0.997300 };

enum whichScaleExtrema{ S, slopx, slopy, none };

enum priorTypes { CUSTOM, GAUSSIAN, CONSTRAINED, MIXED};

enum tuningAlgo {SIMPLEX, AM};

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
        double pivot;
		std::vector < std::vector < std::vector <double> > > bestFit_123Sigmas;
		std::vector < std::vector <double> > slopX_123Sigmas;
		std::vector < std::vector <double> > slopY_123Sigmas;
		std::vector < std::vector < std::vector <double> > > paramDistributionHistograms; // vector: {bins, edges}
    
        // asymmetric stuff:
        double slop_x_minus;
        double slop_y_minus;
        std::vector < std::vector <double> > slopX_minus_123Sigmas;
        std::vector < std::vector <double> > slopY_minus_123Sigmas;

};

class TRK
{
	public:
		//constructors
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess);
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess);

		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, Priors priorsObject);
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, Priors priorsObject);
    

		//default constructor:
		TRK();

		//dataset 
		unsigned long N, M;
		std::vector <double> SigXVec, SigYVec;
		std::vector <double> x, y, sx, sy, w; //datapoints; errorbars
		double datawidth, x_min, x_max;
		void getDataWidth();
        void checkZeroErrorBars();

		//core algorithms
		void performTRKFit(); //finds optimum scale AND calculates uncertainties
		void performTRKFit(double scale); //perform fit on some provided scale (for example, if they already know optimum scale, they can just start with this) and calculates uncertainties
		void performSimpleTRKFit(); //finds optimum scale and and performs TRK fit but without finding uncertainties
        void performSimpleTRKFit(double scale); //given some scale, performs TRK fit without finding uncertainties.

		//results
		Results results;
    
        //asymmetric distribution tools
        double slop_x_minus_guess = -1.0;  // negative asymmetric slop
        double slop_y_minus_guess = -1.0;
        std::vector <double> sx_minus, sy_minus; // negative asymmetric error bars
        double cumulNorm(double z);
        bool hasAsymEB = false;
        bool hasAsymSlop = false;
        void checkAsym();
        double modifiedChiSquaredAsym(std::vector <double> allparams, double s);
        double dunDxAsym(double mtn, std::vector <double> Sigs2, int quadSig_xn2Ind, int quadSig_yn2Ind, double s);
        double cmNorm(double z);
        double zAsym(double x, double quadSig_xn2, double quadSig_yn2, double xn_shifted, double yn_shifted, std::vector <double> shifts, double x_tn, double y_tn, double m_tn);
        double pnAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, double x_tn, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s);
        double singlePointLnLAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, double x_tn, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s);
        double findBestTangentAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, std::vector <double> x_tn_vec, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s);
        std::vector <double> getAsymShifts(std::vector <double> allparams, int n);
        std::vector <double> getAsymSigs2(std::vector <double> allparams, int n);
        std::vector <double> tangentParallelAsym(std::vector<double> allparams, int n, double s);
        
        // Asym MCMC
        double likelihoodAsym(std::vector <double> allparams);
        std::vector <double> tangentParallelLikelihoodAsym(std::vector<double> allparams, int n);

		//simplex tools

		std::vector <double> avoidNegativeSlop(std::vector <double> vertex, unsigned long n);
		std::vector <double> pegToZeroSlop(std::vector <double> vertex);
		double evalWPriors(double(TRK::*f)(std::vector <double>, double), std::vector <double> vertex, double s);
		double simplex_size = 0.1;

		//scales

		double s, s_sx, s_sy;
		double a, b;
		whichScaleExtrema whichExtrema = none;			
		whichScaleExtrema whichExtremaX = none;
		whichScaleExtrema whichExtremaY = none;

		//tolerances

		double pegToZeroTol = 0.004;
		double slopTol = 1e-4;
		double simplexTol = 1e-5;

		//tangent points, parameters, and guesses
		std::vector <double> x_t_slopx, x_t_slopy, x_t_a, x_t_b, x_t_s;
		std::vector <double> params_slopx, params_slopy, params_a, params_b, params_s, allparams_s, iterative_allparams_guess;
		std::vector <double> params_guess, params_sigmas_guess;
        double slop_x_guess, slop_y_guess, slop_x_sigma_guess, slop_y_sigma_guess, slop_x_minus_sigma_guess, slop_y_minus_sigma_guess;
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
        void getBetterGuess();

		double innerR2_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
		std::vector <double> findCentroid(std::vector <std::vector <double> > vertices);

		double R2TRK_prime_as();
		double R2TRK_prime_sb();
		double R2TRK_prime_as0(double s0, std::vector <double> x_t_s0, std::vector <double> params_s0);
		double R2TRK_prime_s0b(double s0, std::vector <double> x_t_s0, std::vector <double> params_s0);


		//tangent finding
		double findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec, double s);
		std::vector <double> approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xr1);
		std::vector <double> tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg);
		double twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1);
		double newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess);
		std::vector <double> tangentCubicSolver(double A, double B, double C, double D);
		double root_bound = 10;
		std::vector <double> tangentParallel(std::vector<double> params, double slop_x, double slop_y, int n, double s);

		//statistics
		bool hasPriors;
		Priors priorsObject;
		double regularChiSquared(std::vector <double> params);
        double regularChiSquaredWSlop(std::vector <double> allparams, double s);
		double modifiedChiSquared(std::vector <double> allparams, double s);
		double normal(double x, double mu, double sig);
		double singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn, double s);
		double likelihood(std::vector <double> allparams);
		double priors(std::vector <double> allparams);
		double posterior(std::vector <double> allparams, std::vector <double> allparams_trial);
		double stDevUnweighted(std::vector <double> x);
		std::vector <double> tangentParallelLikelihood(std::vector<double> params, double slop_x, double slop_y, int n);
		double getMedian(std::vector<double> y);
		double getMedian(int trueCount, std::vector<double> w, std::vector<double> y);
        double getAverage(std::vector <double> x, std::vector <double> w);
        double getAverage(std::vector <double> x);
        double min(double a, double b);
        double max(double a, double b);
        bool isEqual(double x, double y, double maxRelativeError, double maxAbsoluteError);
        double getMode(int trueCount, std::vector<double> w, std::vector<double> y);
        std::vector <std::vector <double> > getHistogram(std::vector <double> data);
        std::vector <std::vector <double> > getHistogram(std::vector <double> data, std::vector <double> weights);
        
    
    

		//function pointers
		double (*yc)(double, std::vector <double>);
		double (*dyc)(double, std::vector <double>);
		double (*ddyc)(double, std::vector <double>);
        
        double (TRK::*selectedChiSq)(std::vector <double>, double) = &TRK::modifiedChiSquared;
        double (TRK::*selectedLikelihood)(std::vector <double>) = &TRK::likelihood;


		// MCMC/uncertainty calculation
        tuningAlgo thisTuningAlgo = AM;
    
		std::vector <std::vector <double >> methastPosterior(int R, int burncount, std::vector <double> sigmas_guess);
		std::vector <std::vector <double >> checkSlopSignMCMC(std::vector <std::vector <double >> result_final);
		std::vector <double> pegToNonZeroDelta(std::vector <double> vertex, std::vector <double> lastvertex);
		std::vector <double> optimizeMetHastDeltas(int burncount, std::vector <double> delta_guess);
		double innerMetHastSimplex(int burncount, std::vector <double> delta, double best_ratio);
        double metHastRatio(std::vector <double> X_trial, std::vector <double> X_i);
		double rnorm(double mu, double sig);
		double runiform(double a, double b);
		std::vector <std::vector <std::vector <double> > > lowerBar(std::vector <std::vector <double> > allparam_samples);
		void calculateUncertainties();
		std::vector <double> allParamsFinalDeltas;
		bool goodDeltasFound = false;
        void guessMCMCDeltas();

		int R = 100000; //adjustable; could make accessible by users later
		int burncount = 10000;
        double best_ratio = 0.325;
        double simplexSuperShrink = 1e-3;
    
        bool useLogPosterior = false;
        bool currentlyOptimizingProposal = false;
    

		//pivot points
		std::vector < std::vector <std::vector <double > > > NDcombos;
		std::vector < std::vector <double> > NDcombination;
		void getCombos(std::vector <std::vector <double> > total, int k, int offset);
		void findPivots();
		static double pivot;
		//double(*pf)(std::vector <double>, std::vector <double>); //pointer to pivotFunction: arguments of std::vector <double> params1, std::vector <double> params2
		double pivotTol = 1e-3;
		double weightPivot(std::vector <double> params1, std::vector <double> params2, std::vector <double> oldPivots, double newPivot);
		double pivotFunc(std::vector <double> params1, std::vector <double> params2);
        std::vector < std::vector <std::vector <double > > > directCombos(std::vector < std::vector <double> > params_sample, int comboCount);
        std::vector <double> removeOutlierPivots(std::vector <double> pivots);
    
        double (*linearizedIntercept)(std::vector <double>);
        double (*linearizedSlope)(std::vector <double>);
    

		bool getCombosFromSampleDirectly = true;
		bool weightPivots = true;
    
		bool writePivots = false;
    
        bool pivotMedian = false;
        bool pivotMean = false;
        bool pruneOutlierPivots = true;
        double pruneWidth = 10.0;
        int pivotR = 10000; //1000 too low
        int randomSampleCount = 450;
        int maxCombos = 100000;
        int maxPivotIter = 10;
        int pivotBurnIn = 1000;
        bool pivotHalfSampleMode= false;
        bool modeInterceptGuess = false;
        bool averageIntercepts = false;
    
        void getPivotGuess();
    
        bool pivotPointActive = false;
        std::vector <double> pivotPointParamsGuess;
    
		// OTHER TOOLS
        double getPeakCoord(std::vector <double> x, std::vector <double> w);
        std::vector < std::vector <double> > transpose(std::vector < std::vector <double> > array);

		// SETTINGS
		bool outputDistributionToFile = false;
		bool cpp17MultiThread = false;
		bool cpp11MultiThread = true;
		bool openMPMultiThread = false;
		bool findPivotPoints = false;
        bool showSimplexSteps = false;
		int maxThreads = 8;
};

//global functions

template <class vec>
vec slice(vec v, int begin, int end) {
	vec x;
	int len = end - begin;

	for (int i = 0; i < len; i++) {
		x.push_back(v[begin + i]);
	}

	return x;
};

template <class vec>
vec concat(vec v, vec u ) {
    vec x;
    for (int i = 0; i < v.size(); i++){
        x.push_back(v[i]);
    }
    
    for (int i = 0; i < u.size(); i++){
        x.push_back(u[i]);
    }
    
    return x;
};

template <class BidiIter >
BidiIter random_unique(BidiIter begin, BidiIter end, size_t num_random) {
	size_t left = std::distance(begin, end);
	while (num_random--) {
		BidiIter r = begin;
		std::advance(r, rand() % left);
		std::swap(*begin, *r);
		++begin;
		--left;
	}
	return begin;
};

double noPrior(double param);

std::vector <double> minMax(std::vector <double> vec);
std::vector <int> argMinMax(std::vector <double> x);

double twoPointNR(double(*y)(double, std::vector <double>), double(*dy)(double, std::vector <double>), double(*ddy)(double, std::vector <double>), std::vector <double> params, double xguess, double xguessp1);

std::vector <double> cubicSolver(double A, double B, double C, double D);

//for testing only

clock_t startTimer();

double secElapsed(clock_t t_i);

void writeResults(TRK TRKobj, double t_sec, std::string filename);

double toRad(double deg);

double toDeg(double rad);

std::vector <std::vector <double > > getData(std::string fileName, int dataSize);

std::vector <std::vector <double > > getData(std::string fileName);
