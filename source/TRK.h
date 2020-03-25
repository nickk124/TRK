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
        double pivot2;
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
    private:
        // dataset
        unsigned long bigM; //bigM is total number of model AND slop params
        double datawidth, x_min, x_max;
        std::vector <double> SigXVec, SigYVec;
    
        void getDataWidth();
        void checkZeroErrorBars();
    
    
        // downhill simplex (likelihood maximization)
        double simplex_size = 0.1;
    
        double evalWPriors(double(TRK::*f)(std::vector <double>, double), std::vector <double> vertex, double s);
        std::vector <double> avoidNegativeSlop(std::vector <double> vertex, unsigned long n);
        std::vector <double> pegToZeroSlop(std::vector <double> vertex);
        std::vector <double> downhillSimplex(double(TRK::*f)(std::vector <double>, double), std::vector <double> allparams_guess, double s);
        std::vector <double> findCentroid(std::vector <std::vector <double> > vertices);
    
    
        // asymmetric uncertainties
        bool hasAsymEB = false;
        bool hasAsymSlop = false;
    
        void checkAsym();
        double modifiedChiSquaredAsym(std::vector <double> allparams, double s);
        double likelihoodAsym(std::vector <double> allparams);
        double singlePointLnLAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, double x_tn, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn);
        double pnAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, double x_tn, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn);
        double dunDxAsym(double mtn, std::vector <double> Sigs2, int quadSig_xn2Ind, int quadSig_yn2Ind, double s);
        double zAsym(double x, double quadSig_xn2, double quadSig_yn2, double xn_shifted, double yn_shifted, std::vector <double> shifts, double x_tn, double y_tn, double m_tn);
        double findBestTangentAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, std::vector <double> x_tn_vec, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn);
        std::vector <double> getAsymShifts(std::vector <double> allparams, int n);
        std::vector <double> getAsymSigs2(std::vector <double> allparams, int n);
        std::vector <double> tangentParallelAsym(std::vector<double> allparams, int n, double s);
        std::vector <double> tangentParallelLikelihoodAsym(std::vector<double> allparams, int n);
    
    
        // fitting scales
        double s, s_sx, s_sy;
        double a, b;
        whichScaleExtrema whichExtrema = none;
        whichScaleExtrema whichExtremaX = none;
        whichScaleExtrema whichExtremaY = none;

    
        // iterative tolerances
        double pegToZeroTol = 0.004;
    
    
        // tangent points, parameters, and guesses
        double slop_x_guess, slop_y_guess, slop_x_sigma_guess, slop_y_sigma_guess, slop_x_minus_sigma_guess, slop_y_minus_sigma_guess;
        std::vector <double> x_t_slopx, x_t_slopy, x_t_a, x_t_b, x_t_s;
        std::vector <double> params_slopx, params_slopy, params_a, params_b, params_s, allparams_s, iterative_allparams_guess;
        std::vector <double> params_guess, params_sigmas_guess;
        std::vector <double> allparams_guess, allparams_sigmas_guess;

    
        // scale optimization
        bool firstGuess = true;
        double slopYGuess;
        double slopYScaleGuess = 1.0;
        std::vector <double (TRK::*)()> optimizeList = {&TRK::optimize_s_SlopX, &TRK::optimize_s_SlopY};
    
        void getBetterSlopYGuess(double slop_y, double s);
        void optimizeScale();
        void getBetterGuess();
        double innerSlopX_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
        double innerSlopY_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
        double innerR2_iter_Simplex(std::vector <double> ss, std::vector <double> allparams_guess, double s0);
        double iterateR2_OptimumScale(double s0);
        double optimize_s_SlopX();
        double optimize_s_SlopY();
        double optimize_s0_R2();
        double optimize_s_prime_R2(double s0);
        double innerR2_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
        double R2TRK_prime_as();
        double R2TRK_prime_sb();
        double R2TRK_prime_as0(double s0, std::vector <double> x_t_s0, std::vector <double> params_s0);
        double R2TRK_prime_s0b(double s0, std::vector <double> x_t_s0, std::vector <double> params_s0);
    
    
        // tangent point finding
        double root_bound = 10;
    
        double findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec, double s);
        double twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1);
        double newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess);
        std::vector <double> approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xr1);
        std::vector <double> tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg);
        std::vector <double> tangentCubicSolver(double A, double B, double C, double D);
        std::vector <double> tangentParallel(std::vector<double> params, double slop_x, double slop_y, int n, double s);

    
        //statistics
        bool hasPriors;
        Priors priorsObject;
        
        bool isEqual(double x, double y, double maxRelativeError, double maxAbsoluteError);
        double regularChiSquared(std::vector <double> params);
        double regularChiSquaredWSlop(std::vector <double> allparams, double s);
        double modifiedChiSquared(std::vector <double> allparams, double s);
        double normal(double x, double mu, double sig);
        double cmNorm(double z); //cumulative unit normal
        double singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn, double s);
        double likelihood(std::vector <double> allparams);
        double likelihood1D(std::vector <double> allparams);
        double priors(std::vector <double> allparams);
        double posterior(std::vector <double> allparams, std::vector <double> allparams_trial);
        double stDevUnweighted(std::vector <double> x);
        double getMedian(std::vector<double> y);
        double getMedian(int trueCount, std::vector<double> w, std::vector<double> y);
        double getAverage(std::vector <double> x, std::vector <double> w);
        double getAverage(std::vector <double> x);
        double min(double a, double b);
        double max(double a, double b);
        double getMode(int trueCount, std::vector<double> w, std::vector<double> y);
        std::vector <std::vector <double> > getHistogram(std::vector <double> data);
        std::vector <std::vector <double> > getHistogram(std::vector <double> data, std::vector <double> weights);
        std::vector <double> tangentParallelLikelihood(std::vector<double> params, double slop_x, double slop_y, int n);
    
        
        //function pointers
        double (*yc)(double, std::vector <double>);
        double (*dyc)(double, std::vector <double>);
        double (*ddyc)(double, std::vector <double>);
        double (TRK::*selectedChiSq)(std::vector <double>, double) = &TRK::modifiedChiSquared;
        double (TRK::*selectedLikelihood)(std::vector <double>) = &TRK::likelihood;
    
    
        // MCMC/uncertainty calculation and RNG
        bool useLogPosterior = false;
        bool currentlyOptimizingProposal = false;
        double best_ratio = 0.325;
        double simplexSuperShrink = 1e-3;
        std::vector <double> allParamsFinalDeltas;
        tuningAlgo thisTuningAlgo = AM;
    
        void calculateUncertainties();
        void guessMCMCDeltas();
        double innerMetHastSimplex(int burncount, std::vector <double> delta, double best_ratio);
        double metHastRatio(std::vector <double> X_trial, std::vector <double> X_i);
        double rnorm(double mu, double sig);
        double runiform(double a, double b);
        std::vector <double> pegToNonZeroDelta(std::vector <double> vertex, std::vector <double> lastvertex);
        std::vector <double> optimizeMetHastDeltas(int burncount, std::vector <double> delta_guess);
        std::vector <std::vector <double >> methastPosterior(int R, int burncount, std::vector <double> sigmas_guess);
        std::vector <std::vector <double >> checkSlopSignMCMC(std::vector <std::vector <double >> result_final);
        std::vector <std::vector <std::vector <double> > > lowerBar(std::vector <std::vector <double> > allparam_samples);
    
        
        // pivot points
        bool getCombosFromSampleDirectly = true;
        bool weightPivots = true;
        bool pivotMedian = false;
        bool pivotMean = false;
        bool pruneOutlierPivots = true;
        bool pivotPointActive = false;
        bool pivotHalfSampleMode= false;
        bool modeInterceptGuess = false;
        bool averageIntercepts = false;
        int pivotR = 5000; //1000 too low; 5000 seems sufficient, but 10,000 works for sure
        int randomSampleCount = 450;
        int maxCombos = 50000; // 50,000 seems sufficient, but 100,000 works for sure
        int maxPivotIter = 1; // 1 is usually sufficient, as successive iterations seem to only jump around (may not be true for linearIZED models, not just linear, however)
        int pivotBurnIn = 1000;
        double pivotTol = 1e-1;
        double pruneWidth = 10.0;
        std::vector <double> pivotPointParamsGuess;
        std::vector < std::vector <double> > NDcombination;
        std::vector < std::vector <std::vector <double > > > NDcombos;

        void findPivots();
        void getPivotGuess();
        void getCombos(std::vector <std::vector <double> > total, int k, int offset);
        double weightPivot(std::vector <double> params1, std::vector <double> params2, std::vector <double> oldPivots, double newPivot);
        double pivotFunc(std::vector <double> params1, std::vector <double> params2);
        double pivotFunc2(std::vector <double> params1, std::vector <double> params2); // for two-pivot models
        std::vector < std::vector <std::vector <double > > > directCombos(std::vector < std::vector <double> > params_sample, int comboCount);
        std::vector <double> removeOutlierPivots(std::vector <double> pivots);


        // OTHER TOOLS
        double getPeakCoord(std::vector <double> x, std::vector <double> w);
        std::vector < std::vector <double> > transpose(std::vector < std::vector <double> > array);
        void showResults(bool didScaleOp, bool didMCMC);
        void checkVerbose();
    
    
        // GENERAL SETTINGS
        bool showSimplexSteps = false;
        bool verboseAsymmetric = false; // show info/steps about asymmetric uncertainty fitting
        int maxThreads = 8;
        bool do1DFit = false; // does not need to be changed by the user, do 1D fits simply using the correct TRK constructor. Only used for testing
    
    
	public:
		// CONSTRUCTORS
        // weighted
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess);
        // unweighted
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess);
        // with priors, weighted
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, Priors priorsObject);
        // with priors, unweighted
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, Priors priorsObject);
    
    
        // 1D statistic (weighted)
        TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess);
        // equal weights/unweighted
        TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess);
        // with priors (weighted)
        TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess, Priors priorsObject);
        // equal weights/unweighted with priors
        TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess, Priors priorsObject);
    

		// default constructor:
		TRK();

    
		// dataset
        unsigned long N, M; // number of data points and number of model parameters, respectively (not including slop)
		std::vector <double> x, y, sx, sy, w; //datapoints; errorbars

    
		// core algorithms
		void performTRKFit(); //finds optimum scale AND calculates uncertainties
		void performTRKFit(double scale); //perform fit on some provided scale (for example, if they already know optimum scale, they can just start with this) and calculates uncertainties
		void performSimpleTRKFit(); //finds optimum scale and and performs TRK fit but without finding uncertainties
        void performSimpleTRKFit(double scale); //given some scale, performs TRK fit without finding uncertainties.

    
		// results
		Results results;
    
    
        // asymmetric uncertainties
        double slop_x_minus_guess = -1.0;  // negative asymmetric slop
        double slop_y_minus_guess = -1.0;
        std::vector <double> sx_minus, sy_minus; // negative asymmetric error bars

    
		// iterative tolerances
		double simplexTol = 1e-5;


		// MCMC/uncertainty calculation and RNG
		int R = 100000; // MCMC sample size (excluding burn-in)
		int burncount = 10000; // MCMC "burn in" count
    

		// pivot points / linearized model parameter correlation removal
        bool twoPivots = false; // two pivot points in the model?
        static double pivot;
        static double pivot2; // for two-pivot models
        double (*linearizedIntercept)(std::vector <double>);
        double (*linearizedSlope)(std::vector <double>);
        double (*linearizedIntercept2)(std::vector <double>); // for two-pivot point models, e.g. broken-linear (experimental)
        double (*linearizedSlope2)(std::vector <double>);
    
        
		// SETTINGS
        bool findPivotPoints = false;
		bool cpp17MultiThread = false; // should only have one of these multithread
		bool cpp11MultiThread = true;
		bool openMPMultiThread = false;
    
        // OUTPUT SETTINGS
        bool printResults = true;
        bool outputDistributionToFile = false; // outputs MCMC-sampled parameter distribution to a text file (must specify path below)
        std::string outputPath;
        bool verbose = false; // sets the below three to true automatically
        bool verboseMCMC = false; // show MCMC steps
        bool verbosePivotPoints = false; // show pivot point finding steps
    
        // TESTING
        bool writePivots = false;
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
