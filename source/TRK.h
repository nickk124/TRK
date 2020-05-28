/*
Trotter Reichart Konz (TRK) Official Codebase
Author: Nick C. Konz
See license at https://github.com/nickk124/TRK
*/

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

enum whichScaleExtrema{ S, slopx, slopy, NONE };

enum priorTypes { CUSTOM, CUSTOM_JOINT, GAUSSIAN, CONSTRAINED, MIXED};

enum samplingMethod {ARWMH, AIES};   // Adaptive Random Walk Metropolis Hastings, and Affine Invariant Ensemble Sampler

enum pivotPointFindingMethod {DIST, REGRESSION, PEARSON_GSS, PEARSON_};  // Distribution generation, vs. regression of confidence ellipse, vs. minimizing pearson correlation of said ellipse

enum ParallelizationBackEnd {CPP11, OPENMP};

class Priors
{
public:
	//constructors:
	Priors(priorTypes priorType, std::vector < std::vector <double> > params); //Only Gaussian or only bounded/constrained
	Priors(priorTypes priorType, std::vector < std::vector <double> > gaussianParams, std::vector < std::vector <double> > paramBounds); //mixed
	Priors(priorTypes priorType, std::vector <double(*)(double)> priorsPDFs); // custom
    Priors(priorTypes priorType, double(*jointPriorsPDF)(std::vector <double>));// custom joint
	Priors();

	priorTypes priorType;
	std::vector < std::vector <double> > gaussianParams; // a vector that contains a vector of mu and sigma for the guassian prior of each param. If no prior, then just use NANs for one or both mu and sigma
	std::vector < std::vector <double> > paramBounds; // a vector that contains vectors of the bounds of each param. If not bounded, use NANs, and if there's only one bound, use NAN for the other "bound".
	std::vector <double(*)(double)> priorsPDFs; // a vector that contains (pointers to) the custom prior probability distribution functions for each parameter. If no prior (uninformative) for a parameter, use the function noPrior()
    double (*jointPriorsPDF)(std::vector <double>); // a pointer to a function that is the joint prior probability, i.e. it takes an argument of a vector of the model params (including slop, last), and returns the joint prior.
};

struct Results
{
    public:
        double slop_x;
        double slop_y;
        double optimumScale, minimumScale, maximumScale;
        double chisquared; // not exactly chi-squared; really -2ln L
        std::vector <double> pivots;
        std::vector <double> bestFitParams;
        std::vector < std::vector <double> > slopX_123Sigmas;
        std::vector < std::vector <double> > slopY_123Sigmas;
        std::vector < std::vector < std::vector <double> > > bestFit_123Sigmas;
        std::vector < std::vector < std::vector <double> > > paramDistributionHistograms; // vector: {bins, edges}
    
        // asymmetric uncertainties
        double slop_x_minus;
        double slop_y_minus;
        std::vector < std::vector <double> > slopX_minus_123Sigmas;
        std::vector < std::vector <double> > slopY_minus_123Sigmas;

};

class TRK // main class
{
    public:
        // NESTED CLASSES
        class Statistics
        {
            friend TRK;
            TRK &trk; // parent class
        
            public:
                // constructors/destructor
                Statistics(TRK &trk);
                ~Statistics();
            
            private:
                // function pointers
                double (Statistics::*selectedChiSq)(std::vector <double>, double) = &Statistics::modifiedChiSquared;
                double (Statistics::*selectedLikelihood)(std::vector <double>) = &Statistics::likelihood;
            
                typedef double (TRK::Statistics::* statsfunc)(std::vector<double>);
            
                // priors
                bool hasPriors;
                Priors priorsObject;
                double priors(std::vector <double> allparams);
            
                // regression
                std::vector <double> simpleLinearRegression(std::vector <double> x, std::vector <double> y);
            
                // correlation
                double pearsonCorrelation(std::vector <double> x, std::vector <double> y);
            
            
                // goodness of fit metrics
                double regularChiSquared(std::vector <double> params);
                double regularChiSquaredWSlop(std::vector <double> allparams, double s);
                double modifiedChiSquared(std::vector <double> allparams, double s);
            
                // probability distributions
                double normal(double x, double mu, double sig);
                double cmNorm(double z); //cumulative unit normal
                double stretch_pdf(double z);
                
                // likelihoods and posteriors
                double likelihood1D(std::vector <double> allparams);
                double singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn, double s);
                double likelihood(std::vector <double> allparams);
            
                // likelihoods and posteriors with asymmetric uncertainties
                double modifiedChiSquaredAsym(std::vector <double> allparams, double s);
                double likelihoodAsym(std::vector <double> allparams);
            
                // statistics (in the literal sense)
                double stDevUnweighted(std::vector <double> x);
                double getAverage(std::vector <double> x, std::vector <double> w);
                double getAverage(std::vector <double> x);
                double getMode(int trueCount, std::vector<double> w, std::vector<double> y);
                double getPeakCoord(std::vector <double> x, std::vector <double> w);
                double min(double a, double b);
                double max(double a, double b);
            
                // histograms
                std::vector <std::vector <double> > getHistogram(std::vector <double> data);
                std::vector <std::vector <double> > getHistogram(std::vector <double> data, std::vector <double> weights);
        };
    
        class Optimization
        {
            friend TRK; // need TRK instances to be able to access private members of this
            TRK &trk; // parent class
            
            public:
                // constructors/destructor
                Optimization(TRK &trk);
                ~Optimization();
            
                // settings
                double simplexTol = 1e-5;  // downhill simplex fitting tolerance
                int max_simplex_iters = 10000;   // maximum number of iterations for downhill simplex
                bool showSimplexSteps = false;
            
            
            private:
            
                // downhill simplex method for finding the minimum of some N-D function
                std::vector <double> downhillSimplex(double(TRK::Statistics::*f)(std::vector <double>, double), std::vector <double> allparams_guess, double s);
                
                // downhill simplex tools
                double evalWPriors(double(TRK::Statistics::*f)(std::vector <double>, double), std::vector <double> vertex, double s);
                std::vector <double> avoidNegativeSlop(std::vector <double> vertex, unsigned long n);
                std::vector <double> pegToZeroSlop(std::vector <double> vertex);
                std::vector <double> findCentroid(std::vector <std::vector <double> > vertices);
                void getBetterGuess();
            
                // settings
                double simplex_size = 0.1;
        };

        class ScaleOptimization
        {
            friend TRK;
            TRK &trk; // parent class
            
            public:
                // constructors/destructor
                ScaleOptimization(TRK &trk);
                ~ScaleOptimization();
            
                // fitting scales
                double s, a, b;
            
                // settings
                bool verbose = false;
            
            private:
                // fitting scales
                whichScaleExtrema whichExtrema = NONE;
                whichScaleExtrema whichExtremaX = NONE;
                whichScaleExtrema whichExtremaY = NONE;
            
                // iterative tolerances
                double pegToZeroTol = 0.004;
            
                // guesses
                bool firstGuess = true;
                double slopYGuess;
                double slopYScaleGuess = 1.0;
                std::vector <double> iterative_allparams_guess;
            
                void getBetterSlopYGuess(double slop_y, double s);
            
                // core optimization / slop vs. scale
                std::vector <double (ScaleOptimization::*)()> optimizeList = {&ScaleOptimization::optimize_s_SlopX, &ScaleOptimization::optimize_s_SlopY};
                
                void optimizeScale();
                double optimize_s_SlopX();
                double optimize_s_SlopY();
                double innerSlopX_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
                double innerSlopY_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
            
                // TRK correlation coefficient
                double optimize_s0_R2();
                double optimize_s_prime_R2(double s0);
                double iterateR2_OptimumScale(double s0);
                double innerR2_iter_Simplex(std::vector <double> ss, std::vector <double> allparams_guess, double s0);
                double innerR2_Simplex(std::vector <double> ss, std::vector <double> allparams_guess);
                double R2TRK_prime_as();
                double R2TRK_prime_sb();
                double R2TRK_prime_as0(double s0, std::vector <double> x_t_s0, std::vector <double> params_s0);
                double R2TRK_prime_s0b(double s0, std::vector <double> x_t_s0, std::vector <double> params_s0);
        };

        class TangentPointMethods
        {
            friend TRK;
            TRK &trk; // parent class
            
            public:
                // constructors/destructor
                TangentPointMethods(TRK &trk);
                ~TangentPointMethods();
            
            private:
                // settings
                double root_bound = 10;
                
                // root finders
                double twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1);
                double newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess);
                
                // tangent point choosing
                double findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec, double s);
                std::vector <double> tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg);
            
                // approximations
                std::vector <double> approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xr1);
                std::vector <double> tangentCubicSolver(double A, double B, double C, double D);
            
                // parallel computing methods
                std::vector <double> tangentParallel(std::vector<double> params, double slop_x, double slop_y, int n, double s);
                std::vector <double> tangentParallelLikelihood(std::vector<double> params, double slop_x, double slop_y, int n);
        };

        class MCMC
        {
            friend TRK;
            TRK &trk; // parent class
            
            public:
                // constructors/destructor
                MCMC(TRK &trk);
                ~MCMC();
            
                // settings
                int R = 50000; // MCMC sample size (excluding burn-in)
                int burncount = 5000; // MCMC "burn in" count
                bool verbose = false; // print MCMC steps
            
            private:
                // settings
                samplingMethod thisSamplingMethod = AIES;
                bool useLogPosterior = false;
                bool currentlyOptimizingProposal = false;
                bool parallelizeAIES = false;
                double best_ratio = 0.325;
                double simplexSuperShrink = 1e-3;
                
                // uncertainty estimation
                void calculateUncertainties();
                std::vector <double> allParamsFinalDeltas;
                std::vector <std::vector <std::vector <double> > > lowerBar(std::vector <std::vector <double> > allparam_samples);
                
                // sampling (general)
                double metHastRatio(std::vector <double> X_trial, std::vector <double> X_i);
                std::vector <std::vector <double >> samplePosterior(int R, int burncount, std::vector <double> sigmas_guess);
                
                // Affine Invariant Ensemble Sampler (AIES)
                bool printAIESWalkerEvolution = false;
                int amt_walkers = 2; // number of walkers is amt_walkers * (number of parameters sampled)
                double AIES_a = 2.0;
                double AIES_initial_scaling = 0.1;
                std::vector <double> updateAIESWalker(std::vector <double> X, std::vector <std::vector <double> > YY);
                std::vector <double> parallelUpdateAIESWalkers(std::vector <std::vector <double> > XX, std::vector <std::vector <double> > YY, int k);
                
                // Adaptive Random Walk Metropolis Hastings (ARWMH) sampler
                void guessARWMHDeltas();
                std::vector <double> params_sigmas_guess, allparams_sigmas_guess;
                double slop_x_sigma_guess, slop_y_sigma_guess, slop_x_minus_sigma_guess, slop_y_minus_sigma_guess;
            
                // RNG
                double rnorm(double mu, double sig);
                double runiform(double a, double b);
                double rstretch(double a);

                // tools
                std::vector <double> pegToNonZeroDelta(std::vector <double> vertex, std::vector <double> lastvertex);
                std::vector <std::vector <double >> checkSlopSignMCMC(std::vector <std::vector <double >> result_final);
        };

        class CorrelationRemoval
        {
            friend TRK;
            TRK &trk; // parent class
            
            public:
                // constructors/destructor
                CorrelationRemoval(TRK &trk);
                ~CorrelationRemoval();
            
            
                // general
                static std::vector <double> pivots; // pivot point(s) themselves; to be fixed pre-fitting if desired
                std::vector <double(*)(std::vector <double>)> pivot_intercept_functions; // vector of functions that return intercept parameters(s) of corresponding pivot points
                std::vector <double(*)(std::vector <double>)> pivot_slope_functions; // same, but for slope parameter(s)
                
            
                // settings
                int maxPivotIter = 100; // Maximum number of iterations for the pivot point finder. 5 is usually sufficient, as successive iterations seem to only jump around (may not be true for linearIZED models, not just linear, however)
                bool findPivotPoints = false;
                bool verbose = false; // show pivot point finding steps
            
            
                // testing/experimental
                bool writePivots = false; //outputs pivot point sampling results for single-pivot case
                bool get_pivot_guess = false;
                bool verbose_refit = false;
                bool verbose_pearson = false;
                bool parallelize = false;
            
            private:
                // core
                void findPivots();
            
                // settings
                pivotPointFindingMethod thisPivotMethod = PEARSON;
            
                bool refit_newPivot = true;
                bool refit_with_simplex = true;
            
                bool getCombosFromSampleDirectly = true;
                bool weightPivots = true;
                bool pivotMedian = false;
                bool pivotMean = false;
                bool pivotHalfSampleMode= false;
                bool pruneOutlierPivots = true;
                bool pivotPointActive = false;
            
                int P = 0; // number of pivot points
                int sample_R = 5000; //1000 too low; 5000 seems sufficient, but 10,000 works for sure
                int randomSampleCount = 450;
                int maxCombos = 10000; // 50,000 seems sufficient, but 100,000 works for sure
                int sample_burnIn = 5000;
            
                double tol = 1e-3;
                double pruneWidth = 10.0;
            
                // combinations
                std::vector < std::vector <double> > NDcombination;
                std::vector < std::vector <std::vector <double > > > NDcombos;
                void getCombos(std::vector <std::vector <double> > total, int k, int offset);
                std::vector < std::vector <std::vector <double > > > directCombos(std::vector < std::vector <double> > params_sample, int comboCount);

            
                // guesses
                std::vector <int> intercept_indices; // indices of intercepts and slopes in vector of parameters describing model
                std::vector <int> slope_indices;
                void getPivotGuess();
                void findLinearParamIndices();
                void refitWithNewPivots(std::vector <double> new_pivots);
                std::vector <double> refitAnalytic(std::vector <double> new_pivots);
            
            
                // find pivots by generating a distribution
                void optimizePivots_Distribution();
                double weightPivot(std::vector <double> params1, std::vector <double> params2, std::vector <double> oldPivots, double newPivot, int p);
                double pivotFunc(std::vector <double> params1, std::vector <double> params2, int p);
                std::vector <double> removeOutlierPivots(std::vector <double> pivots);
            
            
                // find pivots using the slope of the correlation ellipse between intercepts and slopes
                void optimizePivots_Regression();
            
            
                // find pivots using the correlation of intercepts and slopes and the golden section search (GSS) method
                std::vector <double> max_pivots_brackets, min_pivots_brackets;
            
                void optimizePivots_Pearson();
                void optimizePivots_Pearson_Loop(); // both serial and parallel
                void writePearsonOptimizationSampling(std::vector <double> b_samples, std::vector <double> m_samples, int iter, int p);
                // std::vector <double> parallelGoldenSectionSearch(bool c_or_d, double pivot, int p, int iter);
                bool checkPearsonOptimizationTolerance(std::vector <double> previous_pivots);
                double getAbsPearsonCorrFromNewPivot(double new_pivot, int p, int iter);
            
            
                // tools
                std::vector <double> findNTiles(int Q);
                void findPivotBrackets();
        };
    
        class Asymmetric    // Asymmetric Uncertainties
        {
            friend TRK;
            TRK &trk; // parent class
            
            public:
                // constructors/destructor
                Asymmetric(TRK &trk);
                ~Asymmetric();
            
                // dataset
                std::vector <double> sx_minus, sy_minus; // negative asymmetric error bars
                
                // guesses
                double slop_x_minus_guess = -1.0;  // negative asymmetric slop
                double slop_y_minus_guess = -1.0;
            
            
            private:
                // settings
                bool hasAsymEB = false;
                bool hasAsymSlop = false;
                bool verbose = false; // show info/steps about asymmetric uncertainty fitting
            
                // tools
                void checkAsym();
            
                // likelihoods/posteriors
                double singlePointLnLAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, double x_tn, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn);
                double pnAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, double x_tn, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn);
                double dunDxAsym(double mtn, std::vector <double> Sigs2, int quadSig_xn2Ind, int quadSig_yn2Ind, double s);
                double zAsym(double x, double quadSig_xn2, double quadSig_yn2, double xn_shifted, double yn_shifted, std::vector <double> shifts, double x_tn, double y_tn, double m_tn);
            
                // tangent points
                std::vector <double> tangentParallelAsym(std::vector<double> allparams, int n, double s);
                std::vector <double> tangentParallelLikelihoodAsym(std::vector<double> allparams, int n);
                double findBestTangentAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, std::vector <double> x_tn_vec, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn);
                std::vector <double> getAsymShifts(std::vector <double> allparams, int n);
                std::vector <double> getAsymSigs2(std::vector <double> allparams, int n);
        };

        class Settings
        {
            friend TRK;
            TRK &trk; // parent class
            
            public:
                // constructors/destructor
                Settings(TRK &trk);
                ~Settings();
            
                // general settings
                ParallelizationBackEnd ParallelizationBackEnd = CPP11;
                int maxThreads = 8; // maximum threads for parallel processing
            
                // output settings
                bool printResults = true;
                bool outputDistributionToFile = false; // outputs MCMC-sampled parameter distribution to a text file (must specify path below)
                std::string outputPath;
                bool verbose = false; // turns on maximum verbosity for all portions of code (see e.g. MCMC.verbose and CorrelationRemoval.verbose for more specific settings)
            
            
            private:
                // general
                bool do1DFit = false; // does not need to be changed by the user, do 1D fits simply using the correct TRK constructor. Only used for testing
        };

        class COVID19 // for covid19 fits
        {
            friend TRK;
            TRK &trk; // parent class
            
            public:
                // constructors/destructor
                COVID19(TRK &trk);
                ~COVID19();
            
                // general
                bool covid19fit = false;
                static double y12;
                static bool logModel;
                static double s;
                static double t_split;
                static double tmed;
                static std::vector <double> fixed_params;
                static std::vector <double> fixed_pivots;
        };
    
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
    

		// default constructor and destructor;
		TRK();
        ~TRK();
    
    
		// dataset
        unsigned long N, M; // number of data points and number of model parameters, respectively (not including slop)
		std::vector <double> x, y, sx, sy, w; //datapoints; errorbars
    
        // modeling
        std::vector <double> fixed_allparams;

		// core algorithms
		void performTRKFit(); //finds optimum scale AND calculates uncertainties
		void performTRKFit(double scale); //perform fit on some provided scale (for example, if they already know optimum scale, they can just start with this) and calculates uncertainties
		void performSimpleTRKFit(); //finds optimum scale and and performs TRK fit but without finding uncertainties
        void performSimpleTRKFit(double scale); //given some scale, performs TRK fit without finding uncertainties.

		// results
		Results results;
        
        // instances of nested classes
        Statistics statistics;
        Optimization optimization;
        ScaleOptimization scaleOptimization;
        TangentPointMethods tangentPointMethods;
        MCMC mcmc;
        CorrelationRemoval correlationRemoval;
        Asymmetric asymmetric;
        Settings settings;
        COVID19 covid19;

    
    private:
    
        // dataset
        void getDataWidth();
        void checkZeroErrorBars();
    
        unsigned long bigM; //bigM is total number of model AND slop params
        double datawidth, x_min, x_max;
        std::vector <double> SigXVec, SigYVec;
    
    
        // guesses
        double slop_x_guess, slop_y_guess;
        std::vector <double> x_t_slopx, x_t_slopy, x_t_a, x_t_b, x_t_s;
        std::vector <double> params_slopx, params_slopy, params_a, params_b, params_s, allparams_s, iterative_allparams_guess;
        std::vector <double> params_guess;
        std::vector <double> allparams_guess;
    
        
        //function pointers
        double (*yc)(double, std::vector <double>);
        double (*dyc)(double, std::vector <double>);
        double (*ddyc)(double, std::vector <double>);
    

        // OTHER TOOLS
        void showResults(bool didScaleOp, bool didMCMC);
        void checkVerbose();
        bool isEqual(double x, double y, double maxRelativeError, double maxAbsoluteError);
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

// non-object/global tools:

// mcmc/sampling
double noPrior(double param);

// statistics
std::vector <double> minMax(std::vector <double> vec);
std::vector <int> argMinMax(std::vector <double> x);
std::vector <int> getSortedIndices(std::vector <double> x);
double getMedian(std::vector<double> y);
double getMedian(int trueCount, std::vector<double> w, std::vector<double> y);

// numerical methods/optimization
double twoPointNR(double(*y)(double, std::vector <double>), double(*dy)(double, std::vector <double>), double(*ddy)(double, std::vector <double>), std::vector <double> params, double xguess, double xguessp1);
std::vector <double> cubicSolver(double A, double B, double C, double D);

// for testing only

clock_t startTimer();

double secElapsed(clock_t t_i);

void writeResults(TRK TRKobj, double t_sec, std::string filename);

double toRad(double deg);

double toDeg(double rad);

std::vector <std::vector <double > > getData(std::string fileName, int dataSize);

std::vector <std::vector <double > > getData(std::string fileName);
