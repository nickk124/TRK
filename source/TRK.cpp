/*
 Trotter Reichart Konz (TRK) Official Codebase
 Author: Nick C. Konz
 See license at https://github.com/nickk124/TRK
 */

#include "TRK.h"

// STATIC VARIABLES ##########################################################################################################

std::vector <double> TRK::CorrelationRemoval::pivots = {};
double TRK::COVID19::y12 = 2.55;
bool TRK::COVID19::logModel = false;
double TRK::COVID19::s = 1.0;
double TRK::COVID19::t_split = 45.5;
double TRK::COVID19::tmed = 1.0;
std::vector <double> TRK::COVID19::fixed_params = {};
std::vector <double> TRK::COVID19::fixed_pivots = {};

// ###########################################################################################################################



// PRIORS ####################################################################################################################

// Constructors
Priors::Priors(priorTypes priorType, std::vector < std::vector <double> > params) { //Gaussian or bounded priors only
	this->priorType = priorType;
	if (priorType == GAUSSIAN) {
		this->gaussianParams = params;
	}
	else if (priorType == CONSTRAINED) {
		this->paramBounds = params;
	}
};

Priors::Priors(priorTypes priorType, std::vector < std::vector <double> > gaussianParams, std::vector < std::vector <double> > paramBounds) { //mixed
	this->priorType = priorType;
	this->paramBounds = paramBounds;
	this->gaussianParams = gaussianParams;
};

Priors::Priors(priorTypes priorType, std::vector <double(*)(double)> priorsPDFs) { //custom priors
	this->priorType = priorType;
	this->priorsPDFs = priorsPDFs;
};


//default constructor
Priors::Priors() {

};

// ###########################################################################################################################



// CONSTRUCTORS ##############################################################################################################

TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess) : optimization(*this), asymmetric(*this), scaleOptimization(*this), tangentPointMethods(*this), statistics(*this), mcmc(*this), correlationRemoval(*this), settings(*this), covid19(*this) {
	this->yc = (*yc);
	this->dyc = (*dyc);
	this->ddyc = (*ddyc);

	this->x = x;
	this->y = y;
	this->w = w;
	this->sx = sx;
	this->sy = sy;

	this->params_guess = params_guess;
	this->slop_x_guess = slop_x_guess;
	this->slop_y_guess = slop_y_guess;

	this->N = x.size();
	this->M = params_guess.size();

	std::vector <double> allparams_guess = params_guess;

	allparams_guess.push_back(slop_x_guess);
	allparams_guess.push_back(slop_y_guess);

	this->allparams_guess = allparams_guess;
    
    this->bigM = allparams_guess.size();

    this->scaleOptimization.whichExtrema = ANY;

    this->scaleOptimization.s = 1.0;

	getDataWidth();

    this->statistics.hasPriors = false;
    
    checkZeroErrorBars();
    
    this->fixed_allparams = std::vector <double>(bigM, NAN);
    
    mcmc.guessARWMHDeltas();
}


//equal weights/unweighted
TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess) : optimization(*this), asymmetric(*this), scaleOptimization(*this), tangentPointMethods(*this), statistics(*this), mcmc(*this), correlationRemoval(*this), settings(*this), covid19(*this) {
	this->yc = (*yc);
	this->dyc = (*dyc);
	this->ddyc = (*ddyc);

	this->x = x;
	this->y = y;
	this->sx = sx;
	this->sy = sy;

	this->params_guess = params_guess;
	this->slop_x_guess = slop_x_guess;
	this->slop_y_guess = slop_y_guess;

	this->N = x.size();
	this->M = params_guess.size();

	std::vector <double> w(N, 1.0); //no weights === equal weights
	this->w = w;

	std::vector <double> allparams_guess = params_guess;

	allparams_guess.push_back(slop_x_guess);
	allparams_guess.push_back(slop_y_guess);

	this->allparams_guess = allparams_guess;
    
    this->bigM = allparams_guess.size();

	this->scaleOptimization.whichExtrema = ANY;

	this->scaleOptimization.s = 1.0;

	getDataWidth();

    this->statistics.hasPriors = false;
    
    checkZeroErrorBars();
    
    this->fixed_allparams = std::vector <double>(bigM, NAN);
    
    mcmc.guessARWMHDeltas();
}


//priors:
//weighted
TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, Priors priorsObject) : optimization(*this), asymmetric(*this), scaleOptimization(*this), tangentPointMethods(*this), statistics(*this), mcmc(*this), correlationRemoval(*this), settings(*this), covid19(*this) {
	this->yc = (*yc);
	this->dyc = (*dyc);
	this->ddyc = (*ddyc);

	this->x = x;
	this->y = y;
	this->w = w;
	this->sx = sx;
	this->sy = sy;

	this->params_guess = params_guess;
	this->slop_x_guess = slop_x_guess;
	this->slop_y_guess = slop_y_guess;

	this->N = x.size();
	this->M = params_guess.size();

	std::vector <double> allparams_guess = params_guess;

	allparams_guess.push_back(slop_x_guess);
	allparams_guess.push_back(slop_y_guess);

	this->allparams_guess = allparams_guess;
    
    this->bigM = allparams_guess.size();

	this->scaleOptimization.whichExtrema = ANY;

	this->scaleOptimization.s = 1.0;

	getDataWidth();

    this->statistics.priorsObject = priorsObject;

    this->statistics.hasPriors = true;
    
    checkZeroErrorBars();
    
    this->fixed_allparams = std::vector <double>(bigM, NAN);
    
    mcmc.guessARWMHDeltas();
}


//equal weights/unweighted
TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, Priors priorsObject) : optimization(*this), asymmetric(*this), scaleOptimization(*this), tangentPointMethods(*this), statistics(*this), mcmc(*this), correlationRemoval(*this), settings(*this), covid19(*this) {
	this->yc = (*yc);
	this->dyc = (*dyc);
	this->ddyc = (*ddyc);

	this->x = x;
	this->y = y;
	this->sx = sx;
	this->sy = sy;

	this->params_guess = params_guess;
	this->slop_x_guess = slop_x_guess;
	this->slop_y_guess = slop_y_guess;

	this->N = x.size();
	this->M = params_guess.size();

	std::vector <double> w(N, 1.0); //no weights === equal weights
	this->w = w;

	std::vector <double> allparams_guess = params_guess;

	allparams_guess.push_back(slop_x_guess);
	allparams_guess.push_back(slop_y_guess);

	this->allparams_guess = allparams_guess;
    
    this->bigM = allparams_guess.size();

	this->scaleOptimization.whichExtrema = ANY;

	this->scaleOptimization.s = 1.0;

	getDataWidth();

    this->statistics.priorsObject = priorsObject;

    this->statistics.hasPriors = true;
    
    checkZeroErrorBars();
    
    this->fixed_allparams = std::vector <double>(bigM, NAN);
    
    mcmc.guessARWMHDeltas();
}


// 1D statistic
TRK::TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess) : optimization(*this), asymmetric(*this), scaleOptimization(*this), tangentPointMethods(*this), statistics(*this), mcmc(*this), correlationRemoval(*this), settings(*this), covid19(*this) {
    this->settings.do1DFit = true;
    
    this->yc = (*yc);

    this->x = x;
    this->y = y;
    this->w = w;
    this->sy = sy;

    this->params_guess = params_guess;
    this->slop_y_guess = slop_y_guess;

    this->N = x.size();
    this->M = params_guess.size();

    std::vector <double> allparams_guess = params_guess;
    
    allparams_guess.push_back(slop_y_guess);

    this->allparams_guess = allparams_guess;
    
    this->bigM = allparams_guess.size();

    this->scaleOptimization.whichExtrema = ANY;

    this->scaleOptimization.s = 1.0;

    getDataWidth();

    this->statistics.hasPriors = false;
    
    checkZeroErrorBars();
    
    this->fixed_allparams = std::vector <double>(bigM, NAN);
    
    mcmc.guessARWMHDeltas();
}


//equal weights/unweighted
TRK::TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess) : optimization(*this), asymmetric(*this), scaleOptimization(*this), tangentPointMethods(*this), statistics(*this), mcmc(*this), correlationRemoval(*this), settings(*this), covid19(*this) {
    this->settings.do1DFit = true;
    
    this->yc = (*yc);

    this->x = x;
    this->y = y;
    this->sy = sy;

    this->params_guess = params_guess;
    this->slop_y_guess = slop_y_guess;

    this->N = x.size();
    this->M = params_guess.size();

    std::vector <double> w(N, 1.0); //no weights === equal weights
    this->w = w;

    std::vector <double> allparams_guess = params_guess;

    allparams_guess.push_back(slop_y_guess);

    this->allparams_guess = allparams_guess;
    
    this->bigM = allparams_guess.size();

    this->scaleOptimization.whichExtrema = ANY;

    this->scaleOptimization.s = 1.0;

    getDataWidth();

    this->statistics.hasPriors = false;
    
    checkZeroErrorBars();
    
    this->fixed_allparams = std::vector <double>(bigM, NAN);
    
    mcmc.guessARWMHDeltas();
}


//priors:
//weighted
TRK::TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess, Priors priorsObject) : optimization(*this), asymmetric(*this), scaleOptimization(*this), tangentPointMethods(*this), statistics(*this), mcmc(*this), correlationRemoval(*this), settings(*this), covid19(*this) {
    this->settings.do1DFit = true;
    
    this->yc = (*yc);

    this->x = x;
    this->y = y;
    this->w = w;
    this->sy = sy;

    this->params_guess = params_guess;
    this->slop_y_guess = slop_y_guess;

    this->N = x.size();
    this->M = params_guess.size();

    std::vector <double> allparams_guess = params_guess;

    allparams_guess.push_back(slop_y_guess);
    
    this->bigM = allparams_guess.size();

    this->allparams_guess = allparams_guess;

    this->scaleOptimization.whichExtrema = ANY;

    this->scaleOptimization.s = 1.0;

    getDataWidth();

    this->statistics.priorsObject = priorsObject;

    this->statistics.hasPriors = true;
    
    checkZeroErrorBars();
    
    this->fixed_allparams = std::vector <double>(bigM, NAN);
    
    mcmc.guessARWMHDeltas();
}


//equal weights/unweighted
TRK::TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess, Priors priorsObject) : optimization(*this), asymmetric(*this), scaleOptimization(*this), tangentPointMethods(*this), statistics(*this), mcmc(*this), correlationRemoval(*this), settings(*this), covid19(*this) {
    this->settings.do1DFit = true;
    
    this->yc = (*yc);

    this->x = x;
    this->y = y;
    this->sy = sy;

    this->params_guess = params_guess;
    this->slop_y_guess = slop_y_guess;

    this->N = x.size();
    this->M = params_guess.size();

    std::vector <double> w(N, 1.0); //no weights === equal weights
    this->w = w;

    std::vector <double> allparams_guess = params_guess;

    allparams_guess.push_back(slop_y_guess);

    this->allparams_guess = allparams_guess;
    
    this->bigM = allparams_guess.size();

    this->scaleOptimization.whichExtrema = ANY;

    this->scaleOptimization.s = 1.0;

    getDataWidth();

    this->statistics.priorsObject = priorsObject;

    this->statistics.hasPriors = true;
    
    checkZeroErrorBars();
    
    this->fixed_allparams = std::vector <double>(bigM, NAN);
    
    mcmc.guessARWMHDeltas();
}


// default
TRK::TRK() : optimization(*this), asymmetric(*this), scaleOptimization(*this), tangentPointMethods(*this), statistics(*this), mcmc(*this), correlationRemoval(*this), settings(*this), covid19(*this) {

}


// destructor
TRK::~TRK() {

}


// nested class constructors/destructors
TRK::Asymmetric::Asymmetric(TRK &trk) : trk(trk) {};
TRK::Asymmetric::~Asymmetric(){};

TRK::ScaleOptimization::ScaleOptimization(TRK &trk) : trk(trk) {};
TRK::ScaleOptimization::~ScaleOptimization(){};

TRK::TangentPointMethods::TangentPointMethods(TRK &trk) : trk(trk) {};
TRK::TangentPointMethods::~TangentPointMethods(){};

TRK::Statistics::Statistics(TRK &trk) : trk(trk) {};
TRK::Statistics::~Statistics(){};

TRK::Optimization::Optimization(TRK &trk) : trk(trk) {};
TRK::Optimization::~Optimization(){};

TRK::MCMC::MCMC(TRK &trk) : trk(trk) {};
TRK::MCMC::~MCMC(){};

TRK::CorrelationRemoval::CorrelationRemoval(TRK &trk) : trk(trk) {};
TRK::CorrelationRemoval::~CorrelationRemoval(){};

TRK::Settings::Settings(TRK &trk) : trk(trk) {};
TRK::Settings::~Settings(){};

TRK::COVID19::COVID19(TRK &trk) : trk(trk) {};
TRK::COVID19::~COVID19(){};

// ###########################################################################################################################



// CORE ALGORITHMS ###########################################################################################################

void TRK::performTRKFit() {//finds optimum scale AND performs TRK fit + uncertainty
    checkVerbose();
    
    asymmetric.checkAsym();
    
    optimization.getBetterGuess(); // once before pivot point optimization, once after
    
    correlationRemoval.getPivotGuess();
    
    scaleOptimization.optimizeScale();

    correlationRemoval.findPivots();

    optimization.getBetterGuess();
    
    mcmc.calculateUncertainties();
    
    showResults(true, true);
}

void TRK::performTRKFit(double scale) {//perform fit on some provided scale (for example, if they already know optimum scale, they can just start with this) and calculates uncertainties
    checkVerbose();
    
    asymmetric.checkAsym();
    
    scaleOptimization.s = scale;
    results.optimumScale = scale;
    
    optimization.getBetterGuess();
    
    correlationRemoval.getPivotGuess();

    correlationRemoval.findPivots();

    results.bestFitParams.clear();

    scaleOptimization.whichExtrema = S;
    allparams_s = optimization.downhillSimplex_Fit(statistics.selectedChiSq, allparams_guess, scale, optimization.showFittingSteps);
    scaleOptimization.whichExtrema = ANY;

    for (int j = 0; j < M; j++) {
        results.bestFitParams.push_back(allparams_s[j]);
    }

    results.slop_x = allparams_s[M];
    results.slop_y = allparams_s[M + 1];
    
    if (asymmetric.hasAsymSlop){
        results.slop_x_minus = allparams_s[M + 2];
        results.slop_y_minus = allparams_s[M + 3];
    }

    optimization.getBetterGuess();
    mcmc.calculateUncertainties();
    
    showResults(false, true);
}

void TRK::performSimpleTRKFit() {//finds optimum scale and performs TRK fit but without finding uncertainties
    checkVerbose();
    
    asymmetric.checkAsym();
    
    optimization.getBetterGuess();
    
    correlationRemoval.getPivotGuess();
    
    scaleOptimization.optimizeScale(); // (stores results in TRK.results)

    correlationRemoval.findPivots();
    
    results.bestFitParams.clear();
    
    scaleOptimization.whichExtrema = S;
    allparams_s = optimization.downhillSimplex_Fit(statistics.selectedChiSq, allparams_guess, scaleOptimization.s, optimization.showFittingSteps);
    scaleOptimization.whichExtrema = ANY;
    
    for (int j = 0; j < M; j++) {
        results.bestFitParams.push_back(allparams_s[j]);
    }
    
    results.slop_x = allparams_s[M];
    results.slop_y = allparams_s[M + 1];
    
    if (asymmetric.hasAsymSlop){
        results.slop_x_minus = allparams_s[M + 2];
        results.slop_y_minus = allparams_s[M + 3];
    }
    
    showResults(true, false);

    return;
}

void TRK::performSimpleTRKFit(double scale) {//given some provided scale, performs TRK fit but without finding uncertainties
    checkVerbose();
    
    asymmetric.checkAsym();
    
    scaleOptimization.s = scale;
    
    optimization.getBetterGuess();

    correlationRemoval.getPivotGuess();

    correlationRemoval.findPivots();

    results.bestFitParams.clear();

    scaleOptimization.whichExtrema = S;
    allparams_s = optimization.downhillSimplex_Fit(statistics.selectedChiSq, allparams_guess, scaleOptimization.s, optimization.showFittingSteps);
    scaleOptimization.whichExtrema = ANY;
    
    for (int j = 0; j < M; j++) {
        results.bestFitParams.push_back(allparams_s[j]);
    }
    
    results.slop_x = allparams_s[M];
    results.slop_y = allparams_s[M + 1];
    
    if (asymmetric.hasAsymSlop){
        results.slop_x_minus = allparams_s[M + 2];
        results.slop_y_minus = allparams_s[M + 3];
    }
    
    showResults(false, false);
    
    return;
}

// ###########################################################################################################################



// MISCELLANEOUS #############################################################################################################

// dataset
void TRK::getDataWidth() {
    std::vector <double> bounds = minMax(x);

    datawidth = std::abs(bounds[1] - bounds[0]);

    x_min = bounds[0];
    x_max = bounds[1];
    
    return;
}

void TRK::checkZeroErrorBars(){
    std::vector <int> badInds;
    bool check = false;
    if (settings.do1DFit){
        for (int i = 0; i < N; i++){
            if (sy[i] == 0){
                badInds.push_back(i);
                check = true;
            }
        }
        if (check){
            printf("Warning: Error bars of zero found at data points (indexing starting at 1) of:\n");
            for (int j = 0; j < (int)badInds.size(); j++){
                printf(" %i,", badInds[j]);
            }
            printf("\n Could create likelihood overflow errors if slop also zero/close to zero.");
        }
        
    } else {
        for (int i = 0; i < N; i++){
            if (sx[i] == 0 || sy[i] == 0){
                badInds.push_back(i);
                check = true;
            }
        }
        if (check){
            printf("Warning: Error bars of zero found at data points (indexing starting at 1) of:\n");
            for (int j = 0; j < (int)badInds.size(); j++){
                printf(" %i,", badInds[j]);
            }
            printf("\n Could create likelihood overflow errors if slop also zero/close to zero.");
        }
    }
    return;
}


// other tools
bool TRK::isEqual(double x, double y, double maxRelativeError = .00000001, double maxAbsoluteError = DBL_MIN)// .000001; .0000001;.00000001
{
    if (std::abs(x - y) < maxAbsoluteError)
    {
        return true;
    }
    double relativeError = (std::abs(y) > std::abs(x) ? std::abs((x - y) / y) : std::abs((x - y) / x));
    if (relativeError <= maxRelativeError)
    {
        return true;
    }
    return false;
}

void TRK::checkVerbose(){
    if (settings.verbose){
        mcmc.verbose = true;
        correlationRemoval.verbose = true;
        asymmetric.verbose = true;
        scaleOptimization.verbose = true;
    }
    return;
}

void TRK::showResults(bool didScaleOp, bool didMCMC){
    printf("\n\nTRK RESULTS:\n\n");
    
    if (didScaleOp && !settings.do1DFit){
        printf("SCALES:\nOptimum scale: %.3e \n", results.optimumScale);
        printf("Minimum scale: %.3e \n", results.minimumScale);
        printf("Maximum scale: %.3e \n\n", results.maximumScale);
    }
    
    if (correlationRemoval.P > 0){
        printf("PIVOT POINTS:\n");
        for (int p = 0; p < correlationRemoval.P; p++){
            printf("Pivot Point %i = %.3e\n", p, results.pivots[p]);
        }
        printf("\n\n");
    }

    printf("BEST FIT MODEL PARAMETERS (including slop):\n");
    for (int k = 0; k < params_guess.size(); k++) {
        printf("%.3e ", results.bestFitParams[k]);
    }
    printf("\n\nBEST FIT SLOP (EXTRINSIC SCATTER) PARAMETERS:\n");
    if (settings.do1DFit){
        printf("y-slop =  %.3e", results.slop_y);
    } else {
        printf("x-slop =  %.3e y-slop = %.3e", results.slop_x, results.slop_y);
    }
    std::cout << std::endl;
    
    if (didMCMC){
        printf("\nUNCERTAINTIES: (- 1 2 3, + 1 2 3): \n");
        for (int k = 0; k < params_guess.size(); k++) { //kth param
            for (int j = 0; j < 2; j++) { // - and + sigmas
                for (int i = 0; i < 3; i++) { // 1, 2 and 3 sigmas
                    printf("%.3e ", results.bestFit_123Sigmas[k][j][i]);
                }
                printf("\t");
            }
            std::cout << std::endl;
        }

        printf("\nSLOP UNCERTAINTIES: (- 1 2 3, + 1 2 3): \n");
        if (!settings.do1DFit){
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < 3; i++) {
                    printf("%.3e ", results.slopX_123Sigmas[j][i]);
                }
                printf("\t");
            }
            std::cout << std::endl;
        }
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 3; i++) {
                printf("%.3e ", results.slopY_123Sigmas[j][i]);
            }
            printf("\t");
        }
        std::cout << std::endl << std::endl;
    }
    
    printf("\nFITNESS:\nchisquared = %.3e\n\n", results.chisquared);
    
    return;
}

// ###########################################################################################################################



// STATISTICS ################################################################################################################

// priors
double TRK::Statistics::priors(std::vector <double> allparams) {
    double jointPrior = 1.0; //uninformative prior by default

    switch (priorsObject.priorType){
        case CONSTRAINED:
            for (int i = 0; i < trk.M; i++) { //check upper bound
                double ub = priorsObject.paramBounds[i][1];
                double lb = priorsObject.paramBounds[i][0];

                if (allparams[i] >= ub && !std::isnan(ub)) {
                    jointPrior = 0.0;
                    break;
                }

                if (allparams[i] <= lb && !std::isnan(lb)) {
                    jointPrior = 0.0;
                    break;
                }
            }

            break;

        case GAUSSIAN:
            for (int i = 0; i < trk.M; i++) { //check upper bound
                double mu = priorsObject.gaussianParams[i][0];
                double sig = priorsObject.gaussianParams[i][1];

                if (std::isnan(mu) || std::isnan(sig)) { //no prior for this parameter
                    return 1.0;
                } else {
                    jointPrior *= normal(allparams[i], mu, sig);
                }
            }

            break;

        case MIXED:

            for (int i = 0; i < trk.M; i++) { //check upper bound
                double ub = priorsObject.paramBounds[i][1];
                double lb = priorsObject.paramBounds[i][0];
                double mu = priorsObject.gaussianParams[i][0];
                double sig = priorsObject.gaussianParams[i][1];

                if (allparams[i] >= ub && !std::isnan(ub)) {
                    jointPrior = 0.0;
                    break;
                }

                if (allparams[i] <= lb && !std::isnan(lb)) {
                    jointPrior = 0.0;
                    break;
                }

                if (!std::isnan(mu) && !std::isnan(sig)) { //gaussian prior for this parameter
                    jointPrior *= normal(allparams[i], mu, sig);
                }

            }

            break;

        case CUSTOM:
            for (int i = 0; i < trk.M; i++) {
                jointPrior *= priorsObject.priorsPDFs[i](allparams[i]);
            }

            break;
        
        case CUSTOM_JOINT:
            jointPrior = priorsObject.jointPriorsPDF(allparams);
            break;
    }

    return jointPrior;
}


// regression
std::vector <double> TRK::Statistics::simpleLinearRegression(std::vector <double> x, std::vector <double> y){ // least-squares solution for b, m of model y = mx + b
    
    double b, m, xbar, ybar, upper_sum = 0.0, lower_sum = 0.0;
    xbar = getAverage(x); ybar = getAverage(y);
    
    for (int i = 0; i < trk.N; i++){
        upper_sum += (x[i] - xbar) * (y[i] - ybar);
        lower_sum += std::pow(x[i] - xbar, 2.0);
    }
    
    m = upper_sum / lower_sum;
    b = ybar - m * xbar;
    return {b, m};
}


// correlation
double TRK::Statistics::pearsonCorrelation(std::vector<double> x, std::vector<double> y){
    double uppersum = 0.0, lowersum1 = 0.0, lowersum2 = 0.0, x_bar = getAverage(x), y_bar = getAverage(y);
    
    for (int i = 0; i < (int) x.size(); i++){
        uppersum += (x[i] - x_bar) * (y[i] - y_bar);
        lowersum1 += std::pow((x[i] - x_bar), 2.0);
        lowersum2 += std::pow((y[i] - y_bar), 2.0);
    }
    
    return uppersum / (std::sqrt(lowersum1) * std::sqrt(lowersum2));
}

double TRK::Statistics::spearmanCorrelation(std::vector<double> x, std::vector<double> y){
    std::vector <double> u, v; // ranks
    u = rankVector(x);
    v = rankVector(y);
    
//    double num, denom1, denom2;
//    double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0;
//    int n = (int) u.size();
//
//    for (int i = 0; i < n; i++){
//        sum1 += u[i] * v[i];
//        sum2 += u[i];
//        sum3 += v[i];
//        sum4 += std::pow(u[i], 2.0);
//        sum6 += std::pow(v[i], 2.0);
//    }
//    sum5 = sum2;
//    sum7 = sum3;
//
//    num = n*sum1 - sum2*sum2;
//    denom1 = n * sum4 - std::pow(sum5, 2.0);
//    denom2 = n * sum6 - std::pow(sum7, 2.0);
//
//    return num/std::sqrt(denom1 * denom2);
    
    return pearsonCorrelation(u, v);
}


// goodness of fit metrics
double TRK::Statistics::regularChiSquared(std::vector <double> params) {
    double sum = 0.0;

    for (int i = 0; i < trk.N; i++) {
        sum += trk.w[i] * std::pow(trk.y[i] - (*trk.yc)(trk.x[i], params), 2.0);
    }
    return sum;
}

double TRK::Statistics::regularChiSquaredWSlop(std::vector <double> allparams, double s) {
    double sum = 0.0;
    
    double sigma = allparams[(int) allparams.size() - 1];
    
    std::vector <double> params;

    for (int i = 0; i < trk.M; i++) {
        params.push_back(allparams[i]);
    }

    for (int i = 0; i < trk.N; i++) {
        sum += trk.w[i]* std::pow((*trk.yc)(trk.x[i], params) - trk.y[i], 2)/(std::pow(trk.sy[i], 2.0) + std::pow(sigma, 2.0)) + 2.0*std::log(std::pow(std::pow(trk.sy[i], 2.0) + std::pow(sigma, 2.0), trk.w[i]/2.0));
    }
    
    return sum;
}

double TRK::Statistics::modifiedChiSquared(std::vector <double> allparams, double s)
{
    std::vector <double> SigXVec, SigYVec;
    std::vector <double> all_x_t(trk.N, 0.0);

    double sum1 = 0.0;
    double sum2 = 0.0;

    double slop_x = allparams[trk.M];
    double slop_y = allparams[trk.M + 1];

    std::vector <double> params;

    for (int i = 0; i < trk.M; i++) {
        params.push_back(allparams[i]);
    }

    switch (trk.settings.ParallelizationBackEnd){
        case OPENMP:
        {
            #pragma omp parallel for num_threads(maxThreads)
            for (int i = 0; i < trk.N; i++)
            {
                std::vector <double> results;
                results = trk.tangentPointMethods.tangentParallel(params, slop_x, slop_y, i, s); //pointer to fn run through MT, arguments to fn
                all_x_t[i] = results[0];
                sum1 += results[1];
                sum2 += results[2];
            }
            break;
        }
        case CPP11:
        {
            int counter = 0, completedThreads = 0, liveThreads = 0;
            std::vector<double> results;
            std::vector< std::future < std::vector < double > > > futureVec;
            futureVec.resize(trk.N);

            for (int i = 0; i < trk.N; i++)
            {
                futureVec[i] = std::async(std::launch::async, &TRK::TangentPointMethods::tangentParallel, trk.tangentPointMethods, params, slop_x, slop_y, i, s); //pointer to fn run through MT, arguments to fn
                counter++;
                liveThreads++;

                if (liveThreads >= trk.settings.maxThreads)
                {
                    for (int i = completedThreads; i < counter; i++)
                    {
                        results = futureVec[i].get();
                        all_x_t[i] = results[0];
                        sum1 += results[1];
                        sum2 += results[2];
                    }
                    completedThreads += liveThreads;
                    liveThreads = 0;
                }
            }
            for (int i = completedThreads; i < trk.N; i++)
            {
                results = futureVec[i].get();
                all_x_t[i] = results[0];
                sum1 += results[1];
                sum2 += results[2];
            }
            break;
        }
        default:
        {
            for (int i = 0; i < trk.N; i++)
            {
                std::vector <double> results;
                results = trk.tangentPointMethods.tangentParallel(params, slop_x, slop_y, i, trk.scaleOptimization.s); //pointer to fn run through MT, arguments to fn
                all_x_t[i] = results[0];
                sum1 += results[1];
                sum2 += results[2];
            }
            break;
        }
    }

    switch (trk.scaleOptimization.whichExtrema) {
        case ANY:
            break;
        case S:
            trk.x_t_s = all_x_t;
            trk.params_s = params;
            break;
        default:
            break;
    }

    switch (trk.scaleOptimization.whichExtremaX) {
    case ANY:
        break;
    case SLOP_X:
        trk.x_t_slopx = all_x_t;
        trk.params_slopx = params;
        break;
    default:
        break;
    }

    switch (trk.scaleOptimization.whichExtremaY) {
    case ANY:
        break;
    case SLOP_Y:
        trk.x_t_slopy = all_x_t;
        trk.params_slopy = params;
        break;
    default:
        break;
    }

    return sum1 - sum2;
}


// probability distributions
double TRK::Statistics::normal(double x, double mu, double sig) {
    return (std::exp((-0.5) * (std::pow((x - mu), 2.0) / (2.0 * std::pow(sig, 2.0)))) / std::sqrt(2.0*PI*std::pow(sig, 2.0)));
}

double TRK::Statistics::stretch_pdf(double z) {
    double a = trk.mcmc.AIES_a;
    return (z >= 1/a && z <= a) ? 1.0/std::sqrt(z) : 0.0;
}


// likelihoods and posteriors
double TRK::Statistics::singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn, double s) {
    double m_tn = trk.dyc(x_tn, params);
    double y_tn = trk.yc(x_tn, params);
    
    return std::pow(y_n - y_tn - m_tn * (x_n - x_tn), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) - std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(s*Sig_yn2, 2.0)));
}

double TRK::Statistics::likelihood(std::vector <double> allparams) {
    std::vector <double> SigXVec, SigYVec;
    std::vector <double> all_x_t(trk.N, 0.0);
    double L = 1.0;
    
    if (trk.mcmc.useLogPosterior){
        L = 0.0;
    }

    double slop_x = allparams[trk.M];
    double slop_y = allparams[trk.M + 1];

    std::vector <double> params;

    for (int i = 0; i < trk.M; i++) {
        params.push_back(allparams[i]);
    }
    
    std::vector <double> results;

    switch (trk.settings.ParallelizationBackEnd){
        case OPENMP:
        {
            #pragma omp parallel for num_threads(maxThreads)
            for (int i = 0; i < trk.N; i++)
            {
                results = trk.tangentPointMethods.tangentParallelLikelihood(params, slop_x, slop_y, i); //pointer to fn run through MT, arguments to fn
                all_x_t[i] = results[0];
                
                if (trk.mcmc.useLogPosterior){
                    L += std::log(results[1]);
                } else {
                    L *= results[1];
                }
            }
            break;
        }
        case CPP11:
        {
            int counter = 0, completedThreads = 0, liveThreads = 0;
            std::vector<double> results;
            std::vector< std::future < std::vector < double > > > futureVec;
            futureVec.resize(trk.N);

            for (int i = 0; i < trk.N; i++)
            {
                futureVec[i] = std::async(std::launch::async, &TRK::TangentPointMethods::tangentParallelLikelihood, trk.tangentPointMethods, params, slop_x, slop_y, i); //pointer to fn run through MT, arguments to fn
                // function returns 1) tangent points and 2) likelihood component
                counter++;
                liveThreads++;

                if (liveThreads >= trk.settings.maxThreads)
                {
                    for (int i = completedThreads; i < counter; i++)
                    {
                        results = futureVec[i].get();
                        all_x_t.push_back(results[0]);
                        if (trk.mcmc.useLogPosterior){
                            L += std::log(results[1]);
                        } else {
                            L *= results[1];
                        }
                    }
                    completedThreads += liveThreads;
                    liveThreads = 0;
                }
            }
            for (int i = completedThreads; i < trk.N; i++)
            {
                results = futureVec[i].get();
                all_x_t.push_back(results[0]);
                if (trk.mcmc.useLogPosterior){
                    L += std::log(results[1]);
                } else {
                    L *= results[1];
                }
            }
        }
        default:
        {
            for (int i = 0; i < trk.N; i++)
            {
                std::vector <double> results;
                results = trk.tangentPointMethods.tangentParallelLikelihood(params, slop_x, slop_y, i); //pointer to fn run through MT, arguments to fn
                all_x_t[i] = results[0];
                if (trk.mcmc.useLogPosterior){
                    L += std::log(results[1]);
                } else {
                    L *= results[1];
                }
            }
            break;
        }
    }
    return L; // returns log L = logL1 + logL2 + ... given L = L1*L2*L3... if trk.mcmc.useLogPosterior == true
}

double TRK::Statistics::likelihood1D(std::vector <double> allparams) {
//    printf("Notice: 1D likelihood (for testing) not currently configured to work with weights.\n");
    
    double L = 1.0;
    
    if (trk.mcmc.useLogPosterior){
        L = 0.0;
    }
    
    double sigma = allparams[(int) allparams.size() - 1];
    
    std::vector <double> params;

    for (int i = 0; i < trk.M; i++) {
        params.push_back(allparams[i]);
    }
    
    double l;

    for (int i = 0; i < trk.N; i++) {
        l = std::exp(-0.5 * trk.w[i]* std::pow(((*trk.yc)(trk.x[i], params) - trk.y[i])/std::sqrt(std::pow(trk.sy[i], 2.0) + std::pow(sigma, 2.0)), 2.0)) / std::pow(std::pow(trk.sy[i], 2.0) + std::pow(sigma, 2.0), trk.w[i]/2.0);
        if (trk.mcmc.useLogPosterior){
            L += std::log(l);
        } else {
            L *= l;
        }
    }
//    printf("%.3e\n",L);
    return L; // returns log L = logL1 + logL2 + ... given L = L1*L2*L3... if useLogPosterior == true
}


// likelihoods and posteriors with asymmetric uncertainties
double TRK::Statistics::modifiedChiSquaredAsym(std::vector <double> allparams, double s)
{
    std::vector <double> all_x_t(trk.N, 0.0);

    double sum = 0.0;

    std::vector <double> params;

    for (int i = 0; i < trk.M; i++) {
        params.push_back(allparams[i]);
    }
    
    std::vector <double> results;

    switch (trk.settings.ParallelizationBackEnd){
        case OPENMP:
        {
            #pragma omp parallel for num_threads(maxThreads)
            for (int i = 0; i < trk.N; i++)
            {
                results = trk.asymmetric.tangentParallelAsym(allparams, i, s); //pointer to fn run through MT, arguments to fn
                all_x_t[i] = results[0];
                sum += results[1];
            }
            break;
        }
    
        case CPP11:
        {
            int counter = 0, completedThreads = 0, liveThreads = 0;
            std::vector< std::future < std::vector < double > > > futureVec;
            futureVec.resize(trk.N);

            for (int i = 0; i < trk.N; i++)
            {
                futureVec[i] = std::async(std::launch::async, &TRK::Asymmetric::tangentParallelAsym, trk.asymmetric, allparams, i, s); //pointer to fn run through MT, arguments to fn
                counter++;
                liveThreads++;

                if (liveThreads >= trk.settings.maxThreads)
                {
                    for (int i = completedThreads; i < counter; i++)
                    {
                        results = futureVec[i].get();
                        all_x_t[i] = results[0];
                        sum += results[1];
                    }
                    completedThreads += liveThreads;
                    liveThreads = 0;
                }
            }
            for (int i = completedThreads; i < trk.N; i++)
            {
                results = futureVec[i].get();
                all_x_t[i] = results[0];
                sum += results[1];
            }
            break;
        }
        default:
        {
            
            for (int i = 0; i < trk.N; i++)
            {
                results = trk.asymmetric.tangentParallelAsym(allparams, i, s); //pointer to fn run through MT, arguments to fn
                all_x_t[i] = results[0];
                sum += results[1];
            }
            break;
        }
    }

    switch (trk.scaleOptimization.whichExtrema) {
        case ANY:
            break;
        case S:
            trk.x_t_s = all_x_t;
            trk.params_s = params;
            break;
        default:
            break;
    }

    switch (trk.scaleOptimization.whichExtremaX) {
    case ANY:
        break;
    case SLOP_X:
        trk.x_t_slopx = all_x_t;
        trk.params_slopx = params;
        break;
    default:
        break;
    }

    switch (trk.scaleOptimization.whichExtremaY) {
    case ANY:
        break;
    case SLOP_Y:
        trk.x_t_slopy = all_x_t;
        trk.params_slopy = params;
        break;
    default:
        break;
    }

    return -2.0 * sum;
}

double TRK::Statistics::likelihoodAsym(std::vector <double> allparams) {
    std::vector <double> all_x_t(trk.N, 0.0);
    double L = 1.0;
    if (trk.mcmc.useLogPosterior){
        L = 0.0;
    }

    std::vector <double> params;

    for (int i = 0; i < trk.M; i++) {
        params.push_back(allparams[i]);
    }

    std::vector <double> results;
    
    switch (trk.settings.ParallelizationBackEnd){
        case OPENMP:
        {
            #pragma omp parallel for num_threads(maxThreads)
            for (int i = 0; i < trk.N; i++)
            {
                results = trk.asymmetric.tangentParallelLikelihoodAsym(allparams, i); //pointer to fn run through MT, arguments to fn
                all_x_t[i] = results[0];
                if (trk.mcmc.useLogPosterior){
                    L += std::log(results[1]);
                } else {
                    L *= results[1];
                }
            }
            break;
        }
        case CPP11:
        {
        //cpp11 multithreading

            int counter = 0, completedThreads = 0, liveThreads = 0;
            std::vector< std::future < std::vector < double > > > futureVec;
            futureVec.resize(trk.N);

            for (int i = 0; i < trk.N; i++)
            {
                futureVec[i] = std::async(std::launch::async, &TRK::Asymmetric::tangentParallelLikelihoodAsym, trk.asymmetric, allparams, i); //pointer to fn run through MT, arguments to fn
                counter++;
                liveThreads++;

                if (liveThreads >= trk.settings.maxThreads)
                {
                    for (int i = completedThreads; i < counter; i++)
                    {
                        results = futureVec[i].get();
                        all_x_t.push_back(results[0]);
                        L *= results[1];
                    }
                    completedThreads += liveThreads;
                    liveThreads = 0;
                }
            }
            for (int i = completedThreads; i < trk.N; i++)
            {
                results = futureVec[i].get();
                all_x_t.push_back(results[0]);
                if (trk.mcmc.useLogPosterior){
                    L += std::log(results[1]);
                } else {
                    L *= results[1];
                }
            }
            break;
        }
        default:
        {
            for (int i = 0; i < trk.N; i++)
            {
                results = trk.asymmetric.tangentParallelLikelihoodAsym(allparams, i); //pointer to fn run through MT, arguments to fn
                all_x_t[i] = results[0];
                if (trk.mcmc.useLogPosterior){
                    L += std::log(results[1]);
                } else {
                    L *= results[1];
                }
            }
            break;
        }
    }
    return L;
}


// statistics (in the literal sense)
double TRK::Statistics::stDevUnweighted(std::vector <double> x) {
    double uppersum = 0.0;

    for (int i = 0; i < x.size(); i++) {
        uppersum += x[i];
    }

    double xbar = uppersum / x.size();

    double sum = 0.0;

    for (int i = 0; i < x.size(); i++) {
        sum += std::pow(x[i] - xbar, 2.0);
    }

    return std::sqrt(sum / (x.size() - 1.0));
}

double TRK::Statistics::getAverage(std::vector <double> x) {
    double top = 0.0;
    unsigned long N = x.size();

    for (int i = 0; i < N; i++) {
        top += x[i];
    }

    return top / N;
}

double TRK::Statistics::getAverage(std::vector <double> x, std::vector <double> w) {
    double top = 0.0;
    double bottom = 0.0;
    unsigned long N = x.size();
    
    for (int i = 0; i < N; i++) {
        top += x[i]*w[i];
        bottom += w[i];
    }
    
    return top / bottom;
}

double TRK::Statistics::getMode(int trueCount, std::vector<double> w, std::vector<double> y)
{
    int k, lowerLimit = 0, upperLimit = trueCount - 1, lowerLimitIn = -1, upperLimitIn = -1, size;
    int finalLower = 0;
    int finalUpper = 1;
    double halfWeightSum = 0, sSum, total, minDist = 999999;
    std::vector<double> sVec;
    
    while (lowerLimit != lowerLimitIn || upperLimit != upperLimitIn)
    {
        //std::cout<< lowerLimit << "\t" << upperLimit << "\n";
        lowerLimitIn = lowerLimit;
        upperLimitIn = upperLimit;
        size = upperLimit - lowerLimit + 1;
        minDist = 999999;
        halfWeightSum = 0;
        for (int i = lowerLimit; i < upperLimit + 1; i++)
        {
            halfWeightSum += w[i];
        }
        halfWeightSum *= .5;
        
        sVec.resize(size, 0.0);
        sSum = .5 * w[lowerLimit];
        sVec[0] = sSum;
        for (int i = lowerLimit + 1; i < lowerLimit + size; i++)
        {
            sSum += w[i - 1] * .5 + w[i] * .5;
            sVec[i - lowerLimit] = sSum;
        }
        
        for (size_t i = 0; i < sVec.size(); i++)
        {
            if ((sVec[i] < halfWeightSum) || trk.isEqual(sVec[i],halfWeightSum))
            {
                total = sVec[i] + halfWeightSum;
                k = (int)i; // was 0
                while (k < sVec.size() && ((sVec[k] < total) || trk.isEqual(sVec[k], total)))
                {
                    k++;
                }
                k--;
                total = std::abs(y[k + lowerLimit] - y[i + lowerLimit]);
                
                
                if (trk.isEqual(total,minDist))
                {
                    finalLower = (int)(Statistics::min(finalLower, i + lowerLimit));
                    finalUpper = (int)(Statistics::max(finalUpper, k + lowerLimit));
                }
                else if (total < minDist)
                {
                    minDist = total;
                    finalLower = (int)i + lowerLimit;
                    finalUpper = k + lowerLimit;
                }
            }
            if ((sVec[i] > halfWeightSum) || trk.isEqual(sVec[i], halfWeightSum))
            {
                total = sVec[i] - halfWeightSum;
                k = (int)i; // was svec.size() - 1
                while (k > -1 && ((sVec[k] > total) || trk.isEqual(sVec[k], total)))
                {
                    k--;
                }
                k++;
                total = std::abs(y[i + lowerLimit] - y[k + lowerLimit]);
                
                
                if (trk.isEqual(total,minDist))
                {
                    finalLower = (int)(min(finalLower, k + lowerLimit));
                    finalUpper = (int)(max(finalUpper, i + lowerLimit));
                }
                else if (total < minDist)
                {
                    minDist = total;
                    finalLower = k + lowerLimit;
                    finalUpper = (int)i + lowerLimit;
                }
            }
        }
        
        lowerLimit = finalLower;
        upperLimit = finalUpper;
        
        sVec.clear();
    }
    
    std::vector<double> newValues(y.begin() + lowerLimit, y.begin() + upperLimit + 1);
    std::vector<double> newWeights(w.begin() + lowerLimit, w.begin() + upperLimit + 1);
    return getMedian((int)newWeights.size(), newWeights, newValues);
}

double TRK::Statistics::getPeakCoord(std::vector <double> x, std::vector <double> w){
    double xPeak;
    std::vector <double> hist, edges;
    
    std::vector < std::vector <double> > histResults = getHistogram(x, w);
    
    hist = histResults[0];
    edges = histResults[1];
    
    int maxInd = argMinMax(hist)[1];
    xPeak = (edges[maxInd + 1] + edges[maxInd]) / 2.0;
    
    return xPeak;
}

double TRK::Statistics::min(double a, double b)
{
    return (a < b ? a : b);
}

double TRK::Statistics::max(double a, double b)
{
    return (a > b ? a : b);
}


// histograms
std::vector <std::vector <double> > TRK::Statistics::getHistogram(std::vector <double> data) {
    unsigned long dataSize = data.size();

    int bincount = std::round(std::sqrt(data.size()));
    std::vector <double> hist; // each element in this is the number of data points in the associated bin. contains arrays
    std::vector <std::vector <double> > bins; // each array is a bin which contains the datapoints

    std::vector <double> extrema = minMax(data);

    double binwidth = std::abs(extrema[0] - extrema[1]) / bincount;

    std::vector <double> edges = { extrema[0] };

    for (int i = 0; i < bincount; i++) {
        edges.push_back(edges[i] + binwidth);
    }

    //making sure first bin starts at slightly below min (to go against rounding errors);
    edges[0] = extrema[0] - extrema[0] *1e-9;
    //making sure last bin goes to slightly above max
    edges[bincount] = extrema[1] + extrema[1] *1e-9;

    std::vector <double> bintemp;

    for (int i = 0; i < bincount; i++) {
        bintemp.clear();
        for (int j = 0; j < dataSize - 1; j++) {
            if ((data[j] >= edges[i]) && (data[j] < edges[i + 1])) {
                bintemp.push_back(data[j]); //adds data to the bin if it is between the bounds of the bin
            }
        }
        if ((data[dataSize - 1] > edges[i]) && (data[dataSize - 1] <= edges[i + 1])) {
            bintemp.push_back(data[dataSize - 1]); //adds data to the bin if it is between the bounds of the bin
        }
        bins.push_back(bintemp); //adds an array of data for that bin to bins, the 2D array.
    }

    for (int i = 0; i < bincount; i++) {
        hist.push_back((double) bins[i].size());
    }


    return { hist, edges };
}

std::vector <std::vector <double> > TRK::Statistics::getHistogram(std::vector <double> data, std::vector <double> weights) {
    unsigned long dataSize = data.size();
    
    int bincount = std::round(std::sqrt(data.size()));
    std::vector <double> hist; // each element in this is the number of data points in the associated bin. contains arrays
    std::vector <std::vector <double> > bins; // each array is a bin which contains the datapoints
    
    std::vector <double> extrema = minMax(data);
    
    double binwidth = std::abs(extrema[0] - extrema[1]) / bincount;
    
    std::vector <double> edges = { extrema[0] };
    
    for (int i = 0; i < bincount; i++) {
        edges.push_back(edges[i] + binwidth);
    }
    
    //making sure first bin starts at slightly below min (to go against rounding errors);
    edges[0] = extrema[0] - extrema[0] *1e-9;
    //making sure last bin goes to slightly above max
    edges[bincount] = extrema[1] + extrema[1] *1e-9;
    
    double bintemp;
    
    for (int i = 0; i < bincount; i++) {
        bintemp = 0.0;
        for (int j = 0; j < dataSize - 1; j++) {
            if ((data[j] >= edges[i]) && (data[j] < edges[i + 1])) {
                bintemp += weights[j]; //adds data to the bin if it is between the bounds of the bin
            }
        }
        if ((data[dataSize - 1] > edges[i]) && (data[dataSize - 1] <= edges[i + 1])) {
            bintemp += weights[dataSize - 1]; //adds data to the bin if it is between the bounds of the bin
        }
        hist.push_back(bintemp); //adds an array of data for that bin to bins, the 2D array.
    }
    
    return { hist, edges };
}

// ###########################################################################################################################



// OPTIMIZATION ##############################################################################################################

// general downhill simplex/ nelder mead method (ND case)
std::vector <double> TRK::Optimization::downhillSimplex(std::function <double(std::vector <double>)> func, std::vector <double> guess, double tolerance, bool show_steps){
    double tol = tolerance;

    unsigned long n = (int) guess.size(); // dimensionality


    double rho = 1.0; //reflection
    double chi = 2.0; //expansion
    double gamma = 0.5; //contraction
    double sigma = 0.5; //shrinkage
    
    // simplex initialization

    std::vector <double> init_point = guess;

    std::vector <std::vector <double> > vertices(n + 1, init_point);
    std::vector <double> optimum;

    int i = 0;
    for (int j = 1; j < n + 1; j++) { //for each simplex node

        vertices[j][i] = guess[i] != 0 ? guess[i] + simplex_size*guess[i] : 0.1; //add initial "step size"
        i += 1;
    }

    std::vector <double> result;
    int it = 0;
    while (true) {
        while (true) {
            // order

            std::vector <int> orderedindices;
            std::vector <double> unorderedEvals; // ( f(x_1), f(x_2), ... f(x_n+1)
            for (int i = 0; i < n + 1; i++) {
                unorderedEvals.push_back(func(vertices[i]));
            }
            orderedindices = getSortedIndices(unorderedEvals);

            std::vector <std::vector <double> > orderedvertices = { vertices[orderedindices[0]] };
            for (int i = 1; i < n + 1; i++) {
                orderedvertices.push_back(vertices[orderedindices[i]]);
            }

            vertices = orderedvertices;

            // reflect
            std::vector <double> refpoint;
            std::vector <std::vector <double> > nvertices;
            for (int i = 0; i < n; i++) {
                nvertices.push_back(vertices[i]);
            }

            std::vector <double> centroid = findCentroid(nvertices);

            for (int i = 0; i < n; i++) {
                refpoint.push_back(centroid[i] + rho * (centroid[i] - vertices[n][i]));
            }

            double fr = func(refpoint);
            double f1 = func(vertices[0]);
            double fn = func(vertices[n - 1]);

            if (f1 <= fr && fr < fn) {
                result = refpoint;
                break;
            }

            //expand
            if (fr < f1) {
                std::vector <double> exppoint;

                for (int i = 0; i < n; i++) {
                    exppoint.push_back(centroid[i] + chi * (refpoint[i] - centroid[i]));
                }

                double fe = func(exppoint);


                if (fe < fr) {
                    result = exppoint;
                    break;
                }
                else if (fe >= fr) {
                    result = refpoint;
                    break;
                }
            }

            //contract

            if (fr >= fn) {
                //outside
                double fnp1 = func(vertices[n]);

                if (fn <= fr && fr < fnp1) {
                    std::vector <double> cpoint;

                    for (int i = 0; i < n; i++) {
                        cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
                    }

                    double fc = func(cpoint);

                    if (fc <= fr) {
                        result = cpoint;
                        break;
                    }
                    else {
                        //shrink

                        std::vector < std::vector <double> > v = { vertices[0] };

                        for (int i = 1; i < n + 1; i++) {
                            std::vector <double> vi;
                            vi.clear();
                            for (int j = 0; j < n; j++) {
                                vi.push_back(vertices[0][j] + sigma * (vertices[i][j] - vertices[0][j]));
                            }
                            v.push_back(vi);
                        }

                        vertices = v;
                    }
                }
                else if (fr >= fnp1) {
                    std::vector <double> ccpoint;

                    for (int i = 0; i < n; i++) {
                        ccpoint.push_back(centroid[i] - gamma * (centroid[i] - vertices[n][i]));
                    }

                    double fcc = func(ccpoint);

                    if (fcc <= fnp1) {
                        result = ccpoint;
                        break;
                    }
                    else {
                        //shrink

                        std::vector < std::vector <double> > v = { vertices[0] };

                        for (int i = 1; i < n + 1; i++) {
                            std::vector <double> vi;
                            vi.clear();
                            for (int j = 0; j < n; j++) {
                                vi.push_back(vertices[0][j] + sigma * (vertices[i][j] - vertices[0][j]));
                            }
                            v.push_back(vi);
                        }

                        vertices = v;
                    }


                }

            }

            //shrink

            std::vector < std::vector <double> > v = { vertices[0] };

            for (int i = 1; i < n + 1; i++) {
                std::vector <double> vi;
                vi.clear();
                for (int j = 0; j < n; j++) {
                    vi.push_back(vertices[0][j] + sigma * (vertices[i][j] - vertices[0][j]));
                }
                v.push_back(vi);
            }

            vertices = v;
        }

        std::vector <std::vector <double> > bettervertices;
        for (int i = 0; i < n; i++) {
            bettervertices.push_back(vertices[i]);
        }

        //check that new node does not have negative slops (fixes it if it does)

        bettervertices.push_back(result);

        vertices = bettervertices;
        
        if (show_steps){
            std::cout << "simplex vertex = ";
            for (int i = 0; i < result.size(); i++) {
                std::cout << result[i] << " ";
            }
            
            double eval = func(result);
            std::cout << "\tfunc = " << eval << "\n";
            
        }
        

        
        
        //test for termination

        std::vector <double> evals;
        for (int i = 0; i < n + 1; i++) {
            evals.push_back(func(vertices[i]));
        }

        if (trk.statistics.stDevUnweighted(evals) < tol) {
            break;
        }
        
        if (it >= max_simplex_iters){
            printf("Downhill simplex exceeded %i iterations; halting...\n", max_simplex_iters);
            break;
        }
        
        if (it % 100 == 0 && show_steps){printf("simplex iteration: %i\n", it);};

        it++;
    }

    optimum = vertices[n];

    return optimum;
    
}

double TRK::Optimization::downhillSimplex_1DWrapper(std::function <double(std::vector <double>)> func, std::vector <double> guess, double tolerance, bool show_steps){
    
    return downhillSimplex(func, guess, tolerance, show_steps)[0];
}

// downhill simplex/ nelder mead method customized for fitting
std::vector <double> TRK::Optimization::downhillSimplex_Fit(double(TRK::Statistics::*f)(std::vector <double>, double), std::vector <double> allparams_guess, double s, bool show_steps) {

    double tol = simplexTol;

    unsigned long n = trk.bigM; //number of model parameters plus two slop parameters
    
    if (trk.asymmetric.hasAsymSlop){
        n += 2;
    }

    double rho = 1.0; //reflection
    double chi = 2.0; //expansion
    double gamma = 0.5; //contraction
    double sigma = 0.5; //shrinkage
    
    // simplex initialization

    std::vector <double> init_point = allparams_guess;

    std::vector <std::vector <double> > vertices(n + 1, init_point);
    std::vector <double> fitted_params;

    int i = 0;
    for (int j = 1; j < n + 1; j++) { //for each simplex node

        vertices[j][i] = allparams_guess[i] != 0 ? allparams_guess[i] + simplex_size*allparams_guess[i] : 0.1; //add initial "step size"
        i += 1;
    }

    std::vector <double> result;
    int it = 0;
    while (true) {
        while (true) {
            // order

            std::vector <int> orderedindices;
            std::vector <double> unorderedEvals; // ( f(x_1), f(x_2), ... f(x_n+1)
            for (int i = 0; i < n + 1; i++) {
                unorderedEvals.push_back(evalWPriors(f, vertices[i], s));
            }
            orderedindices = getSortedIndices(unorderedEvals);

            std::vector <std::vector <double> > orderedvertices = { vertices[orderedindices[0]] };
            for (int i = 1; i < n + 1; i++) {
                orderedvertices.push_back(vertices[orderedindices[i]]);
            }

            vertices = orderedvertices;

            // reflect
            std::vector <double> refpoint;
            std::vector <std::vector <double> > nvertices;
            for (int i = 0; i < n; i++) {
                nvertices.push_back(vertices[i]);
            }

            std::vector <double> centroid = findCentroid(nvertices);

            for (int i = 0; i < n; i++) {
                refpoint.push_back(centroid[i] + rho * (centroid[i] - vertices[n][i]));
            }

            double fr = evalWPriors(f, refpoint, s);
            double f1 = evalWPriors(f, vertices[0], s);
            double fn = evalWPriors(f, vertices[n - 1], s);

            if (f1 <= fr && fr < fn) {
                result = refpoint;
                break;
            }

            //expand
            if (fr < f1) {
                std::vector <double> exppoint;

                for (int i = 0; i < n; i++) {
                    exppoint.push_back(centroid[i] + chi * (refpoint[i] - centroid[i]));
                }

                double fe = evalWPriors(f, exppoint, s);


                if (fe < fr) {
                    result = exppoint;
                    break;
                }
                else if (fe >= fr) {
                    result = refpoint;
                    break;
                }
            }

            //contract

            if (fr >= fn) {
                //outside
                double fnp1 = evalWPriors(f, vertices[n], s);

                if (fn <= fr && fr < fnp1) {
                    std::vector <double> cpoint;

                    for (int i = 0; i < n; i++) {
                        cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
                    }

                    double fc = evalWPriors(f, cpoint, trk.scaleOptimization.s);

                    if (fc <= fr) {
                        result = cpoint;
                        break;
                    }
                    else {
                        //shrink

                        std::vector < std::vector <double> > v = { vertices[0] };

                        for (int i = 1; i < n + 1; i++) {
                            std::vector <double> vi;
                            vi.clear();
                            for (int j = 0; j < n; j++) {
                                vi.push_back(vertices[0][j] + sigma * (vertices[i][j] - vertices[0][j]));
                            }
                            v.push_back(vi);
                        }

                        vertices = v;
                    }
                }
                else if (fr >= fnp1) {
                    std::vector <double> ccpoint;

                    for (int i = 0; i < n; i++) {
                        ccpoint.push_back(centroid[i] - gamma * (centroid[i] - vertices[n][i]));
                    }

                    double fcc = evalWPriors(f, ccpoint, trk.scaleOptimization.s);

                    if (fcc <= fnp1) {
                        result = ccpoint;
                        break;
                    }
                    else {
                        //shrink

                        std::vector < std::vector <double> > v = { vertices[0] };

                        for (int i = 1; i < n + 1; i++) {
                            std::vector <double> vi;
                            vi.clear();
                            for (int j = 0; j < n; j++) {
                                vi.push_back(vertices[0][j] + sigma * (vertices[i][j] - vertices[0][j]));
                            }
                            v.push_back(vi);
                        }

                        vertices = v;
                    }


                }

            }

            //shrink

            std::vector < std::vector <double> > v = { vertices[0] };

            for (int i = 1; i < n + 1; i++) {
                std::vector <double> vi;
                vi.clear();
                for (int j = 0; j < n; j++) {
                    vi.push_back(vertices[0][j] + sigma * (vertices[i][j] - vertices[0][j]));
                }
                v.push_back(vi);
            }

            vertices = v;
        }

        std::vector <std::vector <double> > bettervertices;
        for (int i = 0; i < n; i++) {
            bettervertices.push_back(vertices[i]);
        }

        //check that new node does not have negative slops (fixes it if it does)

        bettervertices.push_back(result);

        vertices = bettervertices;
        
        if (show_steps){
            std::cout << "chi-square parameters at s = " << s << " ";
            for (int i = 0; i < result.size(); i++) {
                std::cout << result[i] << " ";
            }
            
            double fitness = evalWPriors(f, result, trk.scaleOptimization.s);
            std::cout << "fitness = " << fitness << "\n";
            
        }
        
        if (trk.scaleOptimization.whichExtrema == S or trk.scaleOptimization.whichExtrema == ANY){
            double fitness = evalWPriors(f, result, trk.scaleOptimization.s);
            trk.results.chisquared = fitness;
//            printf("\n\nfinal fitness = %.3e\n\n", fitness);
        }
        
        
        //test for termination

        std::vector <double> evals;
        for (int i = 0; i < n + 1; i++) {
            evals.push_back(evalWPriors(f, vertices[i], s));
        }

        if (trk.statistics.stDevUnweighted(evals) < tol) {
            break;
        }
        
        if (it >= max_simplex_iters){
            printf("Downhill simplex exceeded %i iterations; halting...\n", max_simplex_iters);
            break;
        }
        
        if (it % 100 == 0 && show_steps){printf("simplex iteration: %i\n", it);};
//
        it++;
    }

    fitted_params = vertices[n];

    fitted_params = pegToZeroSlop(fitted_params);
    fitted_params = avoidNegativeSlop(fitted_params, n);

    return fitted_params;
}


// downhill simplex tools
std::vector <double> TRK::Optimization::findCentroid(std::vector <std::vector <double> > nvertices) {
    unsigned long n = nvertices.size();

    std::vector <double> centroid;

    for (int i = 0; i < n; i++) {// for each param
        double sum = 0.0;
        for (int j = 0; j < n; j++) { // for each vertex
            sum += nvertices[j][i];
        }
        centroid.push_back(sum / n);
    }

    return centroid;
}

std::vector <double> TRK::Optimization::pegToZeroSlop(std::vector <double> vertex){
    double pegToZeroTol = trk.scaleOptimization.pegToZeroTol;
    if (trk.settings.do1DFit){
        if (std::abs(vertex[trk.M]) <= pegToZeroTol) {
            vertex[trk.M] = 0;
        }
    } else {
        if (std::abs(vertex[trk.M]) <= pegToZeroTol) {
            vertex[trk.M] = 0;
        }
        if (std::abs(vertex[trk.M+1]) <= pegToZeroTol) {
            vertex[trk.M+1] = 0;
        }
        
        if (trk.asymmetric.hasAsymSlop){
            if (std::abs(vertex[trk.M+2]) <= pegToZeroTol) {
                vertex[trk.M+2] = 0;
            }
            if (std::abs(vertex[trk.M+3]) <= pegToZeroTol) {
                vertex[trk.M+3] = 0;
            }
        }
    }
    
    return vertex;
}

std::vector <double> TRK::Optimization::avoidNegativeSlop(std::vector <double> vertex, unsigned long n) {

    int K = 2;
    
    if (trk.asymmetric.hasAsymSlop){
        K = 4;
    }
    
    for (int k = 0; k < K; k++) {
        if (vertex[n - 1 - k] < 0) {
            vertex[n - 1 - k] = std::abs(vertex[n - 1 - k]);
        }
    }
    

    return vertex;
}

double TRK::Optimization::evalWPriors(double(TRK::Statistics::*f)(std::vector <double>, double), std::vector <double> vertex, double s) {
    if (trk.statistics.hasPriors) {
        switch (trk.statistics.priorsObject.priorType) {
            case CUSTOM:
                break;
            case CUSTOM_JOINT:
                break;
            case GAUSSIAN:
                break;
            case CONSTRAINED:

                for (int i = 0; i < trk.M; i++) { //check upper bound
                    double ub = trk.statistics.priorsObject.paramBounds[i][1];
                    double lb = trk.statistics.priorsObject.paramBounds[i][0];

                    if (vertex[i] >= ub && !std::isnan(ub)) {
                        return DBL_MAX;
                    }

                    if (vertex[i] <= lb && !std::isnan(lb)) {
                        return DBL_MAX;
                    }
                }
                break;

            case MIXED:

                for (int i = 0; i < trk.M; i++) { //check upper bound
                    double ub = trk.statistics.priorsObject.paramBounds[i][1];
                    double lb = trk.statistics.priorsObject.paramBounds[i][0];

                    if (vertex[i] >= ub && !std::isnan(ub)) {
                        return DBL_MAX;
                    }

                    if (vertex[i] <= lb && !std::isnan(lb)) {
                        return DBL_MAX;
                    }
                }
                break;
        }
    }
    
    double ret = (trk.statistics.*f)(vertex, s);
//
    if (std::isnan(ret)){
        return DBL_MAX;
    }
    
    return ret;
}

void TRK::Optimization::getBetterGuess(){
    if (trk.settings.do1DFit){
        trk.results.bestFitParams.clear();

        trk.scaleOptimization.whichExtrema = S;
        trk.allparams_s = trk.optimization.downhillSimplex_Fit(trk.statistics.selectedChiSq, trk.allparams_guess, trk.scaleOptimization.s, false);
        trk.scaleOptimization.whichExtrema = ANY;

        for (int j = 0; j < trk.M; j++) {
            trk.results.bestFitParams.push_back(trk.allparams_s[j]);
        }
        trk.results.slop_y = trk.allparams_s[trk.M];
        
        if (trk.asymmetric.hasAsymSlop){
            printf("NOTE: no current support for asymmetric 1D statistics.\n");
        }
    }
    
    for (int j = 0; j < trk.M; j++){
        trk.allparams_guess[j] = trk.results.bestFitParams[j];
    }
    if (trk.settings.do1DFit){
        trk.allparams_guess[trk.M] = trk.results.slop_y;
    } else {
        trk.allparams_guess[trk.M] = trk.results.slop_x;
        trk.allparams_guess[trk.M+1] = trk.results.slop_y;
        
        if (trk.asymmetric.hasAsymSlop){
            trk.allparams_guess[trk.M+2] = trk.results.slop_x_minus;
            trk.allparams_guess[trk.M+3] = trk.results.slop_y_minus;
        }
    }
    
    return;
}


// golden section search: finds the x_min that minimizes some f(x)
void TRK::Optimization::printGSSProgress(double yc, double yd, double a, double b, double c, double d){
    printf("f(x) = %0.3f for x = %0.3f\n", yc, c);
    printf("f(x) = %0.3f for x = %0.3f\n", yd, d);

    printf("GSS progress:\t");
    double printed_x;
    if (yc < yd){
        printed_x = (a + d) / 2.0;
    } else {
        printed_x = (c + b) / 2.0;
    }
    printf("%0.3f \t", printed_x);
    std::cout << std::endl;
    
    return;
}
    
double TRK::Optimization::goldenSectionSearch(std::function <double(double)> func, double min_bracket, double max_bracket, double tolerance){
    
    double K, h, c, d, yc, yd, x_optimum, tol = tolerance;
    double a = min_bracket, b = max_bracket;
    
    h = b - a;
    int iter = 0;
    
    c = a + invphi2 * h;
    d = a + invphi * h;
    yc = func(c);
    yd = func(d);
    K = (double) std::ceil(std::log(tol / h) / std::log(invphi));  // maximum iterations needed for convergence
        
    while (iter <= K){
        if (yc < yd) {
            b = d;
            d = c;
            yd = yc;
            h = invphi * h;
            c = a + invphi2 * h;
        } else {
            a = c;
            c = d;
            yc = yd;
            h = invphi * h;
            d = a + invphi * h;
        }
        
        if (trk.optimization.verbose_GSS){
            printGSSProgress(yc, yd, a, b, c, d);
        }
    }
    
    if (yc < yd){
        x_optimum = (a + d) / 2.0;
    } else {
        x_optimum = (c + b) / 2.0;
    }
    
    return x_optimum;
}

// ###########################################################################################################################



// SCALE OPTIMIZATION ########################################################################################################

// guesses
void TRK::ScaleOptimization::getBetterSlopYGuess(double slop_y, double s) {
    if (firstGuess) {
        firstGuess = false;
        slopYScaleGuess = s;
        slopYGuess = slop_y;
    }

    if (slop_y < slopYScaleGuess && slop_y > 0) {
        slopYScaleGuess = s;
        slopYGuess = slop_y;
    }
}


// core optimization / slope vs. scale
void TRK::ScaleOptimization::optimizeScale() {
    s = 1.0; //initially begin with s = 1
    
    if (trk.settings.do1DFit){
        printf("1D fit: no need for scale optimization.\n");
        return;
    }

    std::vector <double> scale_extrema;

    std::vector <double> s_slops;
    
    if (trk.asymmetric.hasAsymSlop){
        s_slops = {0.0, 0.0, 0.0, 0.0};
    } else {
        s_slops = {0.0, 0.0};
    }

    //optimize simultaneously
    switch (trk.settings.ParallelizationBackEnd){
        case OPENMP:
        {
            #pragma omp parallel for //num_threads(8)
            for (int i = 0; i < 2; i++)
            {
                s_slops[i] = (this->*optimizeList[i])();
            }
            break;
        }
        case CPP11:
        {
            s_slops.clear();
            
            int counter = 0, completedThreads = 0, liveThreads = 0;
            double result;
            std::vector< std::future < double > > futureVec;
            futureVec.resize(2);
            
            for (int i = 0; i < 2; i++)
            {
                futureVec[i] = std::async(std::launch::async, optimizeList[i], this); //pointer to fn run through MT, arguments to fn
                counter++;
                liveThreads++;
                
                if (liveThreads >= trk.settings.maxThreads)
                {
                    for (int i = completedThreads; i < counter; i++)
                    {
                        result = futureVec[i].get();
                        s_slops.push_back(result);
                    }
                    completedThreads += liveThreads;
                    liveThreads = 0;
                }
            }
            for (int i = completedThreads; i < 2; i++)
            {
                result = futureVec[i].get();
                s_slops.push_back(result);
            }
            
            std::vector <double> mM = minMax(s_slops);
            s_slops[0] = mM[0];
            s_slops[1] = mM[1];
            
            break;
        }
        case NONE:
        {
            for (int i = 0; i < 2; i++)
            {
                s_slops[i] = (this->*optimizeList[i])();
            }
            break;
        }
}

    

    double s_slopx = s_slops[0];
    double s_slopy = s_slops[1];

    printf("%.3e %.3e \t %.3e %.3e\n", trk.params_slopx[0], trk.params_slopx[1], trk.params_slopy[0], trk.params_slopy[1]);

    scale_extrema.push_back(s_slopx); //scale_extrema = {s_slopx, s_slopy}

    scale_extrema.push_back(s_slopy);

    std::vector <int> sortedindices = getSortedIndices(scale_extrema);

    a = scale_extrema[sortedindices[0]]; //figures out which scale is a, and which is b, as well as storing the associated best-fit parameters and their associated tangent points for those two extreme scales
    b = scale_extrema[sortedindices[1]];

    printf("extrema: \t a \t b: \n");
    printf(" \t %.3e \t %.3e \n", a, b);

    if (a == s_slopx) {
        trk.x_t_a = trk.x_t_slopx;
        trk.x_t_b = trk.x_t_slopy;

        trk.params_a = trk.params_slopx;
        trk.params_b = trk.params_slopy;
    }
    else if (a == s_slopy) {
        trk.x_t_a = trk.x_t_slopy;
        trk.x_t_b = trk.x_t_slopx;

        trk.params_a = trk.params_slopy;
        trk.params_b = trk.params_slopx;
    }


    printf("%.3e %.3e %.3e %.3e \n", trk.x_t_a[0], trk.x_t_b[0], trk.params_a[0], trk.params_b[0]);

    //determine best s1 (new s) to satistfy R2TRKp(a,s) = R2TRKp(s,b)

    std::cout << "Finding optimum scale..." << std::endl << std::endl;

    double s0 = optimize_s0_R2();

    s = s0;

    double s_final = iterateR2_OptimumScale(s0);

    s = s_final;

    std::cout << "optimum s = " << s << std::endl;

    trk.results.bestFitParams.clear();

    for (int j = 0; j < trk.M; j++) {
        trk.results.bestFitParams.push_back(trk.allparams_s[j]);
    }

    trk.results.slop_x = trk.allparams_s[trk.M];
    trk.results.slop_y = trk.allparams_s[trk.M + 1];
    
    if (trk.asymmetric.hasAsymSlop){
        trk.results.slop_x_minus = trk.allparams_s[trk.M + 2];
        trk.results.slop_y_minus = trk.allparams_s[trk.M + 3];
    }

    trk.results.optimumScale = s;
    trk.results.minimumScale = a;
    trk.results.maximumScale = b;

    return;
}

double TRK::ScaleOptimization::optimize_s_SlopX() {

    iterative_allparams_guess = trk.allparams_guess;

    // before doing any standard simplex movement, here it checks whether the simplex is within the zero "plateau", and if so, it moves it to the boundary.
    // for slop x: move to right until it hits the boundary


    //bracket finding

    double a = 0.0;
    double b = 1.0;
    double trial_s = 1.0;
    double slop_trial_s = innerSlopX_Simplex({ trial_s }, iterative_allparams_guess);

    double inc = trial_s * 0.5;

    if (slop_trial_s > 0) {
        b = trial_s;

        double trial_a = trial_s;

        while (true) {
            trial_a -= inc;

            double slop_trial_a = innerSlopX_Simplex({ trial_a }, iterative_allparams_guess);

            if (slop_trial_a == 0) {
                a = trial_a;
                break;
            }
            else if (slop_trial_a > 0) {
                inc *= 0.5;

                b = trial_a;
            }
        }
    }
    else if (slop_trial_s == 0) {
        a = trial_s;

        double trial_b = trial_s;

        while (true) {
            trial_b += inc;

            double slop_trial_b = innerSlopX_Simplex({ trial_b }, iterative_allparams_guess);

            if (slop_trial_b > 0) {
                b = trial_b;
                break;
            }
            else if (slop_trial_b == 0) {
                a = trial_b;
            }
        }
    }

    //bisection, now that we have brackets [a,b]

    double c, slop_c;
    double tol_bisect = 2e-3;
    double tol_brackets = 1e-3;

    while (true) {
        c = (a + b) / 2;

        whichExtremaX = SLOP_X;
        slop_c = innerSlopX_Simplex({ c }, iterative_allparams_guess);
        whichExtremaX = ANY;

        if (slop_c <= tol_bisect && slop_c > 0) { //convergence criterion
            break;
        }

        if (slop_c > 0) {
            b = c;
        }
        else if (slop_c == 0) {
            a = c;
        }

        if (std::abs(a - b) <= tol_brackets) { //secondary convergence criterion (bracket width)
            break;
        }
    }

    return c;
}

double TRK::ScaleOptimization::optimize_s_SlopY() {

    iterative_allparams_guess = trk.allparams_guess;

    // before doing any standard simplex movement, here it checks whether the simplex is within the zero "plateau", and if so, it moves it to the boundary.
    // for slop x: move to right until it hits the boundary

        //bracket finding

    double a = 0.0;
    double b = 1.0;
    double trial_s = slopYScaleGuess;
    double slop_trial_s = innerSlopY_Simplex({ trial_s }, iterative_allparams_guess);

    double inc = trial_s * 0.5;
    
    int jumpcheck = 0;

    if (slop_trial_s > 0) {
        a = trial_s;

        double trial_b = trial_s;

        while (true) {
            if (jumpcheck > 1){
                inc *= 2.0;
                jumpcheck = 0;
            }
            trial_b += inc;

            double slop_trial_b = innerSlopY_Simplex({ trial_b }, iterative_allparams_guess);

            if (slop_trial_b == 0) {
                b = trial_b;
                break;
            }
            else if (slop_trial_b > 0) {
                a = trial_b;
                jumpcheck += 1;
            }
        }
    }
    else if (slop_trial_s == 0) {
        b = trial_s;

        double trial_a = trial_s;

        while (true) {
            trial_a -= inc;

            double slop_trial_a = innerSlopY_Simplex({ trial_a }, iterative_allparams_guess);

            if (slop_trial_a > 0) {
                a = trial_a;
                break;
            }
            else if (slop_trial_a == 0) {
                inc *= 0.5;
                b = trial_a;
            }
        }
    }

    //bisection, now that we have brackets [a,b]

    double c, slop_c;
    double tol_bisect = 2e-3;
    double tol_brackets = 1e-3;

    while (true) {
        c = (a + b) / 2;

        whichExtremaY = SLOP_Y;
        slop_c = innerSlopY_Simplex({ c }, iterative_allparams_guess);
        whichExtremaY = ANY;

        if ((slop_c <= tol_bisect && slop_c > 0)) { //convergence criterion
            break;
        }

        if (slop_c > 0) {
            a = c;
        }
        else if (slop_c == 0) {
            b = c;
        }

        if (std::abs(a - b) <= tol_brackets) { //secondary convergence criterion (bracket width)
            break;
        }
    }



    return c;
}

double TRK::ScaleOptimization::innerSlopX_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
    trk.allparams_s = trk.optimization.downhillSimplex_Fit(trk.statistics.selectedChiSq, allparams_guess, ss[0], showSimplexSteps);

    if (verbose){
        if (trk.asymmetric.hasAsymSlop){
            printf("%.3e \t %.3e \t %.3e \t %.3e \t %.3e \t(slop x optimization)\n", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1], trk.allparams_s[trk.M + 2], trk.allparams_s[trk.M + 3]);
        } else {
            printf("s=%.3e \t slop_x=%.3e \t slop_y=%.3e \t(slop x optimization) model params: ", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1]);
            for (int j = 0; j < trk.M; j++){
                printf("%.3e\t", trk.allparams_s[j]);
            }
            printf("\n");
        }
    }
    

    //double sec_elapsed = secElapsed(time);

    //printf("%.3e sec, max threads = %i \n", sec_elapsed, maxThreads);

    getBetterSlopYGuess(trk.allparams_s[trk.M + 1], s);

    return trk.allparams_s[trk.M];
}

double TRK::ScaleOptimization::innerSlopY_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
    //s = ss[0];

    
    trk.allparams_s = trk.optimization.downhillSimplex_Fit(trk.statistics.selectedChiSq, allparams_guess, ss[0], showSimplexSteps);
    
    if (verbose){
        if (trk.asymmetric.hasAsymSlop){
            printf("%.3e \t %.3e \t %.3e \t %.3e \t %.3e \t(slop y optimization)\n", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1], trk.allparams_s[trk.M + 2], trk.allparams_s[trk.M + 3]);
        } else {
            printf("s=%.3e \t slop_x=%.3e \t slop_y=%.3e \t(slop y optimization) model params: ", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1]);
            for (int j = 0; j < trk.M; j++){
                printf("%.3e\t", trk.allparams_s[j]);
            }
            printf("\n");
        }
    }

    return trk.allparams_s[trk.M + 1];
}


// TRK correlation coefficient
double TRK::ScaleOptimization::optimize_s0_R2() {

    iterative_allparams_guess = trk.allparams_guess;

    //bracket finding

    double left, right, f_left;
    
    //bisection, now that we have brackets [left,right]

    left = a;
    right = b;

    f_left = innerR2_Simplex({ a }, iterative_allparams_guess);

    double c, f_c;
    double tol_bisect = 1e-4;
    double tol_brackets = 1e-3;

    while (true) {
        c = (left + right) / 2;

        whichExtrema = S;
        f_c = innerR2_Simplex({ c }, iterative_allparams_guess);
        whichExtrema = ANY;

        //printf("%.3e %.3e \n", f_c, c);

        if (std::abs(f_c) <= tol_bisect) { //convergence criterion
            break;
        }

        if (f_c * f_left > 0) {
            left = c;
            f_left = f_c;
        }
        else if (f_c * f_left < 0) {
            right = c;
        }

        if (std::abs(left - right) <= tol_brackets) { //secondary convergence criterion (bracket width)
            break;
        }
    }

    return c;

}

double TRK::ScaleOptimization::optimize_s_prime_R2(double s0) {

    iterative_allparams_guess = trk.allparams_guess;

    //bracket finding

    double left, right, f_left;

    f_left = innerR2_iter_Simplex({ a }, iterative_allparams_guess, s0);

    //bisection, now that we have brackets [left,right]

    left = a;
    right = b;

    double c, f_c;
    double tol_bisect = 1e-4;
    double tol_brackets = 1e-3;

    while (true) {
        //printf("brackets: %.3e %.3e \n", left, right);
        c = (left + right) / 2;

        whichExtrema = S;
        f_c = innerR2_iter_Simplex({ c }, iterative_allparams_guess, s0);
        whichExtrema = ANY;

        if (std::abs(f_c) <= tol_bisect) { //convergence criterion
            break;
        }

        if (f_c * f_left > 0) {
            left = c;
            f_left = f_c;
        }
        else if (f_c * f_left < 0) {
            right = c;
        }

        if (std::abs(left - right) <= tol_brackets) { //secondary convergence criterion (bracket width)
            break;
        }
    }

    return c;

}

double TRK::ScaleOptimization::iterateR2_OptimumScale(double s0) {
    double tol_scale = 1e-3;

    double s1 = 0.0;

    bool tolcheck = false;

    while (!tolcheck) {
        printf("next s0: %.3e \n", s0);

        s1 = optimize_s_prime_R2(s0);
        if (std::abs(s1-s0) <= tol_scale) {
            tolcheck = true;
        }
        std::printf("new s0: %.3e \n", s1);
        s0 = s1;
    }

    return s1;
}

double TRK::ScaleOptimization::innerR2_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
    //s = ss[0];

    whichExtrema = S;
    trk.allparams_s = trk.optimization.downhillSimplex_Fit(trk.statistics.selectedChiSq, allparams_guess, ss[0], showSimplexSteps);
    whichExtrema = ANY;

    if (verbose){
        if (trk.asymmetric.hasAsymSlop){
            printf("%.3e \t %.3e \t %.3e \t %.3e \t %.3e \t(initial R2 optimization)\n", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1], trk.allparams_s[trk.M + 2], trk.allparams_s[trk.M + 3]);
        } else {
            printf("%.3e \t %.3e \t %.3e \t(initial R2 optimization)\n", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1]);
        }
    }

    double R2as = R2TRK_prime_as();
    double R2sb = R2TRK_prime_sb();

    return R2as - R2sb;
}

double TRK::ScaleOptimization::innerR2_iter_Simplex(std::vector <double> ss, std::vector <double> allparams_guess, double s0) {
    //s = ss[0];

    whichExtrema = S;
    trk.allparams_s = trk.optimization.downhillSimplex_Fit(trk.statistics.selectedChiSq, allparams_guess, ss[0], showSimplexSteps);
    whichExtrema = ANY;

    if (verbose){
        if (trk.asymmetric.hasAsymSlop){
            printf("%.3e \t %.3e \t %.3e \t %.3e \t %.3e \t(additional R2 optimization)\n", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1], trk.allparams_s[trk.M + 2], trk.allparams_s[trk.M + 3]);
        } else {
            printf("%.3e \t %.3e \t %.3e \t(additional R2 optimization)\n", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1]);
        }
    }

    double R2as = R2TRK_prime_as0(s0, trk.x_t_s, trk.params_s);
    double R2sb = R2TRK_prime_s0b(s0, trk.x_t_s, trk.params_s);

    return R2as - R2sb;
}

double TRK::ScaleOptimization::R2TRK_prime_as() {
    double R2 = 1.0 / trk.N;

    double sum = 0.0;

    for (int n = 0; n < trk.N; n++) {
        double m_tn_a = trk.dyc(trk.x_t_a[n], trk.params_a);
        double theta_t_a = std::atan(m_tn_a);

        double m_tn_s = trk.dyc(trk.x_t_s[n], trk.params_s);
        double theta_t_s = std::atan(m_tn_s);

        sum += std::pow(std::tan(PI/4.0 - std::abs(theta_t_a - theta_t_s)/ 2.0), 2.0);
    }

    R2 *= sum;

    return R2;
}

double TRK::ScaleOptimization::R2TRK_prime_sb() {
    double R2 = 1.0 / trk.N;

    double sum = 0.0;

    for (int n = 0; n < trk.N; n++) {
        double m_tn_s = trk.dyc(trk.x_t_s[n], trk.params_s);
        double theta_t_s = std::atan(m_tn_s);

        double m_tn_b = trk.dyc(trk.x_t_b[n], trk.params_b);
        double theta_t_b = std::atan(m_tn_b);

        sum += std::pow(std::tan(PI / 4.0 - std::abs(theta_t_s - theta_t_b) / 2.0), 2.0);
    }

    R2 *= sum;

    return R2;
}

double TRK::ScaleOptimization::R2TRK_prime_as0(double s0, std::vector <double> x_t_s1, std::vector <double> params_s1) {
    double R2 = 1.0 / trk.N;

    double sum = 0.0;

    for (int n = 0; n < trk.N; n++) {
        double m_tn_a = trk.dyc(trk.x_t_a[n], trk.params_a);
        double theta_t_a = std::atan(m_tn_a);

        double m_tn_s1 = trk.dyc(x_t_s1[n], params_s1);
        double theta_t_s1 = std::atan(m_tn_s1);

        sum += std::pow(std::tan(PI / 4.0 - std::abs(std::atan(s0*std::tan(theta_t_a)) - std::atan(s0*std::tan(theta_t_s1))) / 2.0), 2.0);
    }

    R2 *= sum;

    return R2;
}

double TRK::ScaleOptimization::R2TRK_prime_s0b(double s0, std::vector <double> x_t_s1, std::vector <double> params_s1) {
    double R2 = 1.0 / trk.N;

    double sum = 0.0;

    for (int n = 0; n < trk.N; n++) {
        double m_tn_b = trk.dyc(trk.x_t_b[n], trk.params_b);
        double theta_t_b = std::atan(m_tn_b);

        double m_tn_s1 = trk.dyc(x_t_s1[n], params_s1);
        double theta_t_s1 = std::atan(m_tn_s1);

        sum += std::pow(std::tan(PI / 4.0 - std::abs(std::atan(s0*std::tan(theta_t_s1)) - std::atan(s0*std::tan(theta_t_b))) / 2.0), 2.0);
    }

    R2 *= sum;

    return R2;
}

// ###########################################################################################################################



// TANGENT POINT METHODS #####################################################################################################

// root finders
double TRK::TangentPointMethods::newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess) {
    double x0 = xguess;
    int itercount = 0;

    //initial iteration
    double f = (trk.yc(x0, params) - y_n) * trk.dyc(x0, params) * Sig_xn2 + (x0 - x_n) * Sig_yn2;  //function we're finding zero of
    double df = (std::pow(trk.dyc(x0, params), 2.0) + (trk.yc(x0, params) - y_n)*trk.ddyc(x0, params)) * Sig_xn2 + Sig_yn2; //derivative of above

    double x1 = x0 - f / df;

    double tol = 1e-3;

    while (std::abs(x1 - x0) > tol) {
        x0 = x1;

        f = (trk.yc(x0, params) - y_n) * trk.dyc(x0, params) * Sig_xn2 + (x0 - x_n) * Sig_yn2;  //function we're finding zero of
        df = (std::pow(trk.dyc(x0, params), 2.0) + (trk.yc(x0, params) - y_n)*trk.ddyc(x0, params)) * Sig_xn2 + Sig_yn2; //derivative of above

        x1 = x0 - f / df;

        itercount += 1;
    }

    return x1;
}

double TRK::TangentPointMethods::twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1)
{
    double tol = 1e-9;
    double xkm1 = xguess;
    double xk = xguessp1;

    double xkp1;
    double r;
    double ykm1;
    double yk;
    double dyk;

    int itercount = 0;

    while (true) {

        ykm1 = (trk.yc(xkm1, params) - y_n) * trk.dyc(xkm1, params) * Sig_xn2 + (xkm1 - x_n) * Sig_yn2;
        yk = (trk.yc(xk, params) - y_n) * trk.dyc(xk, params) * Sig_xn2 + (xk - x_n) * Sig_yn2;  //function we're finding zero of
        dyk = (std::pow(trk.dyc(xk, params), 2.0) + (trk.yc(xk, params) - y_n)*trk.ddyc(xk, params)) * Sig_xn2 + Sig_yn2; //derivative of above

        r = 1.0 - ((yk / ykm1) * (((yk - ykm1) / (xk - xkm1)) / dyk));

        xkp1 = (1.0 - 1.0 / r)*xkm1 + xk / r;

        if (std::abs(xkp1) > root_bound * trk.datawidth) {
            while (std::abs(xkp1) > root_bound * trk.datawidth) {
                xkp1 = xkp1 / 2.0;
            }
        }

        //bisection?
        
        if (yk * ykm1 < 0) { //checks if xk and xkm1 can be used as bisection brackets

            double c, f_c, f_left;
            double left = 0.0;
            double right = 1.0;
            double tol_brackets = 1e-9;

            if (xk < xkm1) {
                left = xk;
                right = xkm1;
            }
            else if (xkm1 < xk) {
                left = xkm1;
                right = xk;
            }

            while (true) {
                c = (left + right) / 2;

                f_c = (trk.yc(c, params) - y_n) * trk.dyc(c, params) * Sig_xn2 + (c - x_n) * Sig_yn2;

                if (std::abs((left - right) / 2.0) <= tol_brackets) { //secondary convergence criterion (bracket width)
                    break;
                }

                f_left = (trk.yc(left, params) - y_n) * trk.dyc(left, params) * Sig_xn2 + (left - x_n) * Sig_yn2;

                if (f_c * f_left > 0){ //same sign
                    left = c;
                }
                else if (f_c * f_left < 0) {
                    right = c;
                }
            }

            return c;
        }
        
        xkm1 = xk;
        xk = xkp1;

        itercount += 1;


        if (itercount > 5000) {
            return NAN;
        }

        if (std::abs(xk - xkm1) <= tol || std::isnan(xk) || std::abs(yk) <= tol) {
            break;
        }
    }
    return xkp1;
}


// tangent point choosing
std::vector <double> TRK::TangentPointMethods::tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg) {
    
    std::vector <double> result;

    double xg1 = xg;

    std::vector <double> xr1vec;
    double xr1old = 0.0;
    
    bool checkcheck = false;

    int itcount = 0;

    while (true) {
        if (xr1vec.size() > 99) {
            printf("100 iterations of tangent finder loop! \n");
            for (int j = 0; j < params.size(); j++) {
                printf("%.3e ", params[j]);
            }
            printf("%.3e  %.3e  %.3e  %.3e \t", Sig_xn2, Sig_yn2, x_n, y_n);
            printf("s = %.3e\n", trk.scaleOptimization.s);
        }
        result.clear();

        double xr1 = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg1, xg1 + std::sqrt(Sig_xn2)/10.0);

        if (std::isnan(xr1) && xr1vec.size() >= 1) { //is a NAN if RF did not converge
            return { xr1vec[0] };
        }

        xr1vec.push_back(xr1);

        if (checkcheck) { //different guess gives same root;

            if (std::abs(xr1 - xr1old) < 1e-3) {
                result.push_back(xr1);
                break;
            }
            if (xr1vec.size() >= 3){
                if (std::abs(xr1vec[xr1vec.size() - 1] - xr1vec[xr1vec.size() - 3]) < 1e-3) { //root finder oscillating between two roots

                    std::vector <double> unorderedbrackets = { xr1vec[xr1vec.size() - 1], xr1vec[xr1vec.size() - 2] };

                    std::vector <double> orderedbrackets = minMax(unorderedbrackets);

                    double xg_mid = (orderedbrackets[1] - orderedbrackets[0]) / 2.0 + orderedbrackets[0];

                    double xr_mid = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg_mid, xg_mid + std::sqrt(Sig_xn2) / 10.0);

                    result = { orderedbrackets[0], xr_mid, orderedbrackets[1] };
                    break;
                }
            }
        }

        double xg2, xg3;

        //Quadratic Approximation

        std::vector <double> allRoots = approxQuadraticRoots(params, x_n, y_n, Sig_xn2, Sig_yn2, xr1); //get approx. roots from quadratic taylor approximation
        std::vector <double> extraRoots;

        double ident_roots_tol = 1e-8;

        if (allRoots.size() == 3) { //either add the two other roots, or no more roots (depending on discriminant of cubic equation)
            //grab to new roots
            for (int i = 0; i < 3; i++) {
                double root = allRoots[i];
                if (std::abs(root - xr1) <= ident_roots_tol) { //checks if roots are (numerically) identical or not
                    continue;
                }
                else if (std::abs(root - xr1) > ident_roots_tol){
                    extraRoots.push_back(root);
                }
            }
        }

        if (extraRoots.size() == 2) { //if have 2 additional, real roots
            double xr2, xr3;
            xg2 = extraRoots[0];
            xg3 = extraRoots[1];

            if (xg2 < xr1 && xr1 < xg3) {
                xr2 = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg2, xg2 - std::sqrt(Sig_xn2) / 10.0);
                xr3 = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg3, xg3 + std::sqrt(Sig_xn2) / 10.0);

                result.push_back(xr2);
                result.push_back(xr1);
                result.push_back(xr3);

                break;
            }
            else {
                std::vector <double> rootVec = { xr1, xg2, xg3 };
                std::sort(rootVec.begin(), rootVec.end());

                xg1 = rootVec[1];
                xr1old = xr1;

                checkcheck = true;
            }
        } else if (extraRoots.size() == 0 && xr1vec.size() == 1) {//if initial quadratic approximation didn't yield any more guesses, try to find roots with guesses of leftmost and rightmost x values

            
            double xr_left = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, trk.x_min, trk.x_min - std::sqrt(Sig_xn2) / 10.0);
            double xr_right = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, trk.x_max, trk.x_max + std::sqrt(Sig_xn2) / 10.0);

            result.push_back(xr1);
            result.push_back(xr_left);
            result.push_back(xr_right);

            break;
        } else if (extraRoots.size() == 0 && xr1vec.size() == 2) { //found two roots but can't find a third
            result = xr1vec;
            break;
        } else if (extraRoots.size() == 3) {//this can happen if the "root" found is very close to being a root but isn't actually one.
            
            //in this case, try again with a different guess.
            double xr_left = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, trk.x_min, trk.x_min - std::sqrt(Sig_xn2) / 10.0);
            double xr_right = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, trk.x_max, trk.x_max + std::sqrt(Sig_xn2) / 10.0);

            result.push_back(xr_left);
            result.push_back(xr_right);

            break;
        }

        itcount += 1;
    }

    return result;
}

double TRK::TangentPointMethods::findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec, double s) {
    std::vector <double> posts;
    long minindex;

    for (int i = 0; i < x_tn_vec.size(); i++) {
        posts.push_back(trk.statistics.singlePointLnL(params, x_n, y_n, Sig_xn2, Sig_yn2, x_tn_vec[i], trk.scaleOptimization.s));
    }

    std::vector<double>::iterator result = std::min_element(std::begin(posts), std::end(posts));
    minindex = std::distance(std::begin(posts), result);

    return x_tn_vec[minindex];
}


// approximations
std::vector <double> TRK::TangentPointMethods::tangentCubicSolver(double A, double B, double C, double D) {
    //cubic solver for three real and distinct roots
    double a1 = B / A;
    double a2 = C / A;
    double a3 = D / A;

    double Q = (a1 * a1 - 3.0 * a2) / 9.0;
    double R = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
    double Qc = std::pow(Q, 3.0);

    double theta = std::acos(R / sqrt(Qc));

    double r1 =  -2 * std::sqrt(Q) * std::cos(theta / 3) - a1 / 3;
    double r2 = -2 * std::sqrt(Q) * std::cos((theta + 2 * PI) / 3) - a1 / 3;
    double r3 = -2 * std::sqrt(Q) * std::cos((theta + 4 * PI) / 3) - a1 / 3;

    std::vector <double> roots = { r1, r2, r3 };
    std::vector <double> goodroots;

    for (int i = 0; i < 3; i++) {
        if (std::abs(roots[i]) < (root_bound * trk.datawidth)) {
            goodroots.push_back(roots[i]);
        }
    }

    return goodroots;
}

std::vector <double> TRK::TangentPointMethods::approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xr1) {
    //using board derivation notation:
    double b = trk.yc(xr1, params) - trk.dyc(xr1, params) * xr1 + (trk.ddyc(xr1, params) / 2.0) * std::pow(xr1, 2.0); //coefficients of quadratic approximation from taylor expansion
    double m = trk.dyc(xr1, params) - trk.ddyc(xr1, params) * xr1;
    double a = trk.ddyc(xr1, params) / 2.0;


    //DIFFERENT FROM NOTATION OF BOARD DERIVATION!
    double A = 2.0 * std::pow(a, 2.0) * Sig_xn2;                                            // coef of x^3
    double B = 3.0 * a * m * Sig_xn2;                                                        //coef of x^2
    double C = (2.0 * a * (b - y_n) + std::pow(m, 2.0)) * Sig_xn2 + Sig_yn2;    //coef of x
    double D = m * (b - y_n) * Sig_xn2 - x_n * Sig_yn2;                            //coef of 1

    double discriminant = 18.0*A*B*C*D - 4.0*std::pow(B, 3.0)*D + std::pow(B, 2.0)*std::pow(C, 2.0) - 4.0*A*std::pow(C, 3.0) - 27.0*std::pow(A, 2.0)*std::pow(D, 2.0);

    std::vector <double> roots;

    if (discriminant > 0) {
        roots = tangentCubicSolver(A, B, C, D);
        return roots;
    }
    //returns no extra roots (empty vector) if the other two roots are
    return roots;
}


// parallel computing methods
std::vector <double> TRK::TangentPointMethods::tangentParallel(std::vector<double> params, double slop_x, double slop_y, int n, double s) {
    double Sig_xn2 = std::pow(trk.sx[n], 2.0) + std::pow(slop_x, 2.0);
    double Sig_yn2 = std::pow(trk.sy[n], 2.0) + std::pow(slop_y, 2.0);

    std::vector <double> x_tn_vec = trk.tangentPointMethods.tangentsFinder(params, trk.x[n], trk.y[n], Sig_xn2, Sig_yn2, trk.x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points

    double x_t = findBestTangent(params, trk.x[n], trk.y[n], Sig_xn2, Sig_yn2, x_tn_vec, s);

    double m_tn = trk.dyc(x_t, params);
    double y_tn = trk.yc(x_t, params);

    double subsum1 = trk.w[n] * std::pow(trk.y[n] - y_tn - m_tn * (trk.x[n] - x_t), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2);
    double subsum2 = trk.w[n] * std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(s*Sig_yn2, 2.0)));

    return { x_t, subsum1, subsum2 };
}

std::vector <double> TRK::TangentPointMethods::tangentParallelLikelihood(std::vector<double> params, double slop_x, double slop_y, int n) {
    double Sig_xn2 = std::pow(trk.sx[n], 2.0) + std::pow(slop_x, 2.0);
    double Sig_yn2 = std::pow(trk.sy[n], 2.0) + std::pow(slop_y, 2.0);

    std::vector <double> x_tn_vec = trk.tangentPointMethods.tangentsFinder(params, trk.x[n], trk.y[n], Sig_xn2, Sig_yn2, trk.x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points

    double x_t = findBestTangent(params, trk.x[n], trk.y[n], Sig_xn2, Sig_yn2, x_tn_vec, trk.scaleOptimization.s);

    double m_tn = trk.dyc(x_t, params);
    double y_tn = trk.yc(x_t, params);

    double l = std::pow((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(trk.scaleOptimization.s*Sig_yn2, 2.0)), trk.w[n]/2.0);
    l *= std::exp(-0.5 * trk.w[n] * (std::pow(trk.y[n] - y_tn - m_tn * (trk.x[n] - x_t), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2)));
    
//    printf("%.3e\n",l);

    return { x_t, l};
}

// ###########################################################################################################################



// MCMC/SAMPLING #############################################################################################################

// uncertainty estimation
void TRK::MCMC::calculateUncertainties() {
    
//    useLogPosterior = false;
    // only needed the above true for pivot point sampling
//    goodDeltasFound = false;
    // recompute step sizes cause they may be different given new pivot point

    std::vector <std::vector <std::vector <double> > > allparam_uncertainties;
    
    if (trk.mcmc.verbose){
        std::cout << "\nSampling Posterior...\n";
    }

    std::vector <std::vector <double >> allparam_samples = samplePosterior(R, burncount, allparams_sigmas_guess);

    if (trk.settings.outputDistributionToFile) {

        std::string fileName = trk.settings.outputPath + std::string("/TRKMCMC_") + std::to_string(trk.allparams_guess[0]) + std::string("_") + std::to_string(R) + std::string(".txt");

        std::ofstream myfile;
        myfile.open(fileName, std::ofstream::trunc);
        
        int m = 0;
        if (trk.asymmetric.hasAsymSlop){
            m = 2;
        }
        
        if (myfile.is_open())
        {
            // filename    a     b     optimum scale    total computation time (s)
            for (int i = 0; i < allparam_samples.size(); i++) {
                for (int j = 0; j < trk.bigM + m; j++) {
                    (j < trk.bigM + m - 1) ? myfile << allparam_samples[i][j] << " " : myfile << allparam_samples[i][j];
                }
                myfile << std::endl;
            }

            myfile.close();
        }
        else std::cout << "Unable to open file";
    }

    if (trk.mcmc.verbose){
        std::cout << "Computing Uncertainties...\n";
    }

    allparam_uncertainties = lowerBar(allparam_samples); //for each parameter including slope, there is a vector containing 1 vector of -sigmas, 1 vector of +sigmas. This vector contains all of those 2-vectors.

    trk.results.bestFit_123Sigmas.clear();

    for (int j = 0; j < trk.M; j++) {
        trk.results.bestFit_123Sigmas.push_back(allparam_uncertainties[j]);
    }
    if (trk.settings.do1DFit){
        trk.results.slopY_123Sigmas = allparam_uncertainties[trk.M];
    } else {
        trk.results.slopX_123Sigmas = allparam_uncertainties[trk.M];
        trk.results.slopY_123Sigmas = allparam_uncertainties[trk.M + 1];
        
        if (trk.asymmetric.hasAsymSlop){
            trk.results.slopX_minus_123Sigmas = allparam_uncertainties[trk.M];
            trk.results.slopY_minus_123Sigmas = allparam_uncertainties[trk.M + 1];
        }
    }

    return;
}

std::vector <std::vector <std::vector <double> > >  TRK::MCMC::lowerBar(std::vector <std::vector <double> > allparam_samples) { //method used to estimate uncertainty of a sampled distribution
    std::vector <double> data, hist, edges, minusSigmas, plusSigmas;
    std::vector <std::vector <std::vector <double> > > allparam_uncertainties;
    std::vector <std::vector <double> > histResults;
    unsigned long totalCount = allparam_samples.size();
    double tolBar = 1e-6;

    trk.results.paramDistributionHistograms.clear();
    
    int m = 0;
    if (trk.asymmetric.hasAsymSlop){
        m = 2;
    }

    for (int j = 0; j < trk.bigM + m; j++) { //for each model param plus slop
        data.clear();

        for (int i = 0; i < totalCount; i++) {
            data.push_back(allparam_samples[i][j]);
        }

        histResults = trk.statistics.getHistogram(data);

        trk.results.paramDistributionHistograms.push_back(histResults);

        hist = histResults[0]; // each number in this is the number of samples within each bin
        edges = histResults[1];

        double mean = trk.statistics.getAverage(data);

        //bar lowering iteration (essentially a bisection algo)

        double low = 0.0;
        double high = minMax(hist)[1];
        double bar = (low + high) / 2.0;
        double leftBound, rightBound, minusSig, plusSig;
        int aboveCount; // number of SAMPLES within the bins above bar

        double ratioIn = 0.0;

        std::vector <int> indicesIn;

        minusSigmas.clear();
        plusSigmas.clear();

        for (int sigNum = 0; sigNum < 3; sigNum++) { //for all of 1, 2 and 3 sigma
            low = 0.0;
            high = minMax(hist)[1];
            while (std::abs(ratioIn - SIGMAS[sigNum]) > tolBar && high - low > tolBar) {
                indicesIn.clear();

                aboveCount = 0;

                for (int k = 0; k < hist.size(); k++) { //determines how many samples are above bar, and which bins there are that lie within
                    if (hist[k] > bar) {
                        aboveCount += hist[k];
                        indicesIn.push_back(k);
                    }
                }

                ratioIn = (double)aboveCount / (double)totalCount;

                if (ratioIn < SIGMAS[sigNum]) { //bar too high
                    high = bar;
                }
                else if (ratioIn >= SIGMAS[sigNum]) {
                    low = bar;
                }

                bar = (low + high) / 2.0;
            }

            unsigned long K = indicesIn.size();
            leftBound = (edges[indicesIn[0]] + edges[indicesIn[0] + 1]) / 2.0;
            rightBound = (edges[indicesIn[K - 1]] + edges[indicesIn[K - 1] + 1]) / 2.0;

            minusSig = leftBound - mean;
            plusSig = rightBound - mean;

            minusSigmas.push_back(minusSig);
            plusSigmas.push_back(plusSig);
        }

        ratioIn = 0.0;

        allparam_uncertainties.push_back({ minusSigmas, plusSigmas });
    }

    return allparam_uncertainties; //for each parameter including slope, there is a vector containing 1 vector of +sigmas, 1 vector of -sigmas. This vector contains all of those 2-vectors.
}


// sampling (general)
double TRK::MCMC::metHastRatio(std::vector <double> X_trial, std::vector <double> X_i){
    double log_a;

    if (trk.statistics.hasPriors) {
        log_a = (trk.statistics.*trk.statistics.selectedLikelihood)(X_trial) - (trk.statistics.*trk.statistics.selectedLikelihood)(X_i) + std::log(trk.statistics.priors(X_trial)) - std::log(trk.statistics.priors(X_i));
            // these likelihoods return log likelihood given useLogPosterior = true; the computation is done WITHIN the function.
    }
    else {
        log_a = (trk.statistics.*trk.statistics.selectedLikelihood)(X_trial) - (trk.statistics.*trk.statistics.selectedLikelihood)(X_i);
            // this returns the log likelihood given useLogPosterior = true; the computation is done WITHIN the function.
    }
    
    return log_a; // returns log post / log post if useLogPosterior == true
}

std::vector <std::vector <double >> TRK::MCMC::samplePosterior(int R, int burncount, std::vector <double> sigmas_guess) {
    
    useLogPosterior = true;

    unsigned long n = trk.bigM;
    std::vector < std::vector <double > > result, result_final;
    
    switch (thisSamplingMethod) {
        case AIES: {
            // INDEXING:
            // iterations: t
            // walkers: j, k
            // coordinates/parameters: i
            
            int L = amt_walkers * (int) trk.bigM;   // number of walkers
            
            // initialize walkers
            std::vector <std::vector <double> > all_walkers(L, std::vector <double> (trk.bigM, 0.0));
            std::vector <std::vector <double> > YY;
            
            for (int j = 0; j < L; j++){
                for (int i = 0; i < trk.bigM; i++){
                    all_walkers[j][i] = rnorm(trk.allparams_guess[i], trk.allparams_guess[i] != 0 ? trk.allparams_guess[i]*AIES_initial_scaling : 0.1);
                }
            }
            
            std::vector <double> X, Y, X_trial;
            // sample
            int sample_count = 0;
            int accept_count = 0;
            double accept_frac = 0.0;
            int tenth = (int) R/10;
            int prog = 0;
            
            
            if (trk.mcmc.parallelizeAIES) {
                switch (trk.settings.ParallelizationBackEnd) {
                    case CPP11:
                    {
                        while (sample_count < R + burncount){
                            std::vector <std::vector <double> > XX, YY;
                            std::size_t const half_size = L / 2;
                            std::vector <std::vector <double> >  lo(all_walkers.begin(), all_walkers.begin() + half_size);
                            std::vector <std::vector <double> >  hi(all_walkers.begin() + half_size, all_walkers.end());
                            
                            std::vector <std::vector <std::vector <double> > > halves = {lo, hi};
                            for (int q = 0; q < 2; q++){ // split walkers into two sets S^q: S^0 is first half of all, S^1 is second half of all.
                                // split all walkers in half:
                                XX = halves[q];
                                YY = halves[(q + 1) % 2];
                                
                                int counter = 0, completedThreads = 0, liveThreads = 0;
                                std::vector<double> res;
                                std::vector< std::future < std::vector < double > > > futureVec;
                                futureVec.resize((int) trk.bigM);

                                for (int i = 0; i < (int) trk.bigM; i++)
                                {
                                    futureVec[i] = std::async(std::launch::async, &TRK::MCMC::parallelUpdateAIESWalkers, this, XX, YY, i); //pointer to fn run through MT, arguments to fn
                                    counter++;
                                    liveThreads++;

                                    if (liveThreads >= trk.settings.maxThreads)
                                    {
                                        for (int i = completedThreads; i < counter; i++)
                                        {
                                            res = futureVec[i].get();
                                            
                                            if (res.size() > trk.bigM) { // rejected
                                                res.pop_back();
                                                result.push_back(res);
                                            }
                                            else { // accepted
                                                all_walkers[i] = res;
                                                result.push_back(res);
                                                accept_count++;
                                               
                                            }
                                            sample_count += 1;
                                        }
                                        completedThreads += liveThreads;
                                        liveThreads = 0;
                                    }
                                }
                                for (int i = completedThreads; i < trk.bigM; i++)
                                {
                                    res = futureVec[i].get();
                                    
                                    if (res.size() > trk.bigM) { // rejected
                                        res.pop_back();
                                        result.push_back(res);
                                    }
                                    else { // accepted
                                        all_walkers[i] = res;
                                        result.push_back(res);
                                        accept_count++;
                                       
                                    }
                                    sample_count += 1;
                                }
                            }
                            
                            // acceptance fraction progress bar
                            if (sample_count % tenth == 0){
                                accept_frac = (double) accept_count / (double) sample_count;
                                if (trk.mcmc.verbose){
                                    printf("Posterior sampling %i%% complete, acceptance fraction currently %0.3f...\n", prog, accept_frac);
                                }
                                prog += 10;
                            }
                        }
                        break;
                    }
                    case OPENMP:
                    {
                        break;
                    }
                    default:
                    {
                        break;
                    }
                }
            } else {
                while (sample_count < R + burncount){
                    // NON_PARALLELIZED: update each walker one-at-a-time
                    for (int k = 0; k < L; k++){
                        // for each kth walker X,
                        X = all_walkers[k];
                        
                        // all walkers excluding X:
                        YY = all_walkers;
                        YY.erase(YY.begin() + k);
                        
                        std::vector <double> res = updateAIESWalker(X, YY);
                        
                        if (res.size() > trk.bigM) { // rejected; add initial point to sample
                            res.pop_back();
                            result.push_back(X);
                        }
                        else { // accepted; add new point to sample
                            all_walkers[k] = res;
                            result.push_back(res);
                            accept_count++;
                           
                        }
                        sample_count += 1;
                    }
                    
                    if (printAIESWalkerEvolution){
                        for (int k = 0; k < L; k++){
                            printf("%0.3f ", all_walkers[k][0]);
                        }
                        printf("\n");
                    }
                    
                    // acceptance fraction progress bar
                    if (sample_count % tenth == 0){
                        accept_frac = (double) accept_count / (double) sample_count;
                        if (trk.mcmc.verbose){
                            printf("Posterior sampling %i%% complete, acceptance fraction currently %0.3f...\n", prog, accept_frac);
                        }
                        prog += 10;
//                        if (prog >= 70) {
//                            std::cout << std::endl;
//                        }
                    }
                }
            }
            
            if (trk.mcmc.verbose){
                printf("Final AIES acceptance ratio: %0.3f \n", accept_frac);
            }
            
            break;
        }
        case ARWMH:
        {
            std::vector <std::vector <double > > cov_i(n, std::vector <double> (n, 0.0));
            for (int j = 0; j < n; j++){
                cov_i[j][j] = std::pow(sigmas_guess[j], 2.0);
            }
            std::vector <std::vector <double > > cov_i1(n, std::vector <double> (n, 0.0));
            
            std::vector <double> allparams_trial, allparams_0; //allparams_0 is the previous step
            double log_a, rand_unif_log, accept_frac = 0.0;

            int accept_count = 0;
            int delta_count = 0;
            int tenth = (int) R/10;
            int prog = 0;
            
            std::vector <double> mu_i = trk.allparams_guess;
            std::vector <double> mu_i1(n, 0.0);
            std::vector <double> X_i = trk.allparams_guess;
            std::vector <double> X_i1;
            std::vector <double> X_trial(n, 0.0);
            double lamb = std::pow(2.38, 2.0) / n;
            double gam_i1 = 1.0;
            int i = 0;

            while (delta_count < R + burncount) {
                bool loopCheck = true;
                //create trial:
                //sample X_i
                while (loopCheck) {
                    for (int j = 0; j < trk.bigM; j++) {
                        X_trial[j] = rnorm(mu_i[j], std::sqrt(lamb * cov_i[j][j]));
                    }

                    log_a = metHastRatio(X_trial, X_i);

                    rand_unif_log = std::log(runiform(0.0, 1.0));

                    if (rand_unif_log <= log_a) { // accept
                        X_i1 = X_trial;
                        loopCheck = false;
                        delta_count += 1;
                        result.push_back(X_i1);
                        accept_count += 1;
                    }
                    else { // reject
                        delta_count += 1;
                        result.push_back(X_i);
                    }
                }

                //update proposal dist params
                
        //        X_i1 = X_trial;

                gam_i1 = 1.0/((double)(i + 1));

                for (int j = 0; j < n; j++){
                    mu_i1[j] = mu_i[j] + gam_i1*(X_i1[j] - mu_i[j]);
                }

                for (int l = 0; l < n; l++){
                    for (int m = 0; m < n; m++){
                        cov_i1[l][m] = cov_i[l][m] + gam_i1*((X_i1[l] - mu_i[l])*(X_i1[m] - mu_i[m])-cov_i[l][m]);
                    }
                }

                mu_i = mu_i1;
                cov_i = cov_i1;
                X_i = X_i1;

                i += 1;

                accept_frac = (double) accept_count / (double)delta_count;
                
                if (delta_count % tenth == 0){
                    if (trk.mcmc.verbose){
                        printf("Posterior sampling %i%% complete, acceptance fraction currently %.3e...\n", prog, accept_frac);
                    }
                    prog += 10;
                }
            }
            
            if (trk.mcmc.verbose){
                printf("final Adaptive Metropolis Hastings step-sizes/proposal dist. deviations:");
                for (int j = 0; j < trk.bigM; j++) {
                    printf("%.3e ", cov_i[j][j]);
                }
                printf("\t final Adaptive Metropolis Hastings acceptance ratio: %0.3f \n", accept_frac);
            }
            break;
        }
        default:
        {
            break;
        }
    }
            
    //cut off burn-in
    for (int i = 0; i < R; i++) {
        result_final.push_back(result[i + burncount]);
    }

    result_final = checkSlopSignMCMC(result_final);
    
    useLogPosterior = false;

    return result_final;
}


// Affine Invariant Ensemble Sampler (AIES)
std::vector <double> TRK::MCMC::updateAIESWalker(std::vector <double> X, std::vector <std::vector <double> > YY){ // X is the walker to be updated with index k, YY is the set of walkers that the complementary walker for X is randomly chosen from, i.e. X should not be in YY
    double a = 2.0; //stretch variable pdf parameter
    std::vector <double> X_trial, Y, res;
    
    // choose some other jth walker Y:
    int j = rand() % (int)YY.size();
    Y = YY[j];
    
    // make proposal vector
    double Z = rstretch(a);
    for (int i = 0; i < trk.bigM; i++){
        X_trial.push_back(Z * X[i] + (1.0 - Z) * Y[i]);
    }
    
    // accept?
    double rand_unif_log = std::log(runiform(0.0, 1.0));
    double log_a = metHastRatio(X_trial, X);
    
    if (rand_unif_log <= log_a + (trk.bigM - 1) * std::log(Z)) { // accept
        res = X_trial;
    }
    
    else {
        res = {concat(X_trial, {NAN})}; //returns vector one size too big if not accepted
    }
    
    return res;
}

std::vector <double> TRK::MCMC::parallelUpdateAIESWalkers(std::vector <std::vector <double> > XX, std::vector <std::vector <double> > YY, int k){ // XX is set of walkers that contains the kth walker that you want to evolve, YY is set of complementary walkers
    return updateAIESWalker(XX[k], YY);
}


// Adaptive Random Walk Metropolis Hastings (ARWMH) sampler
void TRK::MCMC::guessARWMHDeltas(){
    params_sigmas_guess.clear();
    for (int j = 0; j < trk.M; j++){
        params_sigmas_guess.push_back(10.0 * (double)1/trk.N);
    }
    //guessing slops
    slop_x_sigma_guess = trk.statistics.stDevUnweighted(trk.x) / 100.0;
    slop_y_sigma_guess = trk.statistics.stDevUnweighted(trk.y) / 100.0;
    
    slop_x_minus_sigma_guess = slop_x_sigma_guess;
    slop_y_minus_sigma_guess = slop_y_sigma_guess;
    
    std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

    allparams_sigmas_guess.push_back(slop_x_sigma_guess);
    allparams_sigmas_guess.push_back(slop_y_sigma_guess);

    allparams_sigmas_guess = allparams_sigmas_guess;
    
    return;
}


// random number generation
double TRK::MCMC::rnorm(double mu, double sig) {

    double rand;

    std::random_device dev;
    std::mt19937 generator(dev());

    std::normal_distribution <double> dist(mu, sig);
    rand = dist(generator);
    return rand;
}

double TRK::MCMC::runiform(double a, double b) {
    double rand;

    std::random_device dev;
    std::mt19937 generator(dev());

    std::uniform_real_distribution <double> dist(a, b);
    rand = dist(generator);
    return rand;
}

double TRK::MCMC::rstretch(double a){ // distribution to sample "stretching variable" for AIES
    // using rejection sampling
    double M = 5;
    
    double x = rnorm(1.0, 1.0);   // proposal dist q(x), with Mq(x) >= p(x)
    double u = runiform(0.0, 1.0);

    while (u > (trk.statistics.stretch_pdf(x) / (M * trk.statistics.normal(x, 1.0, 1.0)))){
        x = rnorm(1.0, 1.0);
        u = runiform(0.0, 1.0);
    }
    
    return x;
}


// tools
std::vector <std::vector <double >> TRK::MCMC::checkSlopSignMCMC(std::vector <std::vector <double >> result_final) {

    std::vector <std::vector <double >> result_final_fixed;
    std::vector <double> inner;

    for (int i = 0; i < result_final.size(); i++) {
        inner.clear();

        for (int j = 0; j < trk.M; j++) {
            inner.push_back(result_final[i][j]);
        }
        
        inner.push_back(std::abs(result_final[i][trk.M]));
        
        if (!trk.settings.do1DFit){
            inner.push_back(std::abs(result_final[i][trk.M+1]));
        }
        inner.push_back(std::abs(result_final[i][trk.M+1]));

        result_final_fixed.push_back(inner);
    }

    return result_final_fixed;
}

std::vector <double> TRK::MCMC::pegToNonZeroDelta(std::vector <double> vertex, std::vector <double> lastvertex) {

    std::vector <double> vertexfixed = vertex;

    for (int j = 0; j < trk.M + 2; j++) {
        if (vertex[j] == 0.0) {
            vertexfixed[j] = lastvertex[j] * 0.5;
        }
    }

    return vertexfixed;
}

// ###########################################################################################################################



// CORRELATION REMOVAL/ PIVOT POINT METHODS ##################################################################################


// core
void TRK::CorrelationRemoval::findPivots() {
    if (findPivotPoints) {
        printf("Finding pivot point(s)...\n\n");
        
        switch (thisPivotMethod){
            case DIST:
            { // generate a weighted distribution of possible new pivot points and characterize it to find new ones
                optimizePivots_Distribution();
                break;
            }
            case REGRESSION:
            { // fit line to correlation ellipse between intercept vs slope
                optimizePivots_Regression();
                break;
            }
            case PEARSON_GSS | PEARSON_SIMPLEX:
            { // find pivot(s) that minimize abs(pearson correlation) for each set of intercept and slope
                optimizePivots_Correlation();
                break;
            }
            default:
            {
                break;
            }
        }
        
        if (trk.correlationRemoval.verbose){
            for (int p = 0; p < P; p++){
                printf("final pivot point %i: %.3e \n", p, pivots[p]);
            }
        }
        
        trk.results.pivots = pivots;
        
        pivotPointActive = false;

        return;
    } else {
        // just use whatever guess was provided
        trk.results.pivots = pivots;
        return;
    }
}


// combinations
void TRK::CorrelationRemoval::getCombos(std::vector <std::vector <double> > total, int k, int offset) { //ND case in x

    if (k == trk.M) {
        NDcombos.clear();
    }
    if (k == 0) {
        NDcombos.push_back(NDcombination);
        return;
    }
    for (int i = offset; i <= total.size() - k; ++i) {
        NDcombination.push_back(total[i]);
        getCombos(total, k - 1, i + 1);
        NDcombination.pop_back();
    }
}

std::vector < std::vector <std::vector <double > > > TRK::CorrelationRemoval::directCombos(std::vector < std::vector <double> > params_sample, int comboCount){
    std::vector < std::vector <std::vector <double > > > combos;
    combos.clear();
    
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, (int)params_sample.size() - 1); // define the range
    
    for (int j = 0; j < comboCount; j++){
        int i1 = (int) distr(eng);
        int i2 = (int) distr(eng);

        combos.push_back({params_sample[i1], params_sample[i2]});
        
        if (params_sample[i1].size() == 0 || params_sample[i2].size() == 0 ){
            printf("%i %i \n", i1, i2);
        }
    }
    
    return combos;
}


// guesses
void TRK::CorrelationRemoval::findLinearParamIndices() // finds the indices of intercepts and slopes in vector of parameters describing model
{
    double intercept_check, slope_check;
    std::vector <double> intercepts, slopes, placeholder;
    std::vector <int> intercept_inds, slope_inds;
    for (int i = 0; i < (int) trk.allparams_guess.size(); i++){
        placeholder.push_back((double) i);
    }
    for (int p = 0; p < P; p++){
        intercepts.push_back(pivot_intercept_functions[p](placeholder));
        slopes.push_back(pivot_intercept_functions[p](placeholder));
    }
    for (int p = 0; p < P; p++){
        intercept_check = intercepts[p];
        slope_check = slopes[p];
        for (int i = 0; i < (int) placeholder.size(); i++){
            if (intercept_check == placeholder[i]) {intercept_inds.push_back(i);};
            if (slope_check == placeholder[i]) {slope_inds.push_back(i);};
        }
    }
    intercept_indices = intercept_inds;
    slope_indices = slope_inds;
    return;
}

void TRK::CorrelationRemoval::getPivotGuess(){
    if (findPivotPoints){
        P = (int) pivot_intercept_functions.size();
        
        findLinearParamIndices();
        
        if (pivots.size() == 0){
            if (get_pivot_guess){
                pivots.clear();
                
                // sort data, in case x isn't sorted by default
                std::vector <int> orderedindices = getSortedIndices(trk.x);
                std::vector <double> x_s, w_s;
                for (int l = 0; l < trk.N; l++){
                    x_s.push_back(trk.x[orderedindices[l]]);
                    w_s.push_back(trk.w[orderedindices[l]]);
                }
                
                std::vector <double> divisions = findNTiles(P - 1);
                divisions = concat(concat({x_s[0]}, divisions), {x_s[trk.N-1]}); // include first and last datapoints with divisions
                
                std::vector <double> division_x, division_w;
                double pivot;
                for (int p = 0; p < P; p++){
                    for (int i = 0; i < trk.N; i++){
                        if (x_s[i] >= divisions[p] && x_s[i] < divisions[p+1]){
                            division_x.push_back(x_s[i]);
                            division_w.push_back(w_s[i]);
                        }
                    }
                    pivot = trk.statistics.getAverage(division_x, division_w);
                    pivots.push_back(pivot);
                }
            } else {
                pivots.clear();
                for (int p = 0; p < P; p++){
                    pivots.push_back(0.0);
                }
            }
        }
        
    } else {
        P = (int) pivots.size();
    }
    return;
}

std::vector <double> TRK::CorrelationRemoval::refitAnalytic(std::vector <double> new_pivots){ // refit with new pivot point; only intercepts changed
    std::vector <double> allparams_better = trk.allparams_guess, allparams_old = trk.allparams_guess;

    double intercept_new, intercept_old, slope;

    for (int p = 0; p < P; p++){
        intercept_old = pivot_intercept_functions[p](allparams_old);
        slope = pivot_slope_functions[p](allparams_old);

        intercept_new = intercept_old + slope * (new_pivots[p] - pivots[p]);
        
        allparams_better[intercept_indices[p]] = intercept_new;
    }

    return allparams_better;
}

void TRK::CorrelationRemoval::refitWithNewPivots(std::vector <double> new_pivots){
    std::vector <double> allparams_better = refitAnalytic(new_pivots); // determine the best fit
    
    if (refit_with_simplex){
        allparams_better = trk.optimization.downhillSimplex_Fit(trk.statistics.selectedChiSq, allparams_better, trk.scaleOptimization.s, trk.optimization.showFittingSteps);
    }
    
    if (trk.correlationRemoval.verbose_refit){
        printf("re-fit for new pivot point(s); old / new params:\n");

        for (int j = 0; j < (int)allparams_better.size(); j++){
            printf("%.3e %.3e\n", trk.allparams_guess[j], allparams_better[j]);
        }
    }
    
    trk.allparams_guess = allparams_better;
    return;
}


// DIST method: find pivots by generating a distribution
void TRK::CorrelationRemoval::optimizePivots_Distribution(){
    std::vector <std::vector < std::vector <double > > > drawnCombos;
    std::vector <std::vector <double> > old_pivots_sample, allparam_samples, pivot_samples(P, std::vector<double>()), pivotWeights(P, std::vector<double>());
    std::vector <double> onePivots, oneWeights, finalPivots, allparams_better, oldPivots = pivots;
    
    for (int p = 0; p < P; p++){
        old_pivots_sample.push_back(std::vector <double>((int)(randomSampleCount*(randomSampleCount - 1)) / 2, pivots[p]));
    }
    
    for (int p = 1; p < P + 1; p++){
        finalPivots.push_back((double) p);
    }
    
    int iter = 0;

    while (true)
    {
        if (iter > 0){
            pivotPointActive = true;
        }
        
        std::vector < std::vector <double> > param_samples(sample_R, std::vector<double>());
        
        
        // 1. sample parameter space
        
        if (trk.correlationRemoval.verbose && trk.mcmc.verbose){
            printf("\nSampling for pivot points...\n");
        }

        allparam_samples = trk.mcmc.samplePosterior(sample_R, sample_burnIn, trk.mcmc.allparams_sigmas_guess); //allparam_samples is { {allparams0}, {allparams1}, ... }
    
        pivot_samples = std::vector <std::vector <double> >(P, std::vector<double>());
        pivotWeights = std::vector <std::vector <double> >(P, std::vector<double>());
        
        for (int j = 0; j < allparam_samples.size(); j++) {
            param_samples[j] = slice(allparam_samples[j], 0, (int)trk.M);
        }
        
        // 2. draw slope-intercept combos from parameter space
        if (!getCombosFromSampleDirectly) { //this option takes the ~10,000 MH samples, selects ~200 of them, then generates combos out of this subset
            printf("getCombosFromSampleDirectly=false not currently configured for multiple pivot points.\n");
        }
        else { //this option takes the ~10,000 MH samples, then generates combos directly out of this set.
            drawnCombos = directCombos(param_samples, maxCombos);
        }
        
        // 3. compute new pivot points for each of these combos
        double onePivot;
        for (int j = 0; j < drawnCombos.size(); j++) {
            onePivots.clear(); // size P
            bool has_NaN = false;
            for (int p = 0; p < P; p++){
                onePivot = pivots[p] + pivotFunc(drawnCombos[j][0], drawnCombos[j][1], p);
                has_NaN = std::isnan(onePivot);
                onePivots.push_back(onePivot);
            }
            if (!has_NaN) {
                for (int p = 0; p < P; p++){
                    pivot_samples[p].push_back(onePivots[p]);
                }
            }
        }
        
        if (pruneOutlierPivots){
            for (int p = 0; p < P; p++){
                pivot_samples[p] = removeOutlierPivots(pivot_samples[p]);
            }
        }
        
        if (trk.correlationRemoval.verbose){
            printf("Weighting pivot points...\n");
        }
        
        
        // 4. weight these new pivot points
        std::vector <std::vector <double> >all_finalPivots(P, std::vector<double>());
        std::vector <std::vector <double> >finalWeights(P, std::vector<double>());
                                    
        oneWeights.clear();
        double oneWeight;
        if (weightPivots) {
            for (int p = 0; p < P; p++){
                for (int i = 0; i < pivot_samples[0].size(); i++) {
                    oneWeight = weightPivot(drawnCombos[i][0], drawnCombos[i][1], old_pivots_sample[p], pivot_samples[p][i], p);
                    if (!std::isnan(oneWeight)){
                        all_finalPivots[p].push_back(pivot_samples[p][i]);
                        finalWeights[p].push_back(oneWeight);
                    }
                }
            }
        } else {
            for (int p = 0; p < P; p++){
                finalWeights[p] = std::vector<double>((int) pivot_samples[p].size(), 1.0);
            }
            all_finalPivots = pivot_samples;
        }
        
        pivot_samples = all_finalPivots;
        pivotWeights = finalWeights;
        
        // 5. find best new pivot point using some measure of central tendency
        
        for (int p = 0; p < P; p++){
            if (pivotMedian){
                finalPivots[p] = getMedian((int) pivot_samples[p].size(), pivotWeights[p], pivot_samples[p]);
            } else if (pivotMean){
                finalPivots[p] = trk.statistics.getAverage(pivot_samples[p], pivotWeights[p]);
            } else if (pivotHalfSampleMode){
                finalPivots[p] = trk.statistics.getMode((int) pivot_samples[p].size(), pivotWeights[p], pivot_samples[p]);
            } else { //mode
                finalPivots[p] = trk.statistics.getPeakCoord(pivot_samples[p], pivotWeights[p]);
            }
        }
        
        if (writePivots){
            std::string filename = std::string("/Users/nickk124/research/reichart/TRK/TRKrepo_public/diagnostics/") + std::string("TRKpivots") + (getCombosFromSampleDirectly ? "1" : "0") + (weightPivots ? "1_" : "0_") + std::to_string(iter) + std::string("_") + std::to_string(finalPivots[0]) + std::string(".txt");

            std::ofstream myfile(filename, std::ofstream::trunc);
            if (myfile.is_open())
            {
                for (int i = 0; i < pivot_samples[0].size(); i++) {
                    for (int p = 0; p < P; p++){
                        myfile << pivot_samples[p][i] << " " << pivotWeights[p][i] << " ";
                    }
                    myfile << "\n";
                }

                //myfile << "final pivot: " << finalPivot << "\n\n\n\n\n";

                myfile.close();
            }
            else std::cout << "Unable to open file";
        }
            
        
        
        // prints out results
        if (trk.correlationRemoval.verbose){
            for (int p = 0; p < P; p++){
                printf("new, old pivot %i = %0.3f, %0.3f\t", p + 1, finalPivots[p], pivots[p]);
            }
            std::cout << std::endl;
        }
            
        iter++;
        
        // check for iteration halt
        int tolCheck = 0;
        for (int p = 0; p < P; p++){
            if (std::abs(finalPivots[p] - pivots[p]) <= pivot_tol) {tolCheck++;};
        }
        
        if (iter >= maxPivotIter || tolCheck == P) { //termination check
            break;
        }
        
        if (refit_newPivot){refitWithNewPivots(finalPivots);};
        
        // iterate
        
        pivots = finalPivots;
    }
    
    return;
}

double TRK::CorrelationRemoval::pivotFunc(std::vector <double> params1, std::vector <double> params2, int p) {
    double a01 = pivot_intercept_functions[p](params1);
    double a11 = pivot_slope_functions[p](params1);

    double a02 = pivot_intercept_functions[p](params2);
    double a12 = pivot_slope_functions[p](params2);

    return (a02 - a01) / (a11 - a12);
}

double TRK::CorrelationRemoval::weightPivot(std::vector <double> params1, std::vector <double> params2, std::vector <double> old_pivots_sample, double newPivot, int p) {
    std::vector <double> squares(old_pivots_sample.size(), 0.0);
    double w;

    for (int i = 0; i < old_pivots_sample.size(); i++) {
        squares[i] = std::pow(newPivot - old_pivots_sample[i], 2.0);
    }

    double avg = trk.statistics.getAverage(squares);
    
    double b1 = pivot_intercept_functions[p](params1);
    double m1 = pivot_slope_functions[p](params1);

    double b2 = pivot_intercept_functions[p](params2);
    double m2 = pivot_slope_functions[p](params2);

    w = std::pow(avg / std::pow(m2 - m1, 2.0) + std::pow(b1 - b2, 2.0) / std::pow(m2 - m1, 4.0), -1.0);

    return w;
}

std::vector <double> TRK::CorrelationRemoval::removeOutlierPivots(std::vector <double> pivots){
    std::vector <double> newpivots;
    double pivot;
    
    int outCount = 0;
    
    for (int i = 0; i < pivots.size(); i++){
        pivot = pivots[i];
        if (pivot < trk.x_max + pruneWidth*trk.datawidth && pivot > trk.x_min - pruneWidth*trk.datawidth){
            newpivots.push_back(pivot);
        } else {
            outCount++;
        }
    }
    if (trk.correlationRemoval.verbose){
        printf("%i pivots outside of reasonable region \n", outCount);
    }
    
    return newpivots;
}


// REGRESSION method: find pivots using the slope of the correlatino ellipse between intercepts and slopes
void TRK::CorrelationRemoval::optimizePivots_Regression(){
    std::vector <std::vector <double> > allparam_samples;
    std::vector <double> finalPivots, allparams_better, oldPivots = pivots;
    
    for (int p = 1; p < P + 1; p++){
        finalPivots.push_back((double) p);
    }
    
    int iter = 0;

    while (true)
    {
        if (iter > 0){
            pivotPointActive = true;
        }
        
        std::vector < std::vector <double> > param_samples(sample_R, std::vector<double>());
        
        
        // 1. sample parameter space
        
        if (trk.correlationRemoval.verbose && trk.mcmc.verbose){
            printf("\nSampling for pivot points...\n");
        }

        allparam_samples = trk.mcmc.samplePosterior(sample_R, sample_burnIn, trk.mcmc.allparams_sigmas_guess); //allparam_samples is { {allparams0}, {allparams1}, ... }
        
        for (int j = 0; j < allparam_samples.size(); j++) {
            param_samples[j] = slice(allparam_samples[j], 0, (int)trk.M);
        }
        
        // 2. for each pivot point, find new one by estimating slope of correlation ellipse
        for (int p = 0; p < P; p++){
            std::vector <double> b_samples, m_samples, fit_result;
            double ellipse_slope; // slope along axis of correlation ellipse between slope and intercept
            
            for (int j = 0; j < param_samples.size(); j++) {
                b_samples.push_back(pivot_intercept_functions[p](param_samples[j]));
                m_samples.push_back(pivot_slope_functions[p](param_samples[j]));
            }
            
            fit_result = trk.statistics.simpleLinearRegression(m_samples, b_samples);
            ellipse_slope = fit_result[1];
            
            if (verbose){
                printf("ellipse slope %i = %0.3f\t\t", p, ellipse_slope);
            }
            
            finalPivots[p] = pivots[p] - ellipse_slope;
            
            if (writePivots){
                // write b vs m distribution to file
                std::string fileName = trk.settings.outputPath + std::string("/TRKMCMC_PIVOT") + std::to_string(p) + std::string("_") + std::to_string(iter) + std::string("_") + std::to_string(trk.allparams_guess[0]) + std::string("_") + std::to_string(sample_R) + std::string(".txt");

                std::ofstream myfile;
                myfile.open(fileName, std::ofstream::trunc);
                
                int m = 0;
                if (trk.asymmetric.hasAsymSlop){
                    m = 2;
                }
                
                if (myfile.is_open())
                {
                    // filename    a     b     optimum scale    total computation time (s)
                    for (int i = 0; i < b_samples.size(); i++) {
                        myfile << m_samples[i] << " " << b_samples[i] << std::endl;
                    }

                    myfile.close();
                }
                else std::cout << "Unable to open file";
            }
        }
            
        
        
        // prints out results
        if (trk.correlationRemoval.verbose){
            for (int p = 0; p < P; p++){
                printf("new, old pivot %i = %0.3f, %0.3f\t", p + 1, finalPivots[p], pivots[p]);
            }
            std::cout << std::endl;
        }
            
        iter++;
        
        // check for iteration halt
        int tolCheck = 0;
        for (int p = 0; p < P; p++){
            if (std::abs(finalPivots[p] - pivots[p]) <= pivot_tol) {tolCheck++;};
        }
        
        if (iter >= maxPivotIter || tolCheck == P) { //termination check
            break;
        }
        
        if (refit_newPivot){refitWithNewPivots(finalPivots);};
        
        // iterate
        
        pivots = finalPivots;
    }
    
    return;
}


// PEARSON method: find pivots using the correlation of intercepts and slopes and the golden section search (GSS) method
void TRK::CorrelationRemoval::optimizePivots_Correlation(){
    
    if (thisPivotMethod == PEARSON_GSS){
        // initialization for golden section search
        findPivotBrackets();
    }
    
    switch (trk.settings.ParallelizationBackEnd){ // search for multiple pivot points simultaneously
        case CPP11:
            optimizePivots_Correlation_CPP11();
            break;
        case OPENMP:
            optimizePivots_Correlation_Default();
            break;
        default:
            optimizePivots_Correlation_Default();
            break;
    }
    
    return;
}

void TRK::CorrelationRemoval::optimizePivots_Correlation_CPP11(){
    // initialize parallelization objects, to compute ret_vec in parallel
    int counter = 0, completedThreads = 0, liveThreads = 0;
    double result;
    std::vector < std::future < double > > futureVec;
    futureVec.resize(P);
    
                    
    for (int p = 0; p < P; p++){
        switch (thisPivotMethod){
            case PEARSON_SIMPLEX : {
                // function to be minimized:
                std::function <double(std::vector <double>)> correlation_func_simplex = std::bind(&TRK::CorrelationRemoval::getAbsCorrFromNewPivot_Wrapper, this, std::placeholders::_1, p);
                
                std::vector <double> guess = {pivots[p]};
                
                futureVec[p] = std::async(std::launch::async, &TRK::Optimization::downhillSimplex_1DWrapper, trk.optimization, correlation_func_simplex, guess, correlation_tol, showSimplexSteps);
                break;
            }
            case PEARSON_GSS : {
                // function to be minimized:
                std::function <double(double)> correlation_func_GSS = std::bind(&TRK::CorrelationRemoval::getAbsCorrFromNewPivot, this, std::placeholders::_1, p);
                
                futureVec[p] = std::async(std::launch::async, &TRK::Optimization::goldenSectionSearch, trk.optimization, correlation_func_GSS, min_pivots_brackets[p], max_pivots_brackets[p], correlation_tol);
                break;
            }
            case DIST: {
                break;
            }
            case REGRESSION: {
                break;
            }
            default : {
                break;
            }
        }
        counter++;
        liveThreads++;
        
        if (liveThreads >= trk.settings.maxThreads)
        {
            for (int i = completedThreads; i < counter; i++)
            {
                result = futureVec[i].get();
                pivots[p] = result;
            }
            completedThreads += liveThreads;
            liveThreads = 0;
        }
    }
    for (int p = completedThreads; p < P; p++)
    {
        result = futureVec[p].get();
        pivots[p] = result;
    }
    
    return;
}

void TRK::CorrelationRemoval::optimizePivots_Correlation_Default(){
    #pragma omp parallel for num_threads(maxThreads)
    for (int p = 0; p < P; p++){ // each pivot can be optimized independently (the value of the others pivots shouldn't affect it
        switch (thisPivotMethod){
            case PEARSON_GSS : {
                // function to be minimized:
                std::function <double(double)> correlation_func = std::bind(&TRK::CorrelationRemoval::getAbsCorrFromNewPivot, this, std::placeholders::_1, p);
                
                pivots[p] = trk.optimization.goldenSectionSearch(correlation_func, min_pivots_brackets[p], max_pivots_brackets[p], correlation_tol);
                break;
            }
            case PEARSON_SIMPLEX : {
                // function to be minimized:
                std::function <double(std::vector <double>)> correlation_func = std::bind(&TRK::CorrelationRemoval::getAbsCorrFromNewPivot_Wrapper, this, std::placeholders::_1, p);
                
                pivots[p] = trk.optimization.downhillSimplex(correlation_func, {pivots[p]}, correlation_tol, showSimplexSteps)[0];
                break;
            }
            default : {
                break;
            }
        }
    }
    return;
}

void TRK::CorrelationRemoval::writeCorrelationOptimizationSampling(std::vector <double> b_samples, std::vector <double> m_samples, int p)
{
    if (writePivots){
        // write b vs m distribution to file
        std::string fileName = trk.settings.outputPath + std::string("/TRKMCMC_PIVOT") + std::to_string(p + 1) + std::string("_") + std::to_string(trk.mcmc.runiform(0.0, 1.0)) + std::string("_") + std::to_string(trk.allparams_guess[0]) + std::string("_") + std::to_string(sample_R) + std::string(".txt");

        std::ofstream myfile;
        myfile.open(fileName, std::ofstream::trunc);
        
        int m = 0;
        if (trk.asymmetric.hasAsymSlop){
            m = 2;
        }
        
        if (myfile.is_open())
        {
            // filename    a     b     optimum scale    total computation time (s)
            for (int i = 0; i < b_samples.size(); i++) {
                myfile << m_samples[i] << " " << b_samples[i] << std::endl;
            }

            myfile.close();
        }
        else std::cout << "Unable to open file";
    }
                                                               
    return;
}

double TRK::CorrelationRemoval::getAbsCorrFromNewPivot(double new_pivot, int p){
    std::vector < std::vector <double> > param_samples(sample_R, std::vector<double>());
    std::vector <double> b_samples, m_samples, original_pivots = pivots;
    double rxy, abs_rxy; // correlation between slope and intercept
    
    pivots[p] = new_pivot; // set pivot to new value, and sample parameter space with it
    
    if (refit_newPivot){refitWithNewPivots(pivots);}; // find better starting place for sampling
    
    if (trk.correlationRemoval.verbose && trk.mcmc.verbose){
        printf("\nSampling for pivot points...\n");
    }
    std::vector < std::vector <double> > allparam_samples = trk.mcmc.samplePosterior(sample_R, sample_burnIn, trk.mcmc.allparams_sigmas_guess);
    for (int j = 0; j < allparam_samples.size(); j++) {
        param_samples[j] = slice(allparam_samples[j], 0, (int)trk.M);
    }
    for (int j = 0; j < param_samples.size(); j++) {
        b_samples.push_back(pivot_intercept_functions[p](param_samples[j]));
        m_samples.push_back(pivot_slope_functions[p](param_samples[j]));
    }
    
//    printf("Currently using spearman...\n");
    rxy = trk.statistics.pearsonCorrelation(b_samples, m_samples); // compute pearson correlation
//    rxy = trk.statistics.spearmanCorrelation(b_samples, m_samples); // compute spearman correlation
    
    
    
    double spearmanR = trk.statistics.spearmanCorrelation(b_samples, m_samples);
    double pearsonR = trk.statistics.pearsonCorrelation(b_samples, m_samples);
    
    if (std::abs(spearmanR - pearsonR) > 0.2){
        printf("differing correlations; pearson = %f\t spearman = %f\n", trk.statistics.pearsonCorrelation(b_samples, m_samples), trk.statistics.spearmanCorrelation(b_samples, m_samples));
    }
    
    
    abs_rxy = std::isnan(rxy) ? 1.0 : std::abs(rxy); // returns maximally correlated if NaN
    
    writeCorrelationOptimizationSampling(b_samples, m_samples, p);
    
    pivots = original_pivots; //reset pivots to previous values
    
    return abs_rxy; // returns maximally correlated if NaN
}

double TRK::CorrelationRemoval::getAbsCorrFromNewPivot_Wrapper(std::vector <double> new_pivot, int p){
    return getAbsCorrFromNewPivot(new_pivot[0], p);
}


// tools
std::vector <double> TRK::CorrelationRemoval::findNTiles(int Q){ // Q divisions; Q + 1 regions
    std::vector <double> divisions;
    // sort data, in case x isn't sorted by default
    std::vector <int> orderedindices = getSortedIndices(trk.x);
    std::vector <double> x_s, w_s;
    for (int l = 0; l < trk.N; l++){
        x_s.push_back(trk.x[orderedindices[l]]);
        w_s.push_back(trk.w[orderedindices[l]]);
    }
    
    double total = 0.0;
    for (int i = 0; i < trk.N; i++){
        total += x_s[i] * w_s[i];
    }
    
    double runningSum, division;
    for (int q = 1; q < Q + 1; q++){
        runningSum = 0.0;
        int i = 0;
        while (runningSum <= total * q / (Q + 1)){
            runningSum += x_s[i] * w_s[i];
            i++;
        }
        division = (x_s[i-1] + x_s[i-2]) / 2.0;
        divisions.push_back(division);
    }
    
    return divisions;
}

void TRK::CorrelationRemoval::findPivotBrackets(){
    double b1, b2, m1, m2, xp1, xp2, intersection;
    std::vector <double> b, a = {trk.x_min}; // a is min(s), b is maxe(s)
    
    for (int q = 0; q < P - 1; q++){ // qth intersection between lines
        b1 = pivot_intercept_functions[q](trk.allparams_guess);
        b2 = pivot_intercept_functions[q + 1](trk.allparams_guess);
        m1 = pivot_slope_functions[q](trk.allparams_guess);
        m2 = pivot_slope_functions[q + 1](trk.allparams_guess);
        xp1 = pivots[q];
        xp2 = pivots[q + 1];
        
        intersection = (b2 - b1 + m1 * xp1 - m2 * xp2) / (m1 - m2);
        a.push_back(intersection);
        b.push_back(intersection);
    }
    
    b.push_back(trk.x_max);
    
    min_pivots_brackets = a;
    max_pivots_brackets = b;
    
    return;
}

// ###########################################################################################################################



// ASYMMETRIC UNCERTAINTIES ##################################################################################################

// tools
void TRK::Asymmetric::checkAsym(){ //checks to see whether any or all of the asymmetric error bar and slop parameters were provided.
    
    // SLOP
    
    bool negXSlop = false;
    bool negYSlop = false;
    
    if (slop_x_minus_guess >= 0){
        negXSlop = true;
    }

    if (slop_y_minus_guess >= 0){
        negYSlop = true;
    }
    
    if (negXSlop && !negYSlop){
        hasAsymSlop = true;
        slop_y_minus_guess = trk.slop_y_guess;
    }
    
    else if (!negXSlop && negYSlop){
        hasAsymSlop = true;
        slop_x_minus_guess = trk.slop_x_guess;
    }
    
    else if (negXSlop && negYSlop){
        hasAsymSlop = true;
    }
    
    
    // ERROR BARS
    
    bool negXEB = false;
    bool negYEB = false;
    
    if (sx_minus.size() == trk.N){
        negXEB = true;
    }
    
    if (sy_minus.size() == trk.N){
       negYEB = true;
    }
    
    
    if (negXEB && !negYEB){
        hasAsymEB = true;
        sy_minus = trk.sy;
    }
    
    else if (!negXEB && negYEB){
        hasAsymEB = true;
        sx_minus = trk.sx;
    }
    
    else if (negXEB && negYEB){
        hasAsymEB = true;
    }
    
    if (hasAsymSlop || hasAsymEB){
        trk.statistics.selectedChiSq = &TRK::Statistics::modifiedChiSquaredAsym;
        trk.statistics.selectedLikelihood = &TRK::Statistics::likelihoodAsym;
    }
    
    if (verbose){
        printf("Asymmetries: slop: %s\tError bars: %s\n", hasAsymSlop ? "true" : "false", hasAsymEB ? "true" : "false");
    }
    
    if (hasAsymSlop){
        trk.allparams_guess.push_back(slop_x_minus_guess);
        trk.allparams_guess.push_back(slop_y_minus_guess);
    }
    
    if (trk.settings.do1DFit){
        trk.statistics.selectedChiSq = &TRK::Statistics::regularChiSquaredWSlop;
        trk.statistics.selectedLikelihood = &TRK::Statistics::likelihood1D;
        printf("Running 1D fit.\n");
    }
    
    return;
}


// likelihoods/posteriors
double TRK::Asymmetric::dunDxAsym(double mtn, std::vector <double> Sigs2, int quadSig_xn2Ind, int quadSig_yn2Ind, double s){
    double quadSigX2 = Sigs2[quadSig_xn2Ind]; //the correct Sig2 values for the specific quadrant
    double quadSigY2 = Sigs2[quadSig_yn2Ind];

    return (std::pow(mtn, 2.0)*quadSigX2 + quadSigY2) / (std::sqrt(std::pow(mtn*quadSigX2, 2.0) + std::pow(s*quadSigY2, 2.0)));
}

double TRK::Statistics::cmNorm(double z){
    return (std::sqrt(PI)/2.0) * (1.0 + std::erf(z));
}

double TRK::Asymmetric::zAsym(double x, double quadSig_xn2, double quadSig_yn2, double xn_shifted, double yn_shifted, std::vector <double> shifts, double x_tn, double y_tn, double m_tn){ //equation B-6 of thesis
    
    return (quadSig_yn2 * (x - xn_shifted) + std::pow(m_tn, 2.0) * quadSig_xn2 * (x - x_tn - (yn_shifted - y_tn)/m_tn)) / (std::sqrt(quadSig_xn2) * std::sqrt(quadSig_yn2) * std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2));
}

double TRK::Asymmetric::pnAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, double x_tn, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn){
    
    // INITIALIZATIONS
    
    double m_tn = trk.dyc(x_tn, params);
    double y_tn = trk.yc(x_tn, params);
    
    double x1 = xn_shifted;
    double x2 = ((yn_shifted - y_tn)/m_tn) + x_tn;
    
    double dunDx = dunDxAsym(m_tn, Sigs2, quadSig_xn2Ind, quadSig_yn2Ind, s);
    double norm = 2.0 / (PI * (std::sqrt(Sigs2[0]) + std::sqrt(Sigs2[2])) * (std::sqrt(Sigs2[1]) + std::sqrt(Sigs2[3]))); //normalization factor, same for all three integrals
    
    // INTEGRALS
    
    double I1 = 0;
    double I2 = 0;
    double I3 = 0;
    
    std::vector <int> I1inds = {0, 0};
    std::vector <int> I2inds = {quadSig_xn2Ind, quadSig_yn2Ind};
    std::vector <int> I3inds = {0, 0};
    
    
    if (I2inds[0] == 0 && I2inds[1] == 1){
        I1inds = {2, 1};
        I3inds = {0, 3};
    } else if (I2inds[0] == 2 && I2inds[1] == 1){
        I1inds = {2, 3};
        I3inds = {0, 1};
    } else if (I2inds[0] == 2 && I2inds[1] == 3){
        I1inds = {2, 1};
        I3inds = {0, 3};
    } else if (I2inds[0] == 0 && I2inds[1] == 3){
        I1inds = {2, 3};
        I3inds = {0, 1};
    }
    
    // Do ``Integrals''
    
    double quadSig_xn2 = Sigs2[I1inds[0]];
    double quadSig_yn2 = Sigs2[I1inds[1]];
    
    double z1 = zAsym(x1, quadSig_xn2, quadSig_yn2, xn_shifted, yn_shifted, shifts, x_tn, y_tn, m_tn);
    
    I1 = dunDx * norm * std::sqrt(quadSig_xn2) * std::sqrt(quadSig_yn2) * trk.statistics.cmNorm(z1) * std::exp(-0.5 * std::pow((y_tn - yn_shifted - m_tn*(x_tn - xn_shifted))/std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2), 2.0)) / std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2);
    
    
    
    quadSig_xn2 = Sigs2[I2inds[0]];
    quadSig_yn2 = Sigs2[I2inds[1]];
    
    double z2 = zAsym(x2, quadSig_xn2, quadSig_yn2, xn_shifted, yn_shifted, shifts, x_tn, y_tn, m_tn);
    z1 = zAsym(x1, quadSig_xn2, quadSig_yn2, xn_shifted, yn_shifted, shifts, x_tn, y_tn, m_tn);
    
    I2 = dunDx * norm * std::sqrt(quadSig_xn2) * std::sqrt(quadSig_yn2) * (trk.statistics.cmNorm(z2) - trk.statistics.cmNorm(z1)) * std::exp(-0.5 * std::pow((y_tn - yn_shifted - m_tn*(x_tn - xn_shifted))/std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2), 2.0)) / std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2);
    
    
    
    quadSig_xn2 = Sigs2[I3inds[0]];
    quadSig_yn2 = Sigs2[I3inds[1]];
    
    z2 = zAsym(x2, quadSig_xn2, quadSig_yn2, xn_shifted, yn_shifted, shifts, x_tn, y_tn, m_tn);
    
    I3 = dunDx * norm * std::sqrt(quadSig_xn2) * std::sqrt(quadSig_yn2) * (1.0 - trk.statistics.cmNorm(z2)) * std::exp(-0.5 * std::pow((y_tn - yn_shifted - m_tn*(x_tn - xn_shifted))/std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2), 2.0)) / std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2);
    
    
    
    return I1 + I2 + I3;
    
}

double TRK::Asymmetric::singlePointLnLAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, double x_tn, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn){
    
    return -2.0 * std::log(pnAsym(params, xn_shifted, yn_shifted, Sigs2, x_tn, quadSig_xn2Ind, quadSig_yn2Ind, shifts, s, wn));
}


// tangent points
double TRK::Asymmetric::findBestTangentAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, std::vector <double> x_tn_vec, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn) {
    std::vector <double> posts;
    long minindex;

    for (int i = 0; i < x_tn_vec.size(); i++) {
        posts.push_back(singlePointLnLAsym(params, xn_shifted, yn_shifted, Sigs2, x_tn_vec[i], quadSig_xn2Ind, quadSig_yn2Ind, shifts, s, wn));
    }

    std::vector<double>::iterator result = std::min_element(std::begin(posts), std::end(posts));
    minindex = std::distance(std::begin(posts), result);

    return x_tn_vec[minindex];
}

std::vector <double> TRK::Asymmetric::getAsymShifts(std::vector <double> allparams, int n){
    double deltayn = 0.0;
    double deltaxn = 0.0;
    
    std::vector <double> slops = {allparams[trk.M], allparams[trk.M+1]};
    std::vector <double> EBs = {trk.sx[n], trk.sy[n]};
    
    if (hasAsymSlop && !hasAsymEB){
        slops.push_back(allparams[trk.M+2]);
        slops.push_back(allparams[trk.M+3]);
        
        EBs = concat(EBs, EBs);
        
    } else if (!hasAsymSlop && hasAsymEB){
        slops = concat(slops, slops);
        
        EBs.push_back(sx_minus[n]);
        EBs.push_back(sy_minus[n]);
        
    } else if (hasAsymSlop && hasAsymEB){
        slops.push_back(allparams[trk.M+2]);
        slops.push_back(allparams[trk.M+3]);
        
        EBs.push_back(trk.asymmetric.sx_minus[n]);
        EBs.push_back(trk.asymmetric.sy_minus[n]);
    }
    
    // Y SHIFT
    
    double sigmaL = minMax({slops[1], slops[3]})[1];
    double sigmaS = minMax({slops[1], slops[3]})[0];
    double sigmanL = minMax({EBs[1], EBs[3]})[1];
    double sigmanS = minMax({EBs[1], EBs[3]})[0];
    double sigmaMax = minMax({sigmaL, sigmanL})[1];
    
    
    double xi = sigmaS/sigmaL + sigmanS / sigmanL;
    double eta = sigmanL < sigmaL ? sigmanS/sigmanL - sigmaS/sigmaL : sigmaS/sigmaL - sigmanS/sigmanL;
    double r = minMax({sigmaL, sigmanL})[0] / minMax({sigmaL, sigmanL})[1];
    
    double xip = xi <= 1 ? xi : 2.0 - xi;
    double etap = xip == 0 ? 0 : 2.0 * xip * std::pow(0.5*eta/xip + 1.0, std::pow(r, -0.4087)) - xip;
    
    double Nr = -0.5326*std::pow(r,2.0) + 1.5307*r + 0.0019;
    double fXi = xi <= 1 ? 0.2454*std::pow(xi, -1.1452) : 0.2454*std::pow(xi, -0.5203);
    double gEtaP = std::pow(etap,2.0);
    double hXi = -0.042*std::pow(xi,2.0) - 0.1602 * xi + 0.4884;
    
    double deltastr = sigmaMax * Nr * (fXi * gEtaP + hXi);
    
    int i = 1;
    if (slops[1] == slops[3] || EBs[1] == EBs[3]){ //one of the dists is symmetric
        if (sigmanL == EBs[1] || sigmaL == slops[3]){
            i = 1;
        } else if (sigmaMax == EBs[3] || sigmaMax == slops[1]){
            i = -1;
        }
        
        deltayn = i * deltastr;
        
    } else if ((sigmaL == slops[3] && sigmanL == EBs[1]) || (sigmaL = slops[1] && sigmanL == EBs[3])){ //both asymm first case
        if (sigmaMax == EBs[1] || sigmaMax == slops[3]){
            i = 1;
        } else if (sigmaMax == EBs[3] || sigmaMax == slops[1]){
            i = -1;
        }
        
        deltayn = i * deltastr;
        
    } else if ((sigmaL == slops[1] && sigmanL == EBs[1]) || (sigmaL = slops[3] && sigmanL == EBs[3])){ //both asymm second case
        if (sigmaMax == EBs[1] || sigmaMax == slops[3]){
            i = 1;
        } else if (sigmaMax == EBs[3] || sigmaMax == slops[1]){
            i = -1;
        }
        
        double pwr = eta <= 1 ? 0.7413 : -0.1268;
        deltayn = i * deltastr * std::sin(PI/2.0 * etap/xip) * std::pow(eta, pwr);
    }
    
    // X SHIFT
    
    sigmaL = minMax({slops[0], slops[2]})[0];
    sigmaS = minMax({slops[0], slops[2]})[0];
    sigmanL = minMax({EBs[0], EBs[2]})[0];
    sigmanS = minMax({EBs[0], EBs[2]})[0];
    sigmaMax = minMax({sigmaL, sigmanL})[0];
    
    
    xi = sigmaS/sigmaL + sigmanS / sigmanL;
    eta = sigmanL < sigmaL ? sigmanS/sigmanL - sigmaS/sigmaL : sigmaS/sigmaL - sigmanS/sigmanL;
    r = minMax({sigmaL, sigmanL})[0] / minMax({sigmaL, sigmanL})[0];
    
    xip = xi <= 1 ? xi : 2.0 - xi;
    etap = xip == 0 ? 0 : 2.0 * xip * std::pow(0.5*eta/xip + 1.0, std::pow(r, -0.4087)) - xip;
    
    Nr = -0.5326*std::pow(r,2.0) + 1.5307*r + 0.0019;
    fXi = xi <= 1 ? 0.2454*std::pow(xi, -1.1452) : 0.2454*std::pow(xi, -0.5203);
    gEtaP = std::pow(etap,2.0);
    hXi = -0.042*std::pow(xi,2.0) - 0.1602 * xi + 0.4884;
    
    deltastr = sigmaMax * Nr * (fXi * gEtaP + hXi);
    
    if (slops[0] == slops[2] || EBs[0] == EBs[2]){ //one of the dists is symmetric
        if (sigmanL == EBs[0] || sigmaL == slops[2]){
            i = 1;
        } else if (sigmaMax == EBs[2] || sigmaMax == slops[0]){
            i = -1;
        }
        
        deltaxn = i * deltastr;
        
    } else if ((sigmaL == slops[2] && sigmanL == EBs[0]) || (sigmaL = slops[0] && sigmanL == EBs[2])){ //both asymm first case
        if (sigmaMax == EBs[0] || sigmaMax == slops[2]){
            i = 1;
        } else if (sigmaMax == EBs[2] || sigmaMax == slops[0]){
            i = -1;
        }
        
        deltaxn = i * deltastr;
        
    } else if ((sigmaL == slops[0] && sigmanL == EBs[0]) || (sigmaL = slops[2] && sigmanL == EBs[2])){ //both asymm second case
        if (sigmaMax == EBs[0] || sigmaMax == slops[2]){
            i = 1;
        } else if (sigmaMax == EBs[2] || sigmaMax == slops[0]){
            i = -1;
        }
        
        double pwr = eta <= 1 ? 0.7413 : -0.1268;
        deltaxn = i * deltastr * std::sin(PI/2.0 * etap/xip) * std::pow(eta, pwr);
    }
    
    
    return {deltaxn, deltayn};
}

std::vector <double> TRK::Asymmetric::getAsymSigs2(std::vector <double> allparams, int n){
    std::vector <double> Sigs2(4, 0.0);
    std::vector <double> slops = {allparams[trk.M], allparams[trk.M+1]};
    std::vector <double> EBs = {trk.sx[n], trk.sy[n]};
    
    if (hasAsymSlop && !hasAsymEB){
        slops = concat({allparams[trk.M+2], allparams[trk.M+3]}, slops);
        
        EBs = concat(EBs, EBs);
        
    } else if (!hasAsymSlop && hasAsymEB){
        slops = concat(slops, slops);
        
        EBs.push_back(sx_minus[n]);
        EBs.push_back(sy_minus[n]);
        
    } else if (hasAsymSlop && hasAsymEB){
        slops = concat({allparams[trk.M+2], allparams[trk.M+3]}, slops);
        
        EBs.push_back(sx_minus[n]);
        EBs.push_back(sy_minus[n]);
    }
    
    for (int i = 0; i < 4; i++){
        Sigs2[i] = std::pow(slops[i], 2.0) + std::pow(EBs[i], 2.0);
    }
        
    return Sigs2;
}

std::vector <double> TRK::Asymmetric::tangentParallelAsym(std::vector<double> allparams, int n, double s) {
    std::vector <double> params;
    for (int i = 0; i < trk.M; i++) {
        params.push_back(allparams[i]);
    }
    
    // shift centroid
    std::vector <double> shifts = getAsymShifts(allparams, n);
    double xn_shifted = trk.x[n] + shifts[0];
    double yn_shifted = trk.y[n] + shifts[1];
    
    
    // choose correct Sigmas
    std::vector <double> Sigs2 = getAsymSigs2(allparams, n);
    double quadSig_xn2 = 0.0;
    double quadSig_yn2 = 0.0;
    int quadSig_xn2Ind = -1;
    int quadSig_yn2Ind = -1;
    
    double yCxn = trk.yc(xn_shifted, params);
    double dyCxn = trk.dyc(xn_shifted, params);
    int quadrant = 0;
    
    if (yCxn > yn_shifted){
        if (dyCxn < 0){ // Quadrant I
            quadSig_xn2 = Sigs2[0];
            quadSig_yn2 = Sigs2[1];
            quadrant = 1;
            quadSig_xn2Ind = 0;
            quadSig_yn2Ind = 1;
        } else if (dyCxn >= 0){ // QII
            quadSig_xn2 = Sigs2[2];
            quadSig_yn2 = Sigs2[1];
            quadrant = 2;
            quadSig_xn2Ind = 2;
            quadSig_yn2Ind = 1;
        }
    } else if (yCxn <= yn_shifted){
        if (dyCxn < 0){ // QIII
            quadSig_xn2 = Sigs2[2];
            quadSig_yn2 = Sigs2[3];
            quadrant = 3;
            quadSig_xn2Ind = 2;
            quadSig_yn2Ind = 3;
        } else if (dyCxn >= 0){ // QIV
            quadSig_xn2 = Sigs2[0];
            quadSig_yn2 = Sigs2[3];
            quadrant = 4;
            quadSig_xn2Ind = 0;
            quadSig_yn2Ind = 3;
        }
    }
    
    // find tangent point(s)
    
    std::vector <double> x_tn_vec = trk.tangentPointMethods.tangentsFinder(params, xn_shifted, yn_shifted, quadSig_xn2, quadSig_yn2, xn_shifted); // we use x_n as the initial guess for this. gives the three closest tangest points

    double x_t = findBestTangentAsym(params, xn_shifted, yn_shifted, Sigs2, x_tn_vec, quadSig_xn2Ind, quadSig_yn2Ind, shifts, s, trk.w[n]);

    double subsum = std::log(pnAsym(params, xn_shifted, yn_shifted, Sigs2, x_t, quadSig_xn2Ind, quadSig_yn2Ind, shifts, s, trk.w[n]));

    return {x_t, subsum}; // returns x_t_n and ln(p_n)
}

std::vector <double> TRK::Asymmetric::tangentParallelLikelihoodAsym(std::vector<double> allparams, int n) {
    std::vector <double> params;
    for (int i = 0; i < trk.M; i++) {
        params.push_back(allparams[i]);
    }
    
    // shift centroid
    std::vector <double> shifts = getAsymShifts(allparams, n);
    double xn_shifted = trk.x[n] + shifts[0];
    double yn_shifted = trk.y[n] + shifts[1];
    
    
    // choose correct Sigmas
    std::vector <double> Sigs2 = getAsymSigs2(allparams, n);
    double quadSig_xn2 = 0.0;
    double quadSig_yn2 = 0.0;
    int quadSig_xn2Ind = -1;
    int quadSig_yn2Ind = -1;
    
    double yCxn = trk.yc(xn_shifted, params);
    double dyCxn = trk.dyc(xn_shifted, params);
    int quadrant = 0;
    
    if (yCxn > yn_shifted){
        if (dyCxn < 0){ // Quadrant I
            quadSig_xn2 = Sigs2[0];
            quadSig_yn2 = Sigs2[1];
            quadrant = 1;
            quadSig_xn2Ind = 0;
            quadSig_yn2Ind = 1;
        } else if (dyCxn >= 0){ // QII
            quadSig_xn2 = Sigs2[2];
            quadSig_yn2 = Sigs2[1];
            quadrant = 2;
            quadSig_xn2Ind = 2;
            quadSig_yn2Ind = 1;
        }
    } else if (yCxn <= yn_shifted){
        if (dyCxn < 0){ // QIII
            quadSig_xn2 = Sigs2[2];
            quadSig_yn2 = Sigs2[3];
            quadrant = 3;
            quadSig_xn2Ind = 2;
            quadSig_yn2Ind = 3;
        } else if (dyCxn >= 0){ // QIV
            quadSig_xn2 = Sigs2[0];
            quadSig_yn2 = Sigs2[3];
            quadrant = 4;
            quadSig_xn2Ind = 0;
            quadSig_yn2Ind = 3;
        }
    }
    
    // find tangent point(s)
    
    std::vector <double> x_tn_vec = trk.tangentPointMethods.tangentsFinder(params, xn_shifted, yn_shifted, quadSig_xn2, quadSig_yn2, xn_shifted); // we use x_n as the initial guess for this. gives the three closest tangest points

    double x_t = findBestTangentAsym(params, xn_shifted, yn_shifted, Sigs2, x_tn_vec, quadSig_xn2Ind, quadSig_yn2Ind, shifts, trk.scaleOptimization.s, trk.w[n]);

    double l = pnAsym(params, xn_shifted, yn_shifted, Sigs2, x_t, quadSig_xn2Ind, quadSig_yn2Ind, shifts, trk.scaleOptimization.s, trk.w[n]);

    return { x_t, l};
}


// ###########################################################################################################################



// GLOBAL FUNCTIONS/TOOLS ####################################################################################################

// mcmc/sampling
double noPrior(double param) {
    return 1.0;
}


// statistics
std::vector <double> minMax(std::vector <double> vec) {

    // Finding the smallest of all the numbers
    double min = *std::min_element(std::begin(vec), std::end(vec));
    double max = *std::max_element(std::begin(vec), std::end(vec));

    return { min, max };
}

std::vector <int> argMinMax(std::vector <double> x){
    int argMin = (int)std::distance(x.begin(), std::min_element(x.begin(), x.end()));
    int argMax = (int)std::distance(x.begin(), std::max_element(x.begin(), x.end()));
    return {argMin, argMax};
}

std::vector <int> getSortedIndices(std::vector <double> x)
{
    std::vector<int> y(x.size());
    std::size_t n(0);
    std::generate(std::begin(y), std::end(y), [&] { return n++; });

    std::sort(std::begin(y),
        std::end(y),
        [&](int i1, int i2) { return x[i1] < x[i2]; });

    return y;
}

double getMedian(std::vector<double> y)
{
    int high = (int)(floor(y.size() / 2));
    int low = high - 1;
    double runningSum = 0, median = 0;
    double totalSum = y.size();
    if (y.size() > 1)
    {
        if (y.size() % 2 == 0)
        {
            runningSum = y.size() / 2.0 + .5;
        }
        else
        {
            runningSum = y.size() / 2.0;
        }
        median = y[low] + (.5*totalSum - runningSum + 1.0)* (y[high] - y[low]);
    }

    else
    {
        median = y[0];
    }
    return median;

}

double getMedian(int trueCount, std::vector<double> w, std::vector<double> y)
{
    size_t sumCounter = 0;
    double median = 0, totalSum = 0, runningSum = 0;
    for (int i = 0; i < trueCount; i++)
    {
        totalSum += w[i];
    }
    if (trueCount > 1)
    {
        runningSum = w[sumCounter] * .5;
        while (runningSum < .5*totalSum)
        {
            sumCounter++;
            runningSum += w[sumCounter - 1] * .5 + w[sumCounter] * .5;
        }
        if (sumCounter == 0)
        {
            median = y[0];
//            std::cout << median << std::endl;
        }
        else
        {
            median = y[sumCounter - 1] + (.5*totalSum - (runningSum - (w[sumCounter - 1] * .5 + w[sumCounter] * .5))) / (w[sumCounter - 1] * .5 + w[sumCounter] * .5)*(y[sumCounter] - y[sumCounter - 1]);
//            std::cout << median << std::endl;
        }
    }
    else
    {
        median = y[0];
//        std::cout << median << std::endl;
    }
    return median;
}


// numerical methods/optimization
double twoPointNR(double(*y)(double, std::vector <double>), double(*dy)(double, std::vector <double>), double(*ddy)(double, std::vector <double>), std::vector <double> params, double xguess, double xguessp1)
{
    double tol = 1e-9;
    double xkm1 = xguess;
    double xk = xguessp1;

    double xkp1;
    double r;
    double ykm1;
    double yk;
    double dyk;

    int itercount = 0;

    while (true) {

        ykm1 = y(xkm1, params);
        yk = y(xk, params); //function we're finding zero of
        dyk = dy(xk, params); //derivative of above

        r = 1.0 - ((yk / ykm1) * (((yk - ykm1) / (xk - xkm1)) / dyk));

        xkp1 = (1.0 - 1.0 / r)*xkm1 + xk / r;

        xkm1 = xk;
        xk = xkp1;

        itercount += 1;


        if (itercount > 5000) {
            return NAN;
        }

        if (std::abs(xk - xkm1) <= tol || std::isnan(xk) || std::abs(yk) <= tol) {
            break;
        }
    }
    return xkp1;
}

std::vector <double> cubicSolver(double A, double B, double C, double D) {
    double a1 = B / A;
    double a2 = C / A;
    double a3 = D / A;

    double Q = (a1 * a1 - 3.0 * a2) / 9.0;
    double R = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
    double Qc = std::pow(Q, 3.0);
    //double d = Qc - std::pow(R, 2.0);

    double theta = std::acos(R / sqrt(Qc));

    double r1 = -2 * std::sqrt(Q) * std::cos(theta / 3) - a1 / 3;
    double r2 = -2 * std::sqrt(Q) * std::cos((theta + 2 * PI) / 3) - a1 / 3;
    double r3 = -2 * std::sqrt(Q) * std::cos((theta + 4 * PI) / 3) - a1 / 3;

    std::vector <double> roots = { r1, r2, r3 };

    return roots;
}


// testing purposes only
clock_t startTimer() {
    std::cout << "starting timer... \n";

    clock_t t_i = clock();
    return t_i;
}

double secElapsed(clock_t t_i) {
    clock_t t_f = clock() - t_i;
    double sec_elapsed = ((float)t_f) / CLOCKS_PER_SEC;

    printf("%0.3f seconds elapsed \n", sec_elapsed);
    return sec_elapsed;
}

void writeResults(TRK TRKobj, double t_sec, std::string filename) {

    std::ofstream myfile("/Users/nickk124/research/reichart/TRK/TRKrepo_public/diagnostics/TRKresults.txt", std::ofstream::app);
    if (myfile.is_open())
    {
        // filename    a     b     optimum scale    total computation time (s)
        myfile << filename << "\t" <<  TRKobj.results.minimumScale << "\t" << TRKobj.results.maximumScale << "\t" << TRKobj.results.optimumScale << "\t" << t_sec << "\n";
        myfile.close();
    }
    else std::cout << "Unable to open file";

}

double toRad(double deg) {
    return deg * (PI / 180);
}

double toDeg(double rad) {
    return rad * (180 / PI);
}

std::vector <std::vector <double > > getData(std::string fileName, int dataSize) {
    
    std::ifstream readFile;
    readFile.open(fileName);
    
    int columns = 5;
    int rows = dataSize;
    std::vector <double> innerData(dataSize, 0.0);
    std::vector<std::vector<double>> rawData(5, innerData);
    
    for (int row = 0; row < rows; ++row)
    {
        std::string line;
        std::getline(readFile, line);
//        if (!readFile.good())
//            break;
        
        std::stringstream iss(line);
        
        for (int col = 0; col < columns; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            //            if (!iss.good())
            //                break;
            
            std::stringstream convertor(val);
            convertor >> rawData[col][row];
        }
    }
    
    return rawData;
}

std::vector <std::vector <double > > getData(std::string fileName) {
    
    std::ifstream readFile;
    readFile.open(fileName);
    
    int columns = 5;
    int rows = 441;
    std::vector <double> innerData(441, 0.0);
    std::vector<std::vector<double>> rawData(5, innerData);
    
    for (int row = 0; row < rows; ++row)
    {
        std::string line;
        std::getline(readFile, line);
        if (!readFile.good())
            break;
        
        std::stringstream iss(line);
        
        for (int col = 0; col < columns; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            if (!iss.good())
                break;
            
            std::stringstream convertor(val);
            convertor >> rawData[col][row];
        }
    }
    
    return rawData;
}

// ###########################################################################################################################
