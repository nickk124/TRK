/*
 Trotter Reichart Konz (TRK) Official Codebase
 Author: Nick C. Konz
 See license at https://github.com/nickk124/TRK
 */

#include "TRK.h"

namespace TRKLib {
    // STATIC VARIABLES ##########################################################################################################

    std::vector <double> TRK::CorrelationRemoval::pivots = {};
    bool TRK::COVID19::logModel = false;
    int TRK::COVID19::S = 3;
    double TRK::COVID19::y12 = 2.55;
    double TRK::COVID19::s = 1.0;
    double TRK::COVID19::t_split = 45.5;
    double TRK::COVID19::tmed = 1.0;
    std::vector <double> TRK::COVID19::fixed_params = {};
    std::vector <double> TRK::COVID19::fixed_pivots = {};

    // ###########################################################################################################################



    // PRIORS ####################################################################################################################

    // Constructors
    Priors::Priors(enum priorType priorType, std::vector < std::vector <double> > params) { //Gaussian or bounded priors only
        this->priorType = priorType;
        if (priorType == GAUSSIAN) {
            this->gaussianParams = params;
        }
        else if (priorType == CONSTRAINED) {
            this->paramBounds = params;
        }
    };

    Priors::Priors(enum priorType priorType, std::vector < std::vector <double> > gaussianParams, std::vector < std::vector <double> > paramBounds) { //mixed
        this->priorType = priorType;
        this->paramBounds = paramBounds;
        this->gaussianParams = gaussianParams;
    };

    Priors::Priors(enum priorType priorType, std::vector <std::function <double(double)> > priorsPDFs) { //custom priors
        this->priorType = priorType;
        this->priorsPDFs = priorsPDFs;
    };

    Priors::Priors(enum priorType priorType, std::function <double(std::vector <double>)> jointPriorsPDF) { //custom priors
        this->priorType = priorType;
        this->jointPriorsPDF = jointPriorsPDF;
    };


    //default constructor
    Priors::Priors() {};

    // ###########################################################################################################################



    // CONSTRUCTORS ##############################################################################################################

    TRK::TRK(std::function <double(double, std::vector <double>)> yc, std::function <double(double, std::vector <double>)> dyc, std::function <double(double, std::vector <double>)> ddyc, std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess) : statistics(*this), optimization(*this), scaleOptimization(*this), tangentPointMethods(*this),  mcmc(*this), correlationRemoval(*this), asymmetric(*this), settings(*this), covid19(*this) {
        this->yc = yc;
        this->dyc = dyc;
        this->ddyc = ddyc;

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
    TRK::TRK(std::function <double(double, std::vector <double>)> yc, std::function <double(double, std::vector <double>)> dyc, std::function <double(double, std::vector <double>)> ddyc, std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess) : statistics(*this), optimization(*this), scaleOptimization(*this), tangentPointMethods(*this),  mcmc(*this), correlationRemoval(*this), asymmetric(*this), settings(*this), covid19(*this) {
        this->yc = yc;
        this->dyc = dyc;
        this->ddyc = ddyc;

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
    TRK::TRK(std::function <double(double, std::vector <double>)> yc, std::function <double(double, std::vector <double>)> dyc, std::function <double(double, std::vector <double>)> ddyc, std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, Priors priorsObject) : statistics(*this), optimization(*this), scaleOptimization(*this), tangentPointMethods(*this),  mcmc(*this), correlationRemoval(*this), asymmetric(*this), settings(*this), covid19(*this) {
        this->yc = yc;
        this->dyc = dyc;
        this->ddyc = ddyc;

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
    TRK::TRK(std::function <double(double, std::vector <double>)> yc, std::function <double(double, std::vector <double>)> dyc, std::function <double(double, std::vector <double>)> ddyc, std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, Priors priorsObject) : statistics(*this), optimization(*this), scaleOptimization(*this), tangentPointMethods(*this),  mcmc(*this), correlationRemoval(*this), asymmetric(*this), settings(*this), covid19(*this) {
        this->yc = yc;
        this->dyc = dyc;
        this->ddyc = ddyc;

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
    TRK::TRK(std::function <double(double, std::vector <double>)> yc, std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess) : statistics(*this), optimization(*this), scaleOptimization(*this), tangentPointMethods(*this),  mcmc(*this), correlationRemoval(*this), asymmetric(*this), settings(*this), covid19(*this) {
        this->settings.do1DFit = true;
        
        this->yc = yc;

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
    TRK::TRK(std::function <double(double, std::vector <double>)> yc, std::vector <double> x, std::vector <double> y, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess) : statistics(*this), optimization(*this), scaleOptimization(*this), tangentPointMethods(*this),  mcmc(*this), correlationRemoval(*this), asymmetric(*this), settings(*this), covid19(*this) {
        this->settings.do1DFit = true;
        
        this->yc = yc;

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
    TRK::TRK(std::function <double(double, std::vector <double>)> yc, std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess, Priors priorsObject) : statistics(*this), optimization(*this), scaleOptimization(*this), tangentPointMethods(*this),  mcmc(*this), correlationRemoval(*this), asymmetric(*this), settings(*this), covid19(*this) {
        this->settings.do1DFit = true;
        
        this->yc = yc;

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
    TRK::TRK(std::function <double(double, std::vector <double>)> yc, std::vector <double> x, std::vector <double> y, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess, Priors priorsObject) : statistics(*this), optimization(*this), scaleOptimization(*this), tangentPointMethods(*this),  mcmc(*this), correlationRemoval(*this), asymmetric(*this), settings(*this), covid19(*this) {
        this->settings.do1DFit = true;
        
        this->yc = yc;

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
    TRK::TRK() : statistics(*this), optimization(*this), scaleOptimization(*this), tangentPointMethods(*this),  mcmc(*this), correlationRemoval(*this), asymmetric(*this), settings(*this), covid19(*this) {

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
        mcmc.do_mcmc = true;
        
        checkVerbose();
        
        asymmetric.checkAsym();
        
        mcmc.findAIESstartingWidths(); // only if sampling will be used, either for uncertainty computation or pivot points
        
        optimization.getBetterGuess(); // once before pivot point optimization, once after
        
        correlationRemoval.getPivotGuess(); // guesses pivots if no guess provided
        
        scaleOptimization.optimizeScale();

        correlationRemoval.findPivots();

        optimization.getBetterGuess();
        
        mcmc.calculateUncertainties();
        
        showResults(true, true);
    }

    void TRK::performTRKFit(double scale) {//perform fit on some provided scale (for example, if they already know optimum scale, they can just start with this) and calculates uncertainties
        mcmc.do_mcmc = true;
        
        checkVerbose();
        
        asymmetric.checkAsym();
        
        mcmc.findAIESstartingWidths(); // only if sampling will be used, either for uncertainty computation or pivot points
        
        scaleOptimization.s = scale;
        results.optimumScale = scale;
        
        optimization.getBetterGuess();
        
        correlationRemoval.getPivotGuess(); // guesses pivots if no guess provided

        correlationRemoval.findPivots();

        results.bestFitParams.clear();

        scaleOptimization.whichExtrema = S;
        allparams_s = optimization.downhillSimplex_Fit(statistics.selectedChiSq, allparams_guess, scale, optimization.showFittingSteps);
        scaleOptimization.whichExtrema = ANY;

        for (int j = 0; j < M; j++) {
            results.bestFitParams.push_back(allparams_s[j]);
        }

        if (!settings.do1DFit){
            
            results.slop_x = allparams_s[M];
            results.slop_y = allparams_s[M + 1];
        } else {
            results.slop_y = allparams_s[M];
        }
        
        // note: true def of asym slop signs inconsistent thru code,
        // but trust this:
        if (asymmetric.hasAsymSlop){
            if (settings.do1DFit){
                results.slop_y = allparams_s[M];
                results.slop_y_minus = allparams_s[M + 1];
            } else {
                // but the above note isn't implemented yet for
                // the following 2D asym case, as the likelihood/
                // delta shift code needs to be updated analogously
                // to the 1D case.
                results.slop_x = allparams_s[M];
                results.slop_y = allparams_s[M + 1];
                results.slop_x_minus = allparams_s[M + 2];
                results.slop_y_minus = allparams_s[M + 3];
            }
        }

        optimization.getBetterGuess();
        mcmc.calculateUncertainties();
        
        showResults(false, true);
    }

    void TRK::performSimpleTRKFit() {//finds optimum scale and performs TRK fit but without finding uncertainties
        checkVerbose();
        
        asymmetric.checkAsym();
        
        mcmc.findAIESstartingWidths(); // only if sampling will be used, either for uncertainty computation or pivot points
        
        optimization.getBetterGuess();
        
        correlationRemoval.getPivotGuess(); // guesses pivots if no guess provided
        
        scaleOptimization.optimizeScale(); // (stores results in TRK.results)

        correlationRemoval.findPivots();
        
        results.bestFitParams.clear();
        
        scaleOptimization.whichExtrema = S;
        allparams_s = optimization.downhillSimplex_Fit(statistics.selectedChiSq, allparams_guess, scaleOptimization.s, optimization.showFittingSteps);
        scaleOptimization.whichExtrema = ANY;
        
        for (int j = 0; j < M; j++) {
            results.bestFitParams.push_back(allparams_s[j]);
        }
    
        if (!settings.do1DFit){
            
            results.slop_x = allparams_s[M];
            results.slop_y = allparams_s[M + 1];
        } else {
            results.slop_y = allparams_s[M];
        }
        
        // note: true def of asym slop signs inconsistent thru code,
        // but trust this:
        if (asymmetric.hasAsymSlop){
            if (settings.do1DFit){
                results.slop_y = allparams_s[M];
                results.slop_y_minus = allparams_s[M + 1];
            } else {
                // but the above note isn't implemented yet for
                // the following 2D asym case, as the likelihood/
                // delta shift code needs to be updated analogously
                // to the 1D case.
                results.slop_x = allparams_s[M];
                results.slop_y = allparams_s[M + 1];
                results.slop_x_minus = allparams_s[M + 2];
                results.slop_y_minus = allparams_s[M + 3];
            }
        }
        
        showResults(true, false);

        return;
    }

    void TRK::performSimpleTRKFit(double scale) {//given some provided scale, performs TRK fit but without finding uncertainties
        checkVerbose();
        
        asymmetric.checkAsym();
        
        mcmc.findAIESstartingWidths(); // only if sampling will be used, either for uncertainty computation or pivot points
        
        scaleOptimization.s = scale;
        
        optimization.getBetterGuess();

        correlationRemoval.getPivotGuess(); // guesses pivots if no guess provided

        correlationRemoval.findPivots();

        results.bestFitParams.clear();

        scaleOptimization.whichExtrema = S;
        allparams_s = optimization.downhillSimplex_Fit(statistics.selectedChiSq, allparams_guess, scaleOptimization.s, optimization.showFittingSteps);
        scaleOptimization.whichExtrema = ANY;
        
        for (int j = 0; j < M; j++) {
            results.bestFitParams.push_back(allparams_s[j]);
        }
        
        if (!settings.do1DFit){
            
            results.slop_x = allparams_s[M];
            results.slop_y = allparams_s[M + 1];
        } else {
            results.slop_y = allparams_s[M];
        }
        
        // note: true def of asym slop signs inconsistent thru code,
        // but trust this:
        if (asymmetric.hasAsymSlop){
            if (settings.do1DFit){
                results.slop_y = allparams_s[M];
                results.slop_y_minus = allparams_s[M + 1];
            } else {
                // but the above note isn't implemented yet for
                // the following 2D asym case, as the likelihood/
                // delta shift code needs to be updated analogously
                // to the 1D case.
                results.slop_x = allparams_s[M];
                results.slop_y = allparams_s[M + 1];
                results.slop_x_minus = allparams_s[M + 2];
                results.slop_y_minus = allparams_s[M + 3];
            }
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
        if (settings.printResults){
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

            printf("BEST FIT MODEL PARAMETERS:\n");
            for (int k = 0; k < params_guess.size(); k++) {
                printf("%.5e ", results.bestFitParams[k]);
            }
            printf("\n\nBEST FIT SLOP (EXTRINSIC SCATTER) PARAMETERS:\n");
            if (asymmetric.hasAsymSlop){
                if (settings.do1DFit){
                    printf("y-slop+ =  %.5f y-slop- =  %.5f", results.slop_y, results.slop_y_minus);
                } else {
                    printf("x-slop+ =  %.5f y-slop+ = %.5f x-slop- =  %.5f y-slop- = %.5f", results.slop_x, results.slop_y, results.slop_x_minus, results.slop_y_minus);
                }
            } else {
                if (settings.do1DFit){
                    printf("y-slop =  %.5f", results.slop_y);
                } else {
                    printf("x-slop =  %.5f y-slop = %.5f", results.slop_x, results.slop_y);
                }
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
                
                if (asymmetric.hasAsymSlop){
                    printf("\nASYMMETRIC (NEGATIVE) SLOP UNCERTAINTIES: (- 1 2 3, + 1 2 3): \n");
                    if (!settings.do1DFit){
                        for (int j = 0; j < 2; j++) {
                            for (int i = 0; i < 3; i++) {
                                printf("%.3e ", results.slopX_minus_123Sigmas[j][i]);
                            }
                            printf("\t");
                        }
                        std::cout << std::endl;
                    }
                    for (int j = 0; j < 2; j++) {
                        for (int i = 0; i < 3; i++) {
                            printf("%.3e ", results.slopY_minus_123Sigmas[j][i]);
                        }
                        printf("\t");
                    }
                    std::cout << std::endl << std::endl;
                }
            }
            
            printf("\nFITNESS:\nchisquared = %f\n\n", results.fitness);
        }
        
        
        if (covid19.print_custom_output){
           covid19.printCustomResults();
        }
        
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
    double TRK::Statistics::pearsonCorrelation(std::vector <double> x, std::vector <double> y){
        double uppersum = 0.0, lowersum1 = 0.0, lowersum2 = 0.0, x_bar = getAverage(x), y_bar = getAverage(y);
        
        for (int i = 0; i < (int) x.size(); i++){
            uppersum += (x[i] - x_bar) * (y[i] - y_bar);
            lowersum1 += std::pow((x[i] - x_bar), 2.0);
            lowersum2 += std::pow((y[i] - y_bar), 2.0);
        }
        
        return uppersum / (std::sqrt(lowersum1) * std::sqrt(lowersum2));
    }

    double TRK::Statistics::spearmanCorrelation(std::vector <double> x, std::vector <double> y){
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
            sum += trk.w[i] * std::pow(trk.y[i] - trk.yc(trk.x[i], params), 2.0);
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
            sum += trk.w[i]* std::pow(trk.yc(trk.x[i], params) - trk.y[i], 2)/(std::pow(trk.sy[i], 2.0) + std::pow(sigma, 2.0)) + 2.0*std::log(std::pow(std::pow(trk.sy[i], 2.0) + std::pow(sigma, 2.0), trk.w[i]/2.0));
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
            l = std::exp(-0.5 * trk.w[i]* std::pow((trk.yc(trk.x[i], params) - trk.y[i])/std::sqrt(std::pow(trk.sy[i], 2.0) + std::pow(sigma, 2.0)), 2.0)) / std::pow(std::pow(trk.sy[i], 2.0) + std::pow(sigma, 2.0), trk.w[i]/2.0);
            if (trk.statistics.useLogLikelihood1D){
                L += std::log(l);
            } else {
                L *= l;
            }
        }
    //    printf("%.3e\n",L);
        return L; // returns log L = logL1 + logL2 + ... given L = L1*L2*L3... if useLogPosterior == true
    }


    // likelihoods and posteriors with asymmetric uncertainties
    // 1D
    double TRK::Statistics::likelihood1DAsym(std::vector <double> allparams) {
    //    printf("Notice: 1D likelihood (for testing) not currently configured to work with weights.\n");
        
        double lnL = 0.0;
        
        if (trk.asymmetric.use_analytic_1D_asym_likelihood){
            double sigmap = allparams[trk.M];
            double sigmam = allparams[trk.M + 1];
            
            // switching
//            double sigmapold = sigmap;
//            double sigmamold = sigmam;
//            sigmam = sigmapold;
//            sigmap = sigmamold;
            
//            printf("%f\t%f\n", sigmap, sigmam);
            
            std::vector <double> params;

            for (int i = 0; i < trk.M; i++) {
                params.push_back(allparams[i]);
            }
            
            double l;
            double delta_i; // asymmetric shift
            double Sigmap, Sigmam, ymodel, Sigmapostfactor; // capital sigma +, - (convolved uncertainties)

            for (int i = 0; i < trk.N; i++) {
                delta_i = trk.asymmetric.getAsymShift1D(allparams, i);
    //            Sigmap = std::sqrt(std::pow(sigmap, 2.0) + std::pow(trk.sy[i], 2.0));
    //            Sigmam = std::sqrt(std::pow(sigmam, 2.0) + std::pow(trk.asymmetric.sy_minus[i], 2.0));
                Sigmap = std::sqrt(std::pow(sigmam, 2.0) + std::pow(trk.sy[i], 2.0));
                Sigmam = std::sqrt(std::pow(sigmap, 2.0) + std::pow(trk.asymmetric.sy_minus[i], 2.0));
                ymodel = trk.yc(trk.x[i], params);
                
                Sigmapostfactor = trk.y[i] + delta_i >= ymodel ? Sigmap : Sigmam;
                
                l = -1.0 * (trk.w[i] * (std::log(Sigmap + Sigmam) + 0.5 * std::pow((trk.y[i] + delta_i - ymodel) / Sigmapostfactor, 2.0)) - 0.5 * std::log(2.0/PI));
                
                
    //            printf("%f\n", l);
    //
    //            if (l == 0.0){
    //                printf("test\n");
    //            }
                lnL += l;
            }
            
            
//            printf("%f\t%f\n", sigmap, sigmam);
        //    printf("%.3e\n",L);
        } else {
            // all this for constant model case
            
            double syp,sym;
            double prob;
            double pi=3.1415926535898;
            double xn,yn,sigyplus,sigyminus; // error bars (last 2)
            double b; // model "curve"
            double p1,p2;
            double deltay,weight,sigp,sigm; // asym slop components (last 2)
            
            std::vector <double> params;
            
            for (int i = 0; i < trk.M; i++) {
                params.push_back(allparams[i]);
            }
            
            // slop
            syp = allparams[trk.M];
            sym = syp;
            if (trk.asymmetric.hasAsymSlop){
                sym = allparams[trk.M+1];
            }
            
            // switch signs of slops, due to convention of Capital sigmas
            double sypold = syp;
            double symold = sym;
            sym = sypold;
            syp = symold;
            
//            printf("%f\t%f\n", syp, sym);
            
            for (int i = 0; i < trk.N; i++) {
                b = trk.yc(trk.x[i], params);
                // xn,yn = x,y quoted values for data point n
                xn = trk.x[i];
                yn = trk.y[i];
                
                // errorbars:
                sigyplus = trk.sy[i];
                sigyminus = trk.asymmetric.sy_minus[i];
                weight = trk.w[i];
                
//                deltay = trk.asymmetric.getAsymShift1D(allparams, i);
//                yn = yn + deltay; // add delta shift to yn (incorrect)

                sigp = sqrt(syp*syp+sigyplus*sigyplus);
                sigm = sqrt(sym*sym+sigyminus*sigyminus);
                
                // computation of contribution to likelihood
                double normaldist_factor1, normaldist_factor2, x;
                if (sigyplus != 0)
                {
                    x = syp*(b-yn)/sigyplus/sigp;
                    normaldist_factor1 = (1.0 + std::erf(x/std::sqrt(2.0)))  / 2.0;
                    
                    p1 = 1/sqrt(2.0*pi)/sigp*exp(-(b-yn)*(b-yn)/2.0/sigp/sigp)*2.0*syp/(sym+syp)*normaldist_factor1;
                    
                    x = sym*(b-yn)/sigyplus/sigm;
                    normaldist_factor2 = (1.0 + std::erf(x/std::sqrt(2.0))) / 2.0;
                    
                    p2 = 1/sqrt(2.0*pi)/sigm*exp(-(b-yn)*(b-yn)/2.0/sigm/sigm)*2.0*sym/(sym+syp)*(1.0-normaldist_factor2);
                    prob = p1+p2; // contribution to the likelihood function
                }
                else // if symmetric intrinsic errorbar = 0
                {
                    if (yn <= b)
                    {
                        prob = 1/sqrt(2.0*pi)*exp(-(b-yn)*(b-yn)/sigp/sigp/2.0)*2.0/(sym+syp);
                    }
                    else
                    {
                        prob = 1/sqrt(2.0*pi)*exp(-(b-yn)*(b-yn)/sigm/sigm/2.0)*2.0/(sym+syp);
                    }
                }
                
    //            if (yn <= b)
    //            {
    //                z = sqrt((b-yn)*(b-yn)/sigp/sigp);
    //            }
    //            else
    //            {
    //                z = sqrt((b-yn)*(b-yn)/sigm/sigm);
    //            }
                
                lnL += weight*log(prob);
            }
        }
        return lnL;
    }

    double TRK::Statistics::regularChiSquaredWSlopAsym(std::vector <double> allparams, double s) {
        double chi2;
        
        double lnL = likelihood1DAsym(allparams);
        chi2 = -2.0 * lnL;
        
        return chi2;
    }

    // 2D
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

    double TRK::Statistics::getMode(int trueCount, std::vector <double> w, std::vector <double> y)
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
    std::vector <double> TRK::Optimization::downhillSimplex(std::function <double(std::vector <double> )> func, std::vector <double> guess, double tolerance, bool show_steps, int max_iters, bool save_evals){
        double tol = tolerance;
        
        // empty saved evals
        saved_vertices.clear();
        saved_evals.clear();

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
//                    unorderedEvals.push_back(func(vertices[i]));
                    unorderedEvals.push_back(downhillSimplex_func_SaveEvals(func, vertices[i], save_evals));
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

//                double fr = func(refpoint);
//                double f1 = func(vertices[0]);
//                double fn = func(vertices[n - 1]);
                
                double fr = downhillSimplex_func_SaveEvals(func, refpoint, save_evals);
                double f1 = downhillSimplex_func_SaveEvals(func, vertices[0], save_evals);
                double fn = downhillSimplex_func_SaveEvals(func, vertices[n - 1], save_evals);

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

//                    double fe = func(exppoint);
                    double fe = downhillSimplex_func_SaveEvals(func, exppoint, save_evals);


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
//                    double fnp1 = func(vertices[n]);
                    double fnp1 = downhillSimplex_func_SaveEvals(func, vertices[n], save_evals);

                    if (fn <= fr && fr < fnp1) {
                        std::vector <double> cpoint;

                        for (int i = 0; i < n; i++) {
                            cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
                        }

//                        double fc = func(cpoint);
                        double fc = downhillSimplex_func_SaveEvals(func, cpoint, save_evals);

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

//                        double fcc = func(ccpoint);
                        double fcc = downhillSimplex_func_SaveEvals(func, ccpoint, save_evals);

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
                
//                double eval = func(result);
                double eval = downhillSimplex_func_SaveEvals(func, result, save_evals);
                std::cout << "\tfunc = " << eval << "\n";
                
            }
            

            
            
            //test for termination

            std::vector <double> evals;
            for (int i = 0; i < n + 1; i++) {
                evals.push_back(downhillSimplex_func_SaveEvals(func, vertices[i], save_evals));
            }
            
            if (show_steps){
                printf("stdev(fitnesses) = %f\n", trk.statistics.stDevUnweighted(evals));
            }

            if (n > 1 && trk.statistics.stDevUnweighted(evals) < tol) { // n-dimensional
                break;
            }
            
            if (n == 1 && std::abs(vertices[1][0] - vertices[0][0]) < tol) { // 1 dimensional
                break;
            }
            
            if (it >= max_iters){ // iteration check
                if (trk.settings.verbose){
                    printf("Downhill simplex exceeded %i iterations; halting...\n", max_iters);
                }
                break;
            }
            
            if (it % 100 == 0 && show_steps){printf("simplex iteration: %i\n", it);};

            it++;
        }

        optimum = vertices[n];

        return optimum;
        
    }

    double TRK::Optimization::downhillSimplex_func_SaveEvals(std::function <double(std::vector <double> )> func, std::vector <double> vertex, bool save_evals){
        double eval = NAN;
        
        if (save_evals){
            // check to see if stored eval
            for (int i = 0; i < saved_evals.size(); i++){
                bool isEqual = true;
                for (int j = 0; j < vertex.size(); j++){
                    if (vertex[j] != saved_vertices[i][j]){
                        isEqual = false;
                    }
                }
                if (isEqual){
                    eval = saved_evals[i];
                    break;
                }
            }
            
            if (std::isnan(eval)){ // not stored; find eval and store it
                eval = func(vertex);
                
                saved_evals.push_back(eval);
                saved_vertices.push_back(vertex);
            }
            
        } else {
            eval = func(vertex);
        }
        
        return eval;
    }

    double TRK::Optimization::downhillSimplex_1DWrapper(std::function <double(std::vector <double> )> func, std::vector <double> guess, double tolerance, bool show_steps, int max_iters){
        
        return downhillSimplex(func, guess, tolerance, show_steps, max_iters, trk.correlationRemoval.save_simplex_correlation_evals)[0];
    }

    // downhill simplex/ nelder mead method customized for fitting
    std::vector <double> TRK::Optimization::downhillSimplex_Fit(double(TRK::Statistics::*f)(std::vector <double> , double), std::vector <double> allparams_guess, double s, bool show_steps) {
        
//        if (s == 0.0625){
//            printf("test\n");
//        }
        
        double prev_fitness = DBL_MAX;

        unsigned long n = trk.bigM; //number of model parameters plus two symmetric slop parameters
        
        if (trk.asymmetric.hasAsymSlop){
            if (trk.settings.do1DFit){
                n++;
            } else {
                n += 2;
            }
        }

        double rho = 1.0; //reflection
        double chi = 2.0; //expansion
        double gamma = 0.5; //contraction
        double sigma = 0.5; //shrinkage
        
        double tol = simplexTol;// * std::sqrt(n);
        
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
                std::vector <std::vector <double> > nvertices(n, std::vector<double>());
                for (int i = 0; i < n; i++) {
                    nvertices[i] = vertices[i];
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
            
            double fitness = evalWPriors(f, result, trk.scaleOptimization.s);
            
            if (show_steps){
                std::cout << "chi-square parameters at s = " << s << " ";
                for (int i = 0; i < result.size(); i++) {
                    std::cout << result[i] << " ";
                }
                
                std::cout << "fitness = " << fitness << "\n";
                
            }
            
            if (trk.scaleOptimization.whichExtrema == S or trk.scaleOptimization.whichExtrema == ANY){
                double fitness_res = evalWPriors(f, result, trk.scaleOptimization.s);
                trk.results.fitness = fitness_res;
    //            printf("\n\nfinal fitness = %.3e\n\n", fitness);
            }
            
            
            //test for termination

            std::vector <double> evals;
            for (int i = 0; i < n + 1; i++) {
                evals.push_back(evalWPriors(f, vertices[i], s));
            }

            if (simplex_terminate_stdev){
                if (trk.statistics.stDevUnweighted(evals) < tol) {
                    break;
                }
            } else {
                if (std::abs(fitness - prev_fitness) < tol){
                    break;
                }
            }
            
            prev_fitness = fitness;
            
            if (it >= max_simplex_iters){
                if (trk.settings.verbose){
                    printf("Downhill simplex exceeded %i iterations; halting...\n", max_simplex_iters);
                }
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

    std::vector <double> TRK::Optimization::getFullallparams(std::vector <double> allparams_free, std::vector <double> fixed_param_vals, std::vector <bool> fixed_allparams_flags){
        std::vector <double> allparams_full;
            // takes vector of free params and adds values from fixed params
            
        unsigned int count_free = 0, count_fixed = 0;
        for (int j = 0; j < fixed_allparams_flags.size(); j++){
                if (!fixed_allparams_flags[j]){  // free
                    allparams_full.push_back(allparams_free[count_free]);
                    count_free++;
                } else { // fixed
                    allparams_full.push_back(fixed_param_vals[count_fixed]);
                    count_fixed++;
                }
            }
            return allparams_full;
        }

    std::vector <double> TRK::Optimization::downhillSimplex_Fit_Fixed(double(TRK::Statistics::*f)(std::vector <double>, double), std::vector <double> allparams_free_guess, double s, bool show_steps, std::vector <bool> fixed_allparam_flags, std::vector <double> fixed_allparam_vals){
        double prev_fitness = DBL_MAX;

            unsigned long n = trk.bigM; //number of model parameters plus two slop parameters
            
            if (trk.asymmetric.hasAsymSlop){
                if (trk.settings.do1DFit){
                    n++;
                } else {
                    n += 2;
                }
            }
        
        
            // parse through fixed params
            n -= fixed_allparam_vals.size();

            double rho = 1.0; //reflection
            double chi = 2.0; //expansion
            double gamma = 0.5; //contraction
            double sigma = 0.5; //shrinkage
            
            double tol = simplexTol;// * std::sqrt(n);
            
            // simplex initialization

            std::vector <double> init_point = allparams_free_guess;

            std::vector <std::vector <double> > vertices(n + 1, init_point);
            std::vector <double> fitted_params, full_vertex;

            int i = 0;
            for (int j = 1; j < n + 1; j++) { //for each simplex node

                vertices[j][i] = init_point[i] != 0 ? init_point[i] + simplex_size*init_point[i] : 0.1; //add initial "step size"
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
                        full_vertex = getFullallparams(vertices[i], fixed_allparam_vals, fixed_allparam_flags);
                        unorderedEvals.push_back(evalWPriors(f, full_vertex, s));
                    }
                    orderedindices = getSortedIndices(unorderedEvals);

                    std::vector <std::vector <double> > orderedvertices = { vertices[orderedindices[0]] };
                    for (int i = 1; i < n + 1; i++) {
                        orderedvertices.push_back(vertices[orderedindices[i]]);
                    }

                    vertices = orderedvertices;

                    // reflect
                    std::vector <double> refpoint;
                    std::vector <std::vector <double> > nvertices(n, std::vector<double>());
                    for (int i = 0; i < n; i++) {
                        nvertices[i] = vertices[i];
                    }

                    std::vector <double> centroid = findCentroid(nvertices);

                    for (int i = 0; i < n; i++) {
                        refpoint.push_back(centroid[i] + rho * (centroid[i] - vertices[n][i]));
                    }

                    full_vertex = getFullallparams(refpoint, fixed_allparam_vals, fixed_allparam_flags);
                    double fr = evalWPriors(f, full_vertex, s);
                    
                    full_vertex = getFullallparams(vertices[0], fixed_allparam_vals, fixed_allparam_flags);
                    double f1 = evalWPriors(f, full_vertex, s);
                    
                    full_vertex = getFullallparams(vertices[n - 1], fixed_allparam_vals, fixed_allparam_flags);
                    double fn = evalWPriors(f, full_vertex, s);

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

                        full_vertex = getFullallparams(exppoint, fixed_allparam_vals, fixed_allparam_flags);
                        double fe = evalWPriors(f, full_vertex, s);


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
                        full_vertex = getFullallparams(vertices[n], fixed_allparam_vals, fixed_allparam_flags);
                        double fnp1 = evalWPriors(f, full_vertex, s);

                        if (fn <= fr && fr < fnp1) {
                            std::vector <double> cpoint;

                            for (int i = 0; i < n; i++) {
                                cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
                            }
                            
                            full_vertex = getFullallparams(cpoint, fixed_allparam_vals, fixed_allparam_flags);
                            double fc = evalWPriors(f, full_vertex, trk.scaleOptimization.s);

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

                            full_vertex = getFullallparams(ccpoint, fixed_allparam_vals, fixed_allparam_flags);
                            double fcc = evalWPriors(f, full_vertex, trk.scaleOptimization.s);

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
                
                full_vertex = getFullallparams(result, fixed_allparam_vals, fixed_allparam_flags);
                double fitness = evalWPriors(f, full_vertex, trk.scaleOptimization.s);
                
                if (show_steps){
                    std::cout << "chi-square parameters at s = " << s << " ";
                    for (int i = 0; i < result.size(); i++) {
                        std::cout << result[i] << " ";
                    }
                    
                    std::cout << "fitness = " << fitness << "\n";
                    
                }
                
                if (trk.scaleOptimization.whichExtrema == S or trk.scaleOptimization.whichExtrema == ANY){
                    full_vertex = getFullallparams(result, fixed_allparam_vals, fixed_allparam_flags);
                    double fitness_res = evalWPriors(f, full_vertex, trk.scaleOptimization.s);
                    trk.results.fitness = fitness_res;
        //            printf("\n\nfinal fitness = %.3e\n\n", fitness);
                }
                
                
                //test for termination

                std::vector <double> evals;
                for (int i = 0; i < n + 1; i++) {
                    full_vertex = getFullallparams(vertices[i], fixed_allparam_vals, fixed_allparam_flags);
                    evals.push_back(evalWPriors(f, full_vertex, s));
                }

                if (simplex_terminate_stdev){
                    if (trk.statistics.stDevUnweighted(evals) < tol) {
                        break;
                    }
                } else {
                    if (std::abs(fitness - prev_fitness) < tol){
                        break;
                    }
                }
                
                prev_fitness = fitness;
                
                if (it >= max_simplex_iters){
                    if (trk.settings.verbose){
                        printf("Downhill simplex exceeded %i iterations; halting...\n", max_simplex_iters);
                    }
                    break;
                }
                
                if (it % 100 == 0 && show_steps){printf("simplex iteration: %i\n", it);};
        //
                it++;
            }

            fitted_params = getFullallparams(vertices[n], fixed_allparam_vals, fixed_allparam_flags);

            fitted_params = pegToZeroSlop(fitted_params);
            fitted_params = avoidNegativeSlop(fitted_params, n);

            return fitted_params;
    }
    
    // downhill simplex tools
    std::vector <double> TRK::Optimization::findCentroid(std::vector < std::vector <double> > nvertices) {
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
        if (!trk.settings.do1DFit){ // no scale optimization on 1D fits
            if (std::abs(vertex[trk.M]) <= pegToZeroTol){//} && std::abs(vertex[trk.M+1]) > pegToZeroTol) {
                vertex[trk.M] = 0;
            }
            if (std::abs(vertex[trk.M+1]) <= pegToZeroTol){//} && std::abs(vertex[trk.M]) > pegToZeroTol) {
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
//        
//        if (vertex[vertex.size() - 1] < 0 || vertex[vertex.size() - 2] < 0){
//            printf("bandaid fix for negative asym slop issue active\n");
//            return DBL_MAX;
//        }
        
        double ret = (trk.statistics.*f)(vertex, s);
    //
        if (std::isnan(ret)){
            return DBL_MAX;
        }
        
        return ret;
    }

    void TRK::Optimization::getBetterGuess(){
        // re-fitting is done at beginning of scale optimization for 2D case
        
        if (trk.settings.do1DFit){ // re-fitting is done at beginning of scale optimization for 2D case
            trk.results.bestFitParams.clear();

            trk.scaleOptimization.whichExtrema = S;
            trk.allparams_s = trk.optimization.downhillSimplex_Fit(trk.statistics.selectedChiSq, trk.allparams_guess, trk.scaleOptimization.s, showFittingSteps);
            trk.scaleOptimization.whichExtrema = ANY;

            for (int j = 0; j < trk.M; j++) {
                trk.results.bestFitParams.push_back(trk.allparams_s[j]);
            }
            trk.results.slop_y = trk.allparams_s[trk.M];
            
            if (trk.asymmetric.hasAsymSlop){
                trk.results.slop_y_minus = trk.allparams_s[trk.M+1];
            }
        }
        
        if (trk.results.bestFitParams.size() != 0){
            for (int j = 0; j < trk.M; j++){
                trk.allparams_guess[j] = trk.results.bestFitParams[j];
            }
        }
        
        if (trk.settings.do1DFit){
            trk.allparams_guess[trk.M] = trk.results.slop_y;
            
            if (trk.asymmetric.hasAsymSlop){
                trk.allparams_guess[trk.M+1] = trk.results.slop_y_minus;
            }
            
            if (trk.settings.verbose){
                printf("Better parameter guess found to be:\n");
                printVector(trk.allparams_guess);
            }
        }

        //        else {  // re-fitting is done at beginning of scale optimization for 2D case
//            trk.allparams_guess[trk.M] = trk.results.slop_x;
//            trk.allparams_guess[trk.M+1] = trk.results.slop_y;
//
//            if (trk.asymmetric.hasAsymSlop){
//                trk.allparams_guess[trk.M+2] = trk.results.slop_x_minus;
//                trk.allparams_guess[trk.M+3] = trk.results.slop_y_minus;
//            }
//        }
        
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
            if (trk.settings.verbose){
                printf("1D fit: no need for scale optimization.\n");
            }
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
        
        iterative_allparams_guess_slopx = trk.allparams_guess;
        
        if (fix_a){
            whichExtremaX = SLOP_X;
            std::vector <double> a_vec = {a};
            innerSlopX_Simplex(a_vec, iterative_allparams_guess_slopx); // to get x_t_slopx
            whichExtremaX = ANY;
            return a;
        }

        // before doing any standard simplex movement, here it checks whether the simplex is within the zero "plateau", and if so, it moves it to the boundary.
        // for slop x: move to right until it hits the boundary


        //bracket finding

        double a = 0.0;
        double b = 1.0;
        double trial_s = 1.0;
        std::vector <double> trial_s_vec = {trial_s};
        double slop_trial_s = innerSlopX_Simplex(trial_s_vec, iterative_allparams_guess_slopx);

        double inc = trial_s * 0.5;

        if (slop_trial_s > 0) {
            b = trial_s;

            double trial_a = trial_s;

            while (true) {
                trial_a -= inc;
                std::vector <double> trial_a_vec = {trial_a};
                double slop_trial_a = innerSlopX_Simplex(trial_a_vec, iterative_allparams_guess_slopx);

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
                std::vector <double> trial_b_vec = {trial_b};
                double slop_trial_b = innerSlopX_Simplex(trial_b_vec, iterative_allparams_guess_slopx);

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
            std::vector <double> c_vec = {c};
            slop_c = innerSlopX_Simplex(c_vec, iterative_allparams_guess_slopx);
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
        
        iterative_allparams_guess_slopy = trk.allparams_guess;
        
        if (fix_b){
            whichExtremaY = SLOP_Y;
            std::vector <double> b_vec = {b};
            innerSlopX_Simplex(b_vec, iterative_allparams_guess_slopy); // to get x_t_slopy
            whichExtremaY = ANY;
            return b;
        }

        // before doing any standard simplex movement, here it checks whether the simplex is within the zero "plateau", and if so, it moves it to the boundary.
        // for slop x: move to right until it hits the boundary

            //bracket finding

        double a = 0.0;
        double b = 1.0;
        double trial_s = slopYScaleGuess;
        std::vector <double> trial_s_vec = {trial_s};
        double slop_trial_s = innerSlopY_Simplex(trial_s_vec, iterative_allparams_guess_slopy);

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
                std::vector <double> trial_b_vec = {trial_b};
                double slop_trial_b = innerSlopY_Simplex(trial_b_vec, iterative_allparams_guess_slopy);

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
                std::vector <double> trial_a_vec = {trial_a};
                double slop_trial_a = innerSlopY_Simplex(trial_a_vec, iterative_allparams_guess_slopy);

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
            std::vector <double> c_vec = {c};
            slop_c = innerSlopY_Simplex(c_vec, iterative_allparams_guess_slopy);
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
                printf("s=%f \t slop_x=%.3e \t slop_y=%.3e \t(slop x optimization) model params: ", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1]);
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
                printf("s=%f \t slop_x=%.3e \t slop_y=%.3e \t(slop y optimization) model params: ", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1]);
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
        std::vector <double> a_vec = {a};
        f_left = innerR2_Simplex(a_vec, iterative_allparams_guess_slopx);

        double c, f_c;
        double tol_bisect = 1e-4;
        double tol_brackets = 1e-3;

        while (true) {
            c = (left + right) / 2;

            whichExtrema = S;
            std::vector <double> c_vec = {c};
            f_c = innerR2_Simplex(c_vec, iterative_allparams_guess);
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
        std::vector <double> a_vec = {a};
        f_left = innerR2_iter_Simplex(a_vec, iterative_allparams_guess_slopx, s0);

        //bisection, now that we have brackets [left,right]

        left = a;
        right = b;

        double f_c, c = -1.0;
        double tol_bisect = 1e-4;
        double tol_brackets = 1e-3;
        
        bool tolcheck = false;
        
        int max_iter = 10;
        int iter = 0;

        while (!tolcheck) {
            printf("brackets: %.3e %.3e \n", left, right);
            c = (left + right) / 2;

            whichExtrema = S;
            std::vector <double> c_vec = {c};
            f_c = innerR2_iter_Simplex(c_vec, iterative_allparams_guess, s0);
            whichExtrema = ANY;

            if (std::abs(f_c) <= tol_bisect) { //convergence criterion
                tolcheck = true;
            }

            if (f_c * f_left > 0) {
                left = c;
                f_left = f_c;
            }
            else if (f_c * f_left < 0) {
                right = c;
            }
            
            iter++;

            if (std::abs(left - right) <= tol_brackets || iter >= max_iter) { //secondary convergence criterion (bracket width)
                tolcheck = true;
            }
        }

        return c;

    }

    double TRK::ScaleOptimization::iterateR2_OptimumScale(double s0) {
        double tol_scale = 1e-3;

        double s1 = 0.0;

        bool tolcheck = false;
        
        int max_iter = 10;
        int iter = 0;

        while (!tolcheck) {
            printf("next s0: %.3e \n", s0);

            s1 = optimize_s_prime_R2(s0);
            if (std::abs(s1-s0) <= tol_scale) {
                tolcheck = true;
            }
            std::printf("new s0: %.3e \n", s1);
            s0 = s1;
            
            iter++;
            
            if (iter >= max_iter){
                tolcheck = true;
            }
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

        double R2as = R2TRK_prime_as0(s0, trk.x_t_s, trk.params_s);
        double R2sb = R2TRK_prime_s0b(s0, trk.x_t_s, trk.params_s);
        
        if (verbose){
            if (trk.asymmetric.hasAsymSlop){
                printf("%.3e \t %.3e \t %.3e \t %.3e \t %.3e \t(additional R2 optimization)\tR2as = %f\tR2sb = %f\n", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1], trk.allparams_s[trk.M + 2], trk.allparams_s[trk.M + 3], R2as, R2sb);
            } else {
                printf("%.3e \t %.3e \t %.3e \t(additional R2 optimization)\tR2as = %f\tR2sb = %f\n", ss[0], trk.allparams_s[trk.M], trk.allparams_s[trk.M + 1], R2as, R2sb);
            }
        }

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

        double summand, sum = 0.0;
        int nan_count = 0;

        for (int n = 0; n < trk.N; n++) {
            double m_tn_a = trk.dyc(trk.x_t_a[n], trk.params_a);
            double theta_t_a = std::atan(m_tn_a);

            double m_tn_s1 = trk.dyc(x_t_s1[n], params_s1);
            double theta_t_s1 = std::atan(m_tn_s1);
            
            summand = std::pow(std::tan(PI / 4.0 - std::abs(std::atan(s0*std::tan(theta_t_a)) - std::atan(s0*std::tan(theta_t_s1))) / 2.0), 2.0);

            if (std::isnan(summand)){
                nan_count++;
            } else {
                sum += summand;
            }
        }

        R2 *= sum;
        
        if (std::isnan(R2) && verbose){
            printf("Error: NaN R2! Summand NaN count: %i\n", nan_count);
        }

        return R2;
    }

    double TRK::ScaleOptimization::R2TRK_prime_s0b(double s0, std::vector <double> x_t_s1, std::vector <double> params_s1) {
        double R2 = 1.0 / trk.N;

        double summand, sum = 0.0;
        int nan_count = 0;

        for (int n = 0; n < trk.N; n++) {
            double m_tn_b = trk.dyc(trk.x_t_b[n], trk.params_b);
            double theta_t_b = std::atan(m_tn_b);

            double m_tn_s1 = trk.dyc(x_t_s1[n], params_s1);
            double theta_t_s1 = std::atan(m_tn_s1);

            summand = std::pow(std::tan(PI / 4.0 - std::abs(std::atan(s0*std::tan(theta_t_s1)) - std::atan(s0*std::tan(theta_t_b))) / 2.0), 2.0);
            
            if (std::isnan(summand)){
                nan_count++;
            } else {
                sum += summand;
            }
        }

        R2 *= sum;
        
        if (std::isnan(R2) && verbose){
            printf("Error: NaN R2! Summand NaN count: %i\n", nan_count);
        }

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
            
            if (let_use_bisection){
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
    std::vector <double> TRK::TangentPointMethods::tangentParallel(std::vector <double> params, double slop_x, double slop_y, int n, double s) {
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

    std::vector <double> TRK::TangentPointMethods::tangentParallelLikelihood(std::vector <double> params, double slop_x, double slop_y, int n) {
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
    std::vector <std::vector <double> > TRK::MCMC::combineLinearAndNonLinearSamples(std::vector < std::vector < std::vector <double> > > &all_linearparam_samples, std::vector < std::vector <double> > &nonlinear_allparam_samples, std::vector <std::vector <bool> > &fixed_allparams_flags_linear, std::vector <bool> &fixed_allparams_flags_nonlinear){
        
        std::vector <std::vector <double> > allparam_samples(R, std::vector <double>(trk.bigM, 0.0));
        
        // pth vector of fixed_allparams_flags_linear is a vector of bools that are false for the corresponding linear params
        
        // fill up linear param samples
        for (int p = 0; p < trk.correlationRemoval.P; p++){ // indexes pivot point/ pairs of linear params
            std::vector <int> linear_param_indices;
            get_where_boolean(fixed_allparams_flags_linear[p], false, linear_param_indices); // see which indices the linear params have in the total model param vector
            
            for (int i = 0; i < R; i++){ // for each individual sample
                for (int l = 0; l < 2; l++){
                    allparam_samples[i][linear_param_indices[l]] = all_linearparam_samples[p][i][l];
                }
            }
        }
        
        // fill up remaining non-linear samples (elements of fixed_allparams_flags_nonlinear are false if non-linear)
        for (int i = 0; i < R; i++){ // for each individual sample
            std::vector <int> nonlinear_param_indices;
            get_where_boolean(fixed_allparams_flags_nonlinear, false, nonlinear_param_indices); // see which indices the nonlinear params have in the total model param vector
            
            for (int nl = 0; nl < (int) nonlinear_param_indices.size(); nl++){
                allparam_samples[i][nonlinear_param_indices[nl]] = nonlinear_allparam_samples[i][nl];
            }
        }
            
        return allparam_samples;
    }

    std::vector < std::vector <double> > TRK::MCMC::sampleForUncertainties()
    {
        std::vector <std::vector <double> > allparam_samples;
        
        if (trk.mcmc.verbose){
            std::cout << "\nSampling Posterior...\n";
        }
        
        if (trk.correlationRemoval.sampleLinearParams_seperately && trk.correlationRemoval.P > 0) { // sample linear params seperately
//            std::vector <std::vector <double> > allparam_samples(R, std::vector <double>(trk.bigM, 0.0));
            
            bool old_sampleOnlyLinearParams_pivots_setting = trk.correlationRemoval.sampleOnlyLinearParams_pivots;
            trk.correlationRemoval.sampleOnlyLinearParams_pivots = true; // needs to be on for getFixedLinearParams() to work as expected
            
            std::vector < std::vector < std::vector <double> > > all_linearparam_samples;
            
            std::vector <bool> fixed_allparams_flags_nonlinear(trk.bigM, false); // for sampling with only the linear params fixed
            
            std::vector < std::vector <bool> > fixed_allparams_flags_linear;
            
            for (int p = 0; p < trk.correlationRemoval.P; p++){ // sample each set of linear parameters separately
                std::vector <bool> fixed_allparams_flags = trk.correlationRemoval.getFixedLinearParams(p); // sample with only the relavant linear params free
                fixed_allparams_flags_linear.push_back(fixed_allparams_flags); // true if fixed, false if free
                
                std::vector < std::vector <double> > linearparam_samples = trk.mcmc.samplePosterior(R, burncount, allparams_sigmas_guess, fixed_allparams_flags);
                
                all_linearparam_samples.push_back(linearparam_samples);
                
                fixed_allparams_flags_nonlinear[trk.correlationRemoval.intercept_indices[p]] = true; // initializing for sampling only NON-linear params subsequently
                fixed_allparams_flags_nonlinear[trk.correlationRemoval.slope_indices[p]] = true;
            }
            
            trk.correlationRemoval.sampleOnlyLinearParams_pivots = old_sampleOnlyLinearParams_pivots_setting;
            
            // then sample the remaining parameters
            std::vector < std::vector <double> > nonlinear_allparam_samples = samplePosterior(R, burncount, allparams_sigmas_guess, fixed_allparams_flags_nonlinear);
            
            
            // finally, combine them all
            allparam_samples = combineLinearAndNonLinearSamples(all_linearparam_samples, nonlinear_allparam_samples, fixed_allparams_flags_linear, fixed_allparams_flags_nonlinear);
        } else {
            int m = 0;
            if (trk.asymmetric.hasAsymSlop){
                if (trk.settings.do1DFit){
                    m = 1;
                } else {
                    m = 2;
                }
            }
            
            std::vector <bool> fixed_allparams_flags_default(trk.bigM + m, false); // none fixed
            
            allparam_samples = samplePosterior(R, burncount, allparams_sigmas_guess, fixed_allparams_flags_default);
        }
        

        return allparam_samples;
    }

    void TRK::MCMC::calculateUncertainties() {
        
    //    useLogPosterior = false;
        // only needed the above true for pivot point sampling
    //    goodDeltasFound = false;
        // recompute step sizes cause they may be different given new pivot point

        std::vector <std::vector <std::vector <double> > > allparam_uncertainties;
        
        std::vector <std::vector <double> > allparam_samples = sampleForUncertainties();

        if (trk.settings.outputDistributionToFile) {

            std::string fileName = trk.settings.outputPath + std::string("/TRKMCMC_") + std::to_string(trk.allparams_guess[0]) + std::string("_") + std::to_string(R) + std::string(".txt");

            std::ofstream myfile;
            myfile.open(fileName, std::ofstream::trunc);
            
            int m = 0;
            if (trk.asymmetric.hasAsymSlop){
                if (trk.settings.do1DFit){
                    m = 1;
                } else {
                    m = 2;
                }
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

        allparam_uncertainties = lowerBar(allparam_samples);
        //for each parameter including slop, there is a vector containing 1 vector of -sigmas,
        // 1 vector of +sigmas. This vector contains all of those 2-vectors.

        trk.results.bestFit_123Sigmas.clear();
        trk.results.bestFitParamUncertainties.clear();
        
        std::vector <double> simple_uncertainties;

        for (int j = 0; j < trk.M; j++) {
            trk.results.bestFit_123Sigmas.push_back(allparam_uncertainties[j]);

            simple_uncertainties = {allparam_uncertainties[j][0][0], allparam_uncertainties[j][1][0]};
            trk.results.bestFitParamUncertainties.push_back(simple_uncertainties);
        }
        if (trk.settings.do1DFit){
            trk.results.slopY_123Sigmas = allparam_uncertainties[trk.M];
            trk.results.slopYUncertainty = {allparam_uncertainties[trk.M][0][0], allparam_uncertainties[trk.M][1][0]};
            
            if (trk.asymmetric.hasAsymSlop){
                trk.results.slopY_minus_123Sigmas = allparam_uncertainties[trk.M + 1];
            }
        } else {
            trk.results.slopX_123Sigmas = allparam_uncertainties[trk.M];
            trk.results.slopY_123Sigmas = allparam_uncertainties[trk.M + 1];
            
            trk.results.slopXUncertainty = {allparam_uncertainties[trk.M][0][0], allparam_uncertainties[trk.M][1][0]};
            trk.results.slopYUncertainty = {allparam_uncertainties[trk.M + 1][0][0], allparam_uncertainties[trk.M + 1][1][0]};
            
            if (trk.asymmetric.hasAsymSlop){
                trk.results.slopX_minus_123Sigmas = allparam_uncertainties[trk.M + 2];
                trk.results.slopY_minus_123Sigmas = allparam_uncertainties[trk.M + 3];
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
            if (trk.settings.do1DFit){
                m = 1;
            } else {
                m = 2;
            }
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
                if (K > 0){
                    leftBound = (edges[indicesIn[0]] + edges[indicesIn[0] + 1]) / 2.0;
                    rightBound = (edges[indicesIn[K - 1]] + edges[indicesIn[K - 1] + 1]) / 2.0;
                } else {
                    // if the histogram computation messed (due to poorly-defined sampling)
                    leftBound = -DBL_MAX;
                    rightBound = DBL_MAX;
                }

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

    std::vector <double> TRK::MCMC::getFullallparams(std::vector <double> X, std::vector <bool> fixed_allparams_flags){
        // takes vector of params from MCMC sampler and adds values from fixed params
        
        if (X.size() < trk.bigM){
            std::vector <double> X_full = trk.allparams_guess; // default MCMC starting point
            unsigned int count = 0;
            for (int j = 0; j < trk.bigM; j++){
                if (!fixed_allparams_flags[j]){  // true if it is fixed
                    X_full[j] = X[count];
                    count++;
                }
            }
//            X = X_full;
            return X_full;
            
        } else { // all params free
            return X;
        }
    }


    // sampling (general)
    double TRK::MCMC::metHastRatio(std::vector <double> X_trial, std::vector <double> X_i, std::vector <bool> fixed_allparams_flags){
        double log_a;
        
        X_trial = getFullallparams(X_trial, fixed_allparams_flags);
        X_i = getFullallparams(X_i, fixed_allparams_flags);
        
        double logL_trial = (trk.statistics.*trk.statistics.selectedLikelihood)(X_trial);
        double logL_i = (trk.statistics.*trk.statistics.selectedLikelihood)(X_i);

        if (trk.statistics.hasPriors) {
            double logp_trial = std::log(trk.statistics.priors(X_trial));
            double logp_i = std::log(trk.statistics.priors(X_i));
            
            log_a = logL_trial - logL_i + logp_trial - logp_i;
                // these likelihoods return log likelihood given useLogPosterior = true; the computation is done WITHIN the function.
        }
        else {
            log_a = logL_trial - logL_i;
                // this returns the log likelihood given useLogPosterior = true; the computation is done WITHIN the function.
        }
//        std::cout << std::pow(10.0, log_a) << std::endl;
        
        return log_a; // returns log post / log post if useLogPosterior == true
    }

    std::vector <double> TRK::MCMC::getMCMCStartingPoint(std::vector <bool> fixed_allparams_flags){
        std::vector <double> starting_point;
        
        for (int j = 0; j < (int) fixed_allparams_flags.size(); j++){ // if a param is fixed, the correpond element is true
            if (!fixed_allparams_flags[j]) {starting_point.push_back(trk.allparams_guess[j]); };
        }
        
        return starting_point;
    }

    void TRK::MCMC::findAIESstartingWidths(){
        if (trk.correlationRemoval.findPivotPoints || do_mcmc){
            AIES_param_width_estimates.clear();
            
            // default values
            
            int m = 0;
            if (trk.asymmetric.hasAsymSlop){
                if (trk.settings.do1DFit){
                    m = 1;
                } else {
                    m = 2;
                }
            }
        
            for (int j = 0; j < trk.bigM + m; j++){
                AIES_param_width_estimates.push_back(std::abs(trk.allparams_guess[j] != 0 ? trk.allparams_guess[j] * AIES_initial_scaling : 0.1));
            }
            
            if (!initializeAIESWalkersNaively){
        
                std::vector <bool> fixed_allparams_flags_default(trk.bigM, false); // none fixed
                
                std::vector < std::vector <double> > allparam_samples = trk.mcmc.samplePosterior(starting_width_estimate_samplesize, burncount, trk.mcmc.allparams_sigmas_guess, fixed_allparams_flags_default);
                
                for (int j = 0; j < trk.bigM + m; j++){
                    std::vector <double> param_sample;
                    for (int i = 0; i < allparam_samples.size(); i++) {
                        param_sample.push_back(allparam_samples[i][j]);
                    }
                    
                    AIES_param_width_estimates[j] = (trk.statistics.stDevUnweighted(param_sample));
                }
                
                if (verbose){
                    printf("Estimated param uncertainties for AIES:\n");
                    printVector(AIES_param_width_estimates);
                }
            }
        }
        
        return;
    };

    void TRK::MCMC::initializeAIESWalkers(std::vector <std::vector <double> > &all_walkers, std::vector <double> &starting_point, std::vector <bool> &fixed_allparams_flags, int L, int n){
              
        for (int k = 0; k < L; k++){ // k: walkers; j: all params
            int i = 0; // i: free params
            for (int j = 0; j < (int) fixed_allparams_flags.size(); j++){ // if a param is fixed, the correpond element is true
                if (!fixed_allparams_flags[j]) {
//                    double direction = rbool() ? 1.0 : -1.0;
//                    all_walkers[k][i] = starting_point[i] + direction * AIES_param_width_estimates[j];
//                    i++;
                    
                    all_walkers[k][i] = rnorm(starting_point[i], 1e-3 * AIES_param_width_estimates[j]);
                    i++;
                }
            }
        }
        
        
        return;
    }

    std::vector <std::vector <double >> TRK::MCMC::samplePosterior(int R, int burncount, std::vector <double> sigmas_guess, std::vector <bool> fixed_allparams_flags) {
        
        useLogPosterior = true;

        // n is number of free params
        unsigned long n = std::count(fixed_allparams_flags.begin(), fixed_allparams_flags.end(), false); // trk.bigM + asymmetries;
        
        std::vector < std::vector <double > > result, result_final;
        
        std::vector <double> starting_point = getMCMCStartingPoint(fixed_allparams_flags);
        
        switch (thisSamplingMethod) {
            case AIES: {
                // INDEXING:
                // iterations: t
                // walkers: j, k
                // coordinates/parameters: i
                
                int L = amt_walkers * (int) n;   // number of walkers
                
                // initialize walkers
                std::vector <std::vector <double> > all_walkers(L, std::vector <double> (n, 0.0));
                std::vector <std::vector <double> > YY;
                
                initializeAIESWalkers(all_walkers, starting_point, fixed_allparams_flags, L, (int) n);
                
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
                                    futureVec.resize((int) n);

                                    for (int i = 0; i < (int) n; i++)
                                    {
                                        futureVec[i] = std::async(std::launch::async, &TRK::MCMC::parallelUpdateAIESWalkers, this, XX, YY, i, fixed_allparams_flags); //pointer to fn run through MT, arguments to fn
                                        counter++;
                                        liveThreads++;

                                        if (liveThreads >= trk.settings.maxThreads)
                                        {
                                            for (int i = completedThreads; i < counter; i++)
                                            {
                                                res = futureVec[i].get();
                                                
                                                if (res.size() > n) { // rejected
                                                    res.pop_back();
                                                    result.push_back(res);
                                                    
                                                    int m = 0;
                                                    if (trk.asymmetric.hasAsymSlop){
                                                        if (trk.settings.do1DFit){
                                                            m = 1;
                                                        } else {
                                                            m = 2;
                                                        }
                                                    }
                                                    
                                                    for (int j = 0; j < trk.bigM + m; j++) {
                                                        (j < trk.bigM + m - 1) ? std::cout << res[j] << " " : std::cout << res[j];
                                                    }
                                                    std::cout  << std::endl;
                                                }
                                                else { // accepted
                                                    all_walkers[i] = res;
                                                    result.push_back(res);
                                                    
                                                    int m = 0;
                                                    if (trk.asymmetric.hasAsymSlop){
                                                        if (trk.settings.do1DFit){
                                                            m = 1;
                                                        } else {
                                                            m = 2;
                                                        }
                                                    }
                                                    
                                                    for (int j = 0; j < trk.bigM + m; j++) {
                                                        (j < trk.bigM + m - 1) ? std::cout << res[j] << " " : std::cout << res[j];
                                                    }
                                                    std::cout  << std::endl;
                                                    
                                                    accept_count++;
                                                   
                                                }
                                                sample_count += 1;
                                            }
                                            completedThreads += liveThreads;
                                            liveThreads = 0;
                                        }
                                    }
                                    for (int i = completedThreads; i < n; i++)
                                    {
                                        res = futureVec[i].get();
                                        
                                        if (res.size() > n) { // rejected
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
                            
                            
//                            if (sample_count >= 40000){
//                                std::cout << std::endl;
//                            }
                            
                            
                            
                            std::vector <double> res = updateAIESWalker(X, YY, fixed_allparams_flags);
                            
                            if (res.size() > n) { // rejected; add initial point to sample
                                res.pop_back();
                                result.push_back(X);
                                
//                                double logL = (trk.statistics.*trk.statistics.selectedLikelihood)(X);
//                                std::cout << -2*logL << std::endl;
                            }
                            else { // accepted; add new point to sample
                                all_walkers[k] = res;
                                result.push_back(res);
                                
//                                double logL = (trk.statistics.*trk.statistics.selectedLikelihood)(res);
//                                std::cout << -2*logL << std::endl;
                                
                                accept_count++;
                               
                            }
                            sample_count += 1;
                        }
                        
                        if (printAIESWalkerEvolution){
                            for (int k = 0; k < L; k++){
                                printf("%f ", all_walkers[k][0]);
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
                printf("Notice: the ARWMH sampler is mostly deprecated, so it isn't updated for use with fixing parameters with MCMC.\n");
                
                std::vector <std::vector <double > > cov_i(n, std::vector <double> (n, 0.0));
                for (int j = 0; j < n; j++){
                    cov_i[j][j] = std::pow(sigmas_guess[j], 2.0);
                }
                std::vector <std::vector <double > > cov_i1(n, std::vector <double> (n, 0.0));
                
                std::vector <double> allparams_trial, allparams_0; //allparams_0 is the previous step
                double log_a, rand_unif_log, accept_frac = 0.0;

                int accept_count = 0;
                int delta_count = 0;
                int tenth = (int) R / 10;
                int prog = 0;
                
                std::vector <double> mu_i = starting_point;
                std::vector <double> mu_i1(n, 0.0);
                std::vector <double> X_i = starting_point;
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
                        for (int j = 0; j < n; j++) {
                            X_trial[j] = rnorm(mu_i[j], std::sqrt(lamb * cov_i[j][j]));
                        }

                        log_a = metHastRatio(X_trial, X_i, fixed_allparams_flags);

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
                    for (int j = 0; j < n; j++) {
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

        result_final = checkSlopSignMCMC(result_final, (int) n, fixed_allparams_flags);
        
        useLogPosterior = false;

        return result_final;
    }


    // Affine Invariant Ensemble Sampler (AIES)
    std::vector <double> TRK::MCMC::updateAIESWalker(std::vector <double> X, std::vector <std::vector <double> > YY, std::vector <bool> fixed_allparams_flags){ // X is the walker to be updated with index k, YY is the set of walkers that the complementary walker for X is randomly chosen from, i.e. X should not be in YY
        int n = (int) X.size();
        
        double a = 2.0; //stretch variable pdf parameter
        std::vector <double> X_trial, Y, res;
        
        // choose some other jth walker Y:
        int j = rand() % (int)YY.size();
        Y = YY[j];
        
        // make proposal vector
        double Z = rstretch(a);
        for (int i = 0; i < n; i++){
            X_trial.push_back(Z * X[i] + (1.0 - Z) * Y[i]);
        }
        
        // accept?
        double rand_unif_log = std::log(runiform(0.0, 1.0));
        double log_a = metHastRatio(X_trial, X, fixed_allparams_flags);
        
        if (rand_unif_log <= log_a + (n - 1) * std::log(Z)) { // accept
            res = X_trial;
        }
        
        else {
            std::vector <double> nan_vec = {NAN};
            std::vector <double> inner_res = concat(X_trial, nan_vec);
            res = {inner_res}; //returns vector one size too big if not accepted
        }
        
        return res;
    }

    std::vector <double> TRK::MCMC::parallelUpdateAIESWalkers(std::vector <std::vector <double> > XX, std::vector <std::vector <double> > YY, int k, std::vector <bool> fixed_allparams_flags){ // XX is set of walkers that contains the kth walker that you want to evolve, YY is set of complementary walkers
        return updateAIESWalker(XX[k], YY, fixed_allparams_flags);
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

//        allparams_sigmas_guess = allparams_sigmas_guess;
        
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

    bool TRK::MCMC::rbool(){
        auto gen = std::bind(std::uniform_int_distribution<>(0,1),std::default_random_engine());
        bool b = gen();

        return b;
    }


    std::vector <std::vector <double >> TRK::MCMC::checkSlopSignMCMC(std::vector <std::vector <double >> result_final, int n, std::vector <bool> fixed_allparams_flags) { // true if fixed
        int num_slop = 0;
        
        if (!fixed_allparams_flags[fixed_allparams_flags.size() - 1]) {
            num_slop++;
            
        };
        
        if (!trk.settings.do1DFit){
            if (!fixed_allparams_flags[fixed_allparams_flags.size() - 2]) {num_slop++;};
        }

        std::vector <std::vector <double >> result_final_fixed;
        std::vector <double> inner;

        for (int i = 0; i < result_final.size(); i++) {
            inner.clear();

            for (int j = 0; j < n - num_slop; j++) { // free non-slop params
                inner.push_back(result_final[i][j]);
            }
            
//            if (!trk.settings.do1DFit){
//                inner.push_back(std::abs(result_final[i][n ])); // x slop
//            }
//
//            inner.push_back(std::abs(result_final[i][trk.M+1])); // y slop
            
            for (int k = 0; k < num_slop; k++) { // free slop params
                inner.push_back(std::abs(result_final[i][n - num_slop + k]));
            }
            
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
            if (verbose){
                printf("Finding pivot point(s)...\n\n");
            }
            
            if (findPivotsManually){
                optimizePivots_Manual();
                return;
            }
            
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
    void TRK::CorrelationRemoval::getCombos(std::vector < std::vector <double> > total, int k, int offset) { //ND case in x

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
        if (!mustProvideLinearParamIndices){
            double intercept_check, slope_check;
            std::vector <double> intercepts, slopes, placeholder;
            std::vector <int> intercept_inds, slope_inds;
            for (int i = 0; i < (int) trk.allparams_guess.size(); i++){
                placeholder.push_back((double) i);
            }
            for (int p = 0; p < P; p++){
                intercepts.push_back(pivot_intercept_functions[p](placeholder));
                slopes.push_back(pivot_slope_functions[p](placeholder));
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
        }
        
        return;
    }

    void TRK::CorrelationRemoval::getPivotGuess(){
        if (findPivotPoints){
            P = (int) pivot_intercept_functions.size();
            
            final_pivots = std::vector <double>(P, 0.0);
            
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
                    std::vector <double> x_s_vec = {x_s[0]};
                    std::vector <double> x_s_vec_1 = {x_s[trk.N-1]};
                    
                    std::vector <double> cat_1 = concat(x_s_vec, divisions);
                    divisions = concat(cat_1 , x_s_vec_1); // include first and last datapoints with divisions
                    
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
            
            findLinearParamIndices();
        }
        return;
    }

    std::vector <double> TRK::CorrelationRemoval::refitAnalytic(double new_pivot, int p){ // refit with new pivot point; only intercepts changed
        std::vector <double> allparams_better = trk.allparams_guess, allparams_old = trk.allparams_guess;

        double intercept_new, intercept_old, slope;

        if (verbose_refit) {
            printf("\n\n\n(inner refit analytic) old pivot %i, new pivot %i = %f\t%f\n", p + 1, p + 1, pivots[p], new_pivot);
        }
        
        intercept_old = pivot_intercept_functions[p](allparams_old);
        slope = pivot_slope_functions[p](allparams_old);

        intercept_new = intercept_old + slope * (new_pivot - pivots[p]);
        
        allparams_better[intercept_indices[p]] = intercept_new;
        
        return allparams_better;
    }

    std::vector <double> TRK::CorrelationRemoval::refitWithNewPivots(double new_pivot, int p){
        std::vector <double> allparams_better = refitAnalytic(new_pivot, p); // determine the best fit
        

        if (refit_with_simplex){
            allparams_better = trk.optimization.downhillSimplex_Fit(trk.statistics.selectedChiSq, allparams_better, trk.scaleOptimization.s, trk.optimization.showFittingSteps);
        }
        
        if (trk.correlationRemoval.verbose_refit){
            if (refit_with_simplex){
                printf("re-fit for new pivot point(s); old / new params:\n");

                for (int j = 0; j < (int)allparams_better.size(); j++){
                    printf("%.3e %.3e\n", trk.allparams_guess[j], allparams_better[j]);
                }
                
                double fitness_old = (trk.statistics.*trk.statistics.selectedChiSq)(trk.allparams_guess, trk.scaleOptimization.s);
                double fitness_new = (trk.statistics.*trk.statistics.selectedChiSq)(allparams_better, trk.scaleOptimization.s);
                
                printf("old/new fitness: %f / %f\n", fitness_old, fitness_new);
            } else {
                printf("re-fit for new pivot point(s); old / new intercept %i:\n", p + 1);

                printf("%.3e %.3e\n", trk.allparams_guess[intercept_indices[p]], allparams_better[intercept_indices[p]]);
                
                double fitness_old = (trk.statistics.*trk.statistics.selectedChiSq)(trk.allparams_guess, trk.scaleOptimization.s);
                double fitness_new = (trk.statistics.*trk.statistics.selectedChiSq)(allparams_better, trk.scaleOptimization.s);
                
                printf("old/new fitness: %f / %f\n", fitness_old, fitness_new);
            }
        }
        
        return allparams_better;
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
            
            std::vector <bool> fixed_allparams_flags_default(trk.bigM, false); // none fixed
            
            allparam_samples = trk.mcmc.samplePosterior(sample_R, sample_burnIn, trk.mcmc.allparams_sigmas_guess, fixed_allparams_flags_default); //allparam_samples is { {allparams0}, {allparams1}, ... }
        
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
                std::string filename = trk.settings.outputPath + std::string("/TRKpivots") + std::to_string(iter) + std::string("_") + std::to_string(finalPivots[0]) + std::string(".txt");

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
            
            if (refit_newPivot){
                for (int p = 0; p < P; p++){
                    trk.allparams_guess = refitWithNewPivots(finalPivots[p], p);
                }
            };
            
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

            std::vector <bool> fixed_allparams_flags_default(trk.bigM, false); // none fixed
            
            allparam_samples = trk.mcmc.samplePosterior(sample_R, sample_burnIn, trk.mcmc.allparams_sigmas_guess, fixed_allparams_flags_default); //allparam_samples is { {allparams0}, {allparams1}, ... }
            
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
                        if (trk.settings.do1DFit){
                            m = 1;
                        } else {
                            m = 2;
                        }
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
            
            if (refit_newPivot){
                for (int p = 0; p < P; p++){
                    trk.allparams_guess = refitWithNewPivots(finalPivots[p], p);
                }
            };
            
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
        
        if (parallelize){
            switch (trk.settings.ParallelizationBackEnd){ // search for multiple pivot points simultaneously
                case CPP11:
                    optimizePivots_Correlation_CPP11();
                    break;
                case OPENMP:
                    optimizePivots_Correlation_Default(); // calls omp pragma within this
                    break;
                default:
                    optimizePivots_Correlation_Default();
                    break;
            }
        } else {
            optimizePivots_Correlation_Default();
        }
        
        if (refit_newPivot){ // modify the parameters guess (intercepts) based off new pivots
            for (int p = 0; p < P; p++){
                trk.allparams_guess = refitWithNewPivots(final_pivots[p], p);
            }
        };
        
        pivots = final_pivots; // store final pivots
        
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
                    std::function <double(std::vector <double> )> correlation_func_simplex = std::bind(&TRK::CorrelationRemoval::getAbsCorrFromNewPivot_Wrapper, this, std::placeholders::_1, p);
                    
                    std::vector <double> guess = {pivots[p]};
                    
                    futureVec[p] = std::async(std::launch::async, &TRK::Optimization::downhillSimplex_1DWrapper, trk.optimization, correlation_func_simplex, guess, correlation_tol, showSimplexSteps, max_corr_simplex_iters);
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
                    final_pivots[p] = result;
                }
                completedThreads += liveThreads;
                liveThreads = 0;
            }
        }
        for (int p = completedThreads; p < P; p++)
        {
            result = futureVec[p].get();
            final_pivots[p] = result;
        }
        
        return;
    }

    void TRK::CorrelationRemoval::optimizePivots_Correlation_Default(){
        if (trk.settings.ParallelizationBackEnd == OPENMP && parallelize){
            #pragma omp parallel for num_threads(maxThreads)
        }
        for (int p = 0; p < P; p++){ // each pivot can be optimized independently (the value of the others pivots shouldn't affect it
            if (verbose){
                printf("Optimizing pivot %i...\n\n", p + 1);
            }
            
//            if (p != 2){
//                continue;
//            }

//            trk.optimization.getBetterGuess();
            
            
            switch (thisPivotMethod){
                case PEARSON_GSS : {
                    // function to be minimized:
                    std::function <double(double)> correlation_func = std::bind(&TRK::CorrelationRemoval::getAbsCorrFromNewPivot, this, std::placeholders::_1, p);
                    
                    final_pivots[p] = trk.optimization.goldenSectionSearch(correlation_func, min_pivots_brackets[p], max_pivots_brackets[p], correlation_tol);
                    break;
                }
                case PEARSON_SIMPLEX : {
                    // function to be minimized:
                    std::function <double(std::vector <double> )> correlation_func = std::bind(&TRK::CorrelationRemoval::getAbsCorrFromNewPivot_Wrapper, this, std::placeholders::_1, p);
                    
                    std::vector <double> pivot_vec = {pivots[p]};
                    final_pivots[p] = trk.optimization.downhillSimplex(correlation_func, pivot_vec, correlation_tol, showSimplexSteps, max_corr_simplex_iters, save_simplex_correlation_evals)[0];
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
                if (trk.settings.do1DFit){
                    m = 1;
                } else {
                    m = 2;
                }
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

    void TRK::CorrelationRemoval::rejectLinearParamOutliers(std::vector <double> b_samples, std::vector <double> m_samples){
//        if (verbose) {printf("Performing RCR on slope and intercept sample...\n");}
        
        printf("deprecated...\n");
        
//        using namespace RCRLib;
//        RCR rcr_b = RCR(LS_MODE_DL);
//        RCR rcr_m = RCR(LS_MODE_DL);
//
//        rcr_b.performBulkRejection(b_samples);
//        rcr_m.performBulkRejection(m_samples);
//
//        b_samples = rcr_b.result.cleanY;
//        m_samples = rcr_m.result.cleanY;
//
//        b_samples.size() >= m_samples.size() ? b_samples.resize(m_samples.size()) : m_samples.resize(b_samples.size());

//        if (verbose) {printf("fraction of (%f, %f) (b, m) samples rejected by RCR.\n", (double) rcr_b.result.rejectedY.size() / rcr_b.result.originalY.size(), (double) rcr_m.result.rejectedY.size() / rcr_m.result.originalY.size() );}
//
        return;
    }

    std::vector <bool> TRK::CorrelationRemoval::getFixedLinearParams(int p){
        std::vector <bool> fixed_allparams_flags(trk.bigM, false); // everything free by default
        
        if (sampleOnlyLinearParams_pivots) {
            fixed_allparams_flags = std::vector <bool> (trk.bigM, true); // everything but slope and intercept fixed
            
            fixed_allparams_flags[intercept_indices[p]] = false;
            fixed_allparams_flags[slope_indices[p]] = false;
        }
        
        return fixed_allparams_flags; // true if the param is fixed
    }

    void TRK::CorrelationRemoval::getLinearParamSamples(std::vector < std::vector <double> > &allparam_samples, std::vector <double> &b_samples, std::vector <double> &m_samples, int p){
        if (sampleOnlyLinearParams_pivots) { // if so than, the sampled params are just the intercept and slope for this pivot point (need to find which is which, however)
            int intercept_index = 0; // intercept first by default
            if (intercept_indices[p] > slope_indices[p]) { intercept_index = 1; }; // intercept second
            
            for (int j = 0; j < allparam_samples.size(); j++) {
                b_samples.push_back(allparam_samples[j][intercept_index]);
                m_samples.push_back(allparam_samples[j][intercept_index ^ 1]);
            }
            
        } else {
            std::vector < std::vector <double> > param_samples(sample_R, std::vector<double>());
            
            for (int j = 0; j < allparam_samples.size(); j++) {
                param_samples[j] = slice(allparam_samples[j], 0, (int)trk.M);
            }
            for (int j = 0; j < param_samples.size(); j++) {
                b_samples.push_back(pivot_intercept_functions[p](param_samples[j]));
                m_samples.push_back(pivot_slope_functions[p](param_samples[j]));
            }
        }
        
        return;
    }

    double TRK::CorrelationRemoval::getAbsCorrFromNewPivot(double new_pivot, int p){
        std::vector <double> b_samples, m_samples, original_pivots = pivots;
        double rxy, abs_rxy; // correlation between slope and intercept
        
        if (verbose){
            printf("\n\n\nold pivot %i, new pivot %i = %f\t%f\n", p + 1, p + 1, pivots[p], new_pivot);
        }
        
        std::vector <double> initial_allparams_guess = trk.allparams_guess;
        
        if (refit_newPivot){ // find better starting place for sampling given new pivot
            trk.allparams_guess = refitWithNewPivots(new_pivot, p);
        };
        
        pivots[p] = new_pivot; // set pivot to new value, and sample parameter space with it
        
        if (trk.correlationRemoval.verbose && trk.mcmc.verbose){
            printf("\nSampling for pivot points...\n");
        }
        
        std::vector <bool> fixed_allparams_flags = getFixedLinearParams(p); // (potentially) sample with only the relavant linear params free
        
        std::vector < std::vector <double> > allparam_samples = trk.mcmc.samplePosterior(sample_R, sample_burnIn, trk.mcmc.allparams_sigmas_guess, fixed_allparams_flags);
        
        getLinearParamSamples(allparam_samples, b_samples, m_samples, p);
        
        if (RCR_samples){
            rejectLinearParamOutliers(b_samples, m_samples); // (potentially) reject outlier params with RCR
        }

        switch (whichCorrelationUsed){
            case SPEARMAN:
                rxy = trk.statistics.spearmanCorrelation(b_samples, m_samples);
                break;
            case PEARSON:
                rxy = trk.statistics.pearsonCorrelation(b_samples, m_samples);
                break;
        }
        
        if (verbose){
            double spearmanR = trk.statistics.spearmanCorrelation(b_samples, m_samples);
            double pearsonR = trk.statistics.pearsonCorrelation(b_samples, m_samples);
        
            if (std::abs(spearmanR - pearsonR) > 0.2){
                printf("differing correlations; pearson = %f\t spearman = %f\n", trk.statistics.pearsonCorrelation(b_samples, m_samples), trk.statistics.spearmanCorrelation(b_samples, m_samples));
            }
            
            printf("pearson, spearman, pivot %i = %f %f %f\n", p + 1, pearsonR, spearmanR, new_pivot);
        }
        
        
        abs_rxy = std::isnan(rxy) ? 1.0 : std::abs(rxy); // returns maximally correlated if NaN
        
        writeCorrelationOptimizationSampling(b_samples, m_samples, p);
        
        pivots = original_pivots; //reset pivot(s) and intercepts to previous values
        trk.allparams_guess = initial_allparams_guess;
        
        return abs_rxy; // returns maximally correlated if NaN
    }

    double TRK::CorrelationRemoval::getAbsCorrFromNewPivot_Wrapper(std::vector <double> new_pivot, int p){
        return getAbsCorrFromNewPivot(new_pivot[0], p);
    }

    // MANUAL method (for testing only)
    void TRK::CorrelationRemoval::optimizePivots_Manual(){
        printf("FINDING PIVOTS MANUALLY!!!\n");
        
        // fix to best fit
        std::vector <double> fixed_allparam_vals, allparams_new, allparams_best_free_0  = {4.14979e+00, 1.72651e+00, 2.75165e+00, 0.26189, 0.31024};
        double s_best = 8.677e-01;
        trk.scaleOptimization.s = s_best;
        pivots = {
            1.6636E-01, 1.4205
        };
        
        std::vector <bool> fixed_allparam_flags(trk.bigM, false);
        fixed_allparam_flags[3] = true;
        
        
        double b, m, theta = -6.91534e-04, delta_theta = 1.0E-6;
        int total_iter = 20, iter_count = 0;
        
        theta -= (total_iter / 2) * delta_theta;
        
        while (iter_count <= total_iter){
            
            theta += delta_theta;
            
            fixed_allparam_vals = {theta};
            
            //refit with new theta
            trk.scaleOptimization.whichExtrema = S;
            allparams_new = trk.optimization.downhillSimplex_Fit_Fixed(
                                                                     trk.statistics.selectedChiSq,
                                                                     allparams_best_free_0,
                                                                     s_best,
                                                                     trk.optimization.showFittingSteps,
                                                                     fixed_allparam_flags,
                                                                     fixed_allparam_vals
                                                                 );
            trk.scaleOptimization.whichExtrema = ANY;
            
            
            m = std::tan(theta);
            b = allparams_new[2];
            
            printf("%.20f\t%.20f\n", b, m);
            
            iter_count++;
        }
        
        return;
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
        
        if (trk.settings.do1DFit){
            // SLOP
            if (slop_y_minus_guess >= 0){
                hasAsymSlop = true;
            } else {
                slop_y_minus_guess = trk.slop_y_guess;
            }
            
            // ERROR BARS
            if (sy_minus.size() == trk.N){
                hasAsymEB = true;
            } else {
                sy_minus = trk.sy;
            }
        } else {
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
        }
        
        // likelihood selection
        
        if (hasAsymSlop || hasAsymEB){
            if (trk.settings.do1DFit){
                trk.statistics.selectedChiSq = &TRK::Statistics::regularChiSquaredWSlopAsym;
                trk.statistics.selectedLikelihood = &TRK::Statistics::likelihood1DAsym;
                if (trk.settings.verbose){
                    printf("Running 1D asymmetric fit.\n");
                }
            } else {
                trk.statistics.selectedChiSq = &TRK::Statistics::modifiedChiSquaredAsym;
                trk.statistics.selectedLikelihood = &TRK::Statistics::likelihoodAsym;
            }
                
                
            
        } else { // symmetric case
            if (trk.settings.do1DFit){
                trk.statistics.selectedChiSq = &TRK::Statistics::regularChiSquaredWSlop;
                trk.statistics.selectedLikelihood = &TRK::Statistics::likelihood1D;
                if (trk.settings.verbose){
                    printf("Running 1D symmetric fit.\n");
                }
            }
        }
        
        if (verbose){
            printf("Asymmetries: slop: %s\tError bars: %s\n", hasAsymSlop ? "true" : "false", hasAsymEB ? "true" : "false");
        }
        
        // add asymm slops to all params guess
        
        if (hasAsymSlop){
            if (!trk.settings.do1DFit){
                trk.allparams_guess.push_back(slop_x_minus_guess);
            }
            trk.allparams_guess.push_back(slop_y_minus_guess);
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
        
        // NOTE:
        // Make sure to update this function following edits to getAsymShift1D
        
        // TO DO: modify this to use old code (see getAsymShift1D)
        
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
        
        std::vector <double> x_slops = {slops[0], slops[2]};
        std::vector <double> y_slops = {slops[1], slops[3]};
        std::vector <double> x_EBs = {EBs[0], EBs[2]};
        std::vector <double> y_EBs = {EBs[1], EBs[3]};
        
        deltaxn = getAsymShiftSingle(x_slops, x_EBs);
        deltayn = getAsymShiftSingle(y_slops, y_EBs);
        
        return {deltaxn, deltayn};
        
//        // Y SHIFT
//        std::vector <double> sigma_vec = {slops[1], slops[3]};
//        std::vector <double> sigman_vec = {EBs[1], EBs[3]};
//
//        double sigmaL = minMax(sigma_vec)[1];
//        double sigmaS = minMax(sigma_vec)[0];
//        double sigmanL = minMax(sigman_vec)[1];
//        double sigmanS = minMax(sigman_vec)[0];
//
//        std::vector <double> sigmaMax_vec = {sigmaL, sigmanL};
//        double sigmaMax = minMax(sigmaMax_vec)[1];
//
//
//        double xi = sigmaS/sigmaL + sigmanS / sigmanL;
//        double eta = sigmanL < sigmaL ? sigmanS/sigmanL - sigmaS/sigmaL : sigmaS/sigmaL - sigmanS/sigmanL;
//        double r = minMax(sigmaMax_vec)[0] / minMax(sigmaMax_vec)[1];
//
//        double xip = xi <= 1 ? xi : 2.0 - xi;
//        double etap = xip == 0 ? 0 : 2.0 * xip * std::pow(0.5*eta/xip + 1.0, std::pow(r, -0.4087)) - xip;
//
//        double Nr = -0.5326*std::pow(r,2.0) + 1.5307*r + 0.0019;
//        double fXi = xi <= 1 ? 0.2454*std::pow(xi, -1.1452) : 0.2454*std::pow(xi, -0.5203);
//        double gEtaP = std::pow(etap,2.0);
//        double hXi = -0.042*std::pow(xi,2.0) - 0.1602 * xi + 0.4884;
//
//        double deltastr = sigmaMax * Nr * (fXi * gEtaP + hXi);
//
//        int i = 1;
//        if (slops[1] == slops[3] || EBs[1] == EBs[3]){ //one of the dists is symmetric
//            if (sigmanL == EBs[1] || sigmaL == slops[3]){
//                i = 1;
//            } else if (sigmaMax == EBs[3] || sigmaMax == slops[1]){
//                i = -1;
//            }
//
//            deltayn = i * deltastr;
//
//        } else if ((sigmaL == slops[3] && sigmanL == EBs[1]) || (sigmaL = slops[1] && sigmanL == EBs[3])){ //both asymm first case
//            if (sigmaMax == EBs[1] || sigmaMax == slops[3]){
//                i = 1;
//            } else if (sigmaMax == EBs[3] || sigmaMax == slops[1]){
//                i = -1;
//            }
//
//            deltayn = i * deltastr;
//
//        } else if ((sigmaL == slops[1] && sigmanL == EBs[1]) || (sigmaL = slops[3] && sigmanL == EBs[3])){ //both asymm second case
//            if (sigmaMax == EBs[1] || sigmaMax == slops[3]){
//                i = 1;
//            } else if (sigmaMax == EBs[3] || sigmaMax == slops[1]){
//                i = -1;
//            }
//
//            double pwr = eta <= 1 ? 0.7413 : -0.1268;
//            deltayn = i * deltastr * std::sin(PI/2.0 * etap/xip) * std::pow(eta, pwr);
//        }
//
//        // X SHIFT
//
//        sigma_vec = {slops[0], slops[2]};
//        sigman_vec = {EBs[0], EBs[2]};
//
//
//        sigmaL = minMax(sigma_vec)[0];
//        sigmaS = minMax(sigma_vec)[0];
//        sigmanL = minMax(sigman_vec)[0];
//        sigmanS = minMax(sigman_vec)[0];
//
//        sigmaMax_vec = {sigmaL, sigmanL};
//        sigmaMax = minMax(sigmaMax_vec)[0];
//
//
//        xi = sigmaS/sigmaL + sigmanS / sigmanL;
//        eta = sigmanL < sigmaL ? sigmanS/sigmanL - sigmaS/sigmaL : sigmaS/sigmaL - sigmanS/sigmanL;
//        r = minMax(sigmaMax_vec)[0] / minMax(sigmaMax_vec)[0];
//
//        xip = xi <= 1 ? xi : 2.0 - xi;
//        etap = xip == 0 ? 0 : 2.0 * xip * std::pow(0.5*eta/xip + 1.0, std::pow(r, -0.4087)) - xip;
//
//        Nr = -0.5326*std::pow(r,2.0) + 1.5307*r + 0.0019;
//        fXi = xi <= 1 ? 0.2454*std::pow(xi, -1.1452) : 0.2454*std::pow(xi, -0.5203);
//        gEtaP = std::pow(etap,2.0);
//        hXi = -0.042*std::pow(xi,2.0) - 0.1602 * xi + 0.4884;
//
//        deltastr = sigmaMax * Nr * (fXi * gEtaP + hXi);
//
//        if (slops[0] == slops[2] || EBs[0] == EBs[2]){ //one of the dists is symmetric
//            if (sigmanL == EBs[0] || sigmaL == slops[2]){
//                i = 1;
//            } else if (sigmaMax == EBs[2] || sigmaMax == slops[0]){
//                i = -1;
//            }
//
//            deltaxn = i * deltastr;
//
//        } else if ((sigmaL == slops[2] && sigmanL == EBs[0]) || (sigmaL = slops[0] && sigmanL == EBs[2])){ //both asymm first case
//            if (sigmaMax == EBs[0] || sigmaMax == slops[2]){
//                i = 1;
//            } else if (sigmaMax == EBs[2] || sigmaMax == slops[0]){
//                i = -1;
//            }
//
//            deltaxn = i * deltastr;
//
//        } else if ((sigmaL == slops[0] && sigmanL == EBs[0]) || (sigmaL = slops[2] && sigmanL == EBs[2])){ //both asymm second case
//            if (sigmaMax == EBs[0] || sigmaMax == slops[2]){
//                i = 1;
//            } else if (sigmaMax == EBs[2] || sigmaMax == slops[0]){
//                i = -1;
//            }
//
//            double pwr = eta <= 1 ? 0.7413 : -0.1268;
//            deltaxn = i * deltastr * std::sin(PI/2.0 * etap/xip) * std::pow(eta, pwr);
//        }
    }

    double TRK::Asymmetric::getAsymShift1D_old(std::vector <double> allparams, int n){
        double sigL,sigS,signL,signS,sigmax,x,y,f,xi,eta,nf,normf,delplus;
        double sigp = 0.0, sigm = 0.0, signp = 0.0, signm = 0.0;
        double pi = 3.1415926536;
        double delta;
        int idir;
        
        if (hasAsymSlop && !hasAsymEB){
//                slops.push_back(allparams[trk.M+1]);
//                EBs = concat(EBs, EBs);
            
            sigp = allparams[trk.M];
            sigm = allparams[trk.M + 1];
            signp = trk.sy[n];
            signm = trk.sy[n];
        } else if (!hasAsymSlop && hasAsymEB){
//                slops = concat(slops, slops);
//                EBs.push_back(sy_minus[n]);
            
            sigp = allparams[trk.M];
            sigm = allparams[trk.M];
            signp = trk.sy[n];
            signm = sy_minus[n];
        } else if (hasAsymSlop && hasAsymEB){
//                slops.push_back(allparams[trk.M+1]);
//                EBs.push_back(trk.asymmetric.sy_minus[n]);
            
            sigp = allparams[trk.M];
            sigm = allparams[trk.M + 1];
            signp = trk.sy[n];
            signm = sy_minus[n];
        }
        
        // SWITCH SLOP SIGMAS (testing)
//        printf("NOTICE: removed slop switch\n\n\n\n");
        double sigpold = sigp;
        double sigmold = sigm;
        sigm = sigpold;
        sigp = sigmold;
        
//        printf("TESTING SIGMA VALUES FOR 1D ASYM SHIFT.\n\n\n\n");
//        sigp = 10.0 * sigp;
//        sigm = 10.0 * sigm;
        
        
//        printf("old code: sigp = %f\t sigm = %f\t signp = %f\t signm = %f\n\n", sigp, sigm, signp, signm);
        
        //For a given asymmetric Gaussian, "L" indicates the larger sigma and "S" indicates the smaller sigma
        sigL = minMax({sigp,sigm})[1];
        sigS = minMax({sigp,sigm})[0];
        signL = minMax({signp,signm})[1];
        signS = minMax({signp,signm})[0];
        sigmax = minMax({sigL,signL})[1];
        
//        printf("\n\n\n\n%f\t%f\t%f\t%f\t%f\n", sigL, sigS, signL, signS, sigmax);

        //The direction of the shift is determined by the direction of the maximum of all four sigmas
        //(but is not necessarily in the same direction...the model below can give negative delta in certain cases)
        if (sigmax == sigp || sigmax == signp)
        {
            idir = 1;
        }
        else
        {
            idir = -1;
        }
        
        //Model deltas were derived from fits to surfaces in y=sigS/sigL vs x=signS/signL space, where signL=signmax=signp=1.
        //If sigmax=sigL instead of signL, x=sigS/sigL and y=signS/signL.
        //The model gives delta(x,y)/sigmax for different values of f=min(sigL,signL)/max(sigL,signL).
        if (signL >= sigL)
        {
            x = signS/signL;
            y = sigS/sigL;
        }
        else
        {
            x = sigS/sigL;
            y = signS/signL;
        }
        //When sigL and signL are in the same direction, the model is of the form
        //delta(x,y,f)=N(f)*[f(xi)g(eta)+h(xi)], where xi=(y+x) and eta is a *transformed* value of (y-x).
        //The eta transformation was chosen so as to ensure delta->0 at appropriate points in the x-y plane for different values of f,
        //when sigL and signL are in opposite directions.
        //See Adam Trotter or Dan Reichart for details.
        f = minMax({sigL,signL})[0]/sigmax;
        if (f == 0.0)
        {
            nf = 0.0;
        }
        else
        {
            nf = std::pow(f,-0.4087);
        }
        
        xi = x+y;
        
//        printf("%f\n", xi);
        
        if (xi == 0.0 || xi == 2.0)
        {
            eta = 0.0;
        }
        else if (xi <= 1.0)
        {
            eta = 2.0*xi*std::pow(0.5*((y-x)/xi+1.0),nf)-xi;
        }
        else
        {
            if (x == 1.0){
                eta = -(2.0-xi);
            } else
            {
                eta = 2.0*(2.0-xi)*std::pow(0.5*((y-x)/(2.0-xi)+1.0),nf)-(2.0-xi);
            }
        }
        
//        printf("%f\n", eta);
        
        normf = -0.5326*f*f+1.5307*f+0.0019;
        if (xi == 0.0)
        {
            delplus = normf*0.4884;
        }
        else if (xi <= 1.0)
        {
            delplus = normf*(0.2454*std::pow(xi,-1.1452)*eta*eta-0.042*xi*xi-0.1602*xi+0.4884);
        }
        else
        {
            delplus = normf*(0.2454*std::pow(xi,-0.5203)*eta*eta-0.042*xi*xi-0.1602*xi+0.4884);
        }

        //Case of sigL and signL in same direction
        if ((sigL == sigp && signL == signp) || (sigL == sigm && signL == signm))
        {
            delta = delplus;
        }
        else
        //Case of sigL and signL in opposite directions: the above model is multiplied by a sinusoidal (odd) function in eta.
        {
            if (xi == 0.0 || xi == 2.0)
            {
                delta = 0.0;
            }
            else if (xi <= 1.0)
            {
                delta = delplus*std::pow(xi,0.7413)*std::sin(pi/2.0*eta/xi);
            }
            else
            {
                delta = delplus*std::pow(xi,-0.1268)*std::sin(pi/2.0*eta/(2.0-xi));
            }
        }
        //cout << "(sigp,sigm) = (" << sigp << ", " << sigm << "), (signp,signm) = (" << signp << ", " << signm << ")" << endl;
        //cout << "f = " << f << " x = " << x << " y = " << y << " xi = " << xi << " eta = " << eta << " delta = " << delta << endl;
        
        //Scale the model delta by sigmax, and flip its direction if sigmax is in the negative direction
        double delta_n = delta*idir*sigmax;
        
        return delta_n;
    }

double TRK::Asymmetric::getAsymShiftSingle(std::vector <double> slops_1D, std::vector <double> errorbars_1D){ // gets asym delta shift along one direction
    double delta_n = 0;
    
    std::vector <double> sigma_vec = slops_1D;
    std::vector <double> sigman_vec = errorbars_1D;
    
    sigma_vec = {std::abs(sigma_vec[1]), std::abs(sigma_vec[0])};

//        printf("TESTING SIGMA VALUES FOR 1D ASYM SHIFT.\n\n\n\n");
//        sigma_vec = {10.0 * sigma_vec[0], 10.0 * sigma_vec[1]};
    
//        printf("new code: sigp = %f\t sigm = %f\t signp = %f\t signm = %f\n\n", sigma_vec[0], sigma_vec[1], sigman_vec[0], sigman_vec[1]);
    
    // SHIFT
    double sigmaL = minMax(sigma_vec)[1];
    double sigmaS = minMax(sigma_vec)[0];
    double sigmanL = minMax(sigman_vec)[1];
    double sigmanS = minMax(sigman_vec)[0];
    
    std::vector <double> sigmaMax_vec = {sigmaL, sigmanL};
    double sigmaMax = minMax(sigmaMax_vec)[1];
    
    
    int i = -1;
    
    if (sigmaMax == sigma_vec[0] || sigmaMax == sigman_vec[0])
    {
        i = 1;
    }


//            printf("\n\n\n\n%f\t%f\t%f\t%f\t%f\n", sigmaL, sigmaS, sigmanL, sigmanS, sigmaMax);
    
    
    double xi = (sigmaS/sigmaL) + (sigmanS / sigmanL);
    double r = minMax(sigmaMax_vec)[0] / minMax(sigmaMax_vec)[1];
    double eta;

    if (sigmanL < sigmaL){
        eta = (sigmanS/sigmanL) - (sigmaS/sigmaL);
    } else {
        eta = (sigmaS/sigmaL) - (sigmanS/sigmanL);
    }

//            printf("%f\n", xi);
    
    double xip;
    double etap;
    if (xi <= 1){
        xip = xi;
    } else {
        xip = 2.0 - xi;
    }

    if (xip == 0.0){
        etap = 0.0;
    } else {
        etap = 2.0 * xip * std::pow(0.5*((eta/xip) + 1.0), std::pow(r, -0.4087)) - xip;
    }


//            printf("%f\n", etap);
    
    
//        if (xi == 0.0 || xip == 0.0){
//            printf("xi = %f \t xi==0.0: %s\n", xi, xi == 0.0 ? "true" : "false");
//            printf("xip = %f \t xip==0.0: %s\n\n", xip, xip == 0.0 ? "true" : "false");
//        }
    
    double Nr = -0.5326*std::pow(r,2.0) + 1.5307*r + 0.0019;
    
    double fXi;
    if (xi == 0.0){
        fXi = 0.0;
    } else if ( xi <= 1 ){
        fXi = 0.2454*std::pow(xi, -1.1452);
    } else {
        fXi = 0.2454*std::pow(xi, -0.5203);
    }
    
    double gEtaP = std::pow(etap,2.0);
    double hXi = -0.042*std::pow(xi,2.0) - 0.1602 * xi + 0.4884;
    
    double deltastr = sigmaMax * Nr * (fXi * gEtaP + hXi);

    if (sigma_vec[0] == sigma_vec[1] || sigman_vec[0] == sigman_vec[1]){ //one of the dists is symmetric
//            if (sigmanL == sigman_vec[0] || sigmaL == sigma_vec[1]){
//                i = 1;
//            } else if (sigmaMax == sigman_vec[1] || sigmaMax == sigma_vec[0]){
//                i = -1;
//            }
        
        delta_n = i * deltastr;
        
    } else if ((sigmaL == sigma_vec[1] && sigmanL == sigman_vec[0]) || (sigmaL = sigma_vec[0] && sigmanL == sigman_vec[1])){ //both asymm first case
//            if (sigmaMax == sigman_vec[0] || sigmaMax == sigma_vec[1]){
//                i = 1;
//            } else if (sigmaMax == sigman_vec[1] || sigmaMax == sigma_vec[0]){
//                i = -1;
//            }
        
        delta_n = i * deltastr;
        
    } else if ((sigmaL == sigma_vec[0] && sigmanL == sigman_vec[0]) || (sigmaL = sigma_vec[1] && sigmanL == sigman_vec[1])){ //both asymm second case
//            if (sigmaMax == sigman_vec[0] || sigmaMax == sigma_vec[1]){
//                i = 1;
//            } else if (sigmaMax == sigman_vec[1] || sigmaMax == sigma_vec[0]){
//                i = -1;
//            }
        
        double pwr = eta <= 1 ? 0.7413 : -0.1268;
        delta_n = i * deltastr * std::sin((PI/2.0) * (etap/xip)) * std::pow(eta, pwr);
    }
    
    if (std::abs(delta_n) >= 10000){
        printf("Error: asym shift huge\n");
    }
    
    return delta_n;
}

    double TRK::Asymmetric::getAsymShift1D(std::vector <double> allparams, int n){
        double delta_n = 0.0;
        
//        if (use_new_1D_shift_code){
        std::vector <double> slops = {allparams[trk.M]}; // +, -
        std::vector <double> EBs = {trk.sy[n]};
        
        if (hasAsymSlop && !hasAsymEB){
            slops.push_back(allparams[trk.M+1]);
            EBs = concat(EBs, EBs);
        } else if (!hasAsymSlop && hasAsymEB){
            slops = concat(slops, slops);
            EBs.push_back(sy_minus[n]);
        } else if (hasAsymSlop && hasAsymEB){
            slops.push_back(allparams[trk.M+1]);
            EBs.push_back(trk.asymmetric.sy_minus[n]);
        }
        
        delta_n = getAsymShiftSingle(slops, EBs);
    
        // not used below -v
        if (!use_new_1D_shift_code){
            delta_n = getAsymShift1D_old(allparams, n);
        }
            
//        }
//        printf("delta_n = %f\t delta_old = %f\n", delta_n, delta_old);
        
//        delta_n = 1000;
//        printf("%f\n", delta_n);
        
        
        return delta_n;
    }


    std::vector <double> TRK::Asymmetric::getAsymSigs2(std::vector <double> allparams, int n){
        std::vector <double> Sigs2(4, 0.0);
        double slopxplus = allparams[trk.M], slopyplus = allparams[trk.M+1]; // slop x, slop y, both positive
        double EBxplus = trk.sx[n], EByplus = trk.sy[n]; // same for error bars
        
        double slopxminus = slopxplus, slopyminus = slopyplus; // default symmetric case
        double EBxminus = EBxplus, EByminus = EByplus;
        
        if (hasAsymSlop && !hasAsymEB){
            slopxminus = allparams[trk.M + 2];
            slopyminus = allparams[trk.M + 3];
        } else if (!hasAsymSlop && hasAsymEB){
            EBxminus = sx_minus[n];
            EByminus = sy_minus[n];
            
        } else if (hasAsymSlop && hasAsymEB){
            slopxminus = allparams[trk.M + 2];
            slopyminus = allparams[trk.M + 3];
            
            EBxminus = sx_minus[n];
            EByminus = sy_minus[n];
        }
        
        Sigs2[0] = std::pow(slopxminus, 2.0) + std::pow(EBxplus, 2.0);
        Sigs2[1] = std::pow(slopyminus, 2.0) + std::pow(EByplus, 2.0);
        Sigs2[2] = std::pow(slopxplus, 2.0) + std::pow(EBxminus, 2.0);
        Sigs2[3] = std::pow(slopyplus, 2.0) + std::pow(EByminus, 2.0);
        
        return Sigs2; // Sigx+, Sigy+, Sigx-, Sigy-
    }

    std::vector <double> TRK::Asymmetric::tangentParallelAsym(std::vector <double> allparams, int n, double s) {
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

    std::vector <double> TRK::Asymmetric::tangentParallelLikelihoodAsym(std::vector <double> allparams, int n) {
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


    // COVID19 MODELING ##########################################################################################################

    void TRK::COVID19::printCustomResults(){
        
//        printf("\n\n\nFinal COVID-19 custom fit:\n\n\n");
        
        // pivots and median point
        for (int p = 0; p < trk.correlationRemoval.P; p++){
            printf("%.3e ", trk.results.pivots[p]);
        }
        printf("%.3e\n", trk.covid19.tmed);
        
        // linear params
        for (int p = 0; p < trk.correlationRemoval.P; p++){
            printf("%.3e ", trk.results.bestFitParams[trk.correlationRemoval.intercept_indices[p]]);
            printf("%.3e", trk.results.bestFitParams[trk.correlationRemoval.slope_indices[p]]);
            if (p < trk.correlationRemoval.P - 1){
                printf(" ");
            }
        }
        printf("\n");
        
        // smoothing params, slop and chi2
        for (int s = 0; s < S; s++){
            printf("%.3e ", trk.results.bestFitParams[trk.correlationRemoval.slope_indices[trk.correlationRemoval.P - 1] + s + 1]);
        }
        printf("%.3e %.3e\n", trk.results.slop_y, trk.results.fitness);
        
        // polynomial params
        int D = (int) (trk.results.bestFitParams.size() - trk.correlationRemoval.P*2 - S) / 3; // polynomial order
        for (int m = 0; m < 3; m++){
            for (int d = 0; d < D; d++){
                printf("%.3e", trk.results.bestFitParams[trk.correlationRemoval.P*2 + S + d + m * D]);
                if (m < 2 || d < D - 1){
                    printf(" ");
                }
            }
        }
        printf("\n");
        
        
//        for (int p = 0; p < trk.correlationRemoval.P; p++){
//            for (int k = 0; k < 2; k++) { // intercept, and then slope
//                for (int j = 0; j < 2; j++) { // - and + sigmas
//                    for (int i = 0; i < 3; i++) { // 1, 2 and 3 sigmas
//                        printf("%.3e ", trk.results.bestFit_123Sigmas[2*p + k][j][i]);
//                    }
//                    printf("\t");
//                }
//                std::cout << std::endl;
//            }
//        }
        // linear param uncertainties
        
        for (int j = 0; j < 2; j++) { // - and + sigmas
            for (int p = 0; p < trk.correlationRemoval.P; p++){
                printf("%.3e ", trk.results.bestFit_123Sigmas[trk.correlationRemoval.intercept_indices[p]][j][0]);
                printf("%.3e", trk.results.bestFit_123Sigmas[trk.correlationRemoval.slope_indices[p]][j][0]);
                if (p < trk.correlationRemoval.P - 1){
                    printf(" ");
                }
            }
            printf("\n");
        }
        
        for (int j = 0 ; j < 8; j++){
            printf("\n");
        }
    }

    
    // ###########################################################################################################################

    
    // GLOBAL FUNCTIONS/TOOLS ####################################################################################################

    // misc tools
    void get_where_boolean(std::vector <bool> &vec, bool value, std::vector <int> &indices){
        // finds which indices of some boolean vec have the inputted value
        indices.clear();
        
        for (int i = 0; i < vec.size(); i++){
            if (vec[i] == value){
                indices.push_back(i);
            }
        }
        return;
    }


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
//        std::cout << "starting timer... \n";

        clock_t t_i = clock();
        return t_i;
    }

    double secElapsed(clock_t t_i) {
        clock_t t_f = clock() - t_i;
        double sec_elapsed = ((float)t_f) / CLOCKS_PER_SEC;

//        printf("%0.3f seconds elapsed \n", sec_elapsed);
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

    std::vector <std::vector <double > > getData(std::string fileName, int dataSize, int numcols) {
        
        std::ifstream readFile;
        readFile.open(fileName);
        
        int columns = numcols;
        int rows = dataSize;
        std::vector <double> innerData(dataSize, 0.0);
        std::vector<std::vector<double>> rawData(numcols, innerData);
        
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

    std::vector <std::vector <double > > getData(std::string fileName, int numcols) {
        
        std::ifstream readFile;
        readFile.open(fileName);
        
        int columns = numcols;
        int rows = 441;
        std::vector <double> innerData(441, 0.0);
        std::vector<std::vector<double>> rawData(numcols, innerData);
        
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
}
