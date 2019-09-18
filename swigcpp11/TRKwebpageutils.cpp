#include "TRKwebpageutils.h"

Priors getPriors(int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec, int paramCount) {
    int M = paramCount;
    std::vector < std::vector <double> > gaussianParams;
    std::vector < std::vector <double> > paramBounds;
    std::vector <double> temp_params;
    Priors priorsResult;
    
    switch (priorsCheck) {
        case 1: //gaussian
            for (int i = 0; i < M; i++) {
                temp_params.clear();
                
                for (int j = 0; j < 2; j++) {
                    if (hasPriorsVec[i * 2 + j] == 1) {//has a prior
                        temp_params.push_back(priorsParams[i * 2 + j]);
                    }
                    else if (hasPriorsVec[i * 2 + j] == 0) { //no prior
                        temp_params.push_back(NAN);
                    }
                }
                gaussianParams.push_back(temp_params);
            }
            priorsResult = Priors(GAUSSIAN, gaussianParams);
            
            break;
        case 2: //constrained
            for (int i = 0; i < M; i++) {
                temp_params.clear();
                
                for (int j = 0; j < 2; j++) {
                    if (hasPriorsVec[i * 2 + j] == 1) {//has a prior
                        temp_params.push_back(priorsParams[i * 2 + j]);
                    }
                    else if (hasPriorsVec[i * 2 + j] == 0) { //no prior
                        temp_params.push_back(NAN);
                    }
                }
                paramBounds.push_back(temp_params);
            }
            priorsResult = Priors(CONSTRAINED, paramBounds);
            
            break;
        case 3: //mixed
            for (int i = 0; i < M; i++) {
                temp_params.clear();
                
                for (int j = 0; j < 2; j++) {
                    if (hasPriorsVec[i * 2 + j] == 1) {//has a prior
                        temp_params.push_back(priorsParams[i * 2 + j]);
                    }
                    else if (hasPriorsVec[i * 2 + j] == 0) { //no prior
                        temp_params.push_back(NAN);
                    }
                }
                
                gaussianParams.push_back(temp_params);
            }
            
            for (int i = 0; i < M; i++) {
                temp_params.clear();
                
                for (int j = 0; j < 2; j++) {
                    if (hasPriorsVec[2*M + i * 2 + j] == 1) {//has a prior
                        temp_params.push_back(priorsParams[i * 2 + j]);
                    }
                    else if (hasPriorsVec[2*M + i * 2 + j] == 0) { //no prior
                        temp_params.push_back(NAN);
                    }
                }
                paramBounds.push_back(temp_params);
            }
            priorsResult = Priors(MIXED, gaussianParams, paramBounds);
            
            break;
            
        default:
            break;
    }
    
    
    return priorsResult;
}

//if non weighted, in python just use a vector of w's for below.
std::vector <double> requestHandler(int fType, std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> allparamsguess, int dataSize, int pivotCheck, int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec, int opScale, int findUncertainties, double fitScale){
    
    std::vector <double> result;
    //ftype conversion:
//    1: Linear
//    2: Quadratic
//    3: Cubic
//    4: PowerLaw
//    5: Exponential
//    6: Logarithmic
    
    double(*yc)(double, std::vector <double>);
    double(*dyc)(double, std::vector <double>);
    double(*ddyc)(double, std::vector <double>);
    
    switch (fType){
        case 1:
            yc = linear;
            dyc = dLinear;
            ddyc = ddLinear;
            break;
        case 2:
            yc = quadratic;
            dyc = dQuadratic;
            ddyc = ddQuadratic;
            break;
        case 3:
            yc = cubic;
            dyc = dCubic;
            ddyc = ddCubic;
            break;
        case 4:
            yc = powerlaw;
            dyc = dPowerlaw;
            ddyc = ddPowerlaw;
            break;
        case 5:
            yc = exponential;
            dyc = dExponential;
            ddyc = ddExponential;
            break;
        case 6:
            yc = logarithmic;
            dyc = dLogarithmic;
            ddyc = ddLogarithmic;
            break;
        default:
            yc = linear;
            dyc = dLinear;
            ddyc = ddLinear;
            break;
    }
    
    int M = (int)allparamsguess.size() - 2;
    
    std::vector <double> params_guess = slice(allparamsguess, 0, M);
    double slopx_guess = allparamsguess[M];
    double slopy_guess = allparamsguess[M+1];
    
    
    
    TRK trk;
    
    trk.cpp17MultiThread = true;
    if (pivotCheck == 1){
        trk.findPivotPoints = true;
        
        switch (fType){
            case 1:
                trk.linearizedIntercept = linearIntercept;
                trk.linearizedSlope = linearSlope;
                break;
            case 4:
                trk.linearizedIntercept = powerlawIntercept;
                trk.linearizedSlope = powerlawSlope;
                break;
            case 5:
                trk.linearizedIntercept = exponentialIntercept;
                trk.linearizedSlope = exponentialSlope;
                break;
            case 6:
                trk.linearizedIntercept = logarithmicIntercept;
                trk.linearizedSlope = logarithmicSlope;
                break;
            default:
                trk.linearizedIntercept = linearIntercept;
                trk.linearizedSlope = linearSlope;
                break;
        }
        
    }
    
    if (priorsCheck == 1){
        Priors priorsObj = getPriors(priorsCheck, priorsParams, hasPriorsVec, M);
        
        trk = TRK(yc, dyc, ddyc, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess, priorsObj);
    } else {
        trk = TRK(yc, dyc, ddyc, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess);
    }
    
    if (opScale == 1 && findUncertainties == 1){
        trk.performTRKFit();
    } else if (opScale == 0 && findUncertainties == 1){
        trk.performTRKFit(fitScale);
    } else if (opScale == 1 && findUncertainties == 0){
        trk.performSimpleTRKFit();
    } else if (opScale == 0 && findUncertainties == 0){
        trk.performSimpleTRKFit(fitScale);
    }
    
    //    results = {best fit params, slop, - 1 2 3, + 1 2 3 sigmas, s0, a, b, pivot, bincount1, bincount2 ... , hist1, edges1, hist2, edges2 ...
    for (int j = 0; j < M; j++){
        result.push_back(trk.results.bestFitParams[j]);
    }
    
    result.push_back(trk.results.slop_x);
    result.push_back(trk.results.slop_y);
    
    //TRKtest.results.bestFit_123Sigmas[k][j][i]) kth param, jth sign of sigma, ith sigma,
    
    for (int j = 0; j < M; j++){
        std::vector <double> sigmas = trk.results.bestFit_123Sigmas[j][0];
        sigmas = concat(sigmas, trk.results.bestFit_123Sigmas[j][1]);
        result = concat(result, sigmas);
    }
    
    result = concat(result, {trk.results.optimumScale, trk.results.minimumScale, trk.results.maximumScale});
    
    result.push_back(trk.results.pivot);
    
    for (int j = 0; j < M; j++){
        result.push_back((double) trk.results.paramDistributionHistograms[j][0].size()); // paramDistributionHistograms = { {hist1, edges1}, {hist2, edges2}, ... }
    }
    
    for (int j = 0; j < M; j++){
        result = concat(result, trk.results.paramDistributionHistograms[j][0]);
        result = concat(result, trk.results.paramDistributionHistograms[j][1]);
    }
    
    return result;
}
