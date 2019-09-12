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
std::vector <double> requestHandler(int fType, std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> allparamsguess, int dataSize, int pivotCheck, int priorsCheck, std::vector <double> priorsParams, std::vector <int> hasPriorsVec){
    
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
    
    TRK trk = TRK(yc, dyc, ddyc, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess);
    
    trk.cpp17MultiThread = true;
    if (pivotCheck == 1){
        trk.findPivotPoints = true;
    }
    
    return result;
}
