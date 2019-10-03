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
    
//    trk.cpp17MultiThread = true;
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
    
    //printf("%i\n",(int) trk.results.bestFit_123Sigmas.size());
    
    for (int j = 0; j < M; j++){
        if (findUncertainties == 1){
            std::vector <double> sigmas = trk.results.bestFit_123Sigmas[j][0];
            sigmas = concat(sigmas, trk.results.bestFit_123Sigmas[j][1]);
            result = concat(result, sigmas);
        } else if (findUncertainties == 0){
            std::vector <double> sigmas(6, 0.0);
            result = concat(result, sigmas);
        }
    }
    
    if (opScale == 0){
        trk.results.optimumScale = 1.0;
        trk.results.minimumScale = 0.0;
        trk.results.maximumScale = 0.0;
    }
    
    result = concat(result, {trk.results.optimumScale, trk.results.minimumScale, trk.results.maximumScale});
    
    result.push_back(trk.results.pivot);
    
    for (int j = 0; j < M; j++){
        if (findUncertainties == 1){
            result.push_back((double) trk.results.paramDistributionHistograms[j][0].size()); // paramDistributionHistograms = { {hist1, edges1}, {hist2, edges2}, ... }
        } else if (findUncertainties == 0){
            result.push_back(1.0); //each fake hist has bin size of 1.
        }
    }
    
    for (int j = 0; j < M; j++){
        if (findUncertainties == 1){
            result = concat(result, trk.results.paramDistributionHistograms[j][0]);
            result = concat(result, trk.results.paramDistributionHistograms[j][1]);
        } else if (findUncertainties == 0){
            std::vector <double> fakeHist = {1.0};
            std::vector <double> fakeEdge = {0.0, 1.0};
            result = concat(result, fakeHist);
            result = concat(result, fakeEdge);
        }
        
    }
    
    return result;
}

//#include <fstream>
//#include <sstream>
//#include <string>
//
//
////model pivot points
//
//double c2p1BH;
//double c2p2BH;
//
//double c2p1RV;
//double c2p2RV;
//
//double c2pc1;
//
//double yC(double x, std::vector <double> params) {
//    double a0 = params[0];
//    double a1 = params[1];
//
//    return a0 * std::sin(a1 * x);
//}
//
//double dyC(double x, std::vector <double> params) {
//    double a0 = params[0];
//    double a1 = params[1];
//
//    return a0 * a1 * std::cos(a1 * x);
//}
//
//double ddyC(double x, std::vector <double> params) {
//    double a0 = params[0];
//    double a1 = params[1];
//
//    return -1.0 * a0 * a1 * a1 * std::sin(a1 * x);
//}
//
////c1 vs c2
//
//double c1c2(double c2, std::vector <double> params) {
//    double bc1 = params[0];
//    double mc1 = params[1];
//
//    return bc1 + mc1*(c2 - TRK::pivot);
//}
//
//double dc1c2(double c2, std::vector <double> params) {
//    double mc1 = params[1];
//
//    return mc1;
//}
//
//double ddc1c2(double c2, std::vector <double> params) {
//    return 0.0;
//}
//
//double pivotc1c2(std::vector <double> params1, std::vector <double> params2) {
//    double b1 = params1[0];
//    double t1 = params1[1];
//
//    double b2 = params2[0];
//    double t2 = params2[1];
//
//    return (b2 - b1) / (std::tan(t1) - std::tan(t2));
//}
////BH vs c2
//
//double bhc2(double c2, std::vector <double> params) {
//    double b1BH = params[0];
//    double theta1BH = params[1];
//
//    double b2BH = params[2];
//    double theta2BH = params[3];
//
//    return -std::log(std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH)));
//}
//
//double dbhc2(double c2, std::vector <double> params) {
//    double b1BH = params[0];
//    double theta1BH = params[1];
//
//    double b2BH = params[2];
//    double theta2BH = params[3];
//
//    double top = -std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH))*std::tan(theta1BH) - std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH))*std::tan(theta2BH);
//    double bottom = std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH));
//
//    return -top/bottom;
//}
//
//double ddbhc2(double c2, std::vector <double> params) {
//    double b1BH = params[0];
//    double theta1BH = params[1];
//
//    double b2BH = params[2];
//    double theta2BH = params[3];
//
//    double top1 = std::pow(-std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH))*std::tan(theta1BH) - std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH))*std::tan(theta2BH), 2.0);
//    double bottom1 = std::pow(std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH)), 2.0);
//
//    double top2 = std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH))*std::pow(std::tan(theta1BH), 2.0) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH))*std::pow(std::tan(theta2BH), 2.0);
//    double bottom2 = std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH));
//
//    return top1 / bottom1 - top2 / bottom2;
//}
//
////RV vs c2
//
//double rvc2(double c2, std::vector <double> params) {
//    double b1RV = params[0];
//    double theta1RV = params[1];
//
//    double b2RV = params[2];
//    double theta2RV = params[3];
//
//    return std::log(std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV)));
//}
//
//double drvc2(double c2, std::vector <double> params) {
//    double b1RV = params[0];
//    double theta1RV = params[1];
//
//    double b2RV = params[2];
//    double theta2RV = params[3];
//
//    double top = std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV))*std::tan(theta1RV) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV))*std::tan(theta2RV);
//    double bottom = std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV));
//
//    return top / bottom;
//}
//
//double ddrvc2(double c2, std::vector <double> params) {
//    double b1RV = params[0];
//    double theta1RV = params[1];
//
//    double b2RV = params[2];
//    double theta2RV = params[3];
//
//    double top1 = std::pow(std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV))*std::tan(theta1RV) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV))*std::tan(theta2RV), 2.0);
//    double bottom1 = std::pow(std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV)), 2.0);
//
//    double top2 = std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV))*std::pow(std::tan(theta1RV), 2.0) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV))*std::pow(std::tan(theta2RV), 2.0);
//    double bottom2 = std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV));
//
//    return -(top1 / bottom1) + top2 / bottom2;
//}
//
//int main()
//{
//
//
//    std::string filename;
//
//    //filename = "bhc2_data.csv";
//    //filename = "rvc2_data.csv";
//    //    filename = "c1c2_data.csv";                                                                                                    //**********
//    filename = "simplelinear_data.csv";
//
//    filename = "/Users/nickk124/research/reichart/TRK/TRKrepo/testdata/" + filename;
//
//    std::vector <std::vector <double> > data = getData(filename, 9);
//    //    std::vector <std::vector <double> > data = getData(filename, 441);
//
//    std::vector <double> x, y, sx, sy, w;
//
//    for (int i = 0; i < data[0].size(); i++) {
//        x.push_back(data[0][i]);
//        sx.push_back(data[1][i]);
//        y.push_back(data[2][i]);
//        sy.push_back(data[3][i]);
//        w.push_back(data[4][i]);
//    }
//
//    std::vector <double> allparamsguess = {1.0, 1.0, 0.3, 0.3};
//    std::vector <double> priorsParams;
//    std::vector <int> hasPriorsVec;
//
//    std::vector <double> res;
//
//    res = requestHandler(1, x, y, w, sx, sy, allparamsguess, (int)x.size(), 0, 0, priorsParams, hasPriorsVec, 0, 0, 0.1);
//
//    for (int i = 0; i < (int)res.size(); i++){
//        printf("%f\n", res[i]);
//    }
//
//
//    return 0;
//}

