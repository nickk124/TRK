//#include "pch.h"
//#include "exampleModels.h"
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
////    filename = "c1c2_data.csv";                                                                                                    //**********
//    filename = "simplelinear_data.csv";
//
//    filename = "/Users/nickk124/research/reichart/TRK/TRKrepo/testdata/" + filename;
//    std::vector <std::vector <double> > data = getData(filename, 9);
////    std::vector <std::vector <double> > data = getData(filename, 441);
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
//
//
//    typedef double (TRK::*TRKMemFn)(std::vector <double> params, double s); //here, TRKMemFn is the name of the type. A pointer of this type points to a member of a TRK object that has those specific input and output
//
//    //TRKMemFn p = &TRK::modifiedChiSquared;
//
//
//    // pivot points #############################################################################################################
//
//    c2p1BH = -0.01413;
//    c2p2BH = 1.4087;
//
//    c2p1RV = -0.0708;
//    c2p2RV = 1.4953;
//
////    c2pc1 = 1.2403;
//
//    // priors #############################################################################################################
//
//    std::vector <std::vector <double> > bhc2prior_params = { {NAN, NAN}, {PI, 1.5*PI}, {NAN, NAN}, {-0.5*PI, 0} };
//    std::vector <std::vector <double> > rvc2prior_params = { {NAN, NAN}, {PI/2.0, PI}, {NAN, NAN}, {-0.5*PI, 0} };
//    //std::vector <std::vector <double> > testlinprior_params = { {NAN, NAN}, {NAN, NAN}, {0.0, NAN}, {0.0, NAN} }; //including slop for testing
//
//    Priors bhc2Priors = Priors(CONSTRAINED, bhc2prior_params);
//    Priors rvc2Priors = Priors(CONSTRAINED, rvc2prior_params);
//    //Priors testlinPriors = Priors(CONSTRAINED, testlinprior_params);
//
//    // guesses #############################################################################################################
//
//    //std::vector <double> params_guess = { 5.0, 1.7, 2.5, -0.3 };    //rvc2
//    //std::vector <double> params_guess = { 4.0, 4.6, 2.0, -1.1 };                //bhc2                                                                                    //**********
////    std::vector <double> params_guess = { -1.5038, std::tan(toRad(106.953)) };    //c1c2
//    std::vector <double> params_guess = { 0.89, 1.04};     //test lin
//
//    double slopx_guess = 0.3;                                                                                                    //**********
//    double slopy_guess = 0.3;
//
////    std::vector <double> testsigma_guess = { 0.26, 0.13}; //testlin
//    //std::vector <double> testsigma_guess = { 0.01, 0.005, 0.01, 0.01};  //bhc2/rvc2
////    std::vector <double> testsigma_guess = { 0.065, 0.000081};  //c1c2
//
////    double testslop_x_sigma_guess = 0.07700; //testlin
////    double testslop_y_sigma_guess = 0.01500;
//
////    double testslop_x_sigma_guess = 0.0005; //bhc2/rvc2
////    double testslop_y_sigma_guess = 0.0005;
//
////    double testslop_x_sigma_guess = 0.00025; //c1c2
////    double testslop_y_sigma_guess = 0.000813;
//
//
//    std::vector <double> allparams_guess = params_guess;
//
//    allparams_guess.push_back(slopx_guess);
//    allparams_guess.push_back(slopy_guess);
//
//
//    // constructors #############################################################################################################
//
//    TRK TRKtest = TRK(linear, dLinear, ddLinear, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess);
////    TRK TRKtest = TRK(c1c2, dc1c2, ddc1c2, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess);                                                            //**********
//    //TRK TRKtest = TRK(bhc2, dbhc2, ddbhc2, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess, testsigma_guess, testslop_x_sigma_guess, testslop_y_sigma_guess, bhc2Priors);
//    //TRK TRKtest = TRK(rvc2, drvc2, ddrvc2, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess, testsigma_guess, testslop_x_sigma_guess, testslop_y_sigma_guess, rvc2Priors);
//    //TRK TRKtest = TRK(rvc2, drvc2, ddrvc2, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess, testsigma_guess, testslop_x_sigma_guess, testslop_y_sigma_guess);
//
//    //TRKtest.s = 0.257858;     //optimum for test lin
//
//
//    clock_t time = startTimer();
//
//    // test fitness #############################################################################################################
//
//    //TRKtest.s = 1.0;
//    //std::vector <double> adamParams = { 102.965701, toRad(90.120709), 7.425779, toRad(-53.732416), 0.659038, 0.000003 };
//    //double fitness = TRKtest.modifiedChiSquared(adamParams);
//
//    // test simplex #############################################################################################################
//
//    //TRKtest.s = 1.0;
//    /*std::vector <double> fit = TRKtest.downhillSimplex(p, allparams_guess);
//    std::cout << "s = " << TRKtest.s << ":    ";
//    for (int j = 0; j < allparams_guess.size(); j++) {
//        std::cout << fit[j] << "   ";
//    }
//    std::cout << std::endl; */
//
//    // test tangent finder #############################################################################################################
//
//    //TRKtest.s = 1.0;
//    //std::vector <double> x_tn_vec = TRKtest.tangentsFinder({ 2.444524,  4.586552,  2.258013, -1.134365 }, 0.398100, 3.505000, 0.028963, 0.275889, 0.398100);
//
//    // test MCMC/uncertainty #############################################################################################################
//
//    //TRKtest.R = 1000000;
//    //int test_burn = 10000;
//
//    //std::vector <std::vector <double > > test_drawn = TRKtest.methastPosterior(test_R, test_burn, testsigma_guess);
//
//    //TRKtest.calculateUncertainties();
//
//    // test optimum scale finding (last two loops) #############################################################################################################
//    //TRKtest.optimizeScale();
//
//    // test optimum scale finding (last two loops) #############################################################################################################
//
//    //TRKtest.a = 0.129883;
//    /*TRKtest.a = 0.129;
//    TRKtest.b = 0.977539;
//
//    TRKtest.iterative_allparams_guess = allparams_guess;
//
//    TRKtest.whichExtrema = slopx;
//    TRKtest.innerSlopX_Simplex({ TRKtest.a }, TRKtest.iterative_allparams_guess);
//    TRKtest.whichExtrema = slopy;
//    TRKtest.innerSlopY_Simplex({ TRKtest.b }, TRKtest.iterative_allparams_guess);
//    TRKtest.whichExtrema = none;
//
//    TRKtest.x_t_a = TRKtest.x_t_slopx;
//    TRKtest.x_t_b = TRKtest.x_t_slopy;
//
//    TRKtest.params_a = TRKtest.params_slopx;
//    TRKtest.params_b = TRKtest.params_slopy;
//
//    double s0 = TRKtest.optimize_s0_R2();
//
//    TRKtest.s = s0;
//
//    double s_final = TRKtest.iterateR2_OptimumScale(s0);
//
//    TRKtest.s = s_final;
//
//    std::cout << "optimum s = " << TRKtest.s << std::endl;
//
//    TRKtest.results.bestFitParams.clear();
//
//    for (int j = 0; j < TRKtest.M; j++) {
//        TRKtest.results.bestFitParams.push_back(TRKtest.allparams_s[j]);
//    }
//
//    TRKtest.results.slop_x = TRKtest.allparams_s[TRKtest.M];
//    TRKtest.results.slop_y = TRKtest.allparams_s[TRKtest.M + 1];
//
//    TRKtest.results.optimumScale = TRKtest.s;*/
//
//
//
//    // pivot point testing ##############################
//
//    /*std::vector <std::vector <double> > test = {{0,1}, {1,2}, {2,3}, {3,4}, {4,5} };
//
//    TRKtest.getCombos(test, 2, 0);
//
//    TRKtest.NDcombos;*/
//
//    // key algorithm testing #############################################################################################################
//
////    std::vector <std::vector <double> > test = { {0,1}, {1,2}, {2,3}, {3,4}, {4,5} };
////
////    TRKtest.getCombos(test, 2, 0);
////
////    //TRKtest.NDcombos;
////
////    std::vector < std::vector < std::vector <double> > > drawnCombos;
////
////    random_unique(TRKtest.NDcombos.begin(), TRKtest.NDcombos.end(), 3);
//
//    //TRKtest.openMPMultiThread = true;
//
//    TRKtest.findPivotPoints = true;
//    TRKtest.writePivots = true;
//    TRKtest.outputDistributionToFile = true;
//
//    TRKtest.linearizedIntercept = linearIntercept;
//    TRKtest.linearizedSlope = linearSlope;
//
//    TRKtest.performTRKFit();
//
//    printf("Optimum scale: %f \n", TRKtest.results.optimumScale);
//    printf("Minimum scale: %f \n", TRKtest.results.minimumScale);
//    printf("Maximum scale: %f \n", TRKtest.results.maximumScale);
//
//    printf("Fitted parameters (including slop): ");
//    for (int k = 0; k < params_guess.size(); k++) {
//        printf("%f ", TRKtest.results.bestFitParams[k]);
//    }
//    printf(" %f %f", TRKtest.results.slop_x, TRKtest.results.slop_y);
//    std::cout << std::endl;
//
//    printf("Uncertainties: (- 1 2 3, + 1 2 3): \n");
//    for (int k = 0; k < params_guess.size(); k++) { //kth param
//        for (int j = 0; j < 2; j++) { // - and + sigmas
//            for (int i = 0; i < 3; i++) { // 1, 2 and 3 sigmas
//                printf("%f ", TRKtest.results.bestFit_123Sigmas[k][j][i]);
//            }
//            printf("\t");
//        }
//        std::cout << std::endl;
//    }
//
//    printf("Slop Uncertainties: (- 1 2 3, + 1 2 3): \n");
//    for (int j = 0; j < 2; j++) {
//        for (int i = 0; i < 3; i++) {
//            printf("%f ", TRKtest.results.slopX_123Sigmas[j][i]);
//        }
//        printf("\t");
//    }
//    std::cout << std::endl;
//    for (int j = 0; j < 2; j++) {
//        for (int i = 0; i < 3; i++) {
//            printf("%f ", TRKtest.results.slopY_123Sigmas[j][i]);
//        }
//        printf("\t");
//    }
//    std::cout << std::endl;
//
//    double sec_elapsed = secElapsed(time);
//
//    writeResults(TRKtest, sec_elapsed, filename);
//
//    return 0;
//}
