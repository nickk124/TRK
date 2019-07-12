#include "pch.h"
#include "TRK.h"
#include <fstream>
#include <sstream>
#include <string>


//model pivot points

double c2p1BH;
double c2p2BH;

double c2p1RV;
double c2p2RV;

double c2pc1;

double yC(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * std::sin(a1 * x);
}

double dyC(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * a1 * std::cos(a1 * x);
}

double ddyC(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return -1.0 * a0 * a1 * a1 * std::sin(a1 * x);
}

double testFunc(std::vector <double> params) {
	double x = params[0];
	double y = params[1];

	return std::pow(x, 2.0) + std::pow(y, 2.0);
}

double testFunc2(std::vector <double> params) {
	double x = params[0];
	double y = params[1];

	return std::sin(x)*std::sin(y)*(std::pow(x, 2.0) + std::pow(y, 2.0));
}

double linearFunc(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 + a1 * x;
}

double dLin(double x, std::vector <double> params) {
	double a1 = params[1];

	return a1;
}

double ddLin(double x, std::vector <double> params) {

	return 0.0;
}

double bar = 0;

double function_quadratic(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];
	double a2 = params[2];

	return a0 + a1 * (x - bar) + a2 * std::pow((x - bar), 2.0);
}

double function_cubic(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];
	double a2 = params[2];
	double a3 = params[3];

	return a0 + a1 * (x - bar) + a2 * std::pow((x - bar), 2.0) + a3 * std::pow((x - bar), 3.0);
}

double function_powerlaw(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * std::pow((x / std::pow(10, bar)), a1);
}

double function_exponential(double x, std::vector <double> params) {
	double a0 = params[0];
	double a1 = params[1];

	return a0 * std::pow(10, a1*(x - bar));
}

//c1 vs c2

double c1c2(double c2, std::vector <double> params) {
	double bc1 = params[0];
	double tc1 = params[1];


	return bc1 + std::tan(tc1)*(c2 - c2pc1);
}

double dc1c2(double c2, std::vector <double> params) {
	double tc1 = params[1];

	return std::tan(tc1);
}

double ddc1c2(double c2, std::vector <double> params) {
	return 0.0;
}
//BH vs c2

double bhc2(double c2, std::vector <double> params) {
	double b1BH = params[0];
	double theta1BH = params[1];

	double b2BH = params[2];
	double theta2BH = params[3];

	return -std::log(std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH)));
}

double dbhc2(double c2, std::vector <double> params) {
	double b1BH = params[0];
	double theta1BH = params[1];

	double b2BH = params[2];
	double theta2BH = params[3];

	double top = -std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH))*std::tan(theta1BH) - std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH))*std::tan(theta2BH);
	double bottom = std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH));

	return -top/bottom;
}

double ddbhc2(double c2, std::vector <double> params) {
	double b1BH = params[0];
	double theta1BH = params[1];

	double b2BH = params[2];
	double theta2BH = params[3];

	double top1 = std::pow(-std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH))*std::tan(theta1BH) - std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH))*std::tan(theta2BH), 2.0);
	double bottom1 = std::pow(std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH)), 2.0);

	double top2 = std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH))*std::pow(std::tan(theta1BH), 2.0) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH))*std::pow(std::tan(theta2BH), 2.0);
	double bottom2 = std::exp(-b1BH - std::tan(theta1BH)*(c2 - c2p1BH)) + std::exp(-b2BH - std::tan(theta2BH)*(c2 - c2p2BH));

	return top1 / bottom1 - top2 / bottom2;
}

//RV vs c2

double rvc2(double c2, std::vector <double> params) {
	double b1RV = params[0];
	double theta1RV = params[1];

	double b2RV = params[2];
	double theta2RV = params[3];

	return std::log(std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV)));
}

double drvc2(double c2, std::vector <double> params) {
	double b1RV = params[0];
	double theta1RV = params[1];

	double b2RV = params[2];
	double theta2RV = params[3];

	double top = std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV))*std::tan(theta1RV) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV))*std::tan(theta2RV);
	double bottom = std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV));

	return top / bottom;
}

double ddrvc2(double c2, std::vector <double> params) {
	double b1RV = params[0];
	double theta1RV = params[1];

	double b2RV = params[2];
	double theta2RV = params[3];

	double top1 = std::pow(std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV))*std::tan(theta1RV) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV))*std::tan(theta2RV), 2.0);
	double bottom1 = std::pow(std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV)), 2.0);

	double top2 = std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV))*std::pow(std::tan(theta1RV), 2.0) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV))*std::pow(std::tan(theta2RV), 2.0);
	double bottom2 = std::exp(b1RV + std::tan(theta1RV)*(c2 - c2p1RV)) + std::exp(b2RV + std::tan(theta2RV)*(c2 - c2p2RV));

	return -(top1 / bottom1) + top2 / bottom2;
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

int main()
{
	/*
	std::vector <std::vector <double> > data = getData("c1c2_data.csv");

	std::vector <double> x, y, sx, sy, w;

	for (int i = 0; i < data[0].size(); i++) {
		x.push_back(data[0][i]);
		sx.push_back(data[1][i]);
		y.push_back(data[2][i]);
		sy.push_back(data[3][i]);
		w.push_back(data[4][i]);
	}

	std::vector <double> params_guess = { 2.5, -3.5 };
	double slopx_guess = 1.0;
	double slopy_guess = 1.0;

	TRK TRKtest = TRK(linearFunc, dLin, ddLin, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess);

	std::vector <double> allparams_guess = params_guess;

	allparams_guess.push_back(slopx_guess);
	allparams_guess.push_back(slopy_guess);


	typedef double (TRK::*TRKMemFn)(std::vector <double> allparams); //here, TRKMemFn is the name of the type. A pointer of this type points to a member of a TRK object that has those specific input and output

	TRKMemFn p = &TRK::modifiedChiSquared;

	TRKtest.s = 0.25767;

	std::vector <double> fit = TRKtest.downhillSimplex(p, allparams_guess);
	*/

	/*
	std::vector <double> x, y, sx(9, 0.1), sy(9, 0.2), w(9, 1.0);

	x = { 1.7, 2.1, 2.56, 5.65, 4.83, 6.56, 7.71, 8.18, 7.74 };
	y = { 2.64, 3.59, 4.58, 4.53, 4.47, 7.28, 8.57, 8.74, 9.57};



	std::vector <double> params_guess = {0.8, 1.1 };
	double slopx_guess = 1.0;
	double slopy_guess = 1.0;
	*/

	std::string filename;

	//filename = "bhc2_data.csv";
	filename = "rvc2_data.csv";
	//filename = "c1c2_data.csv";																									//**********
	std::vector <std::vector <double> > data = getData(filename, 441);

	//filename = "simplelinear_data.csv";
	//std::vector <std::vector <double> > data = getData(filename, 9);

	std::vector <double> x, y, sx, sy, w;

	for (int i = 0; i < data[0].size(); i++) {
		x.push_back(data[0][i]);
		sx.push_back(data[1][i]);
		y.push_back(data[2][i]);
		sy.push_back(data[3][i]);
		w.push_back(data[4][i]);
	}

	//printf("filename: %s \n", filename);

	/*

	//c1vc2 testing

	std::vector <double> params_guess = { 2.5, -3.5 };
	double slopx_guess = 1.0;
	double slopy_guess = 1.0;


	TRK TRKtest = TRK(linearFunc, dLin, ddLin, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess);

	std::vector <double> allparams_guess = params_guess;

	allparams_guess.push_back(slopx_guess);
	allparams_guess.push_back(slopy_guess);

	//std::vector <double> ap_check = { 2.86933, -3.55862, 0.296807, -0.00100249 };

	//TRKtest.s = 1.0;

	//double r = TRKtest.modifiedChiSquared(ap_check);
	*/

	typedef double (TRK::*TRKMemFn)(std::vector <double> params); //here, TRKMemFn is the name of the type. A pointer of this type points to a member of a TRK object that has those specific input and output

	TRKMemFn p = &TRK::modifiedChiSquared;

	//TRKtest.s = 3.22e-11;

	//pivot points for smoothly broken line models

	c2p1BH = -0.01413;
	c2p2BH = 1.4087;

	//c2p1BH = 0.2;
	//c2p2BH = 1.7;

	c2p1RV = -0.0708;
	c2p2RV = 1.4953;

	c2pc1 = 1.2403;

	//c2pc1 = 1.1;

	//
	//std::vector <double> testS = { 0.1, 0.5, 0.75, 1.0, 1.5, 2.0 };


	std::vector <std::vector <double> > bhc2prior_params = { {NAN, NAN}, {PI, 1.5*PI}, {NAN, NAN}, {-0.5*PI, 0} };
	std::vector <std::vector <double> > rvc2prior_params = { {NAN, NAN}, {PI/2.0, PI}, {NAN, NAN}, {-0.5*PI, 0} };

	//std::vector <std::vector <double> > testlinprior_params = { {NAN, NAN}, {NAN, NAN}, {0.0, NAN}, {0.0, NAN} }; //including slop for testing

	Priors bhc2Priors = Priors(CONSTRAINED, bhc2prior_params);
	Priors rvc2Priors = Priors(CONSTRAINED, rvc2prior_params);
	//Priors testlinPriors = Priors(CONSTRAINED, testlinprior_params);

	std::vector <double> params_guess = { 5.0, 4.6, 2.5, -0.3 };	//rvc2	
	//std::vector <double> params_guess = { 4.0, 4.6, 2.0, -1.1 };				//bhc2																					//**********
	//std::vector <double> params_guess = { -1.5038, toRad(106.953) };    //c1c2
	//std::vector <double> params_guess = { 0.8, 1.0 };     //test lin



	double slopx_guess = 0.3;																									//**********
	double slopy_guess = 0.3;

	//std::vector <double> testsigma_guess = { 0.212500, 0.106250}; //testlin
	std::vector <double> testsigma_guess = { 0.01, 0.005, 0.01, 0.01};  //bhc2/rvc2
	//std::vector <double> testsigma_guess = { 0.01, 0.005};  //c1c2

	//double testslop_x_sigma_guess = 0.006250; //testlin
	//double testslop_y_sigma_guess = 0.002500;

	double testslop_x_sigma_guess = 0.0005; //bhc2/rvc2
	double testslop_y_sigma_guess = 0.0005;



	//TRK TRKtest = TRK(linearFunc, dLin, ddLin, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess, testsigma_guess, testslop_x_sigma_guess, testslop_y_sigma_guess);
	//TRK TRKtest = TRK(c1c2, dc1c2, ddc1c2, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess, testsigma_guess, testslop_x_sigma_guess, testslop_y_sigma_guess);															//**********
	//TRK TRKtest = TRK(bhc2, dbhc2, ddbhc2, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess, testsigma_guess, testslop_x_sigma_guess, testslop_y_sigma_guess, bhc2Priors);
	TRK TRKtest = TRK(rvc2, drvc2, ddrvc2, x, y, w, sx, sy, params_guess, slopx_guess, slopy_guess, testsigma_guess, testslop_x_sigma_guess, testslop_y_sigma_guess, rvc2Priors);

	std::vector <double> allparams_guess = params_guess;

	allparams_guess.push_back(slopx_guess);
	allparams_guess.push_back(slopy_guess);

	//TRKtest.s = 0.418573;     //optimum for test lin
	//TRKtest.s = 0.259519;		//c1c2
	//TRKtest.s = 0.2745;			//bhc2
	TRKtest.s = 0.37568; //rvc2

	//TRKtest.s = 0.228531;

	std::vector <double> fit = TRKtest.downhillSimplex(p, allparams_guess);
	std::cout << "s = " << TRKtest.s << ":    ";
	for (int j = 0; j < allparams_guess.size(); j++) {
		std::cout << fit[j] << "   ";
	}
	std::cout << std::endl;

	clock_t time = startTimer();

	printf("\n scale of %f \n", TRKtest.s);

	
	//int test_R = 200000;
	//int test_burn = 10000;

	//std::vector <std::vector <double > > test_drawn = TRKtest.methastPosterior(test_R, test_burn, testsigma_guess);

	//TRKtest.R = 200000;

	TRKtest.calculateUncertainties();

	printf("Uncertainties: (+ 1 2 3, - 1 2 3): \n");
	for (int k = 0; k < params_guess.size(); k++) {
		for (int j = 0; j < 2; j++) {
			for (int i = 0; i < 3; i++) {
				printf("%f ", TRKtest.results.bestFit_123Sigmas[k][j][i]);
			}
			printf("\t");
		}
		std::cout << std::endl;
	}

	printf(" Slop Uncertainties: (+ 1 2 3, - 1 2 3): \n");
	for (int j = 0; j < 2; j++) {
		for (int i = 0; i < 3; i++) {
			printf("%f ", TRKtest.results.slopX_123Sigmas[j][i]);
		}
		printf("\t");
	}
	std::cout << std::endl;
	for (int j = 0; j < 2; j++) {
		for (int i = 0; i < 3; i++) {
			printf("%f ", TRKtest.results.slopY_123Sigmas[j][i]);
		}
		printf("\t");
	}
	std::cout << std::endl;
	
	//TRKtest.optimizeScale();

	/*
	TRKtest.s = 0.05;

	std::vector <double> fit;

	for (int i = 0; i < 30; i++) {

		fit = TRKtest.downhillSimplex(p, allparams_guess);
		std::cout << "s = " << TRKtest.s << ":    " << fit[0] << "   " << fit[1] << "   " << fit[2] << "   " << fit[3] << "   " << fit[4] << "   " << fit[5] << std::endl;

		TRKtest.s += 0.05;
	}

	*/
	/*
	TRKtest.s = 1.0; //initially begin with s = 1


	TRKtest.a = 0.129272; //figures out which scale is a, and which is b, as well as storing the associated best-fit parameters and their associated tangent points for those two extreme scales
	TRKtest.b = 0.979492;


	TRKtest.s = TRKtest.a;
	TRKtest.whichExtrema = slopx;
	std::vector <double> allparams_a = TRKtest.downhillSimplex(&TRK::modifiedChiSquared, allparams_guess);
	TRKtest.whichExtrema = none;

	TRKtest.s = TRKtest.b;
	TRKtest.whichExtrema = slopy;
	std::vector <double> allparams_b = TRKtest.downhillSimplex(&TRK::modifiedChiSquared, allparams_guess);
	TRKtest.whichExtrema = none;

	for (int i = 0; i < allparams_a.size() - 2; i++) {
		TRKtest.params_a.push_back(allparams_a[i]);
		TRKtest.params_b.push_back(allparams_b[i]);
	}
	for (int i = 0; i < TRKtest.params_a.size(); i++) {
		printf("%f \t %f \n", TRKtest.params_a[i], TRKtest.params_b[i]);
	}

	//determine best s1 (new s) to satistfy R2TRKp(a,s) = R2TRKp(s,b)


	TRKtest.x_t_a = TRKtest.x_t_slopx;
	TRKtest.x_t_b = TRKtest.x_t_slopy;

	std::cout << "finding optimum scale" << std::endl;

	double s0 = TRKtest.optimize_s0_R2();

	printf("s0 = %f \n", s0);

	double sfinal = TRKtest.iterateR2_OptimumScale(s0);

	printf("sfinal = %f \n", s0);

	*/
	double sec_elapsed = secElapsed(time);

	writeResults(TRKtest, sec_elapsed, filename);
	
	/*

	TRKtest.a = 0.175415;
	TRKtest.b = 0.369385;

	TRKtest.iterative_allparams_guess = allparams_guess;

	TRKtest.whichExtrema = slopx;
	TRKtest.innerSlopX_Simplex({ TRKtest.a }, TRKtest.iterative_allparams_guess);
	TRKtest.whichExtrema = slopy;
	TRKtest.innerSlopY_Simplex({ TRKtest.b }, TRKtest.iterative_allparams_guess);
	TRKtest.whichExtrema = none;

	TRKtest.x_t_a = TRKtest.x_t_slopx;
	TRKtest.x_t_b = TRKtest.x_t_slopy;

	TRKtest.params_a = TRKtest.params_slopx;
	TRKtest.params_b = TRKtest.params_slopy;

	TRKtest.optimize_s_R2();
	*/

	return 0;
}