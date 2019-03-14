#include "pch.h"
#include "TRK.h"
#include <fstream>
#include <sstream>
#include <string>


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

	//TRKtest.s = 3.22e-11;

	//
	std::vector <double> testS = { 0.1, 0.5, 0.75, 1.0, 1.5, 2.0 };


	TRKtest.s = 0.01;

	std::vector <double> fit = TRKtest.downhillSimplex(p, allparams_guess);

	std::cout << "s = " << TRKtest.s << ":    " << fit[0] << "   " << fit[1] << "   " << fit[2] << "   " << fit[3] << std::endl;


	return 0;
}