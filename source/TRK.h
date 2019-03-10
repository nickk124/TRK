#include <vector>
#include <iostream>
#include <cmath>

const double PI = 3.1415926535897932384626434;

class TRK
{
	public:
		std::vector <double> tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg);

		double (*yc)(double, std::vector <double>); //pointers to model function and it's first to derivatives WRT x
		double (*dyc)(double, std::vector <double>);
		double (*ddyc)(double, std::vector <double>);

		//constructors
		TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>));

		//default constructor:
		TRK();

	private:
		double newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess);
		double twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1);
		std::vector <double> approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg1);
		std::vector <double> cubicSolver(double A, double B, double C, double D);
};