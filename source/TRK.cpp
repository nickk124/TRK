#include "pch.h"
#include "TRK.h"
#include <algorithm>

// CONSTRUCTORS

TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess) {
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
}

//default
TRK::TRK() {

}

// OTHER ALGORITHMS
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

// FITTING TOOLS
double TRK::newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess) {
	double x0 = xguess;
	int itercount = 0;

	//initial iteration
	double f = (yc(x0, params) - y_n) * dyc(x0, params) * Sig_xn2 + (x0 - x_n) * Sig_yn2;  //function we're finding zero of
	double df = (std::pow(dyc(x0, params), 2.0) + (yc(x0, params) - y_n)*ddyc(x0, params)) * Sig_xn2 + Sig_yn2; //derivative of above

	double x1 = x0 - f / df;

	double tol = 1e-9;

	while (std::abs(x1 - x0) > tol) {
		x0 = x1;

		f = (yc(x0, params) - y_n) * dyc(x0, params) * Sig_xn2 + (x0 - x_n) * Sig_yn2;  //function we're finding zero of
		df = (std::pow(dyc(x0, params), 2.0) + (yc(x0, params) - y_n)*ddyc(x0, params)) * Sig_xn2 + Sig_yn2; //derivative of above

		x1 = x0 - f / df;
		std::cout << x1 << std::endl;

		itercount += 1;
	}

	std::cout << itercount << " iterations." << std::endl;
	return x1;
}

double TRK::twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1)
{
	double tol = 1e-9;
	double xkm1 = xguess;
	//double xk = xguess + xguess / 100.0;
	double xk = xguessp1;

	double xkp1;
	double r;
	double ykm1;
	double yk;
	double dyk;

	int itercount = 0;

	while (std::abs(xk - xkm1) > tol) {
		ykm1 = (yc(xkm1, params) - y_n) * dyc(xkm1, params) * Sig_xn2 + (xkm1 - x_n) * Sig_yn2;
		yk = (yc(xk, params) - y_n) * dyc(xk, params) * Sig_xn2 + (xk - x_n) * Sig_yn2;  //function we're finding zero of
		dyk = (std::pow(dyc(xk, params), 2.0) + (yc(xk, params) - y_n)*ddyc(xk, params)) * Sig_xn2 + Sig_yn2; //derivative of above

		r = 1.0 - ((yk / ykm1) * (((yk - ykm1) / (xk - xkm1)) / dyk));

		xkp1 = (1.0 - 1.0 / r)*xkm1 + xk / r;
		std::cout << xkp1 << std::endl;

		xkm1 = xk;
		xk = xkp1;

		itercount += 1;
	}

	std::cout << itercount << " iterations." << std::endl;
	return xkp1;
}

std::vector <double> TRK::cubicSolver(double A, double B, double C, double D) {
	//cubic solver for three real and distinct roots
	std::vector <double> roots;

	double a1 = B / A;
	double a2 = C / A;
	double a3 = D / A;

	double Q = (a1 * a1 - 3.0 * a2) / 9.0;
	double R = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
	double Qc = std::pow(Q, 3.0);
	double d = Qc - std::pow(R, 2.0);

	double theta = std::acos(R / sqrt(Qc));

	roots.push_back( -2 * std::sqrt(Q) * std::cos(theta / 3) - a1 / 3 );
	roots.push_back(-2 * std::sqrt(Q) * std::cos((theta + 2 * PI) / 3) - a1 / 3 );
	roots.push_back(-2 * std::sqrt(Q) * std::cos((theta + 4 * PI) / 3) - a1 / 3 );

	return roots;
}

// STATISTICS
double TRK::singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn) {
	double m_tn = dyc(x_tn, params);
	double y_tn = yc(x_tn, params);
	
	return std::pow(y_n - y_tn - m_tn * (x_n - x_tn), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) - std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(Sig_yn2, 2.0)));
}

double TRK::modifiedChiSquared(std::vector <double> allparams)
{
	int N = y.size();
	double sum1 = 0.0;
	double sum2 = 0.0;

	double slop_x = allparams[-2];
	double slop_y = allparams[-1];

	std::vector <double> params;

	for (int i = 0; i < M; i++) {
		params.push_back(allparams[i]);
	}

	for (int n = 0; n < N; n++) {
		double Sig_xn2 = std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0);
		double Sig_yn2 = std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0);

		std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], Sig_xn2, Sig_yn2, x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points
		double x_t = findBestTangent(params, x[n], y[n], Sig_xn2, Sig_yn2, x_tn_vec);

		double m_tn = dyc(x_t, params);
		double y_tn = yc(x_t, params);

		sum1 += std::pow(y[n] - y_tn - m_tn * (x[n] - x_t), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2);
		sum2 += std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(Sig_yn2, 2.0)));
	}
	return sum1 - sum2;
}

// TANGENT-POINT ALGORITHMS
std::vector <double> TRK::approxQuadraticRoots(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xr1) {
	//using board derivation notation:
	double b = yc(xr1, params) - dyc(xr1, params) * xr1 + (ddyc(xr1, params) / 2.0) * std::pow(xr1, 2.0); //coefficients of quadratic approximation from taylor expansion
	double m = dyc(xr1, params) - ddyc(xr1, params) * xr1;
	double a = ddyc(xr1, params) / 2.0;


	//DIFFERENT FROM NOTATION OF BOARD DERIVATION!
	double A = 2.0 * std::pow(a, 2.0) * Sig_xn2;											// coef of x^3
	double B = 3.0 * a * m * Sig_xn2;														//coef of x^2
	double C = (2.0 * a * (b - y_n) + std::pow(m, 2.0)) * Sig_xn2 + Sig_yn2;	//coef of x
	double D = m * (b - y_n) * Sig_xn2 - x_n * Sig_yn2;							//coef of 1

	double discriminant = 18.0*A*B*C*D - 4.0*std::pow(B, 3.0)*D + std::pow(B, 2.0)*std::pow(C, 2.0) - 4.0*A*std::pow(C, 3.0) - 27.0*std::pow(A, 2.0)*std::pow(D, 2.0);

	std::vector <double> roots;

	if (discriminant > 0) {
		roots = cubicSolver(A, B, C, D);
		return roots;
	}
	//returns no extra roots (empty vector) if the other two roots are 
	return roots;
}

std::vector <double> TRK::tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg) {
	
	std::vector <double> result;

	double xg1 = xg;

	while (true) {
		result.clear();

		double xr1 = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg1, xg1 + std::sqrt(Sig_xn2));
		double xg2, xg3;

		//Quadratic Approximation

		std::vector <double> allRoots = approxQuadraticRoots(params, x_n, y_n, Sig_xn2, Sig_yn2, xr1); //get approx. roots from quadratic taylor approximation
		std::vector <double> extraRoots;

		if (allRoots.size() == 3) { //either add the two other roots, or no more roots (depending on discriminant of cubic equation)
			//grab to new roots 
			for (int i = 0; i < 3; i++) {
				double root = allRoots[i];
				if (std::abs(root - xr1) >= 1e-9) { //checks if roots are (numerically) identical or not
					extraRoots.push_back(root);
				}
			}
		}

		if (extraRoots.size() == 2) { //if have 2 additional, real roots
			double xr2, xr3;
			xg2 = extraRoots[0];
			xg3 = extraRoots[1];

			if (xg2 < xr1 && xr1 < xg3) {
				xr2 = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg2, xg2 - std::sqrt(Sig_xn2));
				xr3 = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg3, xg3 + std::sqrt(Sig_xn2));

				result.push_back(xr2);
				result.push_back(xr1);
				result.push_back(xr3);

				break;
			}
			else {
				std::vector <double> rootVec = { xr1, xg2, xg3 };
				std::sort(rootVec.begin(), rootVec.end());

				xg1 = rootVec[1];
			}
		} else if (extraRoots.size() == 0) { //if only have one root (the one initially found)
			result.push_back(xr1);
			break;
		}
	}
	return result;
}

double TRK::findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec) {
	std::vector <double> posts;
	int minindex;

	for (int i = 0; i < 3; i++) {
		posts.push_back(singlePointLnL(params, x_n, y_n, Sig_xn2, Sig_yn2, x_tn_vec[i]));
	}

	std::vector<double>::iterator result = std::min_element(std::begin(posts), std::end(posts));
	minindex = std::distance(std::begin(posts), result);

	return x_tn_vec[minindex];
}

// SCALE OPTIMIZATION ALGORITHMS

std::vector <double> TRK::findCentroid(std::vector <std::vector <double> > vertices) {
	int n = M + 2;

	std::vector <double> centroid(n, 0);

	for (int i = 0; i < n; i++) {// for each param plus slop
		double sum = 0.0;
		for (int j = 0; j < n; j++) { // for each vertex
			sum += vertices[j][i];
		}
		centroid[i] = sum/n;
	}

	return centroid;
}

std::vector <double> TRK::downhillSimplex() {

	int n = params_guess.size() + 2; //number of model parameters plus two slop parameters

	double rho = 1.0;
	double chi = 2.0;
	double gamma = 0.5;
	double sigma = 0.5;

	// simplex initialization

	std::vector <double> init_point = params_guess;

	init_point.push_back(slop_x_guess);
	init_point.push_back(slop_y_guess);

	std::vector <std::vector <double> > vertices = { init_point };

	for (int j = 0; j < n; j++) { //for each simplex node
		init_point.clear();
		for (int i = 0; i < M; i++) { //for each parameter (non slop)
			init_point.push_back(params_guess[i] + params_guess[i]); //add initial "step size"
		}
		init_point.push_back(slop_x_guess + slop_x_guess);
		init_point.push_back(slop_y_guess + slop_y_guess);

		vertices.push_back(init_point);
	}

	std::vector <double> result;

	while (true) {
		// order

		std::vector <int> orderedindices;
		std::vector <double> unorderedEvals; // ( f(x_1), f(x_2), ... f(x_n+1)
		for (int i = 0; i < n + 1; i++) {
			unorderedEvals.push_back(modifiedChiSquared(vertices[i]));
		}
		orderedindices = getSortedIndices(unorderedEvals);

		std::vector <std::vector <double> > orderedvertices = { vertices[orderedindices[0]] };
		for (int i = 1; i < n + 1; i++) {
			orderedvertices.push_back(vertices[orderedindices[i]]);
		}

		vertices = orderedvertices;

		// reflect
		std::vector <double> refpoint;
		std::vector <double> centroid = findCentroid(vertices);

		for (int i = 0; i < n; i++) {
			refpoint.push_back(centroid[i] + rho*(centroid[i] - vertices[n + 1][i]));
		}

		double fr = modifiedChiSquared(refpoint);
		double f1 = modifiedChiSquared(vertices[0]);
		double fn = modifiedChiSquared(vertices[n-1]);

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

			double fe = modifiedChiSquared(exppoint);

			if (fe < fr) {
				result = exppoint;
				break;
			}
			else if (fe >= fr) {
				result = refpoint;
				break;
			}

		//contract
		}
	}

	return result;
}