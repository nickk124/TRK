#include "pch.h"
#include "TRK.h"

// PRIORS

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


// CONSTRUCTORS

TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, std::vector <double> params_sigmas_guess, double slop_x_sigma_guess, double slop_y_sigma_guess) {
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

	this->whichExtrema = none;

	this->s = 1.0;

	getDataWidth();

	this->hasPriors = false;

	this->params_sigmas_guess = params_sigmas_guess; 
	this->slop_x_sigma_guess = slop_x_sigma_guess; 
	this->slop_y_sigma_guess = slop_y_sigma_guess;

	std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

	allparams_sigmas_guess.push_back(slop_x_sigma_guess);
	allparams_sigmas_guess.push_back(slop_y_sigma_guess);

	this->allparams_sigmas_guess = allparams_sigmas_guess;
}
//equal weights/unweighted
TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, std::vector <double> params_sigmas_guess, double slop_x_sigma_guess, double slop_y_sigma_guess) {
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

	this->whichExtrema = none;

	this->s = 1.0;

	getDataWidth();

	this->hasPriors = false;

	this->params_sigmas_guess = params_sigmas_guess;
	this->slop_x_sigma_guess = slop_x_sigma_guess;
	this->slop_y_sigma_guess = slop_y_sigma_guess;

	std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

	allparams_sigmas_guess.push_back(slop_x_sigma_guess);
	allparams_sigmas_guess.push_back(slop_y_sigma_guess);

	this->allparams_sigmas_guess = allparams_sigmas_guess;
}


//priors:
//weighted
TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, std::vector <double> params_sigmas_guess, double slop_x_sigma_guess, double slop_y_sigma_guess, Priors priorsObject) {
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

	this->whichExtrema = none;

	this->s = 1.0;

	getDataWidth();

	this->priorsObject = priorsObject;

	this->hasPriors = true;

	this->params_sigmas_guess = params_sigmas_guess;
	this->slop_x_sigma_guess = slop_x_sigma_guess;
	this->slop_y_sigma_guess = slop_y_sigma_guess;

	std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

	allparams_sigmas_guess.push_back(slop_x_sigma_guess);
	allparams_sigmas_guess.push_back(slop_y_sigma_guess);

	this->allparams_sigmas_guess = allparams_sigmas_guess;
}
//equal weights/unweighted
TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, std::vector <double> params_sigmas_guess, double slop_x_sigma_guess, double slop_y_sigma_guess, Priors priorsObject) {
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

	this->whichExtrema = none;

	this->s = 1.0;

	getDataWidth();

	this->priorsObject = priorsObject;

	this->hasPriors = true;

	this->params_sigmas_guess = params_sigmas_guess;
	this->slop_x_sigma_guess = slop_x_sigma_guess;
	this->slop_y_sigma_guess = slop_y_sigma_guess;

	std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

	allparams_sigmas_guess.push_back(slop_x_sigma_guess);
	allparams_sigmas_guess.push_back(slop_y_sigma_guess);

	this->allparams_sigmas_guess = allparams_sigmas_guess;
}

//default
TRK::TRK() {

}

// OTHER ALGORITHMS AND TOOLS
std::vector <double> TRK::minMax(std::vector <double> vec) {

	// Finding the smallest of all the numbers 
	double min = *std::min_element(std::begin(vec), std::end(vec));
	double max = *std::max_element(std::begin(vec), std::end(vec));

	return { min, max };
}

std::vector <double> TRK::slice(std::vector <double> vec, int begin, int end) {
	std::vector <double> x;
	int len = end - begin;

	for (int i = 0; i < len; i++) {
		x.push_back(vec[begin + i]);
	}

	return x;
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

std::vector <double> TRK::findCentroid(std::vector <std::vector <double> > nvertices) {
	int n = nvertices.size();

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

void TRK::getDataWidth() {
	std::vector <double> bounds = minMax(x);

	datawidth = std::abs(bounds[1] - bounds[0]);

	x_min = bounds[0];
	x_max = bounds[1];
}

double TRK::getAverage(std::vector <double> x) {
	double top = 0.0;
	N = x.size();

	for (int i = 0; i < N; i++) {
		top += x[i];
	}

	return top / N;
}

// FITTING TOOLS
double TRK::newtonRaphson(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess) {
	double x0 = xguess;
	int itercount = 0;

	//initial iteration
	double f = (yc(x0, params) - y_n) * dyc(x0, params) * Sig_xn2 + (x0 - x_n) * Sig_yn2;  //function we're finding zero of
	double df = (std::pow(dyc(x0, params), 2.0) + (yc(x0, params) - y_n)*ddyc(x0, params)) * Sig_xn2 + Sig_yn2; //derivative of above

	double x1 = x0 - f / df;

	double tol = 1e-3;

	while (std::abs(x1 - x0) > tol) {
		x0 = x1;

		f = (yc(x0, params) - y_n) * dyc(x0, params) * Sig_xn2 + (x0 - x_n) * Sig_yn2;  //function we're finding zero of
		df = (std::pow(dyc(x0, params), 2.0) + (yc(x0, params) - y_n)*ddyc(x0, params)) * Sig_xn2 + Sig_yn2; //derivative of above

		x1 = x0 - f / df;

		itercount += 1;
	}

	return x1;
}

double TRK::twoPointNR(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xguess, double xguessp1)
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

		ykm1 = (yc(xkm1, params) - y_n) * dyc(xkm1, params) * Sig_xn2 + (xkm1 - x_n) * Sig_yn2;
		yk = (yc(xk, params) - y_n) * dyc(xk, params) * Sig_xn2 + (xk - x_n) * Sig_yn2;  //function we're finding zero of
		dyk = (std::pow(dyc(xk, params), 2.0) + (yc(xk, params) - y_n)*ddyc(xk, params)) * Sig_xn2 + Sig_yn2; //derivative of above

		r = 1.0 - ((yk / ykm1) * (((yk - ykm1) / (xk - xkm1)) / dyk));

		xkp1 = (1.0 - 1.0 / r)*xkm1 + xk / r;

		if (std::abs(xkp1) > root_bound * datawidth) {
			while (std::abs(xkp1) > root_bound * datawidth) {
				xkp1 = xkp1 / 2.0;
			}
		}

		//bisection?
		
		if (yk * ykm1 < 0) { //checks if xk and xkm1 can be used as bisection brackets

			double c, f_c, f_left, left, right;
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

				f_c = (yc(c, params) - y_n) * dyc(c, params) * Sig_xn2 + (c - x_n) * Sig_yn2;

				if (std::abs((left - right) / 2.0) <= tol_brackets) { //secondary convergence criterion (bracket width)
					break;
				}

				f_left = (yc(left, params) - y_n) * dyc(left, params) * Sig_xn2 + (left - x_n) * Sig_yn2;

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

std::vector <double> TRK::pegToZeroSlop(std::vector <double> vertex){

	if (std::abs(vertex[M]) <= pegToZeroTol) {
		vertex[M] = 0;
	}
	if (std::abs(vertex[M+1]) <= pegToZeroTol) {
		vertex[M+1] = 0;
	}
	
	return vertex;
}

std::vector <double> TRK::avoidNegativeSlop(std::vector <double> vertex, int n) {

	for (int k = 0; k < 2; k++) {
		if (vertex[n - 2 + k] < 0) {
			vertex[n - 2 + k] = std::abs(vertex[n - 2 + k]);
		}
	}
	

	return vertex;
}

double TRK::evalWPriors(double(TRK::*f)(std::vector <double>, double), std::vector <double> vertex, double s) {
	if (hasPriors) {
		switch (priorsObject.priorType) {

		case CONSTRAINED:

			for (int i = 0; i < M; i++) { //check upper bound
				double ub = priorsObject.paramBounds[i][1];
				double lb = priorsObject.paramBounds[i][0];

				if (vertex[i] >= ub && !std::isnan(ub)) {
					return DBL_MAX;
				}

				if (vertex[i] <= lb && !std::isnan(lb)) {
					return DBL_MAX;
				}
			}
			break;

		case MIXED:

			for (int i = 0; i < M; i++) { //check upper bound
				double ub = priorsObject.paramBounds[i][1];
				double lb = priorsObject.paramBounds[i][0];

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

	return (this->*f)(vertex, s);
}

std::vector <double> TRK::downhillSimplex(double(TRK::*f)(std::vector <double>, double), std::vector <double> allparams_guess, double s) {

	double tol = simplexTol;

	int n = allparams_guess.size(); //number of model parameters plus two slop parameters

	double rho = 1.0;
	double chi = 2.0;
	double gamma = 0.5;
	double sigma = 0.5;

	// simplex initialization

	std::vector <double> init_point = allparams_guess;

	std::vector <std::vector <double> > vertices(n + 1, init_point);
	std::vector <double> fitted_params;

	int i = 0;
	for (int j = 1; j < n + 1; j++) { //for each simplex node

		vertices[j][i] = allparams_guess[i] + simplex_size*allparams_guess[i]; //add initial "step size"
		i += 1;
	}

	std::vector <double> result;
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

					double fc = evalWPriors(f, cpoint, s);

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

					double fcc = evalWPriors(f, ccpoint, s);

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
		
		
		
		/*std::cout << "chi-square parameters at s = " << s << " ";
		for (int i = 0; i < result.size(); i++) {
			std::cout << result[i] << " ";
		}
		std::cout << "fitness = " << evalWPriors(f, result, s) << "\n";*/
		
		
		//test for termination

		std::vector <double> evals;
		for (int i = 0; i < n + 1; i++) {
			evals.push_back(evalWPriors(f, vertices[i], s));
		}

		if (stDevUnweighted(evals) < tol) {
			break;
		}
	}

	fitted_params = vertices[n];

	fitted_params = pegToZeroSlop(fitted_params);
	fitted_params = avoidNegativeSlop(fitted_params, n);

	return fitted_params;

	std::cout << std::endl << std::endl;
}

std::vector <double> TRK::cubicSolver(double A, double B, double C, double D) {
	//cubic solver for three real and distinct roots
	double a1 = B / A;
	double a2 = C / A;
	double a3 = D / A;

	double Q = (a1 * a1 - 3.0 * a2) / 9.0;
	double R = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
	double Qc = std::pow(Q, 3.0);
	double d = Qc - std::pow(R, 2.0);

	double theta = std::acos(R / sqrt(Qc));

	double r1 =  -2 * std::sqrt(Q) * std::cos(theta / 3) - a1 / 3;
	double r2 = -2 * std::sqrt(Q) * std::cos((theta + 2 * PI) / 3) - a1 / 3;
	double r3 = -2 * std::sqrt(Q) * std::cos((theta + 4 * PI) / 3) - a1 / 3;

	std::vector <double> roots = { r1, r2, r3 };
	std::vector <double> goodroots;

	for (int i = 0; i < 3; i++) {
		if (std::abs(roots[i]) < (root_bound * datawidth)) {
			goodroots.push_back(roots[i]);
		}
	}

	return goodroots;
}

// STATISTICS

double TRK::normal(double x, double mu, double sig) {
	return (std::exp((-0.5) * (std::pow((x - mu), 2.0) / (2.0 * std::pow(sig, 2.0)))) / std::sqrt(2.0*PI*std::pow(sig, 2.0)));
}

double TRK::singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn) {
	double m_tn = dyc(x_tn, params);
	double y_tn = yc(x_tn, params);
	
	return std::pow(y_n - y_tn - m_tn * (x_n - x_tn), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) - std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(Sig_yn2, 2.0)));
}

std::vector <double> TRK::tangentParallel(std::vector<double> params, double slop_x, double slop_y, int n, double s) {
	double Sig_xn2 = std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0);
	double Sig_yn2 = std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0);

	std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], Sig_xn2, Sig_yn2, x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points

	double x_t = findBestTangent(params, x[n], y[n], Sig_xn2, Sig_yn2, x_tn_vec);

	double m_tn = dyc(x_t, params);
	double y_tn = yc(x_t, params);

	double subsum1 = w[n] * std::pow(y[n] - y_tn - m_tn * (x[n] - x_t), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2);
	double subsum2 = w[n] * std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(s*Sig_yn2, 2.0)));

	return { x_t, subsum1, subsum2 };
}

double TRK::modifiedChiSquared(std::vector <double> allparams, double s)
{
	std::vector <double> SigXVec, SigYVec;
	std::vector <double> all_x_t(N, 0.0);

	double sum1 = 0.0;
	double sum2 = 0.0;

	double slop_x = allparams[M];
	double slop_y = allparams[M + 1];

	std::vector <double> params;

	for (int i = 0; i < M; i++) {
		params.push_back(allparams[i]);
	}

	if (cpp17MultiThread) {

		std::vector <int> nn;

		for (int n = 0; n < N; n++) {
			nn.push_back(n);
		}

		for (int n = 0; n < N; n++) {
			SigXVec.push_back(std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0));
			SigYVec.push_back(std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0));
		}

		std::for_each( //parallel tangent point finding
			std::execution::par_unseq,
			nn.begin(),
			nn.end(),
			[&](auto&& n)
		{
			std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], SigXVec[n], SigYVec[n], x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points

			double x_t = findBestTangent(params, x[n], y[n], SigXVec[n], SigYVec[n], x_tn_vec);

			all_x_t[n] = x_t;
		});

		for (int n = 0; n < N; n++) {

			double m_tn = dyc(all_x_t[n], params);
			double y_tn = yc(all_x_t[n], params);

			sum1 += w[n] * std::pow(y[n] - y_tn - m_tn * (x[n] - all_x_t[n]), 2.0) / (std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n]);
			sum2 += w[n] * std::log((std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n]) / (std::pow(m_tn*SigXVec[n], 2.0) + std::pow(s*SigYVec[n], 2.0)));
		}
	} else if (openMPMultiThread && !cpp17MultiThread) {
		//clock_t time = startTimer();

		#pragma omp parallel for num_threads(maxThreads)
		for (int i = 0; i < N; i++)
		{
			std::vector <double> results;
			results = tangentParallel(params, slop_x, slop_y, i, s); //pointer to fn run through MT, arguments to fn
			all_x_t[i] = results[0];
			sum1 += results[1];
			sum2 += results[2];
		}

		//double sec_elapsed = secElapsed(time);

	} else {
		//cpp11 multithreading

		int counter = 0, completedThreads = 0, liveThreads = 0;
		std::vector<double> results;
		std::vector< std::future < std::vector < double > > > futureVec;
		futureVec.resize(N);

		for (int i = 0; i < N; i++)
		{
			futureVec[i] = std::async(std::launch::async, &TRK::tangentParallel, this, params, slop_x, slop_y, i, s); //pointer to fn run through MT, arguments to fn
			counter++;
			liveThreads++;

			if (liveThreads >= maxThreads)
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
		for (int i = completedThreads; i < N; i++)
		{
			results = futureVec[i].get();
			all_x_t[i] = results[0];
			sum1 += results[1];
			sum2 += results[2];
		}

		/*std::vector <std::thread> ths;
		for (int n = 0; n < N; n++) {
			ths.push_back(std::thread(&TRK::tangentParallel, this, params, slop_x, slop_y, n));
		}
		for (auto& th : ths) {
			th.join();
		}*/
		
	}

	switch (whichExtrema) {
		case none:
			break;
		case S:
			x_t_s = all_x_t;
			params_s = params;
			break;
		case slopx:
			x_t_slopx = all_x_t;
			params_slopx = params;
			break;
		case slopy:
			x_t_slopy = all_x_t;
			params_slopy = params;
			break;
		default:
			break;
	}

	return sum1 - sum2;
}

double TRK::regularChiSquared(std::vector <double> params) {
	int N = y.size();

	double sum = 0.0;

	for (int i = 0; i < N; i++) {
		sum += w[i] * std::pow(y[i] - (*yc)(x[i], params), 2.0);
	}
	return sum;
}

std::vector <double> TRK::tangentParallelLikelihood(std::vector<double> params, double slop_x, double slop_y, int n) {
	double Sig_xn2 = std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0);
	double Sig_yn2 = std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0);

	std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], Sig_xn2, Sig_yn2, x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points

	double x_t = findBestTangent(params, x[n], y[n], Sig_xn2, Sig_yn2, x_tn_vec);

	double m_tn = dyc(x_t, params);
	double y_tn = yc(x_t, params);

	double l = w[n] * std::sqrt((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(s*Sig_yn2, 2.0)));
	l *= std::exp(-0.5 * w[n] * (std::pow(y[n] - y_tn - m_tn * (x[n] - x_t), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2)));

	return { x_t, l};
}

double TRK::likelihood(std::vector <double> allparams) {
	std::vector <double> SigXVec, SigYVec;
	std::vector <double> all_x_t(N, 0.0);
	double L = 1.0;

	double slop_x = allparams[M];
	double slop_y = allparams[M + 1];

	std::vector <double> params;

	for (int i = 0; i < M; i++) {
		params.push_back(allparams[i]);
	}

	if (cpp17MultiThread) {

		std::vector <int> nn;

		for (int n = 0; n < N; n++) {
			nn.push_back(n);
		}

		for (int n = 0; n < N; n++) {
			SigXVec.push_back(std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0));
			SigYVec.push_back(std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0));
		}

		std::for_each( //parallel tangent point finding
			std::execution::par_unseq,
			nn.begin(),
			nn.end(),
			[&](auto&& n)
		{
			std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], SigXVec[n], SigYVec[n], x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points

			double x_t = findBestTangent(params, x[n], y[n], SigXVec[n], SigYVec[n], x_tn_vec);

			all_x_t[n] = x_t;
		});

		for (int n = 0; n < N; n++) {
			double m_tn = dyc(all_x_t[n], params);
			double y_tn = yc(all_x_t[n], params);

			L *= w[n] * std::sqrt((std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n]) / (std::pow(m_tn*SigXVec[n], 2.0) + std::pow(s*SigYVec[n], 2.0)));
			L *= std::exp(-0.5 * w[n] * (std::pow(y[n] - y_tn - m_tn * (x[n] - all_x_t[n]), 2.0) / (std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n])));
		}
	} else if (openMPMultiThread && !cpp17MultiThread) {
		//clock_t time = startTimer();

		#pragma omp parallel for num_threads(maxThreads)
		for (int i = 0; i < N; i++)
		{
			std::vector <double> results;
			results = tangentParallelLikelihood(params, slop_x, slop_y, i); //pointer to fn run through MT, arguments to fn
			all_x_t[i] = results[0];
			L *= results[1];
		}

		//double sec_elapsed = secElapsed(time);

	}
	else {
		//cpp11 multithreading

		int counter = 0, completedThreads = 0, liveThreads = 0;
		std::vector<double> results;
		std::vector< std::future < std::vector < double > > > futureVec;
		futureVec.resize(N);

		for (int i = 0; i < N; i++)
		{
			futureVec[i] = std::async(std::launch::async, &TRK::tangentParallelLikelihood, this, params, slop_x, slop_y, i); //pointer to fn run through MT, arguments to fn
			counter++;
			liveThreads++;

			if (liveThreads >= maxThreads)
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
		for (int i = completedThreads; i < N; i++)
		{
			results = futureVec[i].get();
			all_x_t.push_back(results[0]);
			L *= results[1];
		}
	}
	return L;
}

double TRK::priors(std::vector <double> allparams) {
	double jointPrior = 1.0; //uninformative prior by default

	switch (priorsObject.priorType){
		case CONSTRAINED:
			for (int i = 0; i < M; i++) { //check upper bound
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
			for (int i = 0; i < M + 2; i++) { //check upper bound
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

			for (int i = 0; i < M + 2; i++) { //check upper bound
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
			for (int i = 0; i < M + 2; i++) {
				jointPrior *= priorsObject.priorsPDFs[i](allparams[i]);
			}

			break;
	}

	return jointPrior;
}

double TRK::posterior(std::vector <double> allparams) {
	if (hasPriors) {
		return likelihood(allparams) * priors(allparams);
	}
	else {
		return likelihood(allparams);
	}
}

double TRK::stDevUnweighted(std::vector <double> x) {
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

	std::vector <double> xr1vec;
	double xr1old;
	
	bool checkcheck = false;

	int itcount = 0;

	while (true) {
		if (xr1vec.size() > 99) {
			printf("100 iterations of tangent finder loop! \n");
			for (int j = 0; j < params.size(); j++) {
				printf("%f ", params[j]);
			}
			printf("%f  %f  %f  %f \t", Sig_xn2, Sig_yn2, x_n, y_n);
			printf("s = %f\n", s);
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

			
			double xr_left = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, x_min, x_min - std::sqrt(Sig_xn2) / 10.0);
			double xr_right = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, x_max, x_max + std::sqrt(Sig_xn2) / 10.0);

			result.push_back(xr1);
			result.push_back(xr_left);
			result.push_back(xr_right);

			break;
		} else if (extraRoots.size() == 0 && xr1vec.size() == 2) { //found two roots but can't find a third
			result = xr1vec;
			break;
		} else if (extraRoots.size() == 3) {//this can happen if the "root" found is very close to being a root but isn't actually one.
			
			//in this case, try again with a different guess.
			double xr_left = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, x_min, x_min - std::sqrt(Sig_xn2) / 10.0);
			double xr_right = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, x_max, x_max + std::sqrt(Sig_xn2) / 10.0);

			result.push_back(xr_left);
			result.push_back(xr_right);

			break;
		}

		itcount += 1;
	}

	return result;
}

double TRK::findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec) {
	std::vector <double> posts;
	int minindex;

	for (int i = 0; i < x_tn_vec.size(); i++) {
		posts.push_back(singlePointLnL(params, x_n, y_n, Sig_xn2, Sig_yn2, x_tn_vec[i]));
	}

	std::vector<double>::iterator result = std::min_element(std::begin(posts), std::end(posts));
	minindex = std::distance(std::begin(posts), result);

	return x_tn_vec[minindex];
}

// SCALE OPTIMIZATION ALGORITHMS
void TRK::getBetterSlopYGuess(double slop_y, double s) {
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

double TRK::innerSlopX_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
	//s = ss[0];

	//clock_t time = startTimer();

	std::vector <double> allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess, ss[0]);

	printf("%f \t %f \t %f \n", ss[0], allparams_s[M], allparams_s[M + 1]);

	//double sec_elapsed = secElapsed(time);

	//printf("%f sec, max threads = %i \n", sec_elapsed, maxThreads);

	getBetterSlopYGuess(allparams_s[M + 1], s);

	return allparams_s[M];
}

double TRK::innerSlopY_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
	//s = ss[0];

	std::vector <double> allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess, ss[0]);

	printf("%f \t %f \t %f \n", ss[0], allparams_s[M], allparams_s[M + 1]);

	return allparams_s[M + 1];
}

double TRK::innerR2_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
	//s = ss[0];

	whichExtrema = S;
	std::vector <double> allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess, ss[0]);
	whichExtrema = none;

	printf("%f \t %f \t %f \n", ss[0], allparams_s[M], allparams_s[M + 1]);

	double R2as = R2TRK_prime_as();
	double R2sb = R2TRK_prime_sb();

	return R2as - R2sb;
}

double TRK::innerR2_iter_Simplex(std::vector <double> ss, std::vector <double> allparams_guess, double s0) {
	//s = ss[0];

	whichExtrema = S;
	allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess, ss[0]);
	whichExtrema = none;

	printf("%f \t %f \t %f \n", ss[0], allparams_s[M], allparams_s[M + 1]);

	double R2as = R2TRK_prime_as0(s0, x_t_s, params_s);
	double R2sb = R2TRK_prime_s0b(s0, x_t_s, params_s);

	return R2as - R2sb;
}

double TRK::optimize_s_SlopX() {

	iterative_allparams_guess = allparams_guess;

	// before doing any standard simplex movement, here it checks whether the simplex is within the zero "plateau", and if so, it moves it to the boundary.
	// for slop x: move to right until it hits the boundary


	//bracket finding

	double a, b;
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

		whichExtrema = slopx;
		slop_c = innerSlopX_Simplex({ c }, iterative_allparams_guess);
		whichExtrema = none;

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

double TRK::optimize_s_SlopY() {

	iterative_allparams_guess = allparams_guess;

	// before doing any standard simplex movement, here it checks whether the simplex is within the zero "plateau", and if so, it moves it to the boundary.
	// for slop x: move to right until it hits the boundary

		//bracket finding

	double a, b;
	double trial_s = slopYScaleGuess;
	double slop_trial_s = innerSlopY_Simplex({ trial_s }, iterative_allparams_guess);

	double inc = trial_s * 0.5;

	if (slop_trial_s > 0) {
		a = trial_s;

		double trial_b = trial_s;

		while (true) {
			trial_b += inc;

			double slop_trial_b = innerSlopY_Simplex({ trial_b }, iterative_allparams_guess);

			if (slop_trial_b == 0) {
				b = trial_b;
				break;
			}
			else if (slop_trial_b > 0) {
				a = trial_b;
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

		whichExtrema = slopy;
		slop_c = innerSlopY_Simplex({ c }, iterative_allparams_guess);
		whichExtrema = none;

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

double TRK::R2TRK_prime_as() {
	double R2 = 1.0 / N;

	double sum = 0.0;

	for (int n = 0; n < N; n++) {
		double m_tn_a = dyc(x_t_a[n], params_a);
		double theta_t_a = std::atan(m_tn_a);

		double m_tn_s = dyc(x_t_s[n], params_s);
		double theta_t_s = std::atan(m_tn_s);

		sum += std::pow(std::tan(PI/4.0 - std::abs(theta_t_a - theta_t_s)/ 2.0), 2.0);
	}

	R2 *= sum;

	return R2;
}

double TRK::R2TRK_prime_sb() {
	double R2 = 1.0 / N;

	double sum = 0.0;

	for (int n = 0; n < N; n++) {
		double m_tn_s = dyc(x_t_s[n], params_s);
		double theta_t_s = std::atan(m_tn_s);

		double m_tn_b = dyc(x_t_b[n], params_b);
		double theta_t_b = std::atan(m_tn_b);

		sum += std::pow(std::tan(PI / 4.0 - std::abs(theta_t_s - theta_t_b) / 2.0), 2.0);
	}

	R2 *= sum;

	return R2;
}

double TRK::R2TRK_prime_as0(double s0, std::vector <double> x_t_s1, std::vector <double> params_s1) {
	double R2 = 1.0 / N;

	double sum = 0.0;

	for (int n = 0; n < N; n++) {
		double m_tn_a = dyc(x_t_a[n], params_a);
		double theta_t_a = std::atan(m_tn_a);

		double m_tn_s1 = dyc(x_t_s1[n], params_s1);
		double theta_t_s1 = std::atan(m_tn_s1);

		sum += std::pow(std::tan(PI / 4.0 - std::abs(std::atan(s0*std::tan(theta_t_a)) - std::atan(s0*std::tan(theta_t_s1))) / 2.0), 2.0);
	}

	R2 *= sum;

	return R2;
}

double TRK::R2TRK_prime_s0b(double s0, std::vector <double> x_t_s1, std::vector <double> params_s1) {
	double R2 = 1.0 / N;

	double sum = 0.0;

	for (int n = 0; n < N; n++) {
		double m_tn_b = dyc(x_t_b[n], params_b);
		double theta_t_b = std::atan(m_tn_b);

		double m_tn_s1 = dyc(x_t_s1[n], params_s1);
		double theta_t_s1 = std::atan(m_tn_s1);

		sum += std::pow(std::tan(PI / 4.0 - std::abs(std::atan(s0*std::tan(theta_t_s1)) - std::atan(s0*std::tan(theta_t_b))) / 2.0), 2.0);
	}

	R2 *= sum;

	return R2;
}

double TRK::optimize_s0_R2() {

	iterative_allparams_guess = allparams_guess;

	//bracket finding

	double left, right;
	
	//bisection, now that we have brackets [left,right]

	left = a;
	right = b;

	double c, f_c;
	double tol_bisect = 1e-4;
	double tol_brackets = 1e-3;

	while (true) {
		c = (left + right) / 2;

		whichExtrema = S;
		f_c = innerR2_Simplex({ c }, iterative_allparams_guess);
		whichExtrema = none;

		if (std::abs(f_c) <= tol_bisect) { //convergence criterion
			break;
		}

		if (f_c > 0) {
			left = c;
		}
		else if (f_c < 0) {
			right = c;
		}

		if (std::abs(left - right) <= tol_brackets) { //secondary convergence criterion (bracket width)
			break;
		}
	}

	return c;

}

double TRK::optimize_s_prime_R2(double s0) {

	iterative_allparams_guess = allparams_guess;

	//bracket finding

	double left, right;

	//bisection, now that we have brackets [left,right]

	left = a;
	right = b;

	double c, f_c;
	double tol_bisect = 1e-4;
	double tol_brackets = 1e-3;

	while (true) {
		c = (left + right) / 2;

		whichExtrema = S;
		f_c = innerR2_iter_Simplex({ c }, iterative_allparams_guess, s0);
		whichExtrema = none;

		if (std::abs(f_c) <= tol_bisect) { //convergence criterion
			break;
		}

		if (f_c > 0) {
			left = c;
		}
		else if (f_c < 0) {
			right = c;
		}

		if (std::abs(left - right) <= tol_brackets) { //secondary convergence criterion (bracket width)
			break;
		}
	}

	return c;

}

double TRK::iterateR2_OptimumScale(double s0) {
	double tol_scale = 1e-3;

	double s1;

	bool tolcheck = false;

	while (!tolcheck) {
		printf("next s0: %f \n", s0);

		s1 = optimize_s_prime_R2(s0);
		if (std::abs(s1-s0) <= tol_scale) {
			tolcheck = true;
		}
		std::printf("new s0: %f \n", s1);
		s0 = s1;
	}

	return s1;
}

void TRK::optimizeScale() {
	s = 1.0; //initially begin with s = 1

	std::vector <double> scale_extrema;

	std::vector <double> s_slops(2, 0.0);

	//optimize simultaneously
	if (cpp17MultiThread) {
		std::vector <int> nn = { 0, 1 };

		std::for_each( //parallel tangent point finding
			std::execution::par_unseq,
			nn.begin(),
			nn.end(),
			[&](auto&& n)
		{
			s_slops[n] = (this->*optimizeList[n])();
		});

	} else {
		#pragma omp parallel for //num_threads(8)
		for (int i = 0; i < 2; i++)
		{
			s_slops[i] = (this->*optimizeList[i])();
		}
	}

	double s_slopx = s_slops[0];
	double s_slopy = s_slops[1];

	scale_extrema.push_back(s_slopx); //scale_extrema = {s_slopx, s_slopy}

	scale_extrema.push_back(s_slopy);

	std::vector <int> sortedindices = getSortedIndices(scale_extrema);

	a = scale_extrema[sortedindices[0]]; //figures out which scale is a, and which is b, as well as storing the associated best-fit parameters and their associated tangent points for those two extreme scales
	b = scale_extrema[sortedindices[1]];

	printf("extrema: \t a \t b: \n");
	printf(" \t %f \t %f \n", a, b);

	if (a == s_slopx) {
		x_t_a = x_t_slopx;
		x_t_b = x_t_slopy;

		params_a = params_slopx;
		params_b = params_slopy;
	}
	else if (a == s_slopy) {
		x_t_a = x_t_slopy;
		x_t_b = x_t_slopx;

		params_a = params_slopy;
		params_b = params_slopx;
	}

	//determine best s1 (new s) to satistfy R2TRKp(a,s) = R2TRKp(s,b)

	std::cout << "finding optimum scale..." << std::endl;

	double s0 = optimize_s0_R2();

	s = s0;

	double s_final = iterateR2_OptimumScale(s0);

	s = s_final;

	std::cout << "optimum s = " << s << std::endl;

	results.bestFitParams.clear();

	for (int j = 0; j < M; j++) {
		results.bestFitParams.push_back(allparams_s[j]);
	}

	results.slop_x = allparams_s[M];
	results.slop_y = allparams_s[M + 1];

	results.optimumScale = s;
	results.minimumScale = a;
	results.maximumScale = b;

	return;
}


//MCMC

std::vector <std::vector <double >> TRK::checkSlopSignMCMC(std::vector <std::vector <double >> result_final) {

	std::vector <std::vector <double >> result_final_fixed;
	std::vector <double> inner;

	for (int i = 0; i < result_final.size(); i++) {
		inner.clear();

		for (int j = 0; j < M; j++) {
			inner.push_back(result_final[i][j]);
		}
		inner.push_back(std::abs(result_final[i][M]));
		inner.push_back(std::abs(result_final[i][M+1]));

		result_final_fixed.push_back(inner);
	}

	return result_final_fixed;
}

double TRK::innerMetHastSimplex(int burncount, std::vector <double> delta, double best_ratio) { //does burn in + (1000) MCMC samples with given deltas and returns the acceptance ratio
	int sampleCount = 1000;

	std::vector < std::vector <double > > result, result_final;
	std::vector <double> allparams_trial, allparams_0; //allparams_0 is the previous step
	double a, rand_unif, accept_frac;

	int accept_count = 0;
	int delta_count = 0;

	allparams_0 = allparams_guess;


	while (delta_count < sampleCount){// + burncount) {
		//create trial

		allparams_trial.clear();

		for (int j = 0; j < M + 2; j++) {
			allparams_trial.push_back(delta[j] * rnorm(0.0, 1.0) + allparams_0[j]);
		}

		a = posterior(allparams_trial) / posterior(allparams_0);
		rand_unif = runiform(0.0, 1.0);

		if (a >= 1) {
			allparams_0 = allparams_trial;
			delta_count += 1;
			result.push_back(allparams_0);
			accept_count += 1;
		}
		else if (rand_unif <= a) {
			allparams_0 = allparams_trial;
			delta_count += 1;
			result.push_back(allparams_0);
			accept_count += 1;
		}
		else {
			delta_count += 1;
			result.push_back(allparams_0);
		}
	}

	accept_frac = (double)accept_count / (double)delta_count;

	printf("inner simplex acceptance ratio: %f    with deltas: ", accept_frac);

	for (int j = 0; j < M + 2; j++) {
		printf("%f ", delta[j]);
	}
	std::cout << std::endl;

	return std::abs(accept_frac - best_ratio);
}

std::vector <double> TRK::pegToNonZeroDelta(std::vector <double> vertex, std::vector <double> lastvertex) {

	std::vector <double> vertexfixed = vertex;

	for (int j = 0; j < M + 2; j++) {
		if (vertex[j] == 0.0) {
			vertexfixed[j] = lastvertex[j] * 0.5;
		}
	}

	return vertexfixed;
}


std::vector <double> TRK::optimizeMetHastDeltas(int burncount, std::vector <double> delta_guess) {
	double tol = 0.05;
	double optRatio = 0.45;

	int n = delta_guess.size(); //number of model parameters plus two slop parameters

	double rho = 1.0;
	double chi = 2.0;
	double gamma = 0.5;
	double sigma = 0.5;

	// simplex initialization

	std::vector <double> init_point = delta_guess;

	std::vector <std::vector <double> > vertices(n + 1, init_point);
	std::vector <double> best_delta;

	int i = 0;
	for (int j = 1; j < n + 1; j++) { //for each simplex node

		vertices[j][i] = delta_guess[i] +  delta_guess[i]; //add initial "step size"
		i += 1;
	}

	std::vector <double> result;
	while (true) {
		while (true) {
			// order

			std::vector <int> orderedindices;
			std::vector <double> unorderedEvals; // ( f(x_1), f(x_2), ... f(x_n+1)
			for (int i = 0; i < n + 1; i++) {
				unorderedEvals.push_back(innerMetHastSimplex(burncount, vertices[i], optRatio));
				if (unorderedEvals[i] < tol) {
					return vertices[i];
				}
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

			refpoint = pegToNonZeroDelta(refpoint, vertices[n]);

			double fr = innerMetHastSimplex(burncount, refpoint, optRatio);
			if (fr < tol) {
				return refpoint;
			}
			double f1 = innerMetHastSimplex(burncount, vertices[0], optRatio);
			if (f1 < tol) {
				return vertices[0];
			}
			double fn = innerMetHastSimplex(burncount, vertices[n - 1], optRatio);
			if (fn < tol) {
				return vertices[n - 1];
			}

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

				exppoint = pegToNonZeroDelta(exppoint, vertices[n]);

				double fe = innerMetHastSimplex(burncount, exppoint, optRatio);
				if (fe < tol) {
					return exppoint;
				}


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
				double fnp1 = innerMetHastSimplex(burncount, vertices[n], optRatio);
				if (fnp1 < tol) {
					return vertices[n];
				}

				if (fn <= fr && fr < fnp1) {
					std::vector <double> cpoint;

					for (int i = 0; i < n; i++) {
						cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
					}

					cpoint = pegToNonZeroDelta(cpoint, vertices[n]);

					double fc = innerMetHastSimplex(burncount, cpoint, optRatio);
					if (fc < tol) {
						return cpoint;
					}

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

					ccpoint = pegToNonZeroDelta(ccpoint, vertices[n]);

					double fcc = innerMetHastSimplex(burncount, ccpoint, optRatio);
					if (fcc < tol) {
						return ccpoint;
					}

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

		/*

		std::cout << "chi-square parameters at s = " << s << " ";
		for (int i = 0; i < result.size(); i++) {
			std::cout << result[i] << " ";
		}
		std::cout << "fitness = " << evalWPriors(f, vertices[i]) << "\n";
		*/

		//test for termination

		if (innerMetHastSimplex(burncount, vertices[n], optRatio) < tol) {
			break;
		}

	}

	best_delta = vertices[n];

	return best_delta;
}

std::vector <std::vector <double >> TRK::methastPosterior(int R, int burncount, std::vector <double> sigmas_guess) {
	
	//initialization of adaptive delta
	std::vector <double> delta;
	printf("initial delta:");

	for (int j = 0; j < M + 2 ; j++) {
		printf("%f ", sigmas_guess[j]);
	}
	std::cout << std::endl;

	//optimize deltas

	delta = optimizeMetHastDeltas(burncount, sigmas_guess);

	printf("final delta:");

	for (int j = 0; j < M + 2; j++) {
		printf("%f ", delta[j]);
	}
	std::cout << std::endl;

	std::vector < std::vector <double > > result, result_final;
	std::vector <double> allparams_trial, allparams_0; //allparams_0 is the previous step
	double a, rand_unif, accept_frac;

	int accept_count = 0;
	int delta_count = 0;
	double deltafactor = 1.0;
	int down_count = 0;
	int up_count = 0;

	allparams_0 = allparams_guess;

	while (delta_count < R + burncount) {
		//create trial

		allparams_trial.clear();

		for (int j = 0; j < M + 2; j++) {
			allparams_trial.push_back(delta[j] * rnorm(0.0, 1.0) + allparams_0[j]);
		}

		a = posterior(allparams_trial) / posterior(allparams_0);
		rand_unif = runiform(0.0, 1.0);

		if (a >= 1) {
			allparams_0 = allparams_trial;
			delta_count += 1;
			result.push_back(allparams_0);
			accept_count += 1;
		}
		else if (rand_unif <= a) {
			allparams_0 = allparams_trial;
			delta_count += 1;
			result.push_back(allparams_0);
			accept_count += 1;
		}
		else {
			delta_count += 1;
			result.push_back(allparams_0);
		}

		accept_frac = (double) accept_count / (double)delta_count;
	}

	//cut off burn-in

	for (int i = 0; i < R; i++) {
		result_final.push_back(result[i + burncount]);
	}

	printf("number of delta changes: %i \n", up_count + down_count);
	printf("final delta:");
	for (int j = 0; j < M + 2; j++) {
		printf("%f ", delta[j]);
	}
	printf("\t final acceptance ratio: %f \n", accept_frac);

	result_final = checkSlopSignMCMC(result_final);

	return result_final;
}

double TRK::rnorm(double mu, double sig) {

	double rand;

	std::random_device dev;
	std::mt19937 generator(dev());

	std::normal_distribution <double> dist(mu, sig);
	rand = dist(generator);
	return rand;
}

double TRK::runiform(double a, double b) {
	double rand;

	std::random_device dev;
	std::mt19937 generator(dev());

	std::uniform_real_distribution <double> dist(a, b);
	rand = dist(generator);
	return rand;
}

double noPrior(double param) {
	return 1.0;
}

std::vector <std::vector <double> > TRK::getHistogram(std::vector <double> data) {
	int dataSize = data.size();

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

std::vector <std::vector <std::vector <double> > >  TRK::lowerBar(std::vector <std::vector <double> > allparam_samples) { //method used to estimate uncertainty of a sampled distribution
	std::vector <double> data, hist, edges, minusSigmas, plusSigmas;
	std::vector <std::vector <std::vector <double> > > allparam_uncertainties;
	std::vector <std::vector <double> > histResults;
	int totalCount = allparam_samples.size();
	double bar;
	double tolBar = 1e-6;

	results.paramDistributionHistograms.clear();

	for (int j = 0; j < M + 2; j++) { //for each model param plus slop
		data.clear();

		for (int i = 0; i < totalCount; i++) {
			data.push_back(allparam_samples[i][j]);
		}

		histResults = getHistogram(data); 

		results.paramDistributionHistograms.push_back(histResults);

		hist = histResults[0]; // each number in this is the number of samples within each bin
		edges = histResults[1];

		double mean = getAverage(data);

		//bar lowering iteration (essentially a bisection algo)

		double low = 0.0;
		double high = minMax(hist)[1];
		double bar = (low + high) / 2.0;
		double leftBound, rightBound, minusSig, plusSig;
		int aboveCount, K; // number of SAMPLES within the bins above bar

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

			K = indicesIn.size();
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

void TRK::calculateUncertainties() {

	std::vector <std::vector <std::vector <double> > > allparam_uncertainties;

	std::cout << "Sampling Posterior...\n";

	std::vector <std::vector <double >> allparam_samples = methastPosterior(R, burncount, allparams_sigmas_guess);

	if (outputDistributionToFile) {

		std::string fileName = std::string("TRKMCMC_") + std::to_string(allparams_guess[0]) + std::string("_") + std::to_string(R) + std::string(".txt");

		std::ofstream myfile;
		myfile.open(fileName, std::ofstream::trunc);
		if (myfile.is_open())
		{
			// filename    a     b     optimum scale    total computation time (s)
			for (int i = 0; i < allparam_samples.size(); i++) {
				for (int j = 0; j < allparams_guess.size(); j++) {
					myfile << allparam_samples[i][j] << " ";
				}
				myfile << std::endl;
			}

			myfile.close();
		}
		else std::cout << "Unable to open file";
	}

	std::cout << "Computing Uncertainties...\n";

	allparam_uncertainties = lowerBar(allparam_samples); //for each parameter including slope, there is a vector containing 1 vector of -sigmas, 1 vector of +sigmas. This vector contains all of those 2-vectors.

	results.bestFit_123Sigmas.clear();

	for (int j = 0; j < M; j++) {
		results.bestFit_123Sigmas.push_back(allparam_uncertainties[j]);
	}
	results.slopX_123Sigmas = allparam_uncertainties[M];
	results.slopY_123Sigmas = allparam_uncertainties[M + 1];

	return;
}

// PIVOT POINTS

void TRK::findPivots() {
	if (findPivotPoints) {

	}
	else {
		return;
	}
}


// CORE ALGORITHMS/TRK FITS

void TRK::performTRKFit() {//finds optimum scale AND performs TRK fit + uncertainty
	optimizeScale();

	calculateUncertainties();
}
void TRK::performTRKFit(double scale) {//perform fit on some provided scale (for example, if they already know optimum scale, they can just start with this)
	s = scale;

	results.bestFitParams.clear();

	whichExtrema = S;
	allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess, scale);
	whichExtrema = none;

	for (int j = 0; j < M; j++) {
		results.bestFitParams.push_back(allparams_s[j]);
	}

	results.slop_x = allparams_s[M];
	results.slop_y = allparams_s[M + 1];

	calculateUncertainties();
}
void TRK::performSimpleTRKFit() {//finds optimum scale and and performs TRK fit but without finding uncertainties
	optimizeScale(); // (stores results in TRK.results)

	return;
}



//testing purposes only

clock_t startTimer() {
	std::cout << "starting timer... \n";

	clock_t t_i = clock();
	return t_i;
}

double secElapsed(clock_t t_i) {
	clock_t t_f = clock() - t_i;
	double sec_elapsed = ((float)t_f) / CLOCKS_PER_SEC;

	printf("%f seconds elapsed \n", sec_elapsed);
	return sec_elapsed;
}

void writeResults(TRK TRKobj, double t_sec, std::string filename) {

	std::ofstream myfile("TRKresults.txt", std::ofstream::app);
	if (myfile.is_open())
	{	
		// filename    a     b     optimum scale    total computation time (s)
		myfile << filename << "\t" <<  TRKobj.a << "\t" << TRKobj.b << "\t" << TRKobj.s << "\t" << t_sec << "\n";
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