#include "pch.h"
#include "TRK.h"

// PRIORS

// Constructors

/*
Priors::Priors(priorTypes priorType, std::vector <double>(*p)(std::vector <double>, std::vector <double>)) { //custom priors
	this->priorType = priorType;
	this->p = (*p);
};
*/

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

//default constructor
Priors::Priors() {

};


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

	std::vector <double> allparams_guess = params_guess;

	allparams_guess.push_back(slop_x_guess);
	allparams_guess.push_back(slop_y_guess);

	this->allparams_guess = allparams_guess;

	this->whichExtrema = none;

	this->s = 1.0;

	getDataWidth();

	this->hasPriors = false;
}
//equal weights/unweighted
TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess) {
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
}


//priors:
//weighted
TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, Priors priorsObject) {
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
}
//equal weights/unweighted
TRK::TRK(double(*yc)(double, std::vector <double>), double(*dyc)(double, std::vector <double>), double(*ddyc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sx, std::vector <double> sy, std::vector <double> params_guess, double slop_x_guess, double slop_y_guess, Priors priorsObject) {
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
		//std::cout << x1 << std::endl;

		itercount += 1;
	}

	//std::cout << itercount << " iterations." << std::endl;
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

	while (true) {
		/*
		if (std::isnan(xk) || std::isnan(xkm1)) {
			printf("stop");
		}
		*/
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
		//std::cout << xkp1 << std::endl;


		//bisection?
		
		if (yk * ykm1 < 0) { //checks if xk and xkm1 can be used as bisection brackets

			//std::cout << "beginning bisection (tangent point finder) \n";

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
		
		/*
		if (std::isnan(xkp1) || std::isnan(xk) || std::isnan(xkm1)) {
			printf("stop");
		}
		*/

		xkm1 = xk;
		xk = xkp1;

		itercount += 1;

		/*
		if (std::isnan(xkp1) || std::isnan(xk) || std::isnan(xkm1)) {
			printf("stop");
		}
		*/

		if (itercount > 10000) {
			//std::cout << " itercount of >10000 reached; exiting NR" << std::endl;
			return NAN;
		}

		if (std::abs(xk - xkm1) <= tol || std::isnan(xk) || std::abs(yk) <= tol) {
			break;
		}
	}

	//std::cout << itercount << " iterations." << std::endl;
	return xkp1;
}

std::vector <double> TRK::pegToZeroSlop(std::vector <double> vertex){

	double tol = 0.004;

	if (std::abs(vertex[M]) <= tol) {
		vertex[M] = 0;
	}
	if (std::abs(vertex[M+1]) <= tol) {
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

double TRK::evalWPriors(double(TRK::*f)(std::vector <double>), std::vector <double> vertex) {
	double f_test = (this->*f)(vertex);

	if (hasPriors) {
		switch (priorsObject.priorType) {

		case CONSTRAINED:
			bool exclude = false;

			for (int i = 0; i < M; i++) { //check upper bound
				double ub = priorsObject.paramBounds[i][1];
				double lb = priorsObject.paramBounds[i][0];

				if (vertex[i] >= ub && !std::isnan(ub)) {
					exclude = true;
				}

				if (vertex[i] <= lb && !std::isnan(lb)) {
					exclude = true;
				}
			}

			if (exclude) {
				return DBL_MAX;
			}
			else {
				return f_test;
			}
		}
	}
	else {
		return f_test;
	}

}

std::vector <double> TRK::downhillSimplex(double(TRK::*f)(std::vector <double>), std::vector <double> allparams_guess) {

	double tol = 1e-6;

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
				unorderedEvals.push_back(evalWPriors(f, vertices[i]));
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

			//refpoint = pegToZeroSlop(refpoint);

			double fr = evalWPriors(f, refpoint);
			double f1 = evalWPriors(f, vertices[0]);
			double fn = evalWPriors(f, vertices[n - 1]);

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

				//exppoint = pegToZeroSlop(exppoint);

				double fe = evalWPriors(f, exppoint);


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
				double fnp1 = evalWPriors(f, vertices[n]);

				if (fn <= fr && fr < fnp1) {
					std::vector <double> cpoint;

					for (int i = 0; i < n; i++) {
						cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
					}

					//cpoint = pegToZeroSlop(cpoint);

					double fc = evalWPriors(f, cpoint);

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

					//ccpoint = pegToZeroSlop(ccpoint);

					double fcc = evalWPriors(f, ccpoint);

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
		//result = avoidNegativeSlop(result, n);

		bettervertices.push_back(result);

		vertices = bettervertices;
		
		/*
		
		std::cout << "chi-square parameters at s = " << s << " ";
		for (int i = 0; i < result.size(); i++) {
			std::cout << result[i] << " ";
		}
		std::cout << "fitness = " << evalWPriors(f, vertices[i]) << "\n";
		*/

		
		/*
		std::cout << "chi-square minimized parameters at s = " << s << std::endl;
		for (int j = 0; j < result.size() + 1; j++) {
			for (int i = 0; i < result.size(); i++) {
				std::cout << bettervertices[j][i] << " ";
			}
			double X = evalWPriors(f, bettervertices[j]);
			std::cout << "          " << X << std::endl;
		}
		*/
		
		

		//test for termination

		std::vector <double> evals;
		for (int i = 0; i < n + 1; i++) {
			evals.push_back(evalWPriors(f, vertices[i]));
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

	/*
	if (goodroots.size() < 3) {
		std::cout << "some roots out of reasonable boundary";
	}
	*/

	return goodroots;
}

// STATISTICS
double TRK::singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn) {
	double m_tn = dyc(x_tn, params);
	double y_tn = yc(x_tn, params);
	
	return std::pow(y_n - y_tn - m_tn * (x_n - x_tn), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) - std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(Sig_yn2, 2.0)));
}

double TRK::modifiedChiSquared(std::vector <double> allparams)
{

	std::vector <double> all_x_t;

	double sum1 = 0.0;
	double sum2 = 0.0;

	double slop_x = allparams[M];
	double slop_y = allparams[M + 1];

	std::vector <double> params;

	for (int i = 0; i < M; i++) {
		params.push_back(allparams[i]);
	}

	for (int n = 0; n < N; n++) {
		
		/*
		if (n == 43) {
			test += 1;
			printf("stopp \t %i \n", test);
		}
		*/
		
		/*
		
		if (test == 291) {
			printf("STOP \n");
		}
		*/

		double Sig_xn2 = std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0);
		double Sig_yn2 = std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0);

		//std::cout << x[n] << " " << y[n] << std::endl;

		/*
		if (std::abs(x[n] - 0.015) < 1e-4 && std::abs(y[n] - 3.8727) < 1e-4 ){//&& std::abs(params[0] - 2.4855) < 1e-4 && std::abs(params[1] - 4.573) < 1e-4 && std::abs(params[2] - 2.2869) < 1e-4 && std::abs(params[3] - -1.3164) < 1e-4) {
			std::cout << "stop" << std::endl;
		}
		*/

		std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], Sig_xn2, Sig_yn2, x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points

		int ck = 0;

		for (int i = 0; i < x_tn_vec.size(); i++) {
			if (std::isnan(x_tn_vec[i])) {
				ck += 1;
			}
		}

		/*
		if (ck > 0) {
			printf("STOP\n");
		}
		*/
		double x_t = findBestTangent(params, x[n], y[n], Sig_xn2, Sig_yn2, x_tn_vec);
		

		//printf("%f  %f  %f  %f  %f  %f  %f  %f  %f \n", params[0], params[1], params[2], params[3], Sig_xn2, Sig_yn2, x[n], y[n], x_t);

		double m_tn = dyc(x_t, params);
		double y_tn = yc(x_t, params);

		all_x_t.push_back(x_t);

		sum1 += w[n] * std::pow(y[n] - y_tn - m_tn * (x[n] - x_t), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2);
		sum2 += w[n] * std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(s*Sig_yn2, 2.0)));

		//printf("%f \t %f \n", add1, add2);

		/*
		if (add1 > 1000 || add2 > 1000) {
			printf("stop\n");
		}
		*/
		/*
		if (std::isnan(sum1) || std::isnan(sum1)) {
			printf("stop \n");
		}
		*/
		/*
		std::ofstream myfile("C:\\Users\\nickk124\\Documents\\Reichart Research\\TRK\\TRKtangents.txt", std::ofstream::app);
		if (myfile.is_open())
		{
			// filename    a     b     optimum scale    total computation time (s)
			double fit_contrib = w[n] * std::pow(y[n] - y_tn - m_tn * (x[n] - x_t), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) - w[n] * std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(s*Sig_yn2, 2.0)));

			myfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] << " " << Sig_xn2 << " " << Sig_yn2 << " " << x[n] << " " << y[n] << " " << x_t << " " << fit_contrib << " " << x_tn_vec.size() << " ";
			for (int i = 0; i < x_tn_vec.size(); i++) {
				myfile << x_tn_vec[i] << " ";
			}
			myfile << std::endl;
			myfile.close();
		}
		else std::cout << "Unable to open file";
		*/
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

double TRK::likelihood(std::vector <double> allparams) {
	double L = 1.0;

	double slop_x = allparams[M];
	double slop_y = allparams[M + 1];

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

		L *= w[n]*std::sqrt((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(s*Sig_yn2, 2.0)));
		L *= std::exp(-0.5 * w[n] * (std::pow(y[n] - y_tn - m_tn * (x[n] - x_t), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2)));
	}
	return L;
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
		result.clear();

		double xr1 = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, xg1, xg1 + std::sqrt(Sig_xn2)/10.0);

		if (std::isnan(xr1) && xr1vec.size() >= 1) { //is a NAN if RF did not converge
			return { xr1vec[0] };
		}

		xr1vec.push_back(xr1);

		if (checkcheck) { //different guess gives same root; only one root

			if (std::abs(xr1 - xr1old) < 1e-3) {
				result.push_back(xr1);
				break;
			}
			if (xr1vec.size() >= 3){
				if (std::abs(xr1vec[xr1vec.size() - 1] - xr1vec[xr1vec.size() - 3]) < 1e-3) {

					//std::cout << "oscillation!" << std::endl;

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

		//std::cout << stDevUnweighted(allRoots) << "\n";

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
		} else if (extraRoots.size() == 0 && xr1vec.size() == 1) { //if only have one root (the one initially found)

			//if initial quadratic approximation didn't yield any more guesses, try to find roots with guesses of leftmost and rightmost x values
			double xr_left = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, x_min, x_min - std::sqrt(Sig_xn2) / 10.0);
			double xr_right = twoPointNR(params, x_n, y_n, Sig_xn2, Sig_yn2, x_max, x_max + std::sqrt(Sig_xn2) / 10.0);

			result.push_back(xr1);
			result.push_back(xr_left);
			result.push_back(xr_right);

			break;
		}
		else if (extraRoots.size() == 0 && xr1vec.size() == 2) { //found two roots but can't find a third
			result = xr1vec;
			break;
		}

		itcount += 1;
	}

	//std::cout << itcount << std::endl;

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
double TRK::innerSlopX_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
	s = ss[0];

	std::vector <double> allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess);

	//iterative_allparams_guess = allparams_s;

	printf("%f \t %f \t %f \n", s, allparams_s[M], allparams_s[M + 1]);

	return allparams_s[M];
}

double TRK::innerSlopY_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
	s = ss[0];

	std::vector <double> allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess);

	//iterative_allparams_guess = allparams_s;

	printf("%f \t %f \t %f \n", s, allparams_s[M], allparams_s[M + 1]);

	return allparams_s[M + 1];
}

double TRK::innerR2_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
	s = ss[0];

	whichExtrema = S;
	std::vector <double> allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess);
	whichExtrema = none;

	//iterative_allparams_guess = allparams_s;

	printf("%f \t %f \t %f \n", s, allparams_s[M], allparams_s[M + 1]);

	double R2as = R2TRK_prime_as();
	double R2sb = R2TRK_prime_sb();

	return R2as - R2sb;
}

double TRK::innerR2_iter_Simplex(std::vector <double> ss, std::vector <double> allparams_guess, double s0) {
	s = ss[0];

	whichExtrema = S;
	std::vector <double> allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess);
	whichExtrema = none;

	printf("%f \t %f \t %f \n", s, allparams_s[M], allparams_s[M + 1]);

	//iterative_allparams_guess = allparams_s;

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
			//std::cout << "trial a: " << trial_a << std::endl;

			double slop_trial_a = innerSlopX_Simplex({ trial_a }, iterative_allparams_guess);

			if (slop_trial_a == 0) {
				a = trial_a;
				break;
			}
			else if (slop_trial_a > 0) {
				inc *= 0.5;
			}
		}
	}
	else if (slop_trial_s == 0) {
		a = trial_s;

		double trial_b = trial_s;

		while (true) {
			trial_b += inc;
			//std::cout << "trial b: " << trial_b << std::endl;

			double slop_trial_b = innerSlopX_Simplex({ trial_b }, iterative_allparams_guess);

			if (slop_trial_b > 0) {
				b = trial_b;
				break;
			}
			else if (slop_trial_b == 0) {
				inc *= 0.5;
			}
		}
	}

	//bisection, now that we have brackets [a,b]

	//std::cout << "beginning bisection (slop x) \n";

	double c, slop_c;
	double tol_bisect = 2e-3;
	double tol_brackets = 1e-3;

	while (true) {
		//std::cout << "a = " << a << "\t b = " << b << std::endl;
		c = (a + b) / 2;

		whichExtrema = slopx;
		slop_c = innerSlopX_Simplex({ c }, iterative_allparams_guess);
		whichExtrema = none;

		//std::cout << "c = " << c << "\t slop_x_c = " << slop_c << "\n";

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
	double trial_s = 1.0;
	double slop_trial_s = innerSlopY_Simplex({ trial_s }, iterative_allparams_guess);

	double inc = trial_s * 0.5;

	if (slop_trial_s > 0) {
		a = trial_s;

		double trial_b = trial_s;

		while (true) {
			trial_b += inc;
			//std::cout << "trial b: " << trial_b << std::endl;

			double slop_trial_b = innerSlopY_Simplex({ trial_b }, iterative_allparams_guess);

			if (slop_trial_b == 0) {
				b = trial_b;
				break;
			}
			else if (slop_trial_b > 0) {
				inc *= 0.5;
			}
		}
	}
	else if (slop_trial_s == 0) {
		b = trial_s;

		double trial_a = trial_s;

		while (true) {
			trial_a -= inc;
			//std::cout << "trial a: " << trial_a << std::endl;

			double slop_trial_a = innerSlopY_Simplex({ trial_a }, iterative_allparams_guess);

			if (slop_trial_a > 0) {
				a = trial_a;
				break;
			}
			else if (slop_trial_a == 0) {
				inc *= 0.5;
			}
		}
	}

	//bisection, now that we have brackets [a,b]

	//std::cout << "beginning bisection (slop y) \n";

	double c, slop_c;
	double tol_bisect = 2e-3;
	double tol_brackets = 1e-3;

	while (true) {
		//std::cout << "a = " << a << "\t b = " << b << std::endl;
		c = (a + b) / 2;

		whichExtrema = slopy;
		slop_c = innerSlopY_Simplex({ c }, iterative_allparams_guess);
		whichExtrema = none;

		//std::cout << "c = " << c << "\t slop_y_c = " << slop_c << "\n";

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

	//std::cout << "beginning bisection (optimum scale finder) \n";

	double c, f_c;
	double tol_bisect = 1e-4;
	double tol_brackets = 1e-3;

	while (true) {
		//std::cout << "left = " << left << "\t right = " << right << std::endl;
		c = (left + right) / 2;

		whichExtrema = S;
		f_c = innerR2_Simplex({ c }, iterative_allparams_guess);
		whichExtrema = none;

		//std::cout << "c = " << c << "\t R2(a,s) - R2(s,b)  = " << f_c << "\n";

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

	//std::cout << "beginning bisection (optimum scale finder, iterative step) \n";

	double c, f_c;
	double tol_bisect = 1e-4;
	double tol_brackets = 1e-3;

	while (true) {
		//std::cout << "left = " << left << "\t right = " << right << std::endl;
		c = (left + right) / 2;

		whichExtrema = S;
		f_c = innerR2_iter_Simplex({ c }, iterative_allparams_guess, s0);
		whichExtrema = none;

		//std::cout << "s0 = " << s0 << "\t c = " << c << "\t R2(a,s) - R2(s,b)  = " << f_c << "\n";

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
		printf("next s0: %f", s0);

		s1 = optimize_s_prime_R2(s0);
		if (std::abs(s1-s0) <= tol_scale) {
			tolcheck = true;
		}
		std::printf("new s0: %f", s1);
		s0 = s1;
	}

	return s1;
}

void TRK::optimizeScale() {
	s = 1.0; //initially begin with s = 1

	std::vector <double> scale_extrema;

	std::cout << "minimizing slop x" << std::endl;

	double s_slopx = optimize_s_SlopX(); //computes scale which minizes slop_x, and scale which minimizes slop_y

	std::cout << "minimizing slop y" << std::endl;

	double s_slopy = optimize_s_SlopY();


	scale_extrema.push_back(s_slopx); //scale_extrema = {s_slopx, s_slopy}

	scale_extrema.push_back(s_slopy);

	std::vector <int> sortedindices = getSortedIndices(scale_extrema);

	a = scale_extrema[sortedindices[0]]; //figures out which scale is a, and which is b, as well as storing the associated best-fit parameters and their associated tangent points for those two extreme scales
	b = scale_extrema[sortedindices[1]];

	printf("extrema: \t a \t b:");
	printf(" \t %f \t %f", a, b);

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

	std::cout << "finding optimum scale" << std::endl;

	double s0 = optimize_s0_R2();

	s = s0;

	double s_final = iterateR2_OptimumScale(s0);

	s = s_final;

	std::cout << "optimum s = " << s << std::endl;

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

	std::ofstream myfile("C:\\Users\\nickk124\\Documents\\Reichart Research\\TRK\\TRKresults.txt", std::ofstream::app);
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