#include "pch.h"
#include "TRK.h"

// CONSTRUCTORS

//weighted
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
}

//default
TRK::TRK() {

}

// OTHER ALGORITHMS AND TOOLS
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
	double tol = 1e-3;
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
		//std::cout << xkp1 << std::endl;

		xkm1 = xk;
		xk = xkp1;

		itercount += 1;
	}

	//std::cout << itercount << " iterations." << std::endl;
	return xkp1;
}

std::vector <double> TRK::downhillSimplex(double(TRK::*f)(std::vector <double>), std::vector <double> allparams_guess) {

	double tol = 1e-3;

	int n = allparams_guess.size(); //number of model parameters plus two slop parameters

	double rho = 1.0;
	double chi = 2.0;
	double gamma = 0.5;
	double sigma = 0.5;

	// simplex initialization

	std::vector <double> init_point = allparams_guess;

	std::vector <std::vector <double> > vertices(n + 1, init_point);

	int i = 0;
	for (int j = 1; j < n + 1; j++) { //for each simplex node

		vertices[j][i] = allparams_guess[i] + allparams_guess[i]; //add initial "step size"
		i += 1;
	}

	std::vector <double> result;
	while (true) {
		while (true) {
			// order

			std::vector <int> orderedindices;
			std::vector <double> unorderedEvals; // ( f(x_1), f(x_2), ... f(x_n+1)
			for (int i = 0; i < n + 1; i++) {
				unorderedEvals.push_back((this->*f)(vertices[i]));
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

			double fr = (this->*f)(refpoint);
			double f1 = (this->*f)(vertices[0]);
			double fn = (this->*f)(vertices[n - 1]);

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

				double fe = (this->*f)(exppoint);

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
				double fnp1 = (this->*f)(vertices[n]);

				if (fn <= fr && fr < fnp1) {
					std::vector <double> cpoint;

					for (int i = 0; i < n; i++) {
						cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
					}

					double fc = (this->*f)(cpoint);

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

					double fcc = (this->*f)(ccpoint);

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
		bettervertices.push_back(result);

		vertices = bettervertices;

		for (int i = 0; i < result.size(); i++) {
			std::cout << result[i] << " ";
		}
		std::cout << "                        ";

		for (int i = 0; i < 1; i++) {
			std::cout << modifiedChiSquared(result) << " ";
		}

		std::cout << std::endl << std::endl;

		//test for termination

		std::vector <double> evals;
		for (int i = 0; i < n + 1; i++) {
			evals.push_back((this->*f)(vertices[i]));
		}

		if (stDevUnweighted(evals) < tol) {
			break;
		}

	}
	return vertices[n];
}

std::vector <double> TRK::downhillSimplex(double(*f)(std::vector <double>), std::vector <double> allparams_guess) {

	double tol = 1e-3;

	int n = allparams_guess.size(); //number of model parameters plus two slop parameters

	double rho = 1.0;
	double chi = 2.0;
	double gamma = 0.5;
	double sigma = 0.5;

	// simplex initialization

	std::vector <double> init_point = allparams_guess;

	std::vector <std::vector <double> > vertices(n + 1, init_point);

	int i = 0;
	for (int j = 1; j < n + 1; j++) { //for each simplex node

		vertices[j][i] = allparams_guess[i] + allparams_guess[i]; //add initial "step size"
		i += 1;
	}

	std::vector <double> result;
	while (true) {
		while (true) {
			// order

			std::vector <int> orderedindices;
			std::vector <double> unorderedEvals; // ( f(x_1), f(x_2), ... f(x_n+1)
			for (int i = 0; i < n + 1; i++) {
				unorderedEvals.push_back(f(vertices[i]));
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

			double fr = f(refpoint);
			double f1 = f(vertices[0]);
			double fn = f(vertices[n - 1]);

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

				double fe = f(exppoint);

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
				double fnp1 = f(vertices[n]);

				if (fn <= fr && fr < fnp1) {
					std::vector <double> cpoint;

					for (int i = 0; i < n; i++) {
						cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
					}

					double fc = f(cpoint);

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

					double fcc = f(ccpoint);

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
		bettervertices.push_back(result);

		vertices = bettervertices;

		//test for termination

		std::vector <double> evals;
		for (int i = 0; i < n + 1; i++) {
			evals.push_back(f(vertices[i]));
		}

		if (stDevUnweighted(evals) < tol) {
			break;
		}

	}
	return vertices[n];
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
	x_t_slopx.clear();
	x_t_slopy.clear();

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
		double Sig_xn2 = std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0);
		double Sig_yn2 = std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0);

		std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], Sig_xn2, Sig_yn2, x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points
		double x_t = findBestTangent(params, x[n], y[n], Sig_xn2, Sig_yn2, x_tn_vec);

		double m_tn = dyc(x_t, params);
		double y_tn = yc(x_t, params);

		all_x_t.push_back(x_t);

		sum1 += w[n] * std::pow(y[n] - y_tn - m_tn * (x[n] - x_t), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2);
		sum2 += w[n] * std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(s*Sig_yn2, 2.0)));
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

	iterative_allparams_guess = allparams_s;

	return allparams_s[M];
}

double TRK::innerSlopY_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
	s = ss[0];

	std::vector <double> allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess);

	iterative_allparams_guess = allparams_s;

	return allparams_s[M + 1];
}

double TRK::innerR2_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
	s = ss[0];

	whichExtrema = S;
	std::vector <double> allparams_s = downhillSimplex(&TRK::modifiedChiSquared, allparams_guess);
	whichExtrema = none;

	iterative_allparams_guess = allparams_s;

	return R2TRK_prime_as() - R2TRK_prime_sb();
}

double TRK::optimize_s_SlopX() {

	double tol = 1e-3;

	int n = 1; //only parameter is s

	double rho = 1.0;
	double chi = 2.0;
	double gamma = 0.5;
	double sigma = 0.5;

	// simplex initialization

	std::vector < std::vector <double> > vertices = { {s}, {s*1.1} }; //initial points are at s = 1 and s = 1.1

	iterative_allparams_guess = allparams_guess;

	std::vector <double> result;
	while (true) {
		while (true) {
			// order

			std::vector <int> orderedindices;
			std::vector <double> unorderedEvals; // ( f(x_1), f(x_2), ... f(x_n+1)
			for (int i = 0; i < n + 1; i++) {
				unorderedEvals.push_back(innerSlopX_Simplex(vertices[i], iterative_allparams_guess));
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

			double fr = innerSlopX_Simplex(refpoint, iterative_allparams_guess);
			double f1 = innerSlopX_Simplex(vertices[0], iterative_allparams_guess);
			double fn = innerSlopX_Simplex(vertices[n - 1], iterative_allparams_guess);

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

				double fe = innerSlopX_Simplex(exppoint, iterative_allparams_guess);

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
				double fnp1 = innerSlopX_Simplex(vertices[n], iterative_allparams_guess);

				if (fn <= fr && fr < fnp1) {
					std::vector <double> cpoint;

					for (int i = 0; i < n; i++) {
						cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
					}

					double fc = innerSlopX_Simplex(cpoint, iterative_allparams_guess);;

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

					double fcc = innerSlopX_Simplex(ccpoint, iterative_allparams_guess);;

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
		bettervertices.push_back(result);

		vertices = bettervertices;

		for (int i = 0; i < bettervertices.size(); i++){
			std::cout << bettervertices[i][0] << " ";
		}
		std::cout << "                        ";

		for (int i = 0; i < bettervertices.size(); i++) {
			std::cout << innerSlopX_Simplex(bettervertices[i], iterative_allparams_guess) << " ";
		}

		std::cout << std::endl << std::endl;

		//test for termination

		std::vector <double> evals;
		for (int i = 0; i < n + 1; i++) {
			evals.push_back(innerSlopX_Simplex(vertices[i], iterative_allparams_guess));
		}

		if (stDevUnweighted(evals) < tol) {
			break;
		}

	}
	std::vector <double> optimum_s = vertices[n];


	whichExtrema = slopx;
	innerSlopX_Simplex(optimum_s, iterative_allparams_guess); //runs this with optimum s to store tangest points associated with best fit parameters for this optimum s.
	whichExtrema = none;

	return optimum_s[0];

}

double TRK::optimize_s_SlopY() {

	double tol = 1e-3;

	int n = 1; //only parameter is s

	double rho = 1.0;
	double chi = 2.0;
	double gamma = 0.5;
	double sigma = 0.5;

	// simplex initialization

	std::vector < std::vector <double> > vertices = { {s}, {s*1.1} }; //initial points are at s = 1 and s = 1.1

	iterative_allparams_guess = allparams_guess;

	std::vector <double> result;
	while (true) {
		while (true) {
			// order

			std::vector <int> orderedindices;
			std::vector <double> unorderedEvals; // ( f(x_1), f(x_2), ... f(x_n+1)
			for (int i = 0; i < n + 1; i++) {
				unorderedEvals.push_back(innerSlopY_Simplex(vertices[i], iterative_allparams_guess));
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

			double fr = innerSlopY_Simplex(refpoint, iterative_allparams_guess);
			double f1 = innerSlopY_Simplex(vertices[0], iterative_allparams_guess);
			double fn = innerSlopY_Simplex(vertices[n - 1], iterative_allparams_guess);

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

				double fe = innerSlopY_Simplex(exppoint, iterative_allparams_guess);

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
				double fnp1 = innerSlopY_Simplex(vertices[n], iterative_allparams_guess);

				if (fn <= fr && fr < fnp1) {
					std::vector <double> cpoint;

					for (int i = 0; i < n; i++) {
						cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
					}

					double fc = innerSlopY_Simplex(cpoint, iterative_allparams_guess);;

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

					double fcc = innerSlopY_Simplex(ccpoint, iterative_allparams_guess);;

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
		bettervertices.push_back(result);

		vertices = bettervertices;

		//test for termination

		std::vector <double> evals;
		for (int i = 0; i < n + 1; i++) {
			evals.push_back(innerSlopY_Simplex(vertices[i], iterative_allparams_guess));
		}

		if (stDevUnweighted(evals) < tol) {
			break;
		}

	}
	std::vector <double> optimum_s = vertices[n];


	whichExtrema = slopy;
	innerSlopX_Simplex(optimum_s, iterative_allparams_guess); //runs this with optimum s to store tangest points associated with best fit parameters for this optimum s.
	whichExtrema = none;

	return optimum_s[0];

}

double TRK::optimize_s_R2() {
	
	double tol = 1e-3;

	int n = 1; //only parameter is s

	double rho = 1.0;
	double chi = 2.0;
	double gamma = 0.5;
	double sigma = 0.5;

	// simplex initialization

	std::vector < std::vector <double> > vertices = { {s}, {s*1.1} }; //initial points are at s = 1 and s = 1.1

	iterative_allparams_guess = allparams_guess;

	std::vector <double> result;
	while (true) {
		while (true) {
			// order

			std::vector <int> orderedindices;
			std::vector <double> unorderedEvals; // ( f(x_1), f(x_2), ... f(x_n+1)
			for (int i = 0; i < n + 1; i++) {
				unorderedEvals.push_back(innerR2_Simplex(vertices[i], iterative_allparams_guess));
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

			double fr = innerR2_Simplex(refpoint, iterative_allparams_guess);
			double f1 = innerR2_Simplex(vertices[0], iterative_allparams_guess);
			double fn = innerR2_Simplex(vertices[n - 1], iterative_allparams_guess);

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

				double fe = innerR2_Simplex(exppoint, iterative_allparams_guess);

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
				double fnp1 = innerR2_Simplex(vertices[n], iterative_allparams_guess);

				if (fn <= fr && fr < fnp1) {
					std::vector <double> cpoint;

					for (int i = 0; i < n; i++) {
						cpoint.push_back(centroid[i] + gamma * (refpoint[i] - centroid[i]));
					}

					double fc = innerR2_Simplex(cpoint, iterative_allparams_guess);;

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

					double fcc = innerR2_Simplex(ccpoint, iterative_allparams_guess);;

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
		bettervertices.push_back(result);

		vertices = bettervertices;

		//test for termination

		std::vector <double> evals;
		for (int i = 0; i < n + 1; i++) {
			evals.push_back(innerR2_Simplex(vertices[i], iterative_allparams_guess));
		}

		if (stDevUnweighted(evals) < tol) {
			break;
		}

	}
	std::vector <double> optimum_s = vertices[n];


	whichExtrema = S;
	innerSlopX_Simplex(optimum_s, iterative_allparams_guess); //runs this with optimum s to store tangent points associated with best fit parameters for this optimum s.
	whichExtrema = none;

	return optimum_s[0];

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

void TRK::optimizeScale() {
	double tol = 1e-3;

	s = 1; //initially begin with s = 1
	double s1;

	while (true) { //each loop is for a different s
		std::vector <double> scale_extrema;

		double s_slopx = optimize_s_SlopX(); //computes scale which minizes slop_x, and scale which minimizes slop_y
		double s_slopy = optimize_s_SlopY();


		scale_extrema.push_back(s_slopx); //scale_extrema = {s_slopx, s_slopy}

		scale_extrema.push_back(s_slopy);

		std::vector <int> sortedindices = getSortedIndices(scale_extrema);

		double a = scale_extrema[sortedindices[0]]; //figures out which scale is a, and which is b, as well as storing the associated best-fit parameters and their associated tangent points for those two extreme scales
		double b = scale_extrema[sortedindices[1]];

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
		double s1 = optimize_s_R2();

		//makes new s the old one

		if (std::abs(s - s1) < tol) {
			s = s1;
			break;
		}
		s = s1;
	}
	return;
}