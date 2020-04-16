#include "TRK.h"


double TRK::pivot = 0.0;
double TRK::pivot2 = 1.0;
double TRK::covid_y12 = 2.55;
bool TRK::covid_fitInLogSpace = false;
double TRK::covid_s = 1.0;

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
    
    this->bigM = allparams_guess.size();

	this->whichExtrema = none;

	this->s = 1.0;

	getDataWidth();

	this->hasPriors = false;
    
    guessMCMCDeltas();
    
    checkZeroErrorBars();

	std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

	allparams_sigmas_guess.push_back(slop_x_sigma_guess);
	allparams_sigmas_guess.push_back(slop_y_sigma_guess);

	this->allparams_sigmas_guess = allparams_sigmas_guess;
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
    
    this->bigM = allparams_guess.size();

	this->whichExtrema = none;

	this->s = 1.0;

	getDataWidth();

	this->hasPriors = false;
    
    guessMCMCDeltas();
    
    checkZeroErrorBars();

	std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

	allparams_sigmas_guess.push_back(slop_x_sigma_guess);
	allparams_sigmas_guess.push_back(slop_y_sigma_guess);

	this->allparams_sigmas_guess = allparams_sigmas_guess;
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
    
    this->bigM = allparams_guess.size();

	this->whichExtrema = none;

	this->s = 1.0;

	getDataWidth();

	this->priorsObject = priorsObject;

	this->hasPriors = true;
    
    guessMCMCDeltas();
    
    checkZeroErrorBars();

	std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

	allparams_sigmas_guess.push_back(slop_x_sigma_guess);
	allparams_sigmas_guess.push_back(slop_y_sigma_guess);

	this->allparams_sigmas_guess = allparams_sigmas_guess;
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
    
    this->bigM = allparams_guess.size();

	this->whichExtrema = none;

	this->s = 1.0;

	getDataWidth();

	this->priorsObject = priorsObject;

	this->hasPriors = true;
    
    guessMCMCDeltas();
    
    checkZeroErrorBars();

	std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

	allparams_sigmas_guess.push_back(slop_x_sigma_guess);
	allparams_sigmas_guess.push_back(slop_y_sigma_guess);

	this->allparams_sigmas_guess = allparams_sigmas_guess;
}

// 1D statistic

TRK::TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess) {
    this->do1DFit = true;
    
    this->yc = (*yc);

    this->x = x;
    this->y = y;
    this->w = w;
    this->sy = sy;

    this->params_guess = params_guess;
    this->slop_y_guess = slop_y_guess;

    this->N = x.size();
    this->M = params_guess.size();

    std::vector <double> allparams_guess = params_guess;
    
    allparams_guess.push_back(slop_y_guess);

    this->allparams_guess = allparams_guess;
    
    this->bigM = allparams_guess.size();

    this->whichExtrema = none;

    this->s = 1.0;

    getDataWidth();

    this->hasPriors = false;
    
    guessMCMCDeltas();
    
    checkZeroErrorBars();

    std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

    allparams_sigmas_guess.push_back(slop_y_sigma_guess);

    this->allparams_sigmas_guess = allparams_sigmas_guess;
}
//equal weights/unweighted
TRK::TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess) {
    this->do1DFit = true;
    
    this->yc = (*yc);

    this->x = x;
    this->y = y;
    this->sy = sy;

    this->params_guess = params_guess;
    this->slop_y_guess = slop_y_guess;

    this->N = x.size();
    this->M = params_guess.size();

    std::vector <double> w(N, 1.0); //no weights === equal weights
    this->w = w;

    std::vector <double> allparams_guess = params_guess;

    allparams_guess.push_back(slop_y_guess);

    this->allparams_guess = allparams_guess;
    
    this->bigM = allparams_guess.size();

    this->whichExtrema = none;

    this->s = 1.0;

    getDataWidth();

    this->hasPriors = false;
    
    guessMCMCDeltas();
    
    checkZeroErrorBars();

    std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

    allparams_sigmas_guess.push_back(slop_y_sigma_guess);

    this->allparams_sigmas_guess = allparams_sigmas_guess;
}


//priors:
//weighted
TRK::TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> w, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess, Priors priorsObject) {
    this->do1DFit = true;
    
    this->yc = (*yc);

    this->x = x;
    this->y = y;
    this->w = w;
    this->sy = sy;

    this->params_guess = params_guess;
    this->slop_y_guess = slop_y_guess;

    this->N = x.size();
    this->M = params_guess.size();

    std::vector <double> allparams_guess = params_guess;

    allparams_guess.push_back(slop_y_guess);
    
    this->bigM = allparams_guess.size();

    this->allparams_guess = allparams_guess;

    this->whichExtrema = none;

    this->s = 1.0;

    getDataWidth();

    this->priorsObject = priorsObject;

    this->hasPriors = true;
    
    guessMCMCDeltas();
    
    checkZeroErrorBars();

    std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

    allparams_sigmas_guess.push_back(slop_y_sigma_guess);

    this->allparams_sigmas_guess = allparams_sigmas_guess;
}
//equal weights/unweighted
TRK::TRK(double(*yc)(double, std::vector <double>), std::vector <double> x, std::vector <double> y, std::vector <double> sy, std::vector <double> params_guess, double slop_y_guess, Priors priorsObject) {
    this->do1DFit = true;
    
    this->yc = (*yc);

    this->x = x;
    this->y = y;
    this->sy = sy;

    this->params_guess = params_guess;
    this->slop_y_guess = slop_y_guess;

    this->N = x.size();
    this->M = params_guess.size();

    std::vector <double> w(N, 1.0); //no weights === equal weights
    this->w = w;

    std::vector <double> allparams_guess = params_guess;

    allparams_guess.push_back(slop_y_guess);

    this->allparams_guess = allparams_guess;
    
    this->bigM = allparams_guess.size();

    this->whichExtrema = none;

    this->s = 1.0;

    getDataWidth();

    this->priorsObject = priorsObject;

    this->hasPriors = true;
    
    guessMCMCDeltas();
    
    checkZeroErrorBars();

    std::vector <double> allparams_sigmas_guess = params_sigmas_guess;

    allparams_sigmas_guess.push_back(slop_y_sigma_guess);

    this->allparams_sigmas_guess = allparams_sigmas_guess;
}


//default
TRK::TRK() {

}

// OTHER ALGORITHMS AND TOOLS
std::vector <double> minMax(std::vector <double> vec) {

	// Finding the smallest of all the numbers 
	double min = *std::min_element(std::begin(vec), std::end(vec));
	double max = *std::max_element(std::begin(vec), std::end(vec));

	return { min, max };
}

std::vector <int> argMinMax(std::vector <double> x){
    int argMin = (int)std::distance(x.begin(), std::min_element(x.begin(), x.end()));
    int argMax = (int)std::distance(x.begin(), std::max_element(x.begin(), x.end()));
    return {argMin, argMax};
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
	unsigned long n = nvertices.size();

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
    
    return;
}

void TRK::checkZeroErrorBars(){
    std::vector <int> badInds;
    bool check = false;
    if (do1DFit){
        for (int i = 0; i < N; i++){
            if (sy[i] == 0){
                badInds.push_back(i);
                check = true;
            }
        }
        if (check){
            printf("Warning: Error bars of zero found at data points (indexing starting at 1) of:\n");
            for (int j = 0; j < (int)badInds.size(); j++){
                printf(" %i,", badInds[j]);
            }
            printf("\n Could create likelihood overflow errors if slop also zero/close to zero.");
        }
        
    } else {
        for (int i = 0; i < N; i++){
            if (sx[i] == 0 || sy[i] == 0){
                badInds.push_back(i);
                check = true;
            }
        }
        if (check){
            printf("Warning: Error bars of zero found at data points (indexing starting at 1) of:\n");
            for (int j = 0; j < (int)badInds.size(); j++){
                printf(" %i,", badInds[j]);
            }
            printf("\n Could create likelihood overflow errors if slop also zero/close to zero.");
        }
    }
    return;
}

double TRK::getAverage(std::vector <double> x) {
	double top = 0.0;
	unsigned long N = x.size();

	for (int i = 0; i < N; i++) {
		top += x[i];
	}

	return top / N;
}

double TRK::getAverage(std::vector <double> x, std::vector <double> w) {
    double top = 0.0;
    double bottom = 0.0;
    unsigned long N = x.size();
    
    for (int i = 0; i < N; i++) {
        top += x[i]*w[i];
        bottom += w[i];
    }
    
    return top / bottom;
}

double TRK::min(double a, double b)
{
    return (a < b ? a : b);
}
double TRK::max(double a, double b)
{
    return (a > b ? a : b);
}

std::vector < std::vector <double> > TRK::transpose(std::vector < std::vector <double> > array) { //takes the transpose of some input array
    
    int n = (int)array.size();
    int m = (int)array[0].size();
    
    std::vector<double> innertransposedArray(n, 0.0);
    std::vector <std::vector <double> > transposedArray(m, innertransposedArray);
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            transposedArray[j][i] = array[i][j];
    return transposedArray;
};


bool TRK::isEqual(double x, double y, double maxRelativeError = .00000001, double maxAbsoluteError = DBL_MIN)// .000001; .0000001;.00000001
{
    if (std::abs(x - y) < maxAbsoluteError)
    {
        return true;
    }
    double relativeError = (std::abs(y) > std::abs(x) ? std::abs((x - y) / y) : std::abs((x - y) / x));
    if (relativeError <= maxRelativeError)
    {
        return true;
    }
    return false;
}

double TRK::getMode(int trueCount, std::vector<double> w, std::vector<double> y)
{
    int k, lowerLimit = 0, upperLimit = trueCount - 1, lowerLimitIn = -1, upperLimitIn = -1, size;
    int finalLower = 0;
    int finalUpper = 1;
    double halfWeightSum = 0, sSum, total, minDist = 999999;
    std::vector<double> sVec;
    
    while (lowerLimit != lowerLimitIn || upperLimit != upperLimitIn)
    {
        //std::cout<< lowerLimit << "\t" << upperLimit << "\n";
        lowerLimitIn = lowerLimit;
        upperLimitIn = upperLimit;
        size = upperLimit - lowerLimit + 1;
        minDist = 999999;
        halfWeightSum = 0;
        for (int i = lowerLimit; i < upperLimit + 1; i++)
        {
            halfWeightSum += w[i];
        }
        halfWeightSum *= .5;
        
        sVec.resize(size, 0.0);
        sSum = .5 * w[lowerLimit];
        sVec[0] = sSum;
        for (int i = lowerLimit + 1; i < lowerLimit + size; i++)
        {
            sSum += w[i - 1] * .5 + w[i] * .5;
            sVec[i - lowerLimit] = sSum;
        }
        
        for (size_t i = 0; i < sVec.size(); i++)
        {
            if ((sVec[i] < halfWeightSum) || isEqual(sVec[i],halfWeightSum))
            {
                total = sVec[i] + halfWeightSum;
                k = (int)i; // was 0
                while (k < sVec.size() && ((sVec[k] < total) || isEqual(sVec[k], total)))
                {
                    k++;
                }
                k--;
                total = std::abs(y[k + lowerLimit] - y[i + lowerLimit]);
                
                
                if (isEqual(total,minDist))
                {
                    finalLower = (int)(min(finalLower, i + lowerLimit));
                    finalUpper = (int)(max(finalUpper, k + lowerLimit));
                }
                else if (total < minDist)
                {
                    minDist = total;
                    finalLower = (int)i + lowerLimit;
                    finalUpper = k + lowerLimit;
                }
            }
            if ((sVec[i] > halfWeightSum) || isEqual(sVec[i], halfWeightSum))
            {
                total = sVec[i] - halfWeightSum;
                k = (int)i; // was svec.size() - 1
                while (k > -1 && ((sVec[k] > total) || isEqual(sVec[k], total)))
                {
                    k--;
                }
                k++;
                total = std::abs(y[i + lowerLimit] - y[k + lowerLimit]);
                
                
                if (isEqual(total,minDist))
                {
                    finalLower = (int)(min(finalLower, k + lowerLimit));
                    finalUpper = (int)(max(finalUpper, i + lowerLimit));
                }
                else if (total < minDist)
                {
                    minDist = total;
                    finalLower = k + lowerLimit;
                    finalUpper = (int)i + lowerLimit;
                }
            }
        }
        
        lowerLimit = finalLower;
        upperLimit = finalUpper;
        
        sVec.clear();
    }
    
    std::vector<double> newValues(y.begin() + lowerLimit, y.begin() + upperLimit + 1);
    std::vector<double> newWeights(w.begin() + lowerLimit, w.begin() + upperLimit + 1);
    return getMedian((int)newWeights.size(), newWeights, newValues);
}

// ASYMMETRIC DISTRIBUTION TOOLS

void TRK::checkAsym(){ //checks to see whether any or all of the asymmetric error bar and slop parameters were provided.
    
    // SLOP
    
    bool negXSlop = false;
    bool negYSlop = false;
    
    if (slop_x_minus_guess >= 0){
        negXSlop = true;
    }

    if (slop_y_minus_guess >= 0){
        negYSlop = true;
    }
    
    if (negXSlop && !negYSlop){
        hasAsymSlop = true;
        slop_y_minus_guess = slop_y_guess;
    }
    
    else if (!negXSlop && negYSlop){
        hasAsymSlop = true;
        slop_x_minus_guess = slop_x_guess;
    }
    
    else if (negXSlop && negYSlop){
        hasAsymSlop = true;
    }
    
    
    // ERROR BARS
    
    bool negXEB = false;
    bool negYEB = false;
    
    if (sx_minus.size() == N){
        negXEB = true;
    }
    
    if (sy_minus.size() == N){
       negYEB = true;
    }
    
    
    if (negXEB && !negYEB){
        hasAsymEB = true;
        sy_minus = sy;
    }
    
    else if (!negXEB && negYEB){
        hasAsymEB = true;
        sx_minus = sx;
    }
    
    else if (negXEB && negYEB){
        hasAsymEB = true;
    }
    
    if (hasAsymSlop || hasAsymEB){
        selectedChiSq = &TRK::modifiedChiSquaredAsym;
        selectedLikelihood = &TRK::likelihoodAsym;
    }
    
    if (verboseAsymmetric){
        printf("Asymmetries: slop: %s\tError bars: %s\n", hasAsymSlop ? "true" : "false", hasAsymEB ? "true" : "false");
    }
    
    if (hasAsymSlop){
        allparams_guess.push_back(slop_x_minus_guess);
        allparams_guess.push_back(slop_y_minus_guess);
        
        allparams_sigmas_guess.push_back(slop_x_sigma_guess);
        allparams_sigmas_guess.push_back(slop_y_sigma_guess);
    }
    
    if (do1DFit){
        selectedChiSq = &TRK::regularChiSquaredWSlop;
        selectedLikelihood = &TRK::likelihood1D;
        printf("Running 1D fit.\n");
    }
    
    return;
}

double TRK::dunDxAsym(double mtn, std::vector <double> Sigs2, int quadSig_xn2Ind, int quadSig_yn2Ind, double s){
    double quadSigX2 = Sigs2[quadSig_xn2Ind]; //the correct Sig2 values for the specific quadrant
    double quadSigY2 = Sigs2[quadSig_yn2Ind];

    return (std::pow(mtn, 2.0)*quadSigX2 + quadSigY2) / (std::sqrt(std::pow(mtn*quadSigX2, 2.0) + std::pow(s*quadSigY2, 2.0)));
}

double TRK::cmNorm(double z){
    return (std::sqrt(PI)/2.0) * (1.0 + std::erf(z));
}

double TRK::zAsym(double x, double quadSig_xn2, double quadSig_yn2, double xn_shifted, double yn_shifted, std::vector <double> shifts, double x_tn, double y_tn, double m_tn){ //equation B-6 of thesis
    
    return (quadSig_yn2 * (x - xn_shifted) + std::pow(m_tn, 2.0) * quadSig_xn2 * (x - x_tn - (yn_shifted - y_tn)/m_tn)) / (std::sqrt(quadSig_xn2) * std::sqrt(quadSig_yn2) * std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2));
}

double TRK::pnAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, double x_tn, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn){
    
    // INITIALIZATIONS
    
    double m_tn = dyc(x_tn, params);
    double y_tn = yc(x_tn, params);
    
    double x1 = xn_shifted;
    double x2 = ((yn_shifted - y_tn)/m_tn) + x_tn;
    
    double dunDx = dunDxAsym(m_tn, Sigs2, quadSig_xn2Ind, quadSig_yn2Ind, s);
    double norm = 2.0 / (PI * (std::sqrt(Sigs2[0]) + std::sqrt(Sigs2[2])) * (std::sqrt(Sigs2[1]) + std::sqrt(Sigs2[3]))); //normalization factor, same for all three integrals
    
    // INTEGRALS
    
    double I1 = 0;
    double I2 = 0;
    double I3 = 0;
    
    std::vector <int> I1inds = {0, 0};
    std::vector <int> I2inds = {quadSig_xn2Ind, quadSig_yn2Ind};
    std::vector <int> I3inds = {0, 0};
    
    
    if (I2inds[0] == 0 && I2inds[1] == 1){
        I1inds = {2, 1};
        I3inds = {0, 3};
    } else if (I2inds[0] == 2 && I2inds[1] == 1){
        I1inds = {2, 3};
        I3inds = {0, 1};
    } else if (I2inds[0] == 2 && I2inds[1] == 3){
        I1inds = {2, 1};
        I3inds = {0, 3};
    } else if (I2inds[0] == 0 && I2inds[1] == 3){
        I1inds = {2, 3};
        I3inds = {0, 1};
    }
    
    // Do ``Integrals''
    
    double quadSig_xn2 = Sigs2[I1inds[0]];
    double quadSig_yn2 = Sigs2[I1inds[1]];
    
    double z1 = zAsym(x1, quadSig_xn2, quadSig_yn2, xn_shifted, yn_shifted, shifts, x_tn, y_tn, m_tn);
    
    I1 = dunDx * norm * std::sqrt(quadSig_xn2) * std::sqrt(quadSig_yn2) * cmNorm(z1) * std::exp(-0.5 * std::pow((y_tn - yn_shifted - m_tn*(x_tn - xn_shifted))/std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2), 2.0)) / std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2);
    
    
    
    quadSig_xn2 = Sigs2[I2inds[0]];
    quadSig_yn2 = Sigs2[I2inds[1]];
    
    double z2 = zAsym(x2, quadSig_xn2, quadSig_yn2, xn_shifted, yn_shifted, shifts, x_tn, y_tn, m_tn);
    z1 = zAsym(x1, quadSig_xn2, quadSig_yn2, xn_shifted, yn_shifted, shifts, x_tn, y_tn, m_tn);
    
    I2 = dunDx * norm * std::sqrt(quadSig_xn2) * std::sqrt(quadSig_yn2) * (cmNorm(z2) - cmNorm(z1)) * std::exp(-0.5 * std::pow((y_tn - yn_shifted - m_tn*(x_tn - xn_shifted))/std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2), 2.0)) / std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2);
    
    
    
    quadSig_xn2 = Sigs2[I3inds[0]];
    quadSig_yn2 = Sigs2[I3inds[1]];
    
    z2 = zAsym(x2, quadSig_xn2, quadSig_yn2, xn_shifted, yn_shifted, shifts, x_tn, y_tn, m_tn);
    
    I3 = dunDx * norm * std::sqrt(quadSig_xn2) * std::sqrt(quadSig_yn2) * (1.0 - cmNorm(z2)) * std::exp(-0.5 * std::pow((y_tn - yn_shifted - m_tn*(x_tn - xn_shifted))/std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2), 2.0)) / std::sqrt(std::pow(m_tn, 2.0) * quadSig_xn2 + quadSig_yn2);
    
    
    
    return I1 + I2 + I3;
    
}

double TRK::singlePointLnLAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, double x_tn, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn){
    
    return -2.0 * std::log(pnAsym(params, xn_shifted, yn_shifted, Sigs2, x_tn, quadSig_xn2Ind, quadSig_yn2Ind, shifts, s, wn));
}

double TRK::findBestTangentAsym(std::vector <double> params, double xn_shifted, double yn_shifted, std::vector <double> Sigs2, std::vector <double> x_tn_vec, int quadSig_xn2Ind, int quadSig_yn2Ind, std::vector <double> shifts, double s, double wn) {
    std::vector <double> posts;
    long minindex;

    for (int i = 0; i < x_tn_vec.size(); i++) {
        posts.push_back(singlePointLnLAsym(params, xn_shifted, yn_shifted, Sigs2, x_tn_vec[i], quadSig_xn2Ind, quadSig_yn2Ind, shifts, s, wn));
    }

    std::vector<double>::iterator result = std::min_element(std::begin(posts), std::end(posts));
    minindex = std::distance(std::begin(posts), result);

    return x_tn_vec[minindex];
}


std::vector <double> TRK::getAsymShifts(std::vector <double> allparams, int n){
    double deltayn = 0.0;
    double deltaxn = 0.0;
    
    std::vector <double> slops = {allparams[M], allparams[M+1]};
    std::vector <double> EBs = {sx[n], sy[n]};
    
    if (hasAsymSlop && !hasAsymEB){
        slops.push_back(allparams[M+2]);
        slops.push_back(allparams[M+3]);
        
        EBs = concat(EBs, EBs);
        
    } else if (!hasAsymSlop && hasAsymEB){
        slops = concat(slops, slops);
        
        EBs.push_back(sx_minus[n]);
        EBs.push_back(sy_minus[n]);
        
    } else if (hasAsymSlop && hasAsymEB){
        slops.push_back(allparams[M+2]);
        slops.push_back(allparams[M+3]);
        
        EBs.push_back(sx_minus[n]);
        EBs.push_back(sy_minus[n]);
    }
    
    // Y SHIFT
    
    double sigmaL = minMax({slops[1], slops[3]})[1];
    double sigmaS = minMax({slops[1], slops[3]})[0];
    double sigmanL = minMax({EBs[1], EBs[3]})[1];
    double sigmanS = minMax({EBs[1], EBs[3]})[0];
    double sigmaMax = minMax({sigmaL, sigmanL})[1];
    
    
    double xi = sigmaS/sigmaL + sigmanS / sigmanL;
    double eta = sigmanL < sigmaL ? sigmanS/sigmanL - sigmaS/sigmaL : sigmaS/sigmaL - sigmanS/sigmanL;
    double r = minMax({sigmaL, sigmanL})[0] / minMax({sigmaL, sigmanL})[1];
    
    double xip = xi <= 1 ? xi : 2.0 - xi;
    double etap = xip == 0 ? 0 : 2.0 * xip * std::pow(0.5*eta/xip + 1.0, std::pow(r, -0.4087)) - xip;
    
    double Nr = -0.5326*std::pow(r,2.0) + 1.5307*r + 0.0019;
    double fXi = xi <= 1 ? 0.2454*std::pow(xi, -1.1452) : 0.2454*std::pow(xi, -0.5203);
    double gEtaP = std::pow(etap,2.0);
    double hXi = -0.042*std::pow(xi,2.0) - 0.1602 * xi + 0.4884;
    
    double deltastr = sigmaMax * Nr * (fXi * gEtaP + hXi);
    
    int i = 1;
    if (slops[1] == slops[3] || EBs[1] == EBs[3]){ //one of the dists is symmetric
        if (sigmanL == EBs[1] || sigmaL == slops[3]){
            i = 1;
        } else if (sigmaMax == EBs[3] || sigmaMax == slops[1]){
            i = -1;
        }
        
        deltayn = i * deltastr;
        
    } else if ((sigmaL == slops[3] && sigmanL == EBs[1]) || (sigmaL = slops[1] && sigmanL == EBs[3])){ //both asymm first case
        if (sigmaMax == EBs[1] || sigmaMax == slops[3]){
            i = 1;
        } else if (sigmaMax == EBs[3] || sigmaMax == slops[1]){
            i = -1;
        }
        
        deltayn = i * deltastr;
        
    } else if ((sigmaL == slops[1] && sigmanL == EBs[1]) || (sigmaL = slops[3] && sigmanL == EBs[3])){ //both asymm second case
        if (sigmaMax == EBs[1] || sigmaMax == slops[3]){
            i = 1;
        } else if (sigmaMax == EBs[3] || sigmaMax == slops[1]){
            i = -1;
        }
        
        double pwr = eta <= 1 ? 0.7413 : -0.1268;
        deltayn = i * deltastr * std::sin(PI/2.0 * etap/xip) * std::pow(eta, pwr);
        
//        if (EBs[3] == slops[3] && EBs[1] > slops[1] && EBs[1] < EBs[3]){ //footnote case
//            deltayn *= -1;
//        }
    }
    
    // X SHIFT
    
    sigmaL = minMax({slops[0], slops[2]})[0];
    sigmaS = minMax({slops[0], slops[2]})[0];
    sigmanL = minMax({EBs[0], EBs[2]})[0];
    sigmanS = minMax({EBs[0], EBs[2]})[0];
    sigmaMax = minMax({sigmaL, sigmanL})[0];
    
    
    xi = sigmaS/sigmaL + sigmanS / sigmanL;
    eta = sigmanL < sigmaL ? sigmanS/sigmanL - sigmaS/sigmaL : sigmaS/sigmaL - sigmanS/sigmanL;
    r = minMax({sigmaL, sigmanL})[0] / minMax({sigmaL, sigmanL})[0];
    
    xip = xi <= 1 ? xi : 2.0 - xi;
    etap = xip == 0 ? 0 : 2.0 * xip * std::pow(0.5*eta/xip + 1.0, std::pow(r, -0.4087)) - xip;
    
    Nr = -0.5326*std::pow(r,2.0) + 1.5307*r + 0.0019;
    fXi = xi <= 1 ? 0.2454*std::pow(xi, -1.1452) : 0.2454*std::pow(xi, -0.5203);
    gEtaP = std::pow(etap,2.0);
    hXi = -0.042*std::pow(xi,2.0) - 0.1602 * xi + 0.4884;
    
    deltastr = sigmaMax * Nr * (fXi * gEtaP + hXi);
    
    if (slops[0] == slops[2] || EBs[0] == EBs[2]){ //one of the dists is symmetric
        if (sigmanL == EBs[0] || sigmaL == slops[2]){
            i = 1;
        } else if (sigmaMax == EBs[2] || sigmaMax == slops[0]){
            i = -1;
        }
        
        deltaxn = i * deltastr;
        
    } else if ((sigmaL == slops[2] && sigmanL == EBs[0]) || (sigmaL = slops[0] && sigmanL == EBs[2])){ //both asymm first case
        if (sigmaMax == EBs[0] || sigmaMax == slops[2]){
            i = 1;
        } else if (sigmaMax == EBs[2] || sigmaMax == slops[0]){
            i = -1;
        }
        
        deltaxn = i * deltastr;
        
    } else if ((sigmaL == slops[0] && sigmanL == EBs[0]) || (sigmaL = slops[2] && sigmanL == EBs[2])){ //both asymm second case
        if (sigmaMax == EBs[0] || sigmaMax == slops[2]){
            i = 1;
        } else if (sigmaMax == EBs[2] || sigmaMax == slops[0]){
            i = -1;
        }
        
        double pwr = eta <= 1 ? 0.7413 : -0.1268;
        deltaxn = i * deltastr * std::sin(PI/2.0 * etap/xip) * std::pow(eta, pwr);
        
//        if (EBs[2] == slops[2] && EBs[0] > slops[0] && EBs[0] < EBs[2]){ //footnote case
//            deltaxn *= -1;
//        }
    }
    
    
    return {deltaxn, deltayn};
}

std::vector <double> TRK::getAsymSigs2(std::vector <double> allparams, int n){
    std::vector <double> Sigs2(4, 0.0);
    std::vector <double> slops = {allparams[M], allparams[M+1]};
    std::vector <double> EBs = {sx[n], sy[n]};
    
    if (hasAsymSlop && !hasAsymEB){
        slops = concat({allparams[M+2], allparams[M+3]}, slops);
        
        EBs = concat(EBs, EBs);
        
    } else if (!hasAsymSlop && hasAsymEB){
        slops = concat(slops, slops);
        
        EBs.push_back(sx_minus[n]);
        EBs.push_back(sy_minus[n]);
        
    } else if (hasAsymSlop && hasAsymEB){
        slops = concat({allparams[M+2], allparams[M+3]}, slops);
        
        EBs.push_back(sx_minus[n]);
        EBs.push_back(sy_minus[n]);
    }
    
    for (int i = 0; i < 4; i++){
        Sigs2[i] = std::pow(slops[i], 2.0) + std::pow(EBs[i], 2.0);
    }
        
    return Sigs2;
}

std::vector <double> TRK::tangentParallelAsym(std::vector<double> allparams, int n, double s) {
    std::vector <double> params;
    for (int i = 0; i < M; i++) {
        params.push_back(allparams[i]);
    }
    
    // shift centroid
    std::vector <double> shifts = getAsymShifts(allparams, n);
    double xn_shifted = x[n] + shifts[0];
    double yn_shifted = y[n] + shifts[1];
    
    
    // choose correct Sigmas
    std::vector <double> Sigs2 = getAsymSigs2(allparams, n);
    double quadSig_xn2 = 0.0;
    double quadSig_yn2 = 0.0;
    int quadSig_xn2Ind = -1;
    int quadSig_yn2Ind = -1;
    
    double yCxn = yc(xn_shifted, params);
    double dyCxn = dyc(xn_shifted, params);
    int quadrant = 0;
    
    if (yCxn > yn_shifted){
        if (dyCxn < 0){ // Quadrant I
            quadSig_xn2 = Sigs2[0];
            quadSig_yn2 = Sigs2[1];
            quadrant = 1;
            quadSig_xn2Ind = 0;
            quadSig_yn2Ind = 1;
        } else if (dyCxn >= 0){ // QII
            quadSig_xn2 = Sigs2[2];
            quadSig_yn2 = Sigs2[1];
            quadrant = 2;
            quadSig_xn2Ind = 2;
            quadSig_yn2Ind = 1;
        }
    } else if (yCxn <= yn_shifted){
        if (dyCxn < 0){ // QIII
            quadSig_xn2 = Sigs2[2];
            quadSig_yn2 = Sigs2[3];
            quadrant = 3;
            quadSig_xn2Ind = 2;
            quadSig_yn2Ind = 3;
        } else if (dyCxn >= 0){ // QIV
            quadSig_xn2 = Sigs2[0];
            quadSig_yn2 = Sigs2[3];
            quadrant = 4;
            quadSig_xn2Ind = 0;
            quadSig_yn2Ind = 3;
        }
    }
    
    // find tangent point(s)
    
    std::vector <double> x_tn_vec = tangentsFinder(params, xn_shifted, yn_shifted, quadSig_xn2, quadSig_yn2, xn_shifted); // we use x_n as the initial guess for this. gives the three closest tangest points

    double x_t = findBestTangentAsym(params, xn_shifted, yn_shifted, Sigs2, x_tn_vec, quadSig_xn2Ind, quadSig_yn2Ind, shifts, s, w[n]);

    double subsum = std::log(pnAsym(params, xn_shifted, yn_shifted, Sigs2, x_t, quadSig_xn2Ind, quadSig_yn2Ind, shifts, s, w[n]));

    return {x_t, subsum}; // returns x_t_n and ln(p_n)
}

double TRK::modifiedChiSquaredAsym(std::vector <double> allparams, double s)
{
    std::vector <double> all_x_t(N, 0.0);

    double sum = 0.0;

    std::vector <double> params;

    for (int i = 0; i < M; i++) {
        params.push_back(allparams[i]);
    }

    if (cpp17MultiThread) {

//        std::vector <int> nn;
//
//        for (int n = 0; n < N; n++) {
//            nn.push_back(n);
//        }
//
//        for (int n = 0; n < N; n++) {
//            SigXVec.push_back(std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0));
//            SigYVec.push_back(std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0));
//        }
//
//        std::for_each( //parallel tangent point finding
//            std::execution::par_unseq,
//            nn.begin(),
//            nn.end(),
//            [&](auto&& n)
//        {
//            std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], SigXVec[n], SigYVec[n], x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points
//
//            double x_t = findBestTangent(params, x[n], y[n], SigXVec[n], SigYVec[n], x_tn_vec);
//
//            all_x_t[n] = x_t;
//        });
//
//        for (int n = 0; n < N; n++) {
//
//            double m_tn = dyc(all_x_t[n], params);
//            double y_tn = yc(all_x_t[n], params);
//
//            sum1 += w[n] * std::pow(y[n] - y_tn - m_tn * (x[n] - all_x_t[n]), 2.0) / (std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n]);
//            sum2 += w[n] * std::log((std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n]) / (std::pow(m_tn*SigXVec[n], 2.0) + std::pow(s*SigYVec[n], 2.0)));
//        }
    } else if (openMPMultiThread && !cpp17MultiThread) {
        //clock_t time = startTimer();

        #pragma omp parallel for num_threads(maxThreads)
        for (int i = 0; i < N; i++)
        {
            std::vector <double> results;
            results = tangentParallelAsym(allparams, i, s); //pointer to fn run through MT, arguments to fn
            all_x_t[i] = results[0];
            sum += results[1];
        }

        //double sec_elapsed = secElapsed(time);

    } else if (cpp11MultiThread && !cpp17MultiThread) {
        //cpp11 multithreading

        int counter = 0, completedThreads = 0, liveThreads = 0;
        std::vector<double> results;
        std::vector< std::future < std::vector < double > > > futureVec;
        futureVec.resize(N);

        for (int i = 0; i < N; i++)
        {
            futureVec[i] = std::async(std::launch::async, &TRK::tangentParallelAsym, this, allparams, i, s); //pointer to fn run through MT, arguments to fn
            counter++;
            liveThreads++;

            if (liveThreads >= maxThreads)
            {
                for (int i = completedThreads; i < counter; i++)
                {
                    results = futureVec[i].get();
                    all_x_t[i] = results[0];
                    sum += results[1];
                }
                completedThreads += liveThreads;
                liveThreads = 0;
            }
        }
        for (int i = completedThreads; i < N; i++)
        {
            results = futureVec[i].get();
            all_x_t[i] = results[0];
            sum += results[1];
        }

        /*std::vector <std::thread> ths;
        for (int n = 0; n < N; n++) {
            ths.push_back(std::thread(&TRK::tangentParallel, this, params, slop_x, slop_y, n));
        }
        for (auto& th : ths) {
            th.join();
        }*/
        
    }
    else {
        for (int i = 0; i < N; i++)
        {
            std::vector <double> results;
            results = tangentParallelAsym(allparams, i, s); //pointer to fn run through MT, arguments to fn
            all_x_t[i] = results[0];
            sum += results[1];
        }
    }

    switch (whichExtrema) {
        case none:
            break;
        case S:
            x_t_s = all_x_t;
            params_s = params;
            break;
        default:
            break;
    }

    switch (whichExtremaX) {
    case none:
        break;
    case slopx:
        x_t_slopx = all_x_t;
        params_slopx = params;
        break;
    default:
        break;
    }

    switch (whichExtremaY) {
    case none:
        break;
    case slopy:
        x_t_slopy = all_x_t;
        params_slopy = params;
        break;
    default:
        break;
    }

    return -2.0 * sum;
}

// Asymm MCMC tools

double TRK::likelihoodAsym(std::vector <double> allparams) {
    std::vector <double> all_x_t(N, 0.0);
    double L = 1.0;
    if (useLogPosterior){
        L = 0.0;
    }

    std::vector <double> params;

    for (int i = 0; i < M; i++) {
        params.push_back(allparams[i]);
    }

    if (cpp17MultiThread) {

//        std::vector <int> nn;
//
//        for (int n = 0; n < N; n++) {
//            nn.push_back(n);
//        }
//
//        for (int n = 0; n < N; n++) {
//            SigXVec.push_back(std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0));
//            SigYVec.push_back(std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0));
//        }
//
//        std::for_each( //parallel tangent point finding
//            std::execution::par_unseq,
//            nn.begin(),
//            nn.end(),
//            [&](auto&& n)
//        {
//            std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], SigXVec[n], SigYVec[n], x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points
//
//            double x_t = findBestTangent(params, x[n], y[n], SigXVec[n], SigYVec[n], x_tn_vec);
//
//            all_x_t[n] = x_t;
//        });
//
//        for (int n = 0; n < N; n++) {
//            double m_tn = dyc(all_x_t[n], params);
//            double y_tn = yc(all_x_t[n], params);
//
//            L *= w[n] * std::sqrt((std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n]) / (std::pow(m_tn*SigXVec[n], 2.0) + std::pow(s*SigYVec[n], 2.0)));
//            L *= std::exp(-0.5 * w[n] * (std::pow(y[n] - y_tn - m_tn * (x[n] - all_x_t[n]), 2.0) / (std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n])));
//        }
    } else if (openMPMultiThread && !cpp17MultiThread) {
        //clock_t time = startTimer();

        #pragma omp parallel for num_threads(maxThreads)
        for (int i = 0; i < N; i++)
        {
            std::vector <double> results;
            results = tangentParallelLikelihoodAsym(allparams, i); //pointer to fn run through MT, arguments to fn
            all_x_t[i] = results[0];
            if (useLogPosterior){
                L += std::log(results[1]);
            } else {
                L *= results[1];
            }
        }

        //double sec_elapsed = secElapsed(time);

    }
    else if (cpp11MultiThread && !cpp17MultiThread) {
        //cpp11 multithreading

        int counter = 0, completedThreads = 0, liveThreads = 0;
        std::vector<double> results;
        std::vector< std::future < std::vector < double > > > futureVec;
        futureVec.resize(N);

        for (int i = 0; i < N; i++)
        {
            futureVec[i] = std::async(std::launch::async, &TRK::tangentParallelLikelihoodAsym, this, allparams, i); //pointer to fn run through MT, arguments to fn
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
            if (useLogPosterior){
                L += std::log(results[1]);
            } else {
                L *= results[1];
            }
        }
    }
    else {
        for (int i = 0; i < N; i++)
        {
            std::vector <double> results;
            results = tangentParallelLikelihoodAsym(allparams, i); //pointer to fn run through MT, arguments to fn
            all_x_t[i] = results[0];
            if (useLogPosterior){
                L += std::log(results[1]);
            } else {
                L *= results[1];
            }
        }
    }
    return L;
}

std::vector <double> TRK::tangentParallelLikelihoodAsym(std::vector<double> allparams, int n) {
    std::vector <double> params;
    for (int i = 0; i < M; i++) {
        params.push_back(allparams[i]);
    }
    
    // shift centroid
    std::vector <double> shifts = getAsymShifts(allparams, n);
    double xn_shifted = x[n] + shifts[0];
    double yn_shifted = y[n] + shifts[1];
    
    
    // choose correct Sigmas
    std::vector <double> Sigs2 = getAsymSigs2(allparams, n);
    double quadSig_xn2 = 0.0;
    double quadSig_yn2 = 0.0;
    int quadSig_xn2Ind = -1;
    int quadSig_yn2Ind = -1;
    
    double yCxn = yc(xn_shifted, params);
    double dyCxn = dyc(xn_shifted, params);
    int quadrant = 0;
    
    if (yCxn > yn_shifted){
        if (dyCxn < 0){ // Quadrant I
            quadSig_xn2 = Sigs2[0];
            quadSig_yn2 = Sigs2[1];
            quadrant = 1;
            quadSig_xn2Ind = 0;
            quadSig_yn2Ind = 1;
        } else if (dyCxn >= 0){ // QII
            quadSig_xn2 = Sigs2[2];
            quadSig_yn2 = Sigs2[1];
            quadrant = 2;
            quadSig_xn2Ind = 2;
            quadSig_yn2Ind = 1;
        }
    } else if (yCxn <= yn_shifted){
        if (dyCxn < 0){ // QIII
            quadSig_xn2 = Sigs2[2];
            quadSig_yn2 = Sigs2[3];
            quadrant = 3;
            quadSig_xn2Ind = 2;
            quadSig_yn2Ind = 3;
        } else if (dyCxn >= 0){ // QIV
            quadSig_xn2 = Sigs2[0];
            quadSig_yn2 = Sigs2[3];
            quadrant = 4;
            quadSig_xn2Ind = 0;
            quadSig_yn2Ind = 3;
        }
    }
    
    // find tangent point(s)
    
    std::vector <double> x_tn_vec = tangentsFinder(params, xn_shifted, yn_shifted, quadSig_xn2, quadSig_yn2, xn_shifted); // we use x_n as the initial guess for this. gives the three closest tangest points

    double x_t = findBestTangentAsym(params, xn_shifted, yn_shifted, Sigs2, x_tn_vec, quadSig_xn2Ind, quadSig_yn2Ind, shifts, s, w[n]);

    double l = pnAsym(params, xn_shifted, yn_shifted, Sigs2, x_t, quadSig_xn2Ind, quadSig_yn2Ind, shifts, s, w[n]);

    return { x_t, l};
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

            double c, f_c, f_left;
            double left = 0.0;
            double right = 1.0;
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
    
    if (do1DFit){
        if (std::abs(vertex[M]) <= pegToZeroTol) {
            vertex[M] = 0;
        }
    } else {
        if (std::abs(vertex[M]) <= pegToZeroTol) {
            vertex[M] = 0;
        }
        if (std::abs(vertex[M+1]) <= pegToZeroTol) {
            vertex[M+1] = 0;
        }
        
        if (hasAsymSlop){
            if (std::abs(vertex[M+2]) <= pegToZeroTol) {
                vertex[M+2] = 0;
            }
            if (std::abs(vertex[M+3]) <= pegToZeroTol) {
                vertex[M+3] = 0;
            }
        }
    }
	
	return vertex;
}

std::vector <double> TRK::avoidNegativeSlop(std::vector <double> vertex, unsigned long n) {

    int K = 2;
    
    if (hasAsymSlop){
        K = 4;
    }
    
	for (int k = 0; k < K; k++) {
		if (vertex[n - 1 - k] < 0) {
			vertex[n - 1 - k] = std::abs(vertex[n - 1 - k]);
		}
	}
	

	return vertex;
}

double TRK::evalWPriors(double(TRK::*f)(std::vector <double>, double), std::vector <double> vertex, double s) {
	if (hasPriors) {
		switch (priorsObject.priorType) {
            case CUSTOM:
                break;
            case GAUSSIAN:
                break;
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

	unsigned long n = bigM; //number of model parameters plus two slop parameters
    
    if (hasAsymSlop){
        n += 2;
    }

    double rho = 1.0; //reflection
    double chi = 2.0; //expansion
    double gamma = 0.5; //contraction
    double sigma = 0.5; //shrinkage
    
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
		
        if (showSimplexSteps){
            std::cout << "chi-square parameters at s = " << s << " ";
            for (int i = 0; i < result.size(); i++) {
                std::cout << result[i] << " ";
            }
            
            double fitness = evalWPriors(f, result, s);
            std::cout << "fitness = " << fitness << "\n";
            
        }
        
        if (whichExtrema == S or whichExtrema == none){
            double fitness = evalWPriors(f, result, s);
            results.chisquared = fitness;
//            printf("\n\nfinal fitness = %f\n\n", fitness);
        }
		
		
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
}

std::vector <double> TRK::tangentCubicSolver(double A, double B, double C, double D) {
	//cubic solver for three real and distinct roots
	double a1 = B / A;
	double a2 = C / A;
	double a3 = D / A;

	double Q = (a1 * a1 - 3.0 * a2) / 9.0;
	double R = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
	double Qc = std::pow(Q, 3.0);

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

double TRK::stretch_pdf(double z, double a) {
    return z >= 1/a && z <= a ? 1.0/std::sqrt(z) : 0.0;
}

double TRK::singlePointLnL(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double x_tn, double s) {
	double m_tn = dyc(x_tn, params);
	double y_tn = yc(x_tn, params);
	
	return std::pow(y_n - y_tn - m_tn * (x_n - x_tn), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) - std::log((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(s*Sig_yn2, 2.0)));
}

std::vector <double> TRK::tangentParallel(std::vector<double> params, double slop_x, double slop_y, int n, double s) {
	double Sig_xn2 = std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0);
	double Sig_yn2 = std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0);

	std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], Sig_xn2, Sig_yn2, x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points

	double x_t = findBestTangent(params, x[n], y[n], Sig_xn2, Sig_yn2, x_tn_vec, s);

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

//        std::vector <int> nn;
//
//        for (int n = 0; n < N; n++) {
//            nn.push_back(n);
//        }
//
//        for (int n = 0; n < N; n++) {
//            SigXVec.push_back(std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0));
//            SigYVec.push_back(std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0));
//        }
//
//        std::for_each( //parallel tangent point finding
//            std::execution::par_unseq,
//            nn.begin(),
//            nn.end(),
//            [&](auto&& n)
//        {
//            std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], SigXVec[n], SigYVec[n], x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points
//
//            double x_t = findBestTangent(params, x[n], y[n], SigXVec[n], SigYVec[n], x_tn_vec);
//
//            all_x_t[n] = x_t;
//        });
//
//        for (int n = 0; n < N; n++) {
//
//            double m_tn = dyc(all_x_t[n], params);
//            double y_tn = yc(all_x_t[n], params);
//
//            sum1 += w[n] * std::pow(y[n] - y_tn - m_tn * (x[n] - all_x_t[n]), 2.0) / (std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n]);
//            sum2 += w[n] * std::log((std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n]) / (std::pow(m_tn*SigXVec[n], 2.0) + std::pow(s*SigYVec[n], 2.0)));
//        }
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

	} else if (cpp11MultiThread && !cpp17MultiThread) {
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
	else {
		for (int i = 0; i < N; i++)
		{
			std::vector <double> results;
			results = tangentParallel(params, slop_x, slop_y, i, s); //pointer to fn run through MT, arguments to fn
			all_x_t[i] = results[0];
			sum1 += results[1];
			sum2 += results[2];
		}
	}

	switch (whichExtrema) {
		case none:
			break;
		case S:
			x_t_s = all_x_t;
			params_s = params;
			break;
		default:
			break;
	}

	switch (whichExtremaX) {
	case none:
		break;
	case slopx:
		x_t_slopx = all_x_t;
		params_slopx = params;
		break;
	default:
		break;
	}

	switch (whichExtremaY) {
	case none:
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
	unsigned long N = y.size();

	double sum = 0.0;

	for (int i = 0; i < N; i++) {
		sum += w[i] * std::pow(y[i] - (*yc)(x[i], params), 2.0);
	}
	return sum;
}

double TRK::regularChiSquaredWSlop(std::vector <double> allparams, double s) {
    unsigned long N = y.size();

    double sum = 0.0;
    
    double sigma = allparams[(int) allparams.size() - 1];
    
    std::vector <double> params;

    for (int i = 0; i < M; i++) {
        params.push_back(allparams[i]);
    }

    for (int i = 0; i < N; i++) {
        sum += w[i]* std::pow((*yc)(x[i], params) - y[i], 2)/(std::pow(sy[i], 2.0) + std::pow(sigma, 2.0)) + 2.0*std::log(std::pow(std::pow(sy[i], 2.0) + std::pow(sigma, 2.0), w[i]/2.0));
//        L *= std::exp(-0.5 * w[i]* std::pow(((*yc)(x[i], params) - y[i])/std::sqrt(std::pow(sy[i], 2.0) + std::pow(sigma, 2.0)), 2.0)) / std::pow(std::pow(sy[i], 2.0) + std::pow(sigma, 2.0), w[i]/2.0);
    }
    
    return sum;
}

std::vector <double> TRK::tangentParallelLikelihood(std::vector<double> params, double slop_x, double slop_y, int n) {
	double Sig_xn2 = std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0);
	double Sig_yn2 = std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0);

	std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], Sig_xn2, Sig_yn2, x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points

	double x_t = findBestTangent(params, x[n], y[n], Sig_xn2, Sig_yn2, x_tn_vec, s);

	double m_tn = dyc(x_t, params);
	double y_tn = yc(x_t, params);

	double l = std::pow((std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2) / (std::pow(m_tn*Sig_xn2, 2.0) + std::pow(s*Sig_yn2, 2.0)), w[n]/2.0);
	l *= std::exp(-0.5 * w[n] * (std::pow(y[n] - y_tn - m_tn * (x[n] - x_t), 2.0) / (std::pow(m_tn, 2.0)*Sig_xn2 + Sig_yn2)));
    
//    printf("%f\n",l);

	return { x_t, l};
}

double TRK::likelihood(std::vector <double> allparams) {
	std::vector <double> SigXVec, SigYVec;
	std::vector <double> all_x_t(N, 0.0);
	double L = 1.0;
    
    if (useLogPosterior){
        L = 0.0;
    }

	double slop_x = allparams[M];
	double slop_y = allparams[M + 1];

	std::vector <double> params;

	for (int i = 0; i < M; i++) {
		params.push_back(allparams[i]);
	}

	if (cpp17MultiThread) {

//        std::vector <int> nn;
//
//        for (int n = 0; n < N; n++) {
//            nn.push_back(n);
//        }
//
//        for (int n = 0; n < N; n++) {
//            SigXVec.push_back(std::pow(sx[n], 2.0) + std::pow(slop_x, 2.0));
//            SigYVec.push_back(std::pow(sy[n], 2.0) + std::pow(slop_y, 2.0));
//        }
//
//        std::for_each( //parallel tangent point finding
//            std::execution::par_unseq,
//            nn.begin(),
//            nn.end(),
//            [&](auto&& n)
//        {
//            std::vector <double> x_tn_vec = tangentsFinder(params, x[n], y[n], SigXVec[n], SigYVec[n], x[n]); // we use x_n as the initial guess for this. gives the three closest tangest points
//
//            double x_t = findBestTangent(params, x[n], y[n], SigXVec[n], SigYVec[n], x_tn_vec);
//
//            all_x_t[n] = x_t;
//        });
//
//        for (int n = 0; n < N; n++) {
//            double m_tn = dyc(all_x_t[n], params);
//            double y_tn = yc(all_x_t[n], params);
//
//            L *= w[n] * std::sqrt((std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n]) / (std::pow(m_tn*SigXVec[n], 2.0) + std::pow(s*SigYVec[n], 2.0)));
//            L *= std::exp(-0.5 * w[n] * (std::pow(y[n] - y_tn - m_tn * (x[n] - all_x_t[n]), 2.0) / (std::pow(m_tn, 2.0)*SigXVec[n] + SigYVec[n])));
//        }
	} else if (openMPMultiThread && !cpp17MultiThread) {
		//clock_t time = startTimer();

		#pragma omp parallel for num_threads(maxThreads)
		for (int i = 0; i < N; i++)
		{
			std::vector <double> results;
			results = tangentParallelLikelihood(params, slop_x, slop_y, i); //pointer to fn run through MT, arguments to fn
			all_x_t[i] = results[0];
            
            if (useLogPosterior){
                L += std::log(results[1]);
            } else {
                L *= results[1];
            }
		}

		//double sec_elapsed = secElapsed(time);

	}
	else if (cpp11MultiThread && !cpp17MultiThread) {
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
			if (useLogPosterior){
                L += std::log(results[1]);
            } else {
                L *= results[1];
            }
		}
	}
	else {
		for (int i = 0; i < N; i++)
		{
			std::vector <double> results;
			results = tangentParallelLikelihood(params, slop_x, slop_y, i); //pointer to fn run through MT, arguments to fn
			all_x_t[i] = results[0];
			if (useLogPosterior){
                L += std::log(results[1]);
            } else {
                L *= results[1];
            }
            
//            if (isinf(L)){
//                printf("Inf\n");
//            }
//            printf("%f\n",L);
		}
        
	}
	return L; // returns log L = logL1 + logL2 + ... given L = L1*L2*L3... if useLogPosterior == true
}

double TRK::likelihood1D(std::vector <double> allparams) {
//    printf("Notice: 1D likelihood (for testing) not currently configured to work with weights.\n");
    
    double L = 1.0;
    
    if (useLogPosterior){
        L = 0.0;
    }
    
    unsigned long N = y.size();

//    double L = 1.0/(2.0*PI); constant multiple doesn't matter
    
    double sigma = allparams[(int) allparams.size() - 1];
    
    std::vector <double> params;

    for (int i = 0; i < M; i++) {
        params.push_back(allparams[i]);
    }
    
    double l;

    for (int i = 0; i < N; i++) {
        l = std::exp(-0.5 * w[i]* std::pow(((*yc)(x[i], params) - y[i])/std::sqrt(std::pow(sy[i], 2.0) + std::pow(sigma, 2.0)), 2.0)) / std::pow(std::pow(sy[i], 2.0) + std::pow(sigma, 2.0), w[i]/2.0);
        if (useLogPosterior){
            L += std::log(l);
        } else {
            L *= l;
        }
    }
//    printf("%f\n",L);
    return L; // returns log L = logL1 + logL2 + ... given L = L1*L2*L3... if useLogPosterior == true
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
			for (int i = 0; i < M; i++) { //check upper bound
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

			for (int i = 0; i < M; i++) { //check upper bound
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
			for (int i = 0; i < M; i++) {
				jointPrior *= priorsObject.priorsPDFs[i](allparams[i]);
			}

			break;
	}

	return jointPrior;
}

double TRK::posterior(std::vector <double> allparams, std::vector <double> allparams_trial) {
    // DEPRECATED !!!
    double post;
	if (hasPriors) {
        if (useLogPosterior){
            post = std::exp((*this.*selectedLikelihood)(allparams) - (*this.*selectedLikelihood)(allparams_trial)) * (priors(allparams));
            // this returns the log likelihood; the computation is done WITHIN the function.
        }
        else {
            post = (*this.*selectedLikelihood)(allparams) * priors(allparams);
        }

        return post;
	}
	else {
        if (useLogPosterior){
            post = std::exp((*this.*selectedLikelihood)(allparams) - (*this.*selectedLikelihood)(allparams_trial));
            // this returns the log likelihood; the computation is done WITHIN the function.
        }
        else {
            post = (*this.*selectedLikelihood)(allparams);
        }
	}
    
    return post;
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

double TRK::getMedian(std::vector<double> y)
{
	int high = (int)(floor(y.size() / 2));
	int low = high - 1;
	double runningSum = 0, median = 0;
	double totalSum = y.size();
	if (y.size() > 1)
	{
		if (y.size() % 2 == 0)
		{
			runningSum = y.size() / 2.0 + .5;
		}
		else
		{
			runningSum = y.size() / 2.0;
		}
		median = y[low] + (.5*totalSum - runningSum + 1.0)* (y[high] - y[low]);
	}

	else
	{
		median = y[0];
	}
	return median;

}

double TRK::getMedian(int trueCount, std::vector<double> w, std::vector<double> y)
{
	size_t sumCounter = 0;
	double median = 0, totalSum = 0, runningSum = 0;
	for (int i = 0; i < trueCount; i++)
	{
		totalSum += w[i];
	}
	if (trueCount > 1)
	{
		runningSum = w[sumCounter] * .5;
		while (runningSum < .5*totalSum)
		{
			sumCounter++;
			runningSum += w[sumCounter - 1] * .5 + w[sumCounter] * .5;
		}
		if (sumCounter == 0)
		{
			median = y[0];
			std::cout << median << std::endl;
		}
		else
		{
			median = y[sumCounter - 1] + (.5*totalSum - (runningSum - (w[sumCounter - 1] * .5 + w[sumCounter] * .5))) / (w[sumCounter - 1] * .5 + w[sumCounter] * .5)*(y[sumCounter] - y[sumCounter - 1]);
			std::cout << median << std::endl;
		}
	}
	else
	{
		median = y[0];
		std::cout << median << std::endl;
	}
	return median;
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
		roots = tangentCubicSolver(A, B, C, D);
		return roots;
	}
	//returns no extra roots (empty vector) if the other two roots are 
	return roots;
}

std::vector <double> TRK::tangentsFinder(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, double xg) {
	
	std::vector <double> result;

	double xg1 = xg;

	std::vector <double> xr1vec;
    double xr1old = 0.0;
	
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

double TRK::findBestTangent(std::vector <double> params, double x_n, double y_n, double Sig_xn2, double Sig_yn2, std::vector <double> x_tn_vec, double s) {
	std::vector <double> posts;
	long minindex;

	for (int i = 0; i < x_tn_vec.size(); i++) {
		posts.push_back(singlePointLnL(params, x_n, y_n, Sig_xn2, Sig_yn2, x_tn_vec[i], s));
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
    
    
    allparams_s = downhillSimplex(selectedChiSq, allparams_guess, ss[0]);

    if (hasAsymSlop){
        printf("%f \t %f \t %f \t %f \t %f \t(slop x optimization)\n", ss[0], allparams_s[M], allparams_s[M + 1], allparams_s[M + 2], allparams_s[M + 3]);
    } else {
        printf("s=%f \t slop_x=%f \t slop_y=%f \t(slop x optimization) model params: ", ss[0], allparams_s[M], allparams_s[M + 1]);
        for (int j = 0; j < M; j++){
            printf("%f\t", allparams_s[j]);
        }
        printf("\n");
    }
    

	//double sec_elapsed = secElapsed(time);

	//printf("%f sec, max threads = %i \n", sec_elapsed, maxThreads);

	getBetterSlopYGuess(allparams_s[M + 1], s);

	return allparams_s[M];
}

double TRK::innerSlopY_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
	//s = ss[0];

	
    allparams_s = downhillSimplex(selectedChiSq, allparams_guess, ss[0]);
        
	if (hasAsymSlop){
        printf("%f \t %f \t %f \t %f \t %f \t(slop y optimization)\n", ss[0], allparams_s[M], allparams_s[M + 1], allparams_s[M + 2], allparams_s[M + 3]);
    } else {
        printf("s=%f \t slop_x=%f \t slop_y=%f \t(slop y optimization) model params: ", ss[0], allparams_s[M], allparams_s[M + 1]);
        for (int j = 0; j < M; j++){
            printf("%f\t", allparams_s[j]);
        }
        printf("\n");
    }

	return allparams_s[M + 1];
}

double TRK::innerR2_Simplex(std::vector <double> ss, std::vector <double> allparams_guess) {
	//s = ss[0];

	whichExtrema = S;
	allparams_s = downhillSimplex(selectedChiSq, allparams_guess, ss[0]);
	whichExtrema = none;

	if (hasAsymSlop){
        printf("%f \t %f \t %f \t %f \t %f \t(initial R2 optimization)\n", ss[0], allparams_s[M], allparams_s[M + 1], allparams_s[M + 2], allparams_s[M + 3]);
    } else {
        printf("%f \t %f \t %f \t(initial R2 optimization)\n", ss[0], allparams_s[M], allparams_s[M + 1]);
    }

	double R2as = R2TRK_prime_as();
	double R2sb = R2TRK_prime_sb();

	return R2as - R2sb;
}

double TRK::innerR2_iter_Simplex(std::vector <double> ss, std::vector <double> allparams_guess, double s0) {
	//s = ss[0];

	whichExtrema = S;
    allparams_s = downhillSimplex(selectedChiSq, allparams_guess, ss[0]);
	whichExtrema = none;

	if (hasAsymSlop){
        printf("%f \t %f \t %f \t %f \t %f \t(additional R2 optimization)\n", ss[0], allparams_s[M], allparams_s[M + 1], allparams_s[M + 2], allparams_s[M + 3]);
    } else {
        printf("%f \t %f \t %f \t(additional R2 optimization)\n", ss[0], allparams_s[M], allparams_s[M + 1]);
    }

	double R2as = R2TRK_prime_as0(s0, x_t_s, params_s);
	double R2sb = R2TRK_prime_s0b(s0, x_t_s, params_s);

	return R2as - R2sb;
}

double TRK::optimize_s_SlopX() {

	iterative_allparams_guess = allparams_guess;

	// before doing any standard simplex movement, here it checks whether the simplex is within the zero "plateau", and if so, it moves it to the boundary.
	// for slop x: move to right until it hits the boundary


	//bracket finding

    double a = 0.0;
    double b = 1.0;
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

		whichExtremaX = slopx;
		slop_c = innerSlopX_Simplex({ c }, iterative_allparams_guess);
		whichExtremaX = none;

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

    double a = 0.0;
    double b = 1.0;
	double trial_s = slopYScaleGuess;
	double slop_trial_s = innerSlopY_Simplex({ trial_s }, iterative_allparams_guess);

	double inc = trial_s * 0.5;
    
    int jumpcheck = 0;

	if (slop_trial_s > 0) {
		a = trial_s;

		double trial_b = trial_s;

		while (true) {
            if (jumpcheck > 1){
                inc *= 2.0;
                jumpcheck = 0;
            }
			trial_b += inc;

			double slop_trial_b = innerSlopY_Simplex({ trial_b }, iterative_allparams_guess);

			if (slop_trial_b == 0) {
				b = trial_b;
				break;
			}
			else if (slop_trial_b > 0) {
				a = trial_b;
                jumpcheck += 1;
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

		whichExtremaY = slopy;
		slop_c = innerSlopY_Simplex({ c }, iterative_allparams_guess);
		whichExtremaY = none;

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

	double left, right, f_left;
	
	//bisection, now that we have brackets [left,right]

	left = a;
	right = b;

	f_left = innerR2_Simplex({ a }, iterative_allparams_guess);

	double c, f_c;
	double tol_bisect = 1e-4;
	double tol_brackets = 1e-3;

	while (true) {
		c = (left + right) / 2;

		whichExtrema = S;
		f_c = innerR2_Simplex({ c }, iterative_allparams_guess);
		whichExtrema = none;

		//printf("%f %f \n", f_c, c);

		if (std::abs(f_c) <= tol_bisect) { //convergence criterion
			break;
		}

		if (f_c * f_left > 0) {
			left = c;
			f_left = f_c;
		}
		else if (f_c * f_left < 0) {
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

	double left, right, f_left;

	f_left = innerR2_iter_Simplex({ a }, iterative_allparams_guess, s0);

	//bisection, now that we have brackets [left,right]

	left = a;
	right = b;

	double c, f_c;
	double tol_bisect = 1e-4;
	double tol_brackets = 1e-3;

	while (true) {
		//printf("brackets: %f %f \n", left, right);
		c = (left + right) / 2;

		whichExtrema = S;
		f_c = innerR2_iter_Simplex({ c }, iterative_allparams_guess, s0);
		whichExtrema = none;

		if (std::abs(f_c) <= tol_bisect) { //convergence criterion
			break;
		}

		if (f_c * f_left > 0) {
			left = c;
			f_left = f_c;
		}
		else if (f_c * f_left < 0) {
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

    double s1 = 0.0;

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
    
    if (do1DFit){
        printf("1D fit: no need for scale optimization.");
        return;
    }

	std::vector <double> scale_extrema;

    std::vector <double> s_slops;
    
    if (hasAsymSlop){
        s_slops = {0.0, 0.0, 0.0, 0.0};
    } else {
        s_slops = {0.0, 0.0};
    }

	//optimize simultaneously
	if (cpp17MultiThread) {
//        std::vector <int> nn = { 0, 1 };
//
//        std::for_each( //parallel tangent point finding
//            std::execution::par_unseq,
//            nn.begin(),
//            nn.end(),
//            [&](auto&& n)
//        {
//            s_slops[n] = (this->*optimizeList[n])();
//        });
    } else if (cpp11MultiThread && !cpp17MultiThread){
        s_slops.clear();
        
        int counter = 0, completedThreads = 0, liveThreads = 0;
        double result;
        std::vector< std::future < double > > futureVec;
        futureVec.resize(2);
        
        for (int i = 0; i < 2; i++)
        {                                                  //&TRK::tangentParallelLikelihood
            futureVec[i] = std::async(std::launch::async, optimizeList[i], this); //pointer to fn run through MT, arguments to fn
            counter++;
            liveThreads++;
            
            if (liveThreads >= maxThreads)
            {
                for (int i = completedThreads; i < counter; i++)
                {
                    result = futureVec[i].get();
                    s_slops.push_back(result);
                }
                completedThreads += liveThreads;
                liveThreads = 0;
            }
        }
        for (int i = completedThreads; i < 2; i++)
        {
            result = futureVec[i].get();
            s_slops.push_back(result);
        }
        
        std::vector <double> mM = minMax(s_slops);
        s_slops[0] = mM[0];
        s_slops[1] = mM[1];

	} else {
		#pragma omp parallel for //num_threads(8)
		for (int i = 0; i < 2; i++)
		{
			s_slops[i] = (this->*optimizeList[i])();
		}
	}

	double s_slopx = s_slops[0];
	double s_slopy = s_slops[1];

	printf("%f %f \t %f %f\n", params_slopx[0], params_slopx[1], params_slopy[0], params_slopy[1]);

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


	printf("%f %f %f %f \n", x_t_a[0], x_t_b[0], params_a[0], params_b[0]);

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
    
    if (hasAsymSlop){
        results.slop_x_minus = allparams_s[M + 2];
        results.slop_y_minus = allparams_s[M + 3];
    }

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
    
    if (pivotPointActive){
        allparams_0 = pivotPointParamsGuess;
    }
    
    if (posterior(allparams_0, allparams_0) == 0){
        printf("Alert: zero posterior for initial guess of MCMC!\n");
    }
    
	while (delta_count < sampleCount){// + burncount) {
		//create trial

		allparams_trial.clear();

		for (int j = 0; j < M + 2; j++) {
			allparams_trial.push_back(delta[j] * rnorm(0.0, 1.0) + allparams_0[j]);
		}

		a = metHastRatio(allparams_trial, allparams_0);
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

double TRK::metHastRatio(std::vector <double> X_trial, std::vector <double> X_i){
    double log_a;
    if (hasPriors) {
        log_a = (*this.*selectedLikelihood)(X_trial) - (*this.*selectedLikelihood)(X_i) + std::log(priors(X_trial)) - std::log(priors(X_i));
            // these likelihoods return log likelihood given useLogPosterior = true; the computation is done WITHIN the function.
    }
    else {
        log_a = (*this.*selectedLikelihood)(X_trial) - (*this.*selectedLikelihood)(X_i);
            // this returns the log likelihood given useLogPosterior = true; the computation is done WITHIN the function.
    }
    
    return log_a;
    // returns log post / log post if useLogPosterior == true
}

void TRK::guessMCMCDeltas(){
    params_sigmas_guess.clear();
    for (int j = 0; j < M; j++){
        params_sigmas_guess.push_back(10.0 * (double)1/N);
    }
    //guessing slops
    slop_x_sigma_guess = stDevUnweighted(x) / 100.0;
    slop_y_sigma_guess = stDevUnweighted(y) / 100.0;
    
    slop_x_minus_sigma_guess = slop_x_sigma_guess;
    slop_y_minus_sigma_guess = slop_y_sigma_guess;
    
    return;
}

std::vector <std::vector <double >> TRK::samplePosterior(int R, int burncount, std::vector <double> sigmas_guess) {
    
    useLogPosterior = true;

    unsigned long n = bigM;
    std::vector < std::vector <double > > result, result_final;
    
    switch (thisSamplingMethod) {
        case AIES: {
            double a = 2.0; //stretch variable pdf parameter
            // INDEXING:
            // iterations: t
            // walkers: j, k
            // coordinates/parameters: i
            
            int L = 2 * (int) bigM;   // number of walkers
            
            // initialize walkers
            std::vector <std::vector <double> > all_walkers(L, std::vector <double> (bigM, 0.0));
            
            for (int j = 0; j < L; j++){
                for (int i = 0; i < bigM; i++){
                    all_walkers[j][i] = rnorm(allparams_guess[i], allparams_guess[i]/10.0);
                }
            }
            
            double Z, log_a, rand_unif_log;
            int j;
            std::vector <double> X, Y, X_trial;
            // sample
            int sample_count = 0;
            int accept_count = 0;
            double accept_frac = 0.0;
            int tenth = (int) R/10;
            int prog = 0;
            
            while (sample_count < R + burncount) {
                for (int k = 0; k < L; k++){
                    // for each kth walker X,
                    X = all_walkers[k];
                    // choose some other jth walker Y:
                    j = k;
                    while (j == k){
                        j = rand() % L;
                    }
                    Y = all_walkers[j];
                    
                    // make proposal vector
                    X_trial.clear();
                    Z = rstretch(a);
                    for (int i = 0; i < bigM; i++){
                        X_trial.push_back(Z * X[i] + (1.0 - Z) * Y[i]);
                    }
                    
                    // accept?
                    rand_unif_log = std::log(runiform(0.0, 1.0));
                    log_a = metHastRatio(X_trial, X);
                    
                    if (rand_unif_log <= log_a + (bigM - 1) * std::log(Z)) { // accept
                        all_walkers[k] = X_trial;
                        result.push_back(X_trial);
                        accept_count++;
                    }
                    else { // reject
                        result.push_back(X_trial);
                    }
                    sample_count += 1;
                }
                if (sample_count % tenth == 0){
                        accept_frac = (double) accept_count / (double) sample_count;
                        if (verboseMCMC){
                            printf("Posterior sampling %i%% complete, acceptance fraction currently %f...\n", prog, accept_frac);
                        }
                        prog += 10;
                    }
                }
            
            if (verboseMCMC){
                printf("Final AIES acceptance ratio: %f \n", accept_frac);
            }
            
            break;
        }
        case ARWMH: {
            std::vector <std::vector <double > > cov_i(n, std::vector <double> (n, 0.0));
            for (int j = 0; j < n; j++){
                cov_i[j][j] = std::pow(sigmas_guess[j], 2.0);
            }
            std::vector <std::vector <double > > cov_i1(n, std::vector <double> (n, 0.0));
            
            std::vector <double> allparams_trial, allparams_0; //allparams_0 is the previous step
            double log_a, rand_unif_log, accept_frac = 0.0;

            int accept_count = 0;
            int delta_count = 0;
            int tenth = (int) R/10;
            int prog = 0;
            
            std::vector <double> mu_i = allparams_guess;
            std::vector <double> mu_i1(n, 0.0);
            std::vector <double> X_i = allparams_guess;
            std::vector <double> X_i1;
            std::vector <double> X_trial(n, 0.0);
            double lamb = std::pow(2.38, 2.0) / n;
            double gam_i1 = 1.0;
            int i = 0;

            while (delta_count < R + burncount) {
                bool loopCheck = true;
                //create trial:
                //sample X_i
                while (loopCheck) {
                    for (int j = 0; j < bigM; j++) {
                        X_trial[j] = rnorm(mu_i[j], std::sqrt(lamb * cov_i[j][j]));
                    }

                    log_a = metHastRatio(X_trial, X_i);

                    rand_unif_log = std::log(runiform(0.0, 1.0));

                    if (rand_unif_log <= log_a) { // accept
                        X_i1 = X_trial;
                        loopCheck = false;
                        delta_count += 1;
                        result.push_back(X_i1);
                        accept_count += 1;
                    }
                    else { // reject
                        delta_count += 1;
                        result.push_back(X_i);
                    }
                }

                //update proposal dist params
                
        //        X_i1 = X_trial;

                gam_i1 = 1.0/((double)(i + 1));

                for (int j = 0; j < n; j++){
                    mu_i1[j] = mu_i[j] + gam_i1*(X_i1[j] - mu_i[j]);
                }

                for (int l = 0; l < n; l++){
                    for (int m = 0; m < n; m++){
                        cov_i1[l][m] = cov_i[l][m] + gam_i1*((X_i1[l] - mu_i[l])*(X_i1[m] - mu_i[m])-cov_i[l][m]);
                    }
                }

                mu_i = mu_i1;
                cov_i = cov_i1;
                X_i = X_i1;

                i += 1;

                accept_frac = (double) accept_count / (double)delta_count;
                
                if (delta_count % tenth == 0){
                    if (verboseMCMC){
                        printf("Posterior sampling %i%% complete, acceptance fraction currently %f...\n", prog, accept_frac);
                    }
                    prog += 10;
                }
            }
            
            if (verboseMCMC){
                printf("final Adaptive Metropolis Hastings step-sizes/proposal dist. deviations:");
                for (int j = 0; j < bigM; j++) {
                    printf("%f ", cov_i[j][j]);
                }
                printf("\t final Adaptive Metropolis Hastings acceptance ratio: %f \n", accept_frac);
            }
            break;
        }
        default: {
            break;
        }
    }
    //cut off burn-in

    for (int i = 0; i < R; i++) {
        result_final.push_back(result[i + burncount]);
    }

	result_final = checkSlopSignMCMC(result_final);
    
    useLogPosterior = false;

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

double TRK::rstretch(double a){ // distribution to sample "stretching variable" for AIES
    // using rejection sampling
    double M = 5;
    
    double x = rnorm(1.0, 1.0);   // proposal dist q(x), with Mq(x) >= p(x)
    double u = runiform(0.0, 1.0);

    while (u > (stretch_pdf(x, a) / (M * normal(x, 1.0, 1.0)))){
        x = rnorm(1.0, 1.0);
        u = runiform(0.0, 1.0);
    }
    
    return x;
}

double noPrior(double param) {
	return 1.0;
}

std::vector <std::vector <double> > TRK::getHistogram(std::vector <double> data) {
	unsigned long dataSize = data.size();

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

std::vector <std::vector <double> > TRK::getHistogram(std::vector <double> data, std::vector <double> weights) {
    unsigned long dataSize = data.size();
    
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
    
    double bintemp;
    
    for (int i = 0; i < bincount; i++) {
        bintemp = 0.0;
        for (int j = 0; j < dataSize - 1; j++) {
            if ((data[j] >= edges[i]) && (data[j] < edges[i + 1])) {
                bintemp += weights[j]; //adds data to the bin if it is between the bounds of the bin
            }
        }
        if ((data[dataSize - 1] > edges[i]) && (data[dataSize - 1] <= edges[i + 1])) {
            bintemp += weights[dataSize - 1]; //adds data to the bin if it is between the bounds of the bin
        }
        hist.push_back(bintemp); //adds an array of data for that bin to bins, the 2D array.
    }
    
    return { hist, edges };
}

std::vector <std::vector <std::vector <double> > >  TRK::lowerBar(std::vector <std::vector <double> > allparam_samples) { //method used to estimate uncertainty of a sampled distribution
	std::vector <double> data, hist, edges, minusSigmas, plusSigmas;
	std::vector <std::vector <std::vector <double> > > allparam_uncertainties;
	std::vector <std::vector <double> > histResults;
	unsigned long totalCount = allparam_samples.size();
	double tolBar = 1e-6;

	results.paramDistributionHistograms.clear();
    
    int m = 0;
    if (hasAsymSlop){
        m = 2;
    }

	for (int j = 0; j < bigM + m; j++) { //for each model param plus slop
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
		int aboveCount; // number of SAMPLES within the bins above bar

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

			unsigned long K = indicesIn.size();
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
    
//    useLogPosterior = false;
    // only needed the above true for pivot point sampling
//    goodDeltasFound = false;
    // recompute step sizes cause they may be different given new pivot point

	std::vector <std::vector <std::vector <double> > > allparam_uncertainties;
    
    if (verboseMCMC){
        std::cout << "Sampling Posterior...\n";
    }

	std::vector <std::vector <double >> allparam_samples = samplePosterior(R, burncount, allparams_sigmas_guess);

	if (outputDistributionToFile) {

		std::string fileName = outputPath + std::string("/TRKMCMC_") + std::to_string(allparams_guess[0]) + std::string("_") + std::to_string(R) + std::string(".txt");

		std::ofstream myfile;
		myfile.open(fileName, std::ofstream::trunc);
        
        int m = 0;
        if (hasAsymSlop){
            m = 2;
        }
        
		if (myfile.is_open())
		{
			// filename    a     b     optimum scale    total computation time (s)
			for (int i = 0; i < allparam_samples.size(); i++) {
				for (int j = 0; j < bigM + m; j++) {
					myfile << allparam_samples[i][j] << " ";
				}
				myfile << std::endl;
			}

			myfile.close();
		}
		else std::cout << "Unable to open file";
	}

    if (verboseMCMC){
        std::cout << "Computing Uncertainties...\n";
    }

	allparam_uncertainties = lowerBar(allparam_samples); //for each parameter including slope, there is a vector containing 1 vector of -sigmas, 1 vector of +sigmas. This vector contains all of those 2-vectors.

	results.bestFit_123Sigmas.clear();

	for (int j = 0; j < M; j++) {
		results.bestFit_123Sigmas.push_back(allparam_uncertainties[j]);
	}
    if (do1DFit){
        results.slopY_123Sigmas = allparam_uncertainties[M];
    } else {
        results.slopX_123Sigmas = allparam_uncertainties[M];
        results.slopY_123Sigmas = allparam_uncertainties[M + 1];
        
        if (hasAsymSlop){
            results.slopX_minus_123Sigmas = allparam_uncertainties[M];
            results.slopY_minus_123Sigmas = allparam_uncertainties[M + 1];
        }
    }

	return;
}

// PIVOT POINTS
void TRK::getCombos(std::vector <std::vector <double> > total, int k, int offset) { //ND case in x

	if (k == M) {
		NDcombos.clear();
	}
	if (k == 0) {
		NDcombos.push_back(NDcombination);
		return;
	}
	for (int i = offset; i <= total.size() - k; ++i) {
		NDcombination.push_back(total[i]);
		getCombos(total, k - 1, i + 1);
		NDcombination.pop_back();
	}
}

double TRK::pivotFunc(std::vector <double> params1, std::vector <double> params2) {
    double a01 = linearizedIntercept(params1);
	double a11 = linearizedSlope(params1);

	double a02 = linearizedIntercept(params2);
	double a12 = linearizedSlope(params2);

	return (a02 - a01) / (a11 - a12);
}

double TRK::pivotFunc2(std::vector <double> params1, std::vector <double> params2) {
    double a01 = linearizedIntercept2(params1);
    double a11 = linearizedSlope2(params1);

    double a02 = linearizedIntercept2(params2);
    double a12 = linearizedSlope2(params2);

    return (a02 - a01) / (a11 - a12); // (b1 - b2) / (m2 - m1);
}

double TRK::weightPivot(std::vector <double> params1, std::vector <double> params2, std::vector <double> oldPivots, double newPivot) {
	std::vector <double> squares(oldPivots.size(), 0.0);
	double w;

	for (int i = 0; i < oldPivots.size(); i++) {
		squares[i] = std::pow(newPivot - oldPivots[i], 2.0);
	}

	double avg = getAverage(squares);
    
    double b1 = linearizedIntercept(params1);
    double m1 = linearizedSlope(params1);

    double b2 = linearizedIntercept(params2);
    double m2 = linearizedSlope(params2);

	w = std::pow(avg / std::pow(m2 - m1, 2.0) + std::pow(b1 - b2, 2.0) / std::pow(m2 - m1, 4.0), -1.0);

	return w;
}

double TRK::weightPivot2(std::vector <double> params1, std::vector <double> params2, std::vector <double> oldPivots, double newPivot) {
    std::vector <double> squares(oldPivots.size(), 0.0);
    double w;

    for (int i = 0; i < oldPivots.size(); i++) {
        squares[i] = std::pow(newPivot - oldPivots[i], 2.0);
    }

    double avg = getAverage(squares);
    
    double b1 = linearizedIntercept2(params1);
    double m1 = linearizedSlope2(params1);

    double b2 = linearizedIntercept2(params2);
    double m2 = linearizedSlope2(params2);

    w = std::pow(avg / std::pow(m2 - m1, 2.0) + std::pow(b1 - b2, 2.0) / std::pow(m2 - m1, 4.0), -1.0);

    return w;
}

std::vector < std::vector <std::vector <double > > > TRK::directCombos(std::vector < std::vector <double> > params_sample, int comboCount){
    std::vector < std::vector <std::vector <double > > > combos;
    combos.clear();
    
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, (int)params_sample.size() - 1); // define the range
    
    for (int j = 0; j < comboCount; j++){
        int i1 = (int) distr(eng);
        int i2 = (int) distr(eng);

        combos.push_back({params_sample[i1], params_sample[i2]});
        
        if (params_sample[i1].size() == 0 || params_sample[i2].size() == 0 ){
            printf("%i %i \n", i1, i2);
        }
    }
    
    return combos;
}

std::vector <double> TRK::removeOutlierPivots(std::vector <double> pivots){
    std::vector <double> newpivots;
    double pivot;
    
    int outCount = 0;
    
    for (int i = 0; i < pivots.size(); i++){
        pivot = pivots[i];
        if (pivot  < x_max + pruneWidth*datawidth && pivot > x_min - pruneWidth*datawidth){
            newpivots.push_back(pivot);
        } else {
            outCount++;
        }
    }
    if (verbosePivotPoints){
        printf("%i pivots outside of reasonable region \n", outCount);
    }
    
    return newpivots;
}

void TRK::getPivotGuess(){
    if (findPivotPoints){
        if (twoPivots){
            double medX = getMedian((int)N, w, x);
            
            std::vector <double> lowX, highX, lowW, highW;
            
            for (int i = 0; i < N ; i++){
                x[i] >= medX ? highX.push_back(x[i]) : lowX.push_back(x[i]);
                x[i] >= medX ? highW.push_back(w[i]) : lowW.push_back(w[i]);
            }
            pivot = getAverage(lowX, lowW);
            pivot2 = getAverage(highX, highW);
            
            if (verbosePivotPoints){
                printf("p1: %f\t p2: %f \n", pivot, pivot2);
            }
            
        } else {
            pivot = getAverage(x, w);
        }
    }
    return;
}

void TRK::findPivots() {
	if (findPivotPoints) {
        
//        useLogPosterior = true;
        
		std::vector < std::vector <double > > allparam_samples;
		std::vector < std::vector < std::vector <double> > > drawnCombos;
        std::vector <double> pivots, pivotWeights, allPivots, allparams_better;
		std::vector <double> oldPivots((int)(randomSampleCount*(randomSampleCount - 1)) / 2, pivot);
        double onePivot, oneWeight;
        double finalPivot = 1.0;
		int iter = 0;
        
        // for two pivot point special case (e.g. bh v c2 and rv v c2):
        std::vector < std::vector <double > > allparam_samples2;
        std::vector <double> pivots2, pivotWeights2, allPivots2;
        std::vector <double> oldPivots2((int)(randomSampleCount*(randomSampleCount - 1)) / 2, pivot2);
        double onePivot2, oneWeight2;
        double finalPivot2 = 2.0;

		while (true) {
            if (iter > 0){
                pivotPointActive = true;
            }
            
            //refit for better guess for MCMC to avoid zero likelihood
            
            if (modeInterceptGuess){
                printf("modeInterceptGuess=true not currently configured for multiple pivot points.\n");
                allparams_better = allparams_guess;
                
                std::vector <double> preIntercepts;
                
                for (int i = 0; i < N; i++){
                    preIntercepts.push_back(y[i]-allparams_guess[1]*(x[i] - pivot));
                }
                std::vector <double> ones(N, 1.0);
                allparams_better[0] = getPeakCoord(preIntercepts, ones);
                
            } else {
            
                allparams_better = downhillSimplex(selectedChiSq, allparams_guess, s);
                
                if (verbosePivotPoints){
                    printf("re-fit for new pivot point; old / new params:\n");
                
                    for (int j = 0; j < (int)allparams_better.size(); j++){
                        printf("%f %f\n", allparams_guess[j], allparams_better[j]);
                    }
                }
            }
            
            allparams_guess = allparams_better;


			pivots.clear();
			pivotWeights.clear();
            
            pivots2.clear();
            pivotWeights2.clear();
            
            std::vector < std::vector <double> > param_samples(pivotR, { 0.0, 0.0 });
            
            if (verbosePivotPoints){
                printf("Sampling for pivot points...\n");
            }

			allparam_samples = samplePosterior(pivotR, pivotBurnIn, allparams_sigmas_guess); //allparam_samples is { {allparams0}, {allparams1}, ... }
            
            pivotPointParamsGuess = allparams_guess;
            
            if (averageIntercepts){
                std::vector <double> allBs, allBs2;
                for (int i = 0; i < (int)allparam_samples.size(); i++){
                    allBs.push_back(allparam_samples[i][0]);
                    
                    if (twoPivots){
                        allBs2.push_back(allparam_samples[i][2]);
                    }
                }
            
                pivotPointParamsGuess[0] = getAverage(allBs);
                if (twoPivots){
                    pivotPointParamsGuess[2] = getAverage(allBs2);
                }
            }

            for (int j = 0; j < allparam_samples.size(); j++) {
                param_samples[j] = slice(allparam_samples[j], 0, (int)M);
            }
            
			if (!getCombosFromSampleDirectly) { //this option takes the ~10,000 MH samples, selects ~200 of them, then generates combos out of this subset
                printf("getCombosFromSampleDirectly=false not currently configured for multiple pivot points.\n");
				random_unique(param_samples.begin(), param_samples.end(), randomSampleCount);

				param_samples = slice(param_samples, 0, randomSampleCount);
                
                NDcombos.clear();
				getCombos(param_samples, 2, 0); //generates all 2-combos of the parameter space data points

				drawnCombos = NDcombos;
			}
			else { //this option takes the ~10,000 MH samples, then generates combos directly out of this set.
                drawnCombos = directCombos(param_samples, maxCombos);
			}

            if (twoPivots){
                for (int j = 0; j < drawnCombos.size(); j++) {
                    onePivot = pivot + pivotFunc(drawnCombos[j][0], drawnCombos[j][1]);
                    onePivot2 = pivot2 + pivotFunc2(drawnCombos[j][0], drawnCombos[j][1]);
                    if (!std::isnan(onePivot) && !std::isnan(onePivot2)) {
                        pivots.push_back(onePivot);
                        pivots2.push_back(onePivot2);
                        //std::cout << onePivot << std::endl;
                    }
                }
            } else {
                for (int j = 0; j < drawnCombos.size(); j++) {
                    onePivot = pivot + pivotFunc(drawnCombos[j][0], drawnCombos[j][1]);
                    if (!std::isnan(onePivot)) {
                        pivots.push_back(onePivot);
                        //std::cout << onePivot << std::endl;
                    }
                }
            }
            
            if (pruneOutlierPivots){
                if (twoPivots){
                    pivots2 = removeOutlierPivots(pivots2);
                }
                pivots = removeOutlierPivots(pivots);
            }
            
            if (verbosePivotPoints){
                printf("Weighting pivot points...\n");
            }
            
            pivotWeights = std::vector <double>(pivots.size(), 1.0);
            pivotWeights2 = std::vector <double>(pivots2.size(), 1.0);

            std::vector <double> finalPivots, finalWeights;
            std::vector <double> finalPivots2, finalWeights2;
            
			if (weightPivots) {
                if (twoPivots){
                    for (int i = 0; i < pivots.size(); i++) {
                        oneWeight = weightPivot(drawnCombos[i][0], drawnCombos[i][1], oldPivots, pivots[i]);
                        if (!std::isnan(oneWeight)){
                            finalPivots.push_back(pivots[i]);
                            finalWeights.push_back(oneWeight);
                        }
                    }
                    for (int i = 0; i < pivots2.size(); i++) {
                        oneWeight2 = weightPivot2(drawnCombos[i][0], drawnCombos[i][1], oldPivots2, pivots2[i]);
                        if (!std::isnan(oneWeight2)){
                            finalPivots2.push_back(pivots2[i]);
                            finalWeights2.push_back(oneWeight2);
                        }
                    }
                } else {
                    for (int i = 0; i < pivots.size(); i++) {
                        oneWeight = weightPivot(drawnCombos[i][0], drawnCombos[i][1], oldPivots, pivots[i]);
                        if (!std::isnan(oneWeight)){
                            finalPivots.push_back(pivots[i]);
                            finalWeights.push_back(oneWeight);
                        }
                    }
                }
            } else {
                finalWeights = pivotWeights;
                finalPivots = pivots;
                
                finalWeights2 = pivotWeights2;
                finalPivots2 = pivots2;
            }
            
            pivots = finalPivots;
            pivotWeights = finalWeights;
            
            pivots2 = finalPivots2;
            pivotWeights2 = finalWeights2;
            
            if (twoPivots){
                if (pivotMedian){
                    finalPivot2 = getMedian((int) pivots2.size(), pivotWeights2, pivots2);
                } else if (pivotMean){
                    finalPivot2 = getAverage(pivots2, pivotWeights2);
                } else if (pivotHalfSampleMode){
                    finalPivot2 = getMode((int) pivots2.size(), pivotWeights2, pivots2);
                } else { //mode
                    finalPivot2 = getPeakCoord(pivots2, pivotWeights2);
                }
            }
            if (pivotMedian){
                finalPivot = getMedian((int) pivots.size(), pivotWeights, pivots);
            } else if (pivotMean){
                finalPivot = getAverage(pivots, pivotWeights);
            } else if (pivotHalfSampleMode){
                finalPivot = getMode((int) pivots.size(), pivotWeights, pivots);
            } else { //mode
                finalPivot = getPeakCoord(pivots, pivotWeights);
            }

			if (writePivots && !twoPivots) {

                std::string filename = std::string("/Users/nickk124/research/reichart/TRK/TRKrepo_public/diagnostics/") + std::string("TRKpivots") + (getCombosFromSampleDirectly ? "1" : "0") + (weightPivots ? "1_" : "0_") + std::to_string(iter) + std::string("_") + std::to_string(finalPivot) + std::string(".txt");

				std::ofstream myfile(filename, std::ofstream::trunc);
				if (myfile.is_open())
				{
					for (int i = 0; i < pivots.size(); i++) {
						myfile << pivots[i] << " " << pivotWeights[i] << "\n";
					}

					//myfile << "final pivot: " << finalPivot << "\n\n\n\n\n";

					myfile.close();
				}
				else std::cout << "Unable to open file";
			}

            if (twoPivots){
                if (verbosePivotPoints){
                    printf("new pivot1, old pivot1= %f \t, %f \n", finalPivot, pivot);
                    printf("new pivot2, old pivot2= %f \t, %f \n", finalPivot2, pivot2);
                }
                
                allPivots2.push_back(finalPivot2);
            } else {
                if (verbosePivotPoints){
                    printf("new pivot, old pivot= %f \t, %f \n", finalPivot, pivot);
                }
            }
            
            allPivots.push_back(finalPivot);
            iter += 1;
            if (twoPivots){
                if (((stDevUnweighted(allPivots) <= pivotTol && stDevUnweighted(allPivots2) <= pivotTol)|| iter >= maxPivotIter) && maxPivotIter > 1) {
                    if (weightPivots){
                        pivot = getAverage(slice(allPivots, 1, (int)allPivots.size()));
                        pivot2 = getAverage(slice(allPivots2, 1, (int)allPivots2.size()));
                    } else {
                        pivot = getAverage(allPivots);
                        pivot2 = getAverage(allPivots2);
                    }
                    break;
                } else if (maxPivotIter == 1){
                    pivot = finalPivot;
                    pivot2 = finalPivot2;
                    break;
                }
            } else {
                if ((stDevUnweighted(allPivots) <= pivotTol || iter >= maxPivotIter) && maxPivotIter > 1) {
                    if (weightPivots){
                        pivot = getAverage(slice(allPivots, 1, (int)allPivots.size()));
                    } else {
                        pivot = getAverage(allPivots);
                    }
                    break;
                } else if (maxPivotIter == 1){
                    pivot = finalPivot;
                    break;
                }
            }

			oldPivots = pivots;
            if (twoPivots){
                oldPivots2 = pivots2;
            }
		}
        if (twoPivots){
            if (verbosePivotPoints){
                printf("final pivot point 1: %f \n", pivot);
                printf("final pivot point 2: %f \n", pivot2);
            }
            
            results.pivot2 = pivot2;
        } else {
            if (verbosePivotPoints){
                printf("final pivot point: %f \n", pivot);
            }
        }
        
        results.pivot = pivot;
        
        pivotPointActive = false;
        
        allparams_guess =  pivotPointParamsGuess;

		return;
	
	}
	else {
		return;
	}
}

// OTHER TOOLS

double TRK::getPeakCoord(std::vector <double> x, std::vector <double> w){
    double xPeak;
    std::vector <double> hist, edges;
    
    std::vector < std::vector <double> > histResults = getHistogram(x, w);
    
    hist = histResults[0];
    edges = histResults[1];
    
    int maxInd = argMinMax(hist)[1];
    xPeak = (edges[maxInd + 1] + edges[maxInd]) / 2.0;
    
    return xPeak;
}

void TRK::checkVerbose(){
    if (verbose){
        verboseMCMC = true;
        verbosePivotPoints = true;
        verboseAsymmetric = true;
    }
    return;
}

// CORE ALGORITHMS/TRK FITS

void TRK::getBetterGuess(){
    if (do1DFit){
        results.bestFitParams.clear();

        whichExtrema = S;
        allparams_s = downhillSimplex(selectedChiSq, allparams_guess, s);
        whichExtrema = none;

        for (int j = 0; j < M; j++) {
            results.bestFitParams.push_back(allparams_s[j]);
        }
        results.slop_y = allparams_s[M];
        
        if (hasAsymSlop){
            printf("NOTE: no current support for asymmetric 1D statistics.\n");
        }
    }
    
    for (int j = 0; j < M; j++){
        allparams_guess[j] = results.bestFitParams[j];
    }
    if (do1DFit){
        allparams_guess[M] = results.slop_y;
    } else {
        allparams_guess[M] = results.slop_x;
        allparams_guess[M+1] = results.slop_y;
        
        if (hasAsymSlop){
            allparams_guess[M+2] = results.slop_x_minus;
            allparams_guess[M+3] = results.slop_y_minus;
        }
    }
    
    return;
}

void TRK::performTRKFit() {//finds optimum scale AND performs TRK fit + uncertainty
    checkVerbose();
    
    checkAsym();
    
    getPivotGuess();
    
	optimizeScale();

	findPivots();

    getBetterGuess();
	calculateUncertainties();
    
    showResults(true, true);
}

void TRK::performTRKFit(double scale) {//perform fit on some provided scale (for example, if they already know optimum scale, they can just start with this) and calculates uncertainties
    checkVerbose();
    
    checkAsym();
    
	s = scale;
    results.optimumScale = s;
    
//    printf("%f\n", pivot);
    
    getPivotGuess();
    
//    printf("%f\n", pivot);

	findPivots();

	results.bestFitParams.clear();

	whichExtrema = S;
	allparams_s = downhillSimplex(selectedChiSq, allparams_guess, scale);
	whichExtrema = none;

	for (int j = 0; j < M; j++) {
		results.bestFitParams.push_back(allparams_s[j]);
	}

	results.slop_x = allparams_s[M];
	results.slop_y = allparams_s[M + 1];
    
    if (hasAsymSlop){
        results.slop_x_minus = allparams_s[M + 2];
        results.slop_y_minus = allparams_s[M + 3];
    }

    getBetterGuess();
	calculateUncertainties();
    
    showResults(false, true);
}

void TRK::performSimpleTRKFit() {//finds optimum scale and performs TRK fit but without finding uncertainties
    checkVerbose();
    
    checkAsym();
    
    getPivotGuess();
    
	optimizeScale(); // (stores results in TRK.results)

	findPivots();
    
    results.bestFitParams.clear();
    
    whichExtrema = S;
    allparams_s = downhillSimplex(selectedChiSq, allparams_guess, s);
    whichExtrema = none;
    
    for (int j = 0; j < M; j++) {
        results.bestFitParams.push_back(allparams_s[j]);
    }
    
    results.slop_x = allparams_s[M];
    results.slop_y = allparams_s[M + 1];
    
    if (hasAsymSlop){
        results.slop_x_minus = allparams_s[M + 2];
        results.slop_y_minus = allparams_s[M + 3];
    }
    
    showResults(true, false);

	return;
}

void TRK::performSimpleTRKFit(double scale) {//given some provided scale, performs TRK fit but without finding uncertainties
    checkVerbose();
    
    checkAsym();
    
    s = scale;
    
    getPivotGuess();
    
    findPivots();
    
    results.bestFitParams.clear();
    
    whichExtrema = S;
    allparams_s = downhillSimplex(selectedChiSq, allparams_guess, s);
    whichExtrema = none;
    
    for (int j = 0; j < M; j++) {
        results.bestFitParams.push_back(allparams_s[j]);
    }
    
    results.slop_x = allparams_s[M];
    results.slop_y = allparams_s[M + 1];
    
    if (hasAsymSlop){
        results.slop_x_minus = allparams_s[M + 2];
        results.slop_y_minus = allparams_s[M + 3];
    }
    
    showResults(false, false);
    
    return;
}

void TRK::showResults(bool didScaleOp, bool didMCMC){
    printf("\n\nTRK RESULTS:\n\n");
    
    if (didScaleOp && !do1DFit){
        printf("SCALES:\nOptimum scale: %f \n", results.optimumScale);
        printf("Minimum scale: %f \n", results.minimumScale);
        printf("Maximum scale: %f \n\n", results.maximumScale);
    }
    
    if (findPivotPoints){
        printf("PIVOT POINTS:\nPivot Point 1 = %f\n", results.pivot);
        if (twoPivots){
            printf("PIVOT POINTS:\n\nPivot Point 2 = %f\n\n", results.pivot2);
        }
        printf("\n");
    }
    
    if (covid19){
        for (int i = 0; i < (int)params_guess.size()/2; i++){
            printf("intercept %i = %4.3f +/- %4.3f\n", i+1, results.bestFitParams[0+2*i], results.bestFit_123Sigmas[0+2*i][1][0]);
            printf("slope %i = %4.3f +/- %4.3f\n", i+1, results.bestFitParams[1+2*i], results.bestFit_123Sigmas[1+2*i][1][0]);
        }
        
        if (params_guess.size() == 5){
            printf("s = %4.3f +/- %4.3f\n", results.bestFitParams[4], results.bestFit_123Sigmas[4][1][0]);
        }
        
        printf("y_slop = %4.3f +/- %4.3f\n\n", results.slop_y, results.slopY_123Sigmas[1][0]);
        
        std::vector <double> pivs = {results.pivot, results.pivot2};
        
        for (int i = 0; i < (int)params_guess.size()/2; i++){
            printf("pivot %i = %4.3f\n", i+1, pivs[i]);
        }
        
        return;
    }


    printf("FITTED PARAMETERS (including slop):\n");
    for (int k = 0; k < params_guess.size(); k++) {
        printf("%f ", results.bestFitParams[k]);
    }
    printf(" %f %f", results.slop_x, results.slop_y);
    std::cout << std::endl;
    
    if (didMCMC){
        printf("\nUNCERTAINTIES: (- 1 2 3, + 1 2 3): \n");
        for (int k = 0; k < params_guess.size(); k++) { //kth param
            for (int j = 0; j < 2; j++) { // - and + sigmas
                for (int i = 0; i < 3; i++) { // 1, 2 and 3 sigmas
                    printf("%f ", results.bestFit_123Sigmas[k][j][i]);
                }
                printf("\t");
            }
            std::cout << std::endl;
        }

        printf("\nSLOP UNCERTAINTIES: (- 1 2 3, + 1 2 3): \n");
        if (!do1DFit){
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < 3; i++) {
                    printf("%f ", results.slopX_123Sigmas[j][i]);
                }
                printf("\t");
            }
            std::cout << std::endl;
        }
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 3; i++) {
                printf("%f ", results.slopY_123Sigmas[j][i]);
            }
            printf("\t");
        }
        std::cout << std::endl << std::endl;
    }
    
    printf("\nFITNESS\n\nchisquared = %f\n\n", results.chisquared);
    
    return;
}

//global functions

double twoPointNR(double(*y)(double, std::vector <double>), double(*dy)(double, std::vector <double>), double(*ddy)(double, std::vector <double>), std::vector <double> params, double xguess, double xguessp1)
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

		ykm1 = y(xkm1, params);
		yk = y(xk, params); //function we're finding zero of
		dyk = dy(xk, params); //derivative of above

		r = 1.0 - ((yk / ykm1) * (((yk - ykm1) / (xk - xkm1)) / dyk));

		xkp1 = (1.0 - 1.0 / r)*xkm1 + xk / r;

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


std::vector <double> cubicSolver(double A, double B, double C, double D) {
	double a1 = B / A;
	double a2 = C / A;
	double a3 = D / A;

	double Q = (a1 * a1 - 3.0 * a2) / 9.0;
	double R = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
	double Qc = std::pow(Q, 3.0);
	//double d = Qc - std::pow(R, 2.0);

	double theta = std::acos(R / sqrt(Qc));

	double r1 = -2 * std::sqrt(Q) * std::cos(theta / 3) - a1 / 3;
	double r2 = -2 * std::sqrt(Q) * std::cos((theta + 2 * PI) / 3) - a1 / 3;
	double r3 = -2 * std::sqrt(Q) * std::cos((theta + 4 * PI) / 3) - a1 / 3;

	std::vector <double> roots = { r1, r2, r3 };

	return roots;
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

	std::ofstream myfile("/Users/nickk124/research/reichart/TRK/TRKrepo_public/diagnostics/TRKresults.txt", std::ofstream::app);
	if (myfile.is_open())
	{	
		// filename    a     b     optimum scale    total computation time (s)
		myfile << filename << "\t" <<  TRKobj.results.minimumScale << "\t" << TRKobj.results.maximumScale << "\t" << TRKobj.results.optimumScale << "\t" << t_sec << "\n";
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
//        if (!readFile.good())
//            break;
        
        std::stringstream iss(line);
        
        for (int col = 0; col < columns; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            //            if (!iss.good())
            //                break;
            
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
