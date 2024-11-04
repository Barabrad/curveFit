//
//  curveFit.cpp
//  curveFit
//
//  Created by Brad Barakat on 1/2/23.
//  Finished on 1/3/23. Modified on 1/4/23.
//  Modified by Brad Barakat on 10/31/24.

#include "curveFit.hpp"

using namespace std;

/**
 * @brief  This function finds the squared error for a data set and modeling function
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  sqErr : A double for the squared error
 */
double findSqErr(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double>& params) {
    double sqErr = 0.0;
    for (int i = 0; i < data_x.size(); i++) {
        double yi = data_y[i];
        double xi = data_x[i];
        sqErr += pow(yi - fun(xi,params), 2);
    }
    return sqErr;
}

/**
 * @brief  This function finds the derivative of the squared error for a data set and modeling function
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  funDer : A function for a derivative of the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  derErr : A double for the derivative of the squared error
 */
double findDerErr(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), double (*funDer)(double, vector<double>&), vector<double>& params) {
    double derErr = 0.0;
    for (int i = 0; i < data_x.size(); i++) {
        double yi = data_y[i];
        double xi = data_x[i];
        derErr += -2*funDer(xi,params)*(yi - fun(xi,params));
    }
    return derErr;
}

/**
 * @brief  This function finds the index of the parameter that causes the steepest gradient for the squared error
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  funDers : A vector of partial derivatives for the model that each take a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @param  avoidInd : A vector of ints determining which indices to ignore when iterating through the parameter vector
 * @retval  indAndGrad : A vector of doubles containing the index and the steepest gradient, respectively
 */
vector<double> findSteepestGrad(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double (*)(double, vector<double>&)>& funDers, vector<double>& params, vector<int>& avoidInd) {
    int index = 0;
    double maxGrad = 0;
    for (int i = 0; i < funDers.size(); i++) {
        if (avoidInd.at(i) == 0) {
            double derErr_i = findDerErr(data_x, data_y, fun, funDers.at(i), params);
            if (abs(derErr_i) >= abs(maxGrad)) {
                index = i; maxGrad = derErr_i;
            }
        } // If avoidInd.at(i) == 1, it means to avoid that index
    }
    vector<double> indAndGrad = {(double)index, maxGrad};
    return indAndGrad;
}

/**
 * @brief  An overloaded version of findSteepestGrad that has no index to ignore
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  funDers : A vector of partial derivatives for the model that each take a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 */
vector<double> findSteepestGrad(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double (*)(double, vector<double>&)>& funDers, vector<double>& params) {
    vector<int> avoidInd(params.size(), 0); // There are as many parameters as partial derivatives
    return findSteepestGrad(data_x, data_y, fun, funDers, params, avoidInd);
}

/**
 * @brief  This function finds the parameters that minimize the squared error
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  funDers : A vector of partial derivatives for the model that each take a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @param  paramsLims : A vector of vectors of doubles that contains the limits for each parameter
 * @param  derTol : A double for the tolerance in the derivative of the squared error
 * @param  paramTol : A double for the tolerance in the parameters of the model
 * @param  checkOtherHalf : A boolean indicating if the bisection solver should check the other half if no zero was found in the checked half
 * @retval  iter : The number of iterations done
 */
int findFitParams(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double (*)(double, vector<double>&)>& funDers, vector<double>& params, vector<vector<double>>& paramsLims, double derTol, double paramTol, bool checkOtherHalf) {
    // Set up first iteration
    int iter = 0;
    vector<double> maxIndAndGrad = findSteepestGrad(data_x, data_y, fun, funDers, params);
    double maxGrad = maxIndAndGrad.at(1);
    int maxInd = maxIndAndGrad.at(0);
    int numOfAI = 0; // number of avoided indices
    vector<double> lastParams = params;
    vector<int> avoidInd(params.size(), 0);
    
    while (abs(maxGrad) >= derTol) {
        // Continue setting up this iteration
        iter++;
        lastParams = params;
        // Change parameter, and update vectors
        zeroDerFinder(data_x, data_y, fun, funDers.at(maxInd), params, maxInd, paramsLims.at(maxInd), derTol, paramTol, true, checkOtherHalf);
        // What if the program would keep switching back and forth, but to no avail?
        if ((lastParams == params) && (numOfAI == avoidInd.size()-1)) {
            maxGrad = 0; // Violate while loop condition to break from it
            cout << "  The " << params.size() << " parameters did not change in the past " << params.size() << " iterations." << endl;
            cout << "  Broke loop at i = " << iter << endl;
        }
        else {
            // Set up next iteration
            if (lastParams == params) {
                avoidInd.at(maxInd) = 1;
                numOfAI++;
            }
            else {
                for (int i = 0; i < avoidInd.size(); i++) {avoidInd.at(i) = 0;}
                numOfAI = 0;
            }
            maxIndAndGrad = findSteepestGrad(data_x, data_y, fun, funDers, params, avoidInd);
            maxGrad = maxIndAndGrad.at(1);
            maxInd = maxIndAndGrad.at(0);
        }
    }
    
    return iter;
}

/**
 * @brief  A function that finds the zero of the derivative of the squared error
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  funDer : A function for a derivative of the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @param  paramInd : A double for the index of the parameter to change
 * @param  paramLims : A vector of doubles that contains the limits for the parameter to be adusted
 * @param  derTol : A double for the tolerance in the derivative of the squared error
 * @param  paramTol : A double for the tolerance in the parameters of the model
 * @param  firstTime : A boolean indicating if the function was called for the first time (true) or was called recursively (false)
 * @param  checkOtherHalf : A boolean indicating if the bisection solver should check the other half if no zero was found in the checked half
 * @retval  value : The value of the derivate of the squared error at the set of parameters (post-adjusting)
 */
double zeroDerFinder(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), double (*funDer)(double, vector<double>&), vector<double>& params, double paramInd, vector<double>& paramLims, double derTol, double paramTol, bool firstTime, bool checkOtherHalf) {
    // Declare variables
    double lastArg;
    double value;
    // Determine which bound is lower (if this is wrong, the firstTime variable will account for it)
    params.at(paramInd) = paramLims.at(0);
    double neg = findDerErr(data_x, data_y, fun, funDer, params);
    params.at(paramInd) = paramLims.at(1);
    double pos = findDerErr(data_x, data_y, fun, funDer, params);
    if (neg < pos) {neg = paramLims.at(0); pos = paramLims.at(1);}
    else {neg = paramLims.at(1); pos = paramLims.at(0);}
    
    params.at(paramInd) = (paramLims.at(0) + paramLims.at(1))/2; // Start in middle for bisection method
    
    bool firstIter = true;
    do {
        value = findDerErr(data_x, data_y, fun, funDer, params);
        lastArg = params.at(paramInd);
        if (abs(value) > derTol) {
            // Force code to scan other half if this is the recursive call
            if (checkOtherHalf && firstIter && !firstTime) {value = -value; firstIter = false;}
            // Adjust bounds
            if (value < 0) {neg = params.at(paramInd); params.at(paramInd) = (neg + pos)/2;}
            else {pos = params.at(paramInd); params.at(paramInd) = (pos + neg)/2;}
        }
    }
    while ((abs(params.at(paramInd)-lastArg) > paramTol) && (abs(value) > derTol));
    
    // What if there is no zero in the scanned half? We must flip the bounds to look at the other half
    if ((checkOtherHalf && firstTime) && (abs(value) > derTol)) {
        double arg0 = params.at(paramInd);
        double val0 = value;
        vector<double> paramLimsFlip = {paramLims.at(1), paramLims.at(0)};
        value = zeroDerFinder(data_x, data_y, fun, funDer, params, paramInd, paramLims, derTol, paramTol, false, true);
        if (abs(val0) <= abs(value)) {params.at(paramInd) = arg0; value = val0;};
    }
    
    return value;
}
