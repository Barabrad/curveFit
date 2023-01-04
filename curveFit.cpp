//
//  curveFit.cpp
//  curveFit
//
//  Created by Brad Barakat on 1/2/23.
//  Finished on 1/3/23.

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
 * @param  avoidInd : An index to ignore when iterating through the parameter vector
 * @retval  indAndGrad : A vector of doubles containing the index and the steepest gradient, respectively
 */
vector<double> findSteepestGrad(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double (*)(double, vector<double>&)>& funDers, vector<double>& params, int avoidInd) {
    int index = 0;
    double maxGrad = 0;
    for (int i = 0; i < funDers.size(); i++) {
        if (i != avoidInd) {
            double derErr_i = findDerErr(data_x, data_y, fun, funDers.at(i), params);
            if (abs(derErr_i) >= abs(maxGrad)) {
                index = i; maxGrad = derErr_i;
            }
        }
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
    return findSteepestGrad(data_x, data_y, fun, funDers, params, -1);
}

/**
 * @brief  This function finds the parameters that minimize the squared error
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  funDers : A vector of partial derivatives for the model that each take a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @param  derTol : A double for the tolerance in the derivative of the squared error
 * @param  paramTol : A double for the tolerance in the parameters of the model
 * @param  dx : A double for the step away from the parameter value to calculate an approximation of the derivative: f'(x) = (f(x + dx) - f(x - dx)) / (2*dx))
 * @retval  iter : The number of iterations done
 */
int findFitParams(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double (*)(double, vector<double>&)>& funDers, vector<double>& params, vector<vector<double>>& paramLims, double derTol, double paramTol, double dx) {
    // Set up first iteration
    int iter = 0;
    vector<double> maxIndAndGrad = findSteepestGrad(data_x, data_y, fun, funDers, params);
    double maxGrad = maxIndAndGrad.at(1);
    int maxInd = maxIndAndGrad.at(0);
    int lastInd = -1;
    vector<double> lastParams = params;
    int avoidInd = -1;
    
    while ((abs(maxGrad) >= derTol)) {
        // Continue setting up this iteration
        iter++;
        lastParams = params;
        // Change parameter, and update vectors
        minFinderSqErr(data_x, data_y, fun, params, maxInd, paramLims.at(maxInd), derTol, paramTol, 1, dx);
        // What if the program would keep switching back and forth, but to no avail
        if ((lastParams == params) && (lastInd != maxInd) && (lastInd >= 0)) {
            maxGrad = 0; // Violate while loop condition to break from it
            cout << "  Last Params: (" << lastParams.at(0) << ", " << lastParams.at(1) << ")" << endl;
            cout << "  New Params: (" << params.at(0) << ", " << params.at(1) << ")" << endl;
            cout << "  Broke loop at i = " << iter << endl;
        }
        else {
            // Set up next iteration
            lastInd = maxInd;
            maxIndAndGrad = findSteepestGrad(data_x, data_y, fun, funDers, params, avoidInd);
            maxGrad = maxIndAndGrad.at(1);
            maxInd = maxIndAndGrad.at(0);
            if ((lastParams == params) && (lastInd == maxInd)) {
                avoidInd = maxInd;
                maxIndAndGrad = findSteepestGrad(data_x, data_y, fun, funDers, params, avoidInd);
                maxGrad = maxIndAndGrad.at(1);
                maxInd = maxIndAndGrad.at(0);
            }
            else {
                avoidInd = -1;
            }
        }
    }
    
    return iter;
}

/**
 * @brief  An overloaded version of findFitParams() that leaves dx as the default value in the header file
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  funDers : A vector of partial derivatives for the model that each take a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @param  derTol : A double for the tolerance in the derivative of the squared error
 * @param  paramTol : A double for the tolerance in the parameters of the model
 * @retval  iter : The number of iterations done
 */
int findFitParams(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double (*)(double, vector<double>&)>& funDers, vector<double>& params, vector<vector<double>>& paramLims, double derTol, double paramTol) {
    return findFitParams(data_x, data_y, fun, funDers, params, paramLims, derTol, paramTol, DEFAULT_DX);
}

/**
 * @brief  This function attempts to minimize the squared error by changing only one parameter
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @param  paramInd : A double for the index of the parameter to change
 * @param  paramLims : A vector of doubles that contains the limits for the parameter to be adusted
 * @param  derTol : A double for the tolerance in the derivative of the squared error
 * @param  paramTol : A double for the tolerance in the parameters of the model
 * @param  firstTime : An int indicating if the function was called for the first time (1) or was called recursively (0)
 * @param  dx : A double for the step away from the parameter value to calculate an approximation of the derivative: f'(x) = (f(x + dx) - f(x - dx)) / (2*dx))
 * @retval  value : The value of the derivate of the squared error at the set of parameters (post-adjusting)
 */
double minFinderSqErr(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double>& params, double paramInd, vector<double>& paramLims, double derTol, double paramTol, int firstTime, double dx) {
    // Declare variables
    double lastArg;
    double value;
    // Assume the negative bound is the lower one (if this is wrong, the firstTime variable will account for it)
    double neg = paramLims.at(0), pos = paramLims.at(1);
    params.at(paramInd) = (paramLims.at(0) + paramLims.at(1))/2; // Start in middle
    if (neg < pos) {neg = paramLims.at(0); pos = paramLims.at(1);}
    else {neg = paramLims.at(1); pos = paramLims.at(0);}
    
    double value1, value2;
    do {
        params.at(paramInd) -= dx;
        value1 = findSqErr(data_x, data_y, fun, params);
        params.at(paramInd) += 2*dx;
        value2 = findSqErr(data_x, data_y, fun, params);
        value = (value2 - value1)/(2*dx);
        params.at(paramInd) -= dx;
        lastArg = params.at(paramInd);
        if (abs(value) > derTol) {
            if (value < 0) {neg = params.at(paramInd); params.at(paramInd) = (neg + pos)/2;}
            else {pos = params.at(paramInd); params.at(paramInd) = (pos + neg)/2;}
        }
    }
    while ((abs(params.at(paramInd)-lastArg) > paramTol) && (abs(value) > derTol));
    
    // What if there is no zero in the scanned half? We must flip the bounds to look at the other half
    if ((firstTime == 1) && (abs(value) < derTol)) {
        double arg0 = params.at(paramInd);
        double val0 = value;
        vector<double> paramLimsFlip = {paramLims.at(1), paramLims.at(0)};
        value = minFinderSqErr(data_x, data_y, fun, params, paramInd, paramLims, derTol, paramTol, 0, dx);
        if (val0 <= value) {params.at(paramInd) = arg0; value = val0;};
    }
    
    return value;
}

/**
 * @brief  An overloaded version of minFinderSqErr() that leaves dx as the default value in the header file
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @param  paramInd : A double for the index of the parameter to change
 * @param  paramLims : A vector of doubles that contains the limits for the parameter to be adusted
 * @param  derTol : A double for the tolerance in the derivative of the squared error
 * @param  paramTol : A double for the tolerance in the parameters of the model
 * @param  firstTime : An int indicating if the function was called for the first time (1) or was called recursively (0)
 * @retval  value : The value of the derivate of the squared error at the set of parameters (post-adjusting)
 */
double minFinderSqErr(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double>& params, double paramInd, vector<double>& paramLims, double derTol, double paramTol, int firstTime) {
    return minFinderSqErr(data_x, data_y, fun, params, paramInd, paramLims, derTol, paramTol, firstTime, DEFAULT_DX);
}
