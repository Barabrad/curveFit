//
//  curveFit.hpp
//  curveFit
//
//  Created by Brad Barakat on 1/2/23.
//  Finished on 1/3/23. Modified on 1/4/23.

#ifndef curveFit_hpp
#define curveFit_hpp

#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

/**
 * @brief  This function finds the squared error for a data set and modeling function
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  sqErr : A double for the squared error
 */
double findSqErr(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double>& params);

/**
 * @brief  This function finds the derivative of the squared error for a data set and modeling function
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  funDer : A function for a derivative of the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  derErr : A double for the derivative of the squared error
 */
double findDerErr(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), double (*funDer)(double, vector<double>&), vector<double>& params);

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
vector<double> findSteepestGrad(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double (*)(double, vector<double>&)>& funDers, vector<double>& params, vector<int>& avoidInd);

/**
 * @brief  An overloaded version of findSteepestGrad that has no indices to ignore
 * @param  data_x : A vector of doubles for the data's independent variable
 * @param  data_y : A vector of doubles for the data's dependent variable
 * @param  fun : A function for the model that takes a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  funDers : A vector of partial derivatives for the model that each take a double (independent variable) and a vector of doubles (adjustable constants)
 * @param  params : A vector of doubles that contains the constants to adjust
 */
vector<double> findSteepestGrad(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double (*)(double, vector<double>&)>& funDers, vector<double>& params);

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
 * @retval  iter : The number of iterations done
 */
int findFitParams(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), vector<double (*)(double, vector<double>&)>& funDers, vector<double>& params, vector<vector<double>>& paramsLims, double derTol, double paramTol);

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
 * @param  firstTime : An int indicating if the function was called for the first time (1) or was called recursively (0)
 * @retval  value : The value of the derivate of the squared error at the set of parameters (post-adjusting)
 */
double zeroDerFinder(vector<double>& data_x, vector<double>& data_y, double (*fun)(double, vector<double>&), double (*funDer)(double, vector<double>&), vector<double>& params, double paramInd, vector<double>& paramLims, double derTol, double paramTol, int firstTime);

#endif /* curveFit_hpp */
