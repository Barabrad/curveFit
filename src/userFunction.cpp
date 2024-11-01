//
//  userFunction.cpp
//  curveFit
//
//  Created by Brad Barakat on 1/2/23.
//  Modified by Brad Barakat on 10/31/24.
//  This file will be modified for each new function to fit.

#include <ctime> // For time test in main()
#include <cstdio>
#include "curveFit.hpp"

double userFun(double x, vector<double>& params);
double userFun_da(double x, vector<double>& params);
double userFun_db(double x, vector<double>& params);
double userFun_dc(double x, vector<double>& params);
void printVector(vector<double> v);

int main() {
    vector<double> x = {1, 2, 4, 8, 16};
    vector<double> y = {0.7429400422, 1.793159669, 2.075100647, 2.103703704, 2.13058125};
    vector<double> userParams = {2.4, 0.5, -1.65}; // {a, b, c} in y = a*(1 - b*x^c)
    vector<double (*)(double, vector<double>&)> userFunDers = {userFun_da, userFun_db, userFun_dc};
    vector<double> aLims = {1,3};
    vector<double> bLims = {0,1};
    vector<double> cLims = {-3,-1};
    vector<vector<double>> userLims = {aLims, bLims, cLims};
    double derTol = 0.000001; // This is for the gradient. If you make it too small, the code may run for a long time.
    double paramTol = 0.00000005;
    // Ideally, the function is continuous within the parameter limits. If there are gaps or jumps,
    // the binary searcher may be tricked into checking the half that doesn't have the zero.
    // Thus, there will be an option to have the searcher check the other half if the zero is not found.
    bool checkOtherHalf = false; // If the fit is messy, setting this to true will likely almost double the runtime.
    
    clock_t t = clock();
    
    cout << "Old parameters: "; printVector(userParams);
    cout << "Old sqErr = " << findSqErr(x, y, userFun, userParams) << endl;
    int iterations = findFitParams(x, y, userFun, userFunDers, userParams, userLims, derTol, paramTol, checkOtherHalf);
    cout << "Iterations: " << iterations << endl;
    cout << "New parameters: "; printVector(userParams);
    cout << "New sqErr = " << findSqErr(x, y, userFun, userParams) << endl;
    
    t = clock() - t;
    cout << "Finished in " << ((double)t)/CLOCKS_PER_SEC << " seconds.\n";
    return 0;
}

/**
 * @brief  This function is the model that the user wants to fit to data
 * @param  x : A double for the independent variable
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  y : A double for the function output
 */
double userFun(double x, vector<double>& params) {
    // params = [a, b, c]
    // y = a*(1 - b*x^c)
    return params[0]*(1 - params[1]*pow(x,params[2]));
}

/**
 * @brief  This function is a partial derivative of the model that the user wants to fit to data
 * @param  x : A double for the independent variable
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  partial : A double for the partial derivative output
 */
double userFun_da(double x, vector<double>& params) {
    // params = [a, b, c]
    // y = a*(1 - b*x^c)
    // dy/da = 1 - b*x^c
    return 1 - params[1]*pow(x,params[2]);
}

/**
 * @brief  This function is a partial derivative of the model that the user wants to fit to data
 * @param  x : A double for the independent variable
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  partial : A double for the partial derivative output
 */
double userFun_db(double x, vector<double>& params) {
    // params = [a, b, c]
    // y = a*(1 - b*x^c)
    // dy/db = -a*x^c
    return -1*params[0]*pow(x,params[2]);
}

/**
 * @brief  This function is a partial derivative of the model that the user wants to fit to data
 * @param  x : A double for the independent variable
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  partial : A double for the partial derivative output
 */
double userFun_dc(double x, vector<double>& params) {
    // params = [a, b, c]
    // y = a*(1 - b*x^c)
    // dy/dc = -a*b*ln(x)*x^c
    return -1*params[0]*params[1]*log(x)*pow(x,params[2]);
}

/**
 * @brief  This function will print out a vector in a line (using std::cout), and then start a new line
 * @param  v : A vector of doubles that will be printed out
 * @retval  void : void
 */
void printVector(vector<double> v) {
    for (double vi : v) {cout << vi << " ";}
    cout << endl;
}
