//
//  userFunction.cpp
//  curveFit
//
//  Created by Brad Barakat on 1/2/23.
//  Finished on 1/3/23.

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
    vector<double> y = {1.076539549, 2.015280133, 2.293402279, 2.381693867, 2.42240895};
    vector<double> userParams = {2.4, 0.5, -1.65}; // {a, b, c} in y = a*(1 - b*x^c)
    vector<double (*)(double, vector<double>&)> userFunDers = {userFun_da, userFun_db, userFun_dc};
    vector<double> aLims = {1,3};
    vector<double> bLims = {0,1};
    vector<double> cLims = {-2,-1};
    vector<vector<double>> userLims = {aLims, bLims, cLims};
    double derTol = 0.000001; // This is for the gradient. If you make it too small, the code will run for a long time.
    double paramTol = 0.00000005;
    double dx = 0.00001;
    
    clock_t t = clock();
    
    cout << "Old parameters: "; printVector(userParams);
    cout << "Old sqErr = " << findSqErr(x, y, userFun, userParams) << endl;
    int iterations = findFitParams(x, y, userFun, userFunDers, userParams, userLims, derTol, paramTol, dx);
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
    double y = params[0]*(1 - params[1]*pow(x,params[2]));
    return y;
}

/**
 * @brief  This function is a partial derivative of the model that the user wants to fit to data
 * @param  x : A double for the independent variable
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  dyda : A double for the partial derivative output
 */
double userFun_da(double x, vector<double>& params) {
    // params = [a, b]
    // y = a*(1 - b*x^c) -> dy/da = 1 - b*x^c
    double dyda = 1 - params[1]*pow(x,params[2]);
    return dyda;
}

/**
 * @brief  This function is a partial derivative of the model that the user wants to fit to data
 * @param  x : A double for the independent variable
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  dydb : A double for the partial derivative output
 */
double userFun_db(double x, vector<double>& params) {
    // params = [a, b]
    // y = a*(1 - b*x^c) -> dy/db = -a*x^c
    double dydb = -1*params[0]*pow(x,params[2]);
    return dydb;
}

/**
 * @brief  This function is a partial derivative of the model that the user wants to fit to data
 * @param  x : A double for the independent variable
 * @param  params : A vector of doubles that contains the constants to adjust
 * @retval  dydb : A double for the partial derivative output
 */
double userFun_dc(double x, vector<double>& params) {
    // params = [a, b]
    // y = a*(1 - b*x^c) -> dy/dc = -a*b*ln(x)*x^c
    double dydc = -1*params[0]*params[1]*log(x)*pow(x,params[2]);
    return dydc;
}

/**
 * @brief  This function is a partial derivative of the model that the user wants to fit to data
 * @param  v : A vector of doubles that will be printed out
 * @retval  void
 */
void printVector(vector<double> v) {
    for (double vi : v) {
        cout << vi << " ";
    }
    cout << endl;
}
