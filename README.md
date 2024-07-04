# curveFit
Given a function and the data it is modeling, and its partial derivatives with respect to its adjustable parameters, the code can approximate the parameter values needed to minimize the squared error.

The curveFit.cpp and curveFit.hpp files should not be changed (unless there is an issue that emerges). The user should put the function and derivatives in the userFunction.cpp file (an example is given for reference).

To compile via the command line, run the following when in the `src` directory: `g++ userFunction.cpp curveFit.cpp -std=c++11`
