#include "DFT.hpp"
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex> 
#include <vector>
#include <sstream>
#include <string>

#include <math.h>

std::vector<double> hannFunction(int N) {

    std::vector<double> hann(N);

    for(int n = 0; n < N; n++) {
        hann[n] = pow(sin((M_PI*n)/(N-1)),2);
    }

    return hann;
}