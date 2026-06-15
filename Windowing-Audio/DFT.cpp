#include "DFT.hpp"
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex> 
#include <vector>
#include <sstream>
#include <string>

#include <math.h>


std::vector<std::complex<double>> dft(const std::vector<double>& x) {

    double threshold = 0.0000001;

    size_t N = x.size(); //sets N to the number of data points
    std::vector<std::complex<double>> X(N); //a vector composed of ordered pairs of real and complex numbers

    for (size_t k = 0; k < N; k++){
        std::complex<double> sum(0.0, 0.0); //k is the freqency, adding the first row of the matrix

        for (size_t n = 0; n < N; n++) {

            double angle = (-2*M_PI*k*n)/N; //calculating the angle of the DFT equation -(2pi/N)*k*n

            std::complex<double> expTerm = exp(std::complex<double>(0,angle)); //calculaing the complex expontial part of the DFT: e^i*angle

            sum += x[n]*expTerm; //adding the DFT of xn to the sum;

        }
        if (abs(sum.real()) < threshold && abs(sum.imag()) < threshold) {
            X[k] = std::complex<double>(0.0, 0.0);
        } else if (abs(sum.real()) < threshold) {
            X[k] = std::complex<double>(0.0, sum.imag());
        } else if (abs(sum.imag()) < threshold) {
            X[k] = std::complex<double>(sum.real(), 0.0);
        } else {
            X[k] = sum; //inputting the sum of frequency k into the resultant vector
        }
            

    }

    return X; //returning the whole resultant vector
}

std::vector<std::complex<double>> idft(const std::vector<std::complex<double>>& x) {

    size_t N = x.size(); //sets N to the number of data points
    std::vector<std::complex<double>> X(N); //a vector composed of ordered pairs of real and complex numbers

    for (size_t n = 0; n < N; n++){
        std::complex<double> sum = {};

        for (size_t k = 0; k < N; k++) {

            double angle = (2*M_PI*k*n)/N; //calculating the angle of the iDFT equation (2pi/N)*k*n

            std::complex<double> expTerm = exp(std::complex<double>(0,angle)); //calculaing the complex expontial part of the iDFT: e^i*angle

            sum += (x[k]*expTerm); 

        }
       
        X[n] = sum / (double)N; //Xn is equal to the sum divided by N
    }

    return X; 
}
