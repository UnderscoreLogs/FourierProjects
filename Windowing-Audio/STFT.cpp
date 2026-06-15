#include "DFT.hpp"
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex> 
#include <vector>
#include <sstream>
#include <string>

#include <math.h>


std::vector<double> applyWindow(std::vector<double> data, std::vector<double> window) {

    size_t N = data.size();
    std::vector<double> fhatAppliedWindow(N);

    for(size_t i = 0; i < N; i++){
        fhatAppliedWindow[i] = data[i]*window[i]; //multiplying fhat by the windowing funnction
    }

    return fhatAppliedWindow;
}

std::vector<std::vector<std::complex<double>>> stfft(std::vector<double> data, int windowSampleSize, float windowOverlap, std::vector<double> window) {
    
    if(window.size() != windowSampleSize) {
        throw std::runtime_error("Window size mismatch");
    }

    int stepSize = round(windowSampleSize * (1 - windowOverlap));
    int dataSize = data.size();
    int numberOfWindows = floor((dataSize - windowSampleSize) / stepSize) + 1;

    std::cout << "STFFT Number of Windows: " << numberOfWindows << std::endl;

   int totalNeededSize = (numberOfWindows - 1) * stepSize + windowSampleSize;

    if (data.size() < totalNeededSize) {
        //data.resize(totalNeededSize, 0.0);
    }

    std::vector<std::vector<std::complex<double>>> FFTWindows;

    for(int i = 0; i < numberOfWindows; i++){

        int startIndex = i * stepSize;
        int endIndex = startIndex + windowSampleSize;
        //cout << "FFTWindows startIndex: " << startIndex << endl;
        //cout << "FFTWindows endIndex: " << endIndex << endl;    

        // Safety check
        if(endIndex > data.size()) break;

        std::vector<double> dataWindow(data.begin() + startIndex, data.begin() + endIndex);

        FFTWindows.push_back(fft(applyWindow(dataWindow, window)));

        //cout << "FFTWindows stepSize: " << stepSize << endl;
        //cout << "FFTWindows windowSampleSize: " << windowSampleSize << endl;      
    }

    return FFTWindows;
}