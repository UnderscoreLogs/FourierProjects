#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex> 
#include <vector>
#include <sstream>
#include <string>
#include <math.h>
#include <random>
#include <algorithm>
#include <stdexcept>

#include <fstream>
#include <cstdio>

#include "DFT.hpp"

using namespace std;

vector<double> cleanData(double dt, double startTime, double endTime, int N) {

    //creates a vector of data points from a function with frequencies defined below
    double freq1 = 20;
    int freq2 = 450;
    int freq3 = 50;
    int freq4 = 120;

    size_t numDataPoints = N;

    vector<double> fclean(numDataPoints); //y values of the original function just in case we need it.

    //making a vector of all the y values dt apart
    for(int n = 0; n < numDataPoints; n++) {
        double t = n*dt;
        double fCleanValue = {};

        if(t < ((endTime + startTime)/2.0)) {
            fCleanValue = sin(2*M_PI*freq1*t) + sin(2*M_PI*freq2*t);
        } else {
            fCleanValue = sin(2*M_PI*freq3*t) + sin(2*M_PI*freq4*t);
        }
        fclean[n] = fCleanValue;
    }
    return fclean;
}

vector<double> noisedData(double dt, int startTime, int endTime, int N) {

    //adds random noise the clean function defined in cleanData()

    int numDataPoints = N;

    vector<double> clean = cleanData(dt, startTime, endTime, N);

    vector<double> f(numDataPoints); //y values of the original function + random noise

    default_random_engine generator; //the random noise generator
    normal_distribution<double> noise(0,1);

    //making a vector of all the y values dt apart
    for(size_t n = 0; n < numDataPoints; n++) {

        //adding noise to the cleanData
        double fnoisy = clean[n] + 2.5 * noise(generator);
        f[n] = fnoisy;
    }
    return f;
}

vector<double> PSD(vector<complex<double>> fhat, double fs, double U) {
    //calculates the Power Spectrum Density
    size_t n = fhat.size();
    vector<double> PSD(n);

    for (size_t i = 0; i < n; i++){
            PSD[i] = norm(fhat[i]) / (fs*n*U); //magnitude squared divided by n?
    }
    return PSD;
}

vector<complex<double>> filterDFT(vector<complex<double>> dftRaw, vector<double> PSD, double threshold) {
    //takes a generated DFT complex vector and PSD vector and filters the data so that all values below a given threshold in PSD equals 0

    size_t n = dftRaw.size();

    vector<complex<double>> filteredDFT = dftRaw; 

    //set all values below the threshold to 0 in the DFT
    for(size_t i = 0; i < n; i++ ) {
        if(PSD[i] < threshold) {
            filteredDFT[i] = 0.0;
        }
    }
    return filteredDFT;
}

int main() {
    /*----------------------------------------------Data Collection Settings----------------------------------------------------*/

    //data collection settings    
    int N = 1 << 13; //number of datapoints (2^12)
    double startTime = 0;
    double endTime = 1;

    float windowOverlap = 0.5; // how much each window overlaps the previous window *****USE 50% FOR NOW*****
    int windowSampleSize = 1 << 10; //how many samples per window
    double windowSampleStep = round(windowSampleSize*(1-windowOverlap));

    int numberOfWindows = floor((N - windowSampleSize) / windowSampleStep) + 1;
    double windowLength = (endTime - startTime)/(1+(numberOfWindows-1)*(1-windowOverlap));
    double windowStep = windowLength*(1-windowOverlap);

    double dt = (endTime - startTime) / (double)N; //calculating how far each sample is (dt)

    double fs = 1.0/dt; //sampling frequency, samples per second

    vector<double> window = hannFunction(windowSampleSize);

    double U = 0.0; //weight of window at a given point?
    for (double w : window) {
        U += w * w;
    }
    U /= window.size();


    /*----------------------------------------------Data collection----------------------------------------------------*/

    vector<double> dataClean = cleanData(dt, startTime, endTime, N);    //computing the clean data for comparison reasons

    vector<double> data = noisedData(dt, startTime, endTime, N);    //computing the noised data

    /*----------------------------------------------STFFT Computations----------------------------------------------------*/

    int pad = windowSampleStep;
    vector<double> paddedData;

    //start padding
    paddedData.insert(paddedData.end(), pad, 0.0);

    //original data
    paddedData.insert(paddedData.end(), data.begin(), data.end());

    //back padding
    paddedData.insert(paddedData.end(), pad, 0.0);


    vector<vector<complex<double>>> fhatWindows = stfft(paddedData, windowSampleSize, windowOverlap, window);

    vector<vector<double>> psdOfWindows(fhatWindows.size());

    for(int i = 0; i < fhatWindows.size(); i++) {
        psdOfWindows[i] = PSD(fhatWindows[i], fs, U);
    }

    /*------------------------------------------------Average PSD of Windows --------------------------------------------*/
    //Welch's Method
    vector<double> averagePSD(psdOfWindows[0].size());

    for(int i = 0; i < psdOfWindows.size(); i++) {
        for(int j = 0; j < averagePSD.size(); j++)
        averagePSD[j] += psdOfWindows[i][j] / psdOfWindows.size();
    }

    /*----------------------------------------------FFT Filtering----------------------------------------------------*/

    size_t n = fhatWindows[0].size(); //computes the size of the dataset

    vector<vector<complex<double>>> filteredFhatWindows(fhatWindows.size());
    
    for(int i = 0; i < fhatWindows.size(); i++) {

        double maxPSD = *max_element(psdOfWindows[i].begin(), psdOfWindows[i].end()); //finds the bigest PSD value

        double threshold = maxPSD * 0.2; //sets threshold to 20% of the largest value

        filteredFhatWindows[i] = filterDFT(fhatWindows[i], psdOfWindows[i], threshold);
    }

    vector<vector<complex<double>>> denoisedWindows(fhatWindows.size());

    for(int i = 0; i < fhatWindows.size(); i++) {
        denoisedWindows[i] = ifft(filteredFhatWindows[i]);
        denoisedWindows[i].resize(windowSampleSize); //resizing to the real length of the windows, because the FFT adds dummy samples because it needs to be a power of 2
    }

    /*----------------------------------------------Data Reconstruction----------------------------------------------------*/

    int totalSize = paddedData.size();
    vector<double> reconstructedCleanedData(totalSize, 0.0); //holds data

    
    vector<double> norm(totalSize, 0.0); //keeps track of "weight"

    int windowSize = denoisedWindows[0].size(); //REDUNDENT I think
    int stepSize = round(windowSampleSize * (1-windowOverlap)); //calculates each step

    cout << "Recon windowSize: " << windowSize << endl; //debugging

    //window.resize(windowSize);

    for(int k = 0; k < denoisedWindows.size(); k++) {
        int startIndex = k * stepSize;

        for(int n = 0; n < windowSampleSize; n++) {

            if(startIndex + n >= totalSize) break;

            double value = real(denoisedWindows[k][n]);

            reconstructedCleanedData[startIndex + n] += value * window[n];
            norm[startIndex + n] += window[n] * window[n];
        }
    }


    //Normalizing the data
    for (int i = 0; i < totalSize; i++){
        if(norm[i] > 0.1) {
            reconstructedCleanedData[i] /= norm[i];
        }
    }

    vector<double> finalReconCleanedData(reconstructedCleanedData.begin() + pad, reconstructedCleanedData.begin() + pad + data.size());


    /*----------------------------------------------Printing PSD----------------------------------------------------*/

    //creating the frequency (x) axis for graphing
    vector<double> freq(n);
    for (size_t i = 0; i < n; i++) {
        freq[i] = i * fs / n;
    }

   //Creating file with data to be graphed in excel, note: DELETE 'Graphs.txt' TO REGENERATE
    ofstream fout;
    fout.open("PSDs.txt");

    // Header
    fout << "Hz";
    for (int w = 0; w < psdOfWindows.size(); w++) {
        fout << '\t' << "t: " << startTime + windowStep*w << " -> " << startTime + windowStep*w + windowLength;
    }
    fout << "\tAverage PSD:";
    fout << '\n';

    // Data
    for (size_t k = 0; k < n/2; k++) {

        fout << freq[k];

        for (int w = 0; w < psdOfWindows.size(); w++) {

            fout << '\t' << psdOfWindows[w][k];
        }

        fout << '\t' << averagePSD[k];
        fout << '\n';
    }

    fout.close();

    /*----------------------------------------------Printing The Denoised Data----------------------------------------------------*/

    ofstream fout2("ReconstructedSignal.txt");

    // Header
    fout2 << "Time\tClean\tNoisy\tReconstructed\n";

    // Data
    for(int i = 0; i < data.size(); i++) {
        double time = i * dt;

        fout2 << time << '\t'
            << dataClean[i] << '\t'
            << data[i] << '\t'
            << finalReconCleanedData[i] << '\n';
    }

    fout2.close();

    /*----------------------------------------------Printing The Denoised Data Windows----------------------------------------------------*/
    
    ofstream fout3("DenoisedWindows.txt");

    // Header
    fout3 << "Time";
    for (int w = 0; w < psdOfWindows.size(); w++) {
        fout3 << '\t' << "t: " << startTime + windowStep*w << " -> " << startTime + windowStep*w + windowLength;
    }
    fout3 << '\n';

    // Data
    for(int i = 0; i < denoisedWindows[0].size(); i++) {
        double time = i * dt;
        fout3 << time;

        for(int j = 0; j < denoisedWindows.size(); j++) {
            fout3 << '\t' << denoisedWindows[j][i].real();
        }

        fout3 << '\n';

    }

    fout3.close();


    cout << "Done." << endl;

    return 0;
}