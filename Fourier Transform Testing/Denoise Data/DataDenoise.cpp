#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex> 
#include <vector>
#include <sstream>
#include <string>
#include <math.h>
#include <random>

#include <fstream>
#include <cstdio>

#include "DFT.hpp"

using namespace std;

vector<double> cleanData(double dt, int startTime, int endTime, int N) {

    //creates a vector of data points from a function with frequencies defined below
    double freq1 = 50;
    int freq2 = 120;
    int freq3 = 20;
    int freq4 = 450;

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

    int numDataPoints = (endTime - startTime)/dt;

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

vector<double> PSD(vector<complex<double>> fhat) {
    //calculates the Power Spectrum Density
    size_t n = fhat.size();
    vector<double> PSD(n);

    for (size_t i = 0; i < n; i++){
            PSD[i] = norm(fhat[i]) / n; //magnitude squared divided by n?
    }
    return PSD;
}

vector<complex<double>> filterDFT(vector<complex<double>> dftRaw, vector<double> PSD, int threshold) {
    //takes a generated DFT complex vector and PSD vector and filters the data so that all values below a given threshold in PSD equals 0

    size_t n = dftRaw.size();

    vector<complex<double>> filteredDFT = dftRaw; 

    //set all values below the threshold to 0 in the DFT
    for(size_t i = 0; i < n; i++ ) {
        if(PSD[i] < threshold) {
            filteredDFT[i] = (0,0);
        }
    }
    return filteredDFT;
}

int main() {

    //data collection settings
    int N = 2 << 12;    
    double dt = 1.0/pow(2,12);
    int startTime = 0;
    int endTime = 1;

    vector<double> dataClean = cleanData(dt, startTime, endTime, N);    //computing the clean data for comparison reasons

    vector<double> data = noisedData(dt, startTime, endTime, N);    //computing the noised data

    vector<complex<double>> fhat = fft(data);    //computing the DFT

    vector<double> psd = PSD(fhat);     //computing the PSD

    vector<complex<double>> filteredDFT = filterDFT(fhat, psd, 100); //filtering the DFT so that all values below the threshold equals 0, to be used to reconstruct the clean data

    vector<complex<double>> denoisedData = ifft(filteredDFT); //Inverse DFT of the filterDFT to reconstruct the orginal dataset minus the random noise

    int n = fhat.size(); //computes the size of the dataset

    //creating the frequency (x) axis for graphing
    vector<double> freq(n);
    for (int i = 0; i < n; i++) {
        freq[i] = i/(dt*n);
    }

   //Creating file with data to be graphed in excel, note: DELETE 'Graphs.txt' TO REGENERATE
    ofstream fout;
    fout.open("Graphs.txt");

    //Time vs Amplitude(clean), Time vs Amplidtude, and Hz vs. Power Spectrum Densisity Graphs Datasets
    fout << "Time" << '\t' << "Amplitude(clean)" << '\t' << '\t' << "Time" << '\t' << "Amplitude" << '\t' << '\t' << "Hz" << '\t' << "PSD" << '\t' << '\t' << "Time" << '\t' << "Amplitude(denoised)" << endl;
    cout << "Time" << '\t' << "Amplitude(clean)" << '\t' << '\t' << "Time" << '\t' << "Amplitude" << '\t' << '\t' << "Hz" << '\t' << "PSD" << '\t' << '\t' << "Time" << '\t' << "Amplitude(denoised)" << endl;

    for (int i = 0; i < n; i++) {
        double t = i*dt;
        // '\t' means tab which excel takes to mean a new collumn
        fout << t << '\t' << dataClean[i] << '\t' << '\t' << t << '\t' << data[i] << '\t' << '\t' << freq[i] << '\t' << psd[i] << '\t' << '\t' << t << '\t' << denoisedData[i].real() << endl; //send to file
        cout << t << '\t' << dataClean[i] << '\t' << '\t' << t << '\t' << data[i] << '\t' << '\t' << freq[i] << '\t' << psd[i] << '\t' << '\t' << t << '\t' << denoisedData[i].real() << endl; //send to console
    }
    fout.close();

    return 0;
}