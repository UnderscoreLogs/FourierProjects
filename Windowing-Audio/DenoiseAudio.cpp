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

#include <utility>

#include <sndfile.h>

#include "DFT.hpp"

using namespace std;

//reading .wav file, extracting data, IDK how this work, AI did it for me
tuple<vector<double>, int, int> readWav(const string& filename) {


    //AI code
    SF_INFO info;
    SNDFILE* file = sf_open(filename.c_str(), SFM_READ, &info);

    if (!file) {
        cerr << "Error opening file\n";
        exit(1);
    }

    vector<double> buffer(info.frames * info.channels);
    sf_read_double(file, buffer.data(), buffer.size());
    sf_close(file);

    vector<double> samples(info.frames);

    // Convert to mono if stereo
    for (int i = 0; i < info.frames; i++) {
        samples[i] = buffer[i * info.channels];
    }

    cout << "Sample rate: " << info.samplerate << "\n";
    cout << "Total frames: " << info.frames << "\n";

    return {samples, info.samplerate, info.frames};
}

//writng the result to .wav, AI also wrote this
void writeWav(const std::string& filename, const std::vector<double>& samples, int sampleRate) {

    SF_INFO info;
    info.channels = 1;                 // mono
    info.samplerate = sampleRate;
    info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;

    SNDFILE* file = sf_open(filename.c_str(), SFM_WRITE, &info);

    if (!file) {
        std::cerr << "Error opening output file\n";
        exit(1);
    }

    sf_write_double(file, samples.data(), samples.size());
    sf_close(file);
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

    float windowOverlap = 0.5; // how much each window overlaps the previous window *****USE 50% FOR NOW*****
    int windowSampleSize = 1 << 15; //how many samples per window
    double windowSampleStep = round(windowSampleSize*(1-windowOverlap));


    vector<double> window = hannFunction(windowSampleSize);

    double U = 0.0; //weight of window at a given point?
    for (double w : window) {
        U += w * w;
    }
    U /= window.size();


    /*----------------------------------------------Data collection----------------------------------------------------*/
    string fileName = {};
    cout << "Enter .wav file name: ";
    cin >> fileName;

    auto [data, fs, N] = readWav(fileName);

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
        if(norm[i] > 1e-8) {
            reconstructedCleanedData[i] /= norm[i];
        }
    }

    vector<double> finalReconCleanedData(reconstructedCleanedData.begin() + pad, reconstructedCleanedData.begin() + pad + data.size());

    /*----------------------------------------------Writing The Denoised Data to the .wav file----------------------------------------------------*/

    //writes the cleaned data to a .wav
    writeWav("clean.wav", finalReconCleanedData, fs);

    double freqResolutionFFT = (double)fs / (double)N;// What would be the regular FFT frequency resolution
    double freqResolutionSTFT = (double)fs / (double)windowSize; // the STFT frequency resolution

    cout << "FFT Δf = " << freqResolutionFFT << " Hz" << endl;
    cout << "STFT Δf = " << freqResolutionSTFT << " Hz" << endl;

    cout << "Denoise completed." << endl;

    return 0;
}