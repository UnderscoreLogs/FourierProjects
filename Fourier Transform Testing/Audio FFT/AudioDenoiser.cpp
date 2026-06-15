#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex> 
#include <vector>
#include <string>
#include <algorithm>

#include <sndfile.h>

using namespace std;

//Data resizer, makes the size of the dataset be a power of 2 (required for FFT), AI a
size_t nextPowerOfTwo(size_t n) {
    size_t p = 1;
    while (p < n) {
        p *= 2;
    }
    return p;
}

// -------------------FFT & IFFT-----------------------------------------------------------------------
std::vector<std::complex<double>> fftRecursive(std::vector<std::complex<double>> samples){
    
    size_t N = samples.size();

    //cout << "FFT called with N = " << N << endl; //debug for checking length of the called vector

    size_t M = N/2; //half the size of samples

    if(N <= 0 || (N & (N-1)) != 0) { //checks dataset if it is a power of 2, if it is not then it gives an errror
        std::cout << "Error: Dataset is not a power of 2." << std::endl;
        return {};
    }

    if(N <= 1) { //preventing something..?
        return samples;
    }

    std::vector<std::complex<double>> even(M);
    std::vector<std::complex<double>> odd(M);

    for(size_t i = 0; i < N/2; i++){
        even[i] = samples[2*i];
        odd[i] = samples[2*i+1];
    }

    even = fftRecursive(even);
    odd = fftRecursive(odd);

    std::vector<std::complex<double>> fftOutput(N);

    for(size_t k = 0; k<M; k++) {

        double angle = (-2*M_PI*k)/N; //calculating the angle of the DFT equation -(2pi/N)*k*n

        std::complex<double> expTerm = exp(std::complex<double>(0,angle));     

        std::complex<double> twiddleOdd = odd[k] * expTerm;

        fftOutput[k] = even[k] + twiddleOdd;
        fftOutput[k+M] = even[k] - twiddleOdd;
    }
    return fftOutput;
}

std::vector<std::complex<double>> fft(std::vector<double> realSamples) {

    //converts the real samples to complex and then returns the actual fft
    size_t N = realSamples.size();

    std::vector<std::complex<double>> complexSamples(N);

    for(size_t i = 0; i<N; i++) {
        complexSamples[i] += realSamples[i];
    }

    return fftRecursive(complexSamples);
}

std::vector<std::complex<double>> ifftRecursive(std::vector<std::complex<double>> samples){
    
    size_t N = samples.size();

    //cout << "FFT called with N = " << N << endl; //debug for checking length of the called vector

    size_t M = N/2; //half the size of samples

    if(N <= 0 || (N & (N-1)) != 0) { //checks dataset if it is a power of 2, if it is not then it gives an errror
        std::cout << "Error: Dataset is not a power of 2." << std::endl;
        return {};
    }

    if(N <= 1) { //preventing something..?
        return samples;
    }

    std::vector<std::complex<double>> even(M);
    std::vector<std::complex<double>> odd(M);

    for(size_t i = 0; i < N/2; i++){
        even[i] = samples[2*i];
        odd[i] = samples[2*i+1];
    }

    even = ifftRecursive(even);
    odd = ifftRecursive(odd);

    std::vector<std::complex<double>> fftOutput(N);

    for(size_t k = 0; k<M; k++) {

        double angle = (2*M_PI*k)/N; //calculating the angle of the DFT equation -(2pi/N)*k*n

        std::complex<double> expTerm = exp(std::complex<double>(0,angle));     

        std::complex<double> twiddleOdd = odd[k] * expTerm;

        fftOutput[k] = (even[k] + twiddleOdd);
        fftOutput[k+M] = (even[k] - twiddleOdd);
    }
    return fftOutput;
}

std::vector<std::complex<double>> ifft(std::vector<std::complex<double>> samples) {

    size_t N = samples.size();

    std::vector<std::complex<double>> x = ifftRecursive(samples);

    for(size_t i = 0; i<N; i++) {
        x[i] /= (double)N;
    }

    return x;
}

//reading .wav file, extracting data, IDK how this work, AI did it for me
vector<double> readWav(const string& filename) {


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

    return samples;
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

//calculating PSD (used for analyis)
vector<double> PSD(vector<complex<double>> fhat) {
    //calculates the Power Spectrum Density
    size_t n = fhat.size();
    vector<double> PSD(n);

    for (size_t i = 0; i < n; i++){
            PSD[i] = norm(fhat[i]) / n; //magnitude squared divided by n?
    }
    return PSD;
}

//takes a generated DFT complex vector and PSD vector and filters the data so that all values below a given threshold in PSD equals 0
vector<complex<double>> filterDFT(vector<complex<double>> dftRaw, vector<double> PSD, int threshold) {

    size_t n = dftRaw.size();

    vector<complex<double>> filteredDFT = dftRaw; 

    //set all values below the threshold to 0 in the DFT
    for(size_t i = 0; i < n; i++ ) {
        if(PSD[i] < threshold) {
            filteredDFT[i] = complex<double>(0.0,0.0);
        }
    }
    return filteredDFT;
}

int main() {

    string fileName = {};
    cout << "Enter .wav file name: ";
    cin >> fileName;

    vector<double> samples = readWav(fileName); //turns .wav data into a vector

    cout << "Loaded " << samples.size() << " samples\n"; // states number of samples loaded

    size_t originalSize = samples.size(); //calculates the original size of the data
    size_t fftSize = nextPowerOfTwo(originalSize); //turns the size of the data into a power of 2

    samples.resize(fftSize, 0.0); //resizes data by adding 0.0s

    vector<complex<double>> FFT = fft(samples); //calculates the fft

    vector<double> psd = PSD(FFT); //calculates the Power Spectrum Density

    double maxPSD = *max_element(psd.begin(), psd.end()); //finds the bigest PSD value
    double threshold = maxPSD * 0.01; //sets threshold to 1% of the largest value

    vector<complex<double>> filteredFFT = filterDFT(FFT, psd, threshold); //filters the FFT based on the PSD and a given threshold

    double Fs = 44100.0; // sample rate
    size_t N = filteredFFT.size();

    cout << "Detected Frequencies: " << endl;
    for(size_t i = 0; i < N/2; i++) { // only first half
        if(abs(filteredFFT[i]) > 0.0) {
            double freq = i * Fs / N;
            cout << freq << " Hz" << endl;
        }
    }

    vector<complex<double>> denoised = ifft(filteredFFT); //takes the IFFT to reconstruct the data without noise

    //turns the denoised data into a real vector of the same size as the original data
    vector<double> realSignal(originalSize); 

    for (size_t i = 0; i < originalSize; i++) {
        realSignal[i] = denoised[i].real();
    }

    //writes the cleaned data to a .wav
    writeWav("clean.wav", realSignal, 44100);

    cout << "Denoise completed." << endl;

    return 0;
}