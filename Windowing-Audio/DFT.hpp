#ifndef DFT_HPP
#define DFT_HPP

#include <vector>
#include <complex>

std::vector<std::complex<double>> dft(const std::vector<double>& x);

std::vector<std::complex<double>> idft(const std::vector<std::complex<double>>& x);

std::vector<std::complex<double>> fft(const std::vector<double> realSamples);

std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>> samples);

std::vector<double> hannFunction(const int N);

std::vector<std::vector<std::complex<double>>> stfft(std::vector<double> data, int windowSampleSize, float windowOverlap, std::vector<double> window);

#endif

