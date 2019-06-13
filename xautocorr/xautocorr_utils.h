//
// xautocorr_utils.h : contains the function declaration for xautocorr.cpp

#pragma once
#define _USE_MATH_DEFINES

#pragma region includes

#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments
#include <fftw3.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <random>
#include <numeric>
#include <tuple>

#pragma endregion

namespace bpo = boost::program_options;

// a dislocation is a (double, double, int) tuple for (posx,posy,type)
using disl = std::tuple<double, double, int>;

using pair = std::pair<double, double>;

// Wigner-Seitz, Box-counting and Gauss-convolving; na: default value
enum method { na, ws, bc, gc };

// fill up the whole vector with some arbitrary content
void randomfill(std::vector<double>&);

// calculate the correlation by "hand"
void autocorr(const std::vector<double>&, std::vector<double>&);

// makes the sum of the vector values 1
void normalize(std::vector<double>&);

// makes the sum of the vector<vector> values 1
void normalize(std::vector<std::vector<double>>&, double = 1);

//absolute value square for complex numbers
void abs_val2(fftw_complex*, int);

// calculates the autocorrelation using fftw
std::vector<double> autoCorrelation1D(const std::vector<double>&);

// calculates the fourier values with fftw, and then calculates the absolute values of the elements
std::vector<double> FourierAbsVal1D(const std::vector<double>&);

std::istream& operator >> (std::istream&, disl&);
std::ostream& operator << (std::ostream&, const disl&);

// distance measured with periodic boundary condition on size [-0.5; 0.5)
double dist(double, double);

// the square of the 2D distance with periodic boundary condition on size [-0.5; 0.5) × [-0.5; 0.5)
double distsq(const pair&, const pair&);

// the 2D distance with periodic boundary condition on size [-0.5; 0.5) × [-0.5; 0.5)
double dist(const pair&, const pair&);

template <typename T> std::ostream& operator << (std::ostream&, const std::vector<std::vector<T>>&);
template <typename T> std::vector<std::vector<T>> gnuplot_flip(const std::vector<std::vector<T>>&);


// sandbox playing area
void test_fourier_and_corr(std::vector<double>&);