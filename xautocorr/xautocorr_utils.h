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

// Wigner-Seitz, Box-counting and Gauss-smoothing; na: default value, wspn: Wigner-Seitz positive and negative, wsts: Wigner-Seitz total and signer, bc: box-counting, gc: Gauss-convolution
enum method { na, wspn, wsts, bc, gs };

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

// separate line elements with tab
template <typename T> std::ostream& operator << (std::ostream&, const std::vector<T>&);

// separate line elements with tab, and different lines with \n
template <typename T> std::ostream& operator << (std::ostream&, const std::vector<std::vector<T>>&);

// measures how many points belong to each dislocation on a mesh with samp × samp number of points
void measure_area(std::vector<disl>&, size_t samp);

// measures the density along a line in samp number of points; if samp > disl.size(), it supersamples and calculates average
// calculates total and signed density to linedensity_a and ~_b: if cid < p it supposes that disl is positive
// deduce the right template <bool> void measure_density based on the number of give arguments (see definition in cpp)
void measure_density(const std::vector<disl>& dislocs, size_t samp, std::vector<std::vector<double>>& map_a, std::vector<std::vector<double>>& map_b);

// measures the density along a line in samp number of points; if samp > disl.size(), it supersamples and calculates average
// deduce the right template <bool> void measure_density based on the number of give arguments (see definition in cpp)
void measure_density(const std::vector<disl>& dislocs, size_t samp, std::vector<std::vector<double>>& map);

//searches the index cid of the dislocation closest to coordinate (posx,posy) which must be closer than lastdist, otherwise cid left untouched
void nearestDislIndex(const std::vector<disl>&, size_t& cid, double& lastdist, double& lastdistsq, double posx, double posy);


// returns a vector with pairwise average unique values from map in increasing order
void gnuplotlevels(const std::vector<std::vector<double>>& map, std::string fname);


// sandbox playing area
void test_fourier_and_corr(std::vector<double>&);