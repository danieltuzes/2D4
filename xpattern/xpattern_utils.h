//
// xpattern_utils.h : contains the function declarations for xpattern.cpp

#pragma once

#define VERSION_xpattern_utils 1.3
/* changelog
# 1.3
Xfname is implemented

# 1.2
* normalize is introduced
* disl reading in normalizes x coordinate
* dist uses normalie

# 1.1
* new sandbox methods to test analytical 1D Fourier transformation algorithm to test how well a Dirac-delta sum can be transformed
* new: c2r_map is implemented so that complex maps can be printed out for gnuplot

# 1.0
* first version working with the initial idea

*/

#define _USE_MATH_DEFINES

#pragma region includes

#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments
#include <complex>
#define M_i std::complex<double>(0,1)
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

// transfer a coordinate or distance into the range of [-0.5; 1.5]
void normalize(double& n);

// a dislocation is a (double, double, int) tuple for (posx,posy,type)
using disl = std::tuple<double, double, int>;

using pair = std::pair<double, double>;

// Wigner-Seitz, Box-counting and Gauss-smoothing; na: default value, wspn: Wigner-Seitz positive and negative, wsts: Wigner-Seitz total and signer, bc: box-counting, gc: Gauss-convolution, df: direct Fourier
enum class method { na, wspn, wsts, bc, gs, df };

#pragma region stream and file operators
std::istream& operator >> (std::istream&, disl&);
std::ostream& operator << (std::ostream&, const disl&);

// separate line elements with tab
template <typename T> std::ostream& operator << (std::ostream&, const std::vector<T>&);

// separate line elements with tab, and different lines with \n
template <typename T> std::ostream& operator << (std::ostream&, const std::vector<std::vector<T>>&);

// cuts away the folder names from the input string: in.substr(in.find_last_of('/') + 1);
std::string Xfname(std::string inputFnameWithPath);
#pragma endregion

// distance measured with periodic boundary condition on size [-0.5; 0.5)
double dist(double, double);

// the square of the 2D distance with periodic boundary condition on size [-0.5; 0.5) × [-0.5; 0.5)
double distsq(const pair&, const pair&);

// the 2D distance with periodic boundary condition on size [-0.5; 0.5) × [-0.5; 0.5)
double dist(const pair&, const pair&);

// measures how many points belong to each dislocation on a mesh with samp × samp number of points
void measure_area(std::vector<disl>&, size_t samp);

// measures how many points belong to each dislocation (in the range of beginIndex till endIndex)
// on a mesh with samp × samp number of points
void measure_area(std::vector<disl>& dislocs, size_t beginIndex, size_t endIndex, size_t samp);

// measures the density along a line in samp number of points; if samp > disl.size(), it supersamples and calculates average
// if pozNneg calculates the positive and negative density stored in linedensity_a and ~b, otherwise,
//      calculates the total and signed density to linedensity_a and ~_b
// if cid < dislocs.size()/2 it supposes that disl is positive
template <bool pozNneg> void measure_density(const std::vector<disl>& dislocs, size_t yid, size_t samp, std::vector<double>& linedensity_a, std::vector<double>& linedensity_b);

// measures the density line by line, see the overloaded function with argument std::vector<double>& linedensity_a and ~b
template <bool pozNneg> void measure_density(const std::vector<disl>& dislocs, size_t samp, std::vector<std::vector<double>>& map_a, std::vector<std::vector<double>>& map_b);

//searches the index cid of the dislocation closest to coordinate (posx,posy) which must be closer than lastdist, otherwise cid left untouched
void nearestDislIndex(const std::vector<disl>&, size_t& cid, double& lastdist, double& lastdistsq, double posx, double posy);

//searches the index cid of the dislocation (in the range of beginIndex till endIndex)
// closest to coordinate (posx,posy) which must be closer than lastdist, otherwise cid left untouched
void nearestDislIndex(const std::vector<disl>&, size_t beginIndex, size_t endIndex, size_t& cid, double& lastdist, double& lastdistsq, double posx, double posy);

// returns a vector with pairwise average unique values from map in increasing order
void gnuplotlevels(const std::vector<std::vector<double>>& map, std::string fname);

// casts a complex map into a real map to an easier use with gnuplot, makes no transformation
std::vector<std::vector<double>> c2r_map(const std::vector<std::vector<std::complex<double>>>& in);

#pragma region Fourier analysis and correlations

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

// calculates the fourier values of linedensity, add it to k
void addFourierNorm1D(std::vector<double>& F_absVal, const std::vector<double>& linedensity);

#pragma endregion

#pragma region sandbox playing area

// random adatokkal feltölt egy 1D-s tömböt, aztán megnézi a Fourier komponenseit és az autokorrelációját, hogy megbizonyosdjak róla, hogy jól tudom használni az FFTW-t, ugyanis az autokorrelációt én is jól ki tudom számolni
void test_fourier_and_corr();

// ha már megy az FFTW módszere, akkor megnézem, hogy dirac delták összegére ki tudom-e kézzel én is számolni a Fourier trafót. A Dirac delták összegét boxcounting segítségével elkenem, azon pedig FFTW már jól használható
void test_dirac(int);

// feltölti a vektort random értékekkel
void randomfill(std::vector<double>&);

#pragma endregion