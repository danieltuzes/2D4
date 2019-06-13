//
// xautocorr_utils.cpp : contains the function definitions for xautocorr.cpp

#include "xautocorr_utils.h"

#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments
#include <fftw3.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <random>
#include <numeric>
#include <cmath>

namespace bpo = boost::program_options;

void randomfill(std::vector<double>& in)
{
    std::mt19937 generator(1000);
    std::uniform_real_distribution<double> distr(0, 2);
    int n = static_cast<int>(in.size());
    for (int i = 0; i < n; ++i)
        in[i] = distr(generator) / 10. + sin(distr(generator) / 10. + i / 50. * (2 * 3.14159));
}

void autocorr(const std::vector<double>& data_in, std::vector<double>& data_out)
{
    int n = static_cast<int>(data_in.size());
    for (int shift = 0; shift < n; ++shift)
    {
        data_out[shift] = 0;
        for (int i = 0; i < n; ++i)
            data_out[shift] += data_in[i] * data_in[(i + shift) % n];
    }
}

void normalize(std::vector<double>& data)
{
    double sum = std::accumulate(data.begin(), data.end(), (double)0);
    std::for_each(data.begin(), data.end(), [sum](double& val) {val /= sum; });
}

void normalize(std::vector<std::vector<double>>& in, double norm)
{
    double sum = 0;
    for (auto const& row : in)
        for (auto const& value : row)
            sum += value;
    sum /= norm;

    for (auto& row : in)
        for (auto& value : row)
            value /= sum;
}

std::vector<double> autoCorrelation1D(const std::vector<double>& data)
{
    static int size = (int)data.size();
    std::vector<double> ret(size, 0);
    fftw_complex* tmp;
    tmp = new fftw_complex[(size / 2 + 1) * 2];

    fftw_plan a = fftw_plan_dft_r2c_1d(size, const_cast<double*>(&data.front()), tmp, FFTW_ESTIMATE);
    fftw_execute(a);
    abs_val2(tmp, size);

    fftw_plan b = fftw_plan_dft_c2r_1d(size, tmp, &ret[0], FFTW_ESTIMATE);
    fftw_execute(b);

    fftw_destroy_plan(b);
    fftw_destroy_plan(a);
    for (int i = 0; i < size; ++i)
        ret[i] /= size * size;

    delete[] tmp;
    tmp = nullptr;
    return ret;
}

std::vector<double> FourierAbsVal1D(const std::vector<double>& data)
{
    static int lsize = (int)data.size();
    static int psize = lsize / 2 + 1;
    std::vector<double> ret(psize, 0);
    fftw_complex* tmp;
    tmp = new fftw_complex[(lsize / 2 + 1) * 2];

    fftw_plan a = fftw_plan_dft_r2c_1d(lsize, const_cast<double*>(&data.front()), tmp, FFTW_ESTIMATE);
    fftw_execute(a);

    fftw_destroy_plan(a);
    for (int i = 0; i < psize; ++i)
        ret[i] = tmp[i][0] * tmp[i][0] + tmp[i][1] * tmp[i][1];

    delete[] tmp;
    tmp = nullptr;
    return ret;
}

void abs_val2(fftw_complex* c, int size) //absolute value square for complex numbers
{
    for (int j = 0; j < size; ++j)
    {
        c[j][0] = c[j][0] * c[j][0] + c[j][1] * c[j][1];
        c[j][1] = 0;
    }
}

std::istream& operator >> (std::istream& in, disl& d)
{
    double posx = 0, posy = 0;
    int type;

    if (in)
        in >> posx;
    else
        return in;

    if (in)
        in >> posy;
    else
        return in;

    if (in)
        in >> type;
    else
        return in;

    std::get<0>(d) = posx;
    std::get<1>(d) = posy;
    std::get<2>(d) = type;
    return in;
}

std::ostream& operator << (std::ostream& o, const disl& d)
{
    o << std::get<0>(d) << "\t" << std::get<1>(d) << "\t" << std::get<2>(d);
    return o;
}

double dist(double coord1, double coord2)
{
    double diff = abs(coord1 - coord2); // distance wo periodic bc
    double diff_a = 1 - diff; // alternative distance
    return diff < diff_a ? diff : diff_a;
}

double distsq(const pair& a, const pair& b)
{
    double xdist = dist(a.first , b.first );
    double ydist = dist(a.second, b.second);
    return xdist * xdist + ydist * ydist;
}

double dist(const pair& a, const pair& b)
{
    return sqrt(distsq(a, b));
}

template <typename T> std::ostream& operator << (std::ostream& o, const std::vector<std::vector<T>>& map)
{
    for (auto const& map_line : map)
    {
        for (auto const& ddensity : map_line)
            o << ddensity << "\t";
        o << "\n";
    }
    return o;
}

template std::ostream& operator << (std::ostream&, const std::vector<std::vector<int>>&);    // template instantiation
template std::ostream& operator << (std::ostream&, const std::vector<std::vector<double>>&); // template instantiation


template <typename T> std::vector<std::vector<T>> gnuplot_flip(const std::vector<std::vector<T>>& map)
{
    auto size = map.size();
    std::vector<std::vector<T>> ret(size, std::vector<T>(size, 0));

    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            ret[size-1-i][j] = map[i][j];

    return ret;
}

template std::vector<std::vector<int>> gnuplot_flip(const std::vector<std::vector<int>>&);       // template instantiation
template std::vector<std::vector<double>> gnuplot_flip(const std::vector<std::vector<double>>&); // template instantiation
