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
    double xdist = dist(a.first, b.first);
    double ydist = dist(a.second, b.second);
    return xdist * xdist + ydist * ydist;
}

double dist(const pair& a, const pair& b)
{
    return sqrt(distsq(a, b));
}

template <typename T> std::ostream& operator << (std::ostream& o, const std::vector<T>& line)
{
    for (auto const& value : line)
        o << value << "\t";

    return o;
}

template <typename T> std::ostream& operator << (std::ostream& o, const std::vector<std::vector<T>>& map)
{
    for (auto const& map_line : map)
        o << map_line << "\n";

    return o;
}

template std::ostream& operator << (std::ostream&, const std::vector<std::vector<int>>&);    // template instantiation
template std::ostream& operator << (std::ostream&, const std::vector<std::vector<double>>&); // template instantiation


void measure_area(std::vector<disl>& dislocs, size_t samp)
{
    size_t size = dislocs.size();
    double lastdist = 1; // the last smallest distance multiplied with sqrt(2), always smaller than sqrt(0.5)
    double lastdistsq = 1; // the last smallest distance square
    size_t cid = size; // the id of the last closest dislocation
    for (size_t y = 0; y < samp; ++y) // samp = res * subs, it defines a more fine mesh
    {
        double meshy = (y + 0.5) / samp - 0.5; // the y coordinate of the investigated mesh point
        for (size_t x = 0; x < samp; ++x)
        {
            double meshx = (x + 0.5) / samp - 0.5; // the x coordinate of the investigated mesh point
            nearestDislIndex(dislocs, cid, lastdist, lastdistsq, meshx, meshy);
            std::get<2>(dislocs[cid]) += 1;
            lastdist += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
            lastdistsq = lastdist * lastdist;
        }
        lastdist += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
        lastdistsq = lastdist * lastdist;
    }
}

// if !is_single, calculates total and signed density to linedensity_a and ~_b: if cid < p it supposes that disl is positive
template <bool is_single> void measure_density(const std::vector<disl>& dislocs, size_t yid, size_t samp, std::vector<double>& linedensity_a, std::vector<double>& linedensity_b)
{
    size_t size = linedensity_a.size(); // number of density points
    size_t subs = samp / size; // subsampling inside a density cell

    double lastdist = 1; // the last smallest distance multiplied with sqrt(2), always smaller than sqrt(0.5)
    double lastdistsq = 1; // the last smallest distance square
    size_t cid = dislocs.size(); // the id of the last closest dislocation

    for (size_t y = yid * subs; y < (yid + 1) * subs; ++y) // for the different subsampling y values
    {
        double meshy = (y + 0.5) / samp - 0.5; // the position in y of the subcell
        for (size_t x = 0; x < samp; ++x) // along a whole x line
        {
            double meshx = (x + 0.5) / samp - 0.5; // the position in x of the subcell
            nearestDislIndex(dislocs, cid, lastdist, lastdistsq, meshx, meshy); // cid is set to the nearest dislocation
            if (is_single)
                linedensity_a[x / subs] += (double)(size * size) / std::get<2>(dislocs[cid]); // the nominator is chosen such a way that the integral of the map density is the particle number
            else
            {
                linedensity_a[x / subs] += (double)(size * size) / std::get<2>(dislocs[cid]); // the nominator is chosen such a way that the integral of the map density is the particle number
                if (cid < dislocs.size() / 2)
                    linedensity_b[x / subs] += (double)(size * size) / std::get<2>(dislocs[cid]); // the nominator is chosen such a way that the integral of the map density is the particle number
                else
                    linedensity_b[x / subs] -= (double)(size * size) / std::get<2>(dislocs[cid]); // the nominator is chosen such a way that the integral of the map density is the particle number
            }
            lastdist += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
            lastdistsq = lastdist * lastdist;
        }
        lastdist += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
        lastdistsq = lastdist * lastdist;
    }
}

void measure_density(const std::vector<disl>& dislocs, size_t yid, size_t samp, std::vector<double>& linedensity_a, std::vector<double>& linedensity_b)
{
    measure_density<false>(dislocs, yid, samp, linedensity_a, linedensity_b);
}

void measure_density(const std::vector<disl>& dislocs, size_t yid, size_t samp, std::vector<double>& linedensity_a)
{
    measure_density<true>(dislocs, yid, samp, linedensity_a, linedensity_a);
}

// if !is_single, calculates total and signed density to map_a and ~_b: if cid < p it supposes that disl is positive
template <bool is_single> void measure_density(const std::vector<disl>& dislocs, size_t samp, std::vector<std::vector<double>>& map_a, std::vector<std::vector<double>>& map_b)
{
    for (size_t yid = 0; yid < map_a.size(); ++yid)
        measure_density<is_single>(dislocs, yid, samp, map_a[yid], map_b[yid]);
}

void measure_density(const std::vector<disl>& dislocs, size_t samp, std::vector<std::vector<double>>& map_a, std::vector<std::vector<double>>& map_b)
{
    measure_density<false>(dislocs, samp, map_a, map_b);
}

void measure_density(const std::vector<disl>& dislocs, size_t samp, std::vector<std::vector<double>>& map)
{
    measure_density<true>(dislocs, samp, map, map);
}

template void measure_density<false>(const std::vector<disl>&, size_t, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
template void measure_density<true>(const std::vector<disl>&, size_t, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);

void nearestDislIndex(const std::vector<disl>& dislocs, size_t& cid, double& lastdist, double& lastdistsq, double posx, double posy)
{
    size_t size = dislocs.size(); // the size of the dislocation array
    for (size_t i = 0; i < size; ++i) // find the nearest dislocation
    {
        double disly = std::get<1>(dislocs[i]); // the y coordinate of the dislocation
        double disty = dist(disly, posy); // the distance of the point from the dislocation in y direction
        if (disty < lastdist) // the distance in y must not be larger than the distance
        {
            double dislx = std::get<0>(dislocs[i]); // the x coordinate of the dislocation
            double distx = dist(dislx, posx); // the distance of the point from the dislocation in x direction
            if (distx < lastdist && // the distance in x must not be larger than the distance
                distx * distx + disty * disty < lastdistsq) // eventually it must be the closest
            {
                lastdistsq = distx * distx + disty * disty;
                lastdist = sqrt(lastdistsq);
                cid = i;
            }
        }
    }
}

void gnuplotlevels(const std::vector<std::vector<double>>& map, std::string fname)
{
    size_t size = map.size(); // the 2D map
    std::vector<double> levels(size * size, 0); // stores the 2D map in 1D line
    auto it = levels.begin();

    for (const auto& line : map)
        it = std::copy_n(line.begin(), size, it); // copy the 2D map into the 1D line
    std::sort(levels.begin(), levels.end());

    it = std::unique(levels.begin(), levels.end());
    levels.resize(std::distance(levels.begin(), it));
    std::ofstream ofile(fname);
    if (!ofile)
    {
        std::cerr << "Cannot write " << fname << ". Program terminates." << std::endl;
        exit(-1);
    }
    ofile << "set cntrparam levels discrete ";
    if (levels.size() == 1)
    {
        ofile << 0;
        std::cerr << "Warning: too few densities are found. cntrparam levels cannot be set." << std::endl;
        return;
    }
    for (it = levels.begin(); it != levels.end() - 2; ++it)
        ofile << (*it + *(it + 1)) / 2 << ", ";
    ofile << (*it + *(it + 1)) / 2;

}

#pragma region sandbox
void test_fourier_and_corr(std::vector<double>& in)
{
    const int res = 1024;
    randomfill(in);                 // fill in the vector data with random values for testing
    std::vector<double> correlation(res, 0);
    auto correlation_fftw = autoCorrelation1D(in); // calculate correlation with fftw
    auto FourierComp_fftw = FourierAbsVal1D(in);   // the fourier component of the function
    autocorr(in, correlation);                     // autocorrelation is calculated by hand

    std::ofstream ofile("test.txt");
    if (!ofile)
    {
        std::cerr << "Cannot write file.";
        exit(0);
    }
    else
        ofile << std::setprecision(16) << "# linedensity[i]\tcorrelation[i]/" << res << "\tcorrelation_fftw[i]\tFourierComp_fftw[i]\n";

    for (size_t i = 0; i < res; ++i)
        ofile << in[i] << "\t" << correlation[i] / res << "\t" << correlation_fftw[i] << "\t" << FourierComp_fftw[i < res - i ? i : res - i] << "\n";

}
#pragma endregion