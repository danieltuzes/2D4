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
    double diff = std::abs(coord1 - coord2); // distance wo periodic bc
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

void measure_area(std::vector<disl>& dislocs, size_t beginIndex, size_t endIndex, size_t samp)
{
    double lastdist = 1; // the last smallest distance multiplied with sqrt(2), always smaller than sqrt(0.5)
    double lastdistsq = 1; // the last smallest distance square
    size_t cid = std::numeric_limits<size_t>::max(); // the id of the last closest dislocation; hopefully, out of range
    for (size_t y = 0; y < samp; ++y) // samp = res * subs, it defines a more fine mesh
    {
        double meshy = (y + 0.5) / samp - 0.5; // the y coordinate of the investigated mesh point
        for (size_t x = 0; x < samp; ++x)
        {
            double meshx = (x + 0.5) / samp - 0.5; // the x coordinate of the investigated mesh point
            nearestDislIndex(dislocs, beginIndex, endIndex, cid, lastdist, lastdistsq, meshx, meshy);
            std::get<2>(dislocs[cid]) += 1;
            lastdist += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
            lastdistsq = lastdist * lastdist;
        }
        lastdist += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
        lastdistsq = lastdist * lastdist;
    }
}

template <bool pozNneg> void measure_density(const std::vector<disl>& dislocs, size_t yid, size_t samp, std::vector<double>& linedensity_a, std::vector<double>& linedensity_b)
{
    double area = pow(linedensity_a.size(), 2); // number of density points
    size_t subs = samp / linedensity_a.size(); // subsampling inside a density cell
    size_t size = dislocs.size(); // the number of dislocations

    double lastdist_a = 1; // the last smallest distance multiplied with sqrt(2), always smaller than sqrt(0.5)
    double lastdist_b = 1; // the last smallest distance multiplied with sqrt(2), always smaller than sqrt(0.5)
    double lastdistsq_a = 1; // the last smallest distance square
    double lastdistsq_b = 1; // the last smallest distance square
    size_t cid_a = std::numeric_limits<size_t>::max(); // the id of the last closest dislocation, hopefully out of scope
    size_t cid_b = std::numeric_limits<size_t>::max(); // the id of the last closest dislocation, hopefully out of scope

    for (size_t y = yid * subs; y < (yid + 1) * subs; ++y) // for the different subsampling y values
    {
        double meshy = (y + 0.5) / samp - 0.5; // the position in y of the subcell
        for (size_t x = 0; x < samp; ++x) // along a whole x line
        {
            double meshx = (x + 0.5) / samp - 0.5; // the position in x of the subcell
            if (pozNneg) // positive and negative map
            {
                nearestDislIndex(dislocs, 0, size / 2, cid_a, lastdist_a, lastdistsq_a, meshx, meshy); // cid is set to the nearest dislocation, only the first size/2 dislocation are positive
                nearestDislIndex(dislocs, size / 2, size, cid_b, lastdist_b, lastdistsq_b, meshx, meshy); // cid is set to the nearest dislocation, only the second size/2 dislocation are negative
                linedensity_a[x / subs] += area / std::get<2>(dislocs[cid_a]); // the nominator is chosen such a way that the integral of the map density is the particle number
                linedensity_b[x / subs] += area / std::get<2>(dislocs[cid_b]); // the nominator is chosen such a way that the integral of the map density is the particle number
                lastdist_b += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
                lastdistsq_b = lastdist_b * lastdist_b;
            }
            else // total and kappa map
            {
                nearestDislIndex(dislocs, cid_a, lastdist_a, lastdistsq_a, meshx, meshy); // cid is set to the nearest dislocation
                double disl_d = area / std::get<2>(dislocs[cid_a]); // the density caused by that dislocation
                linedensity_a[x / subs] += disl_d; // the nominator is chosen such a way that the integral of the map density is the particle number
                if (cid_a < dislocs.size() / 2)
                    linedensity_b[x / subs] += disl_d; // the nominator is chosen such a way that the integral of the map density is the particle number
                else
                    linedensity_b[x / subs] -= disl_d; // the nominator is chosen such a way that the integral of the map density is the particle number
            }
            lastdist_a += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
            lastdistsq_a = lastdist_a * lastdist_a;

        }
        lastdist_a += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
        lastdistsq_a = lastdist_a * lastdist_a;
        if (pozNneg)
        {
            lastdist_b += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
            lastdistsq_b = lastdist_b * lastdist_b;
        }
    }
}

// if !pozNneg, calculates total and signed density to map_a and ~_b: if cid < p it supposes that disl is positive
template <bool pozNneg> void measure_density(const std::vector<disl>& dislocs, size_t samp, std::vector<std::vector<double>>& map_a, std::vector<std::vector<double>>& map_b)
{
    for (size_t yid = 0; yid < map_a.size(); ++yid)
        measure_density<pozNneg>(dislocs, yid, samp, map_a[yid], map_b[yid]);
}

template void measure_density<false>(const std::vector<disl>&, size_t, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);

template void measure_density<true>(const std::vector<disl>&, size_t, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);

void nearestDislIndex(const std::vector<disl>& dislocs, size_t& cid, double& lastdist, double& lastdistsq, double posx, double posy)
{
    nearestDislIndex(dislocs, 0, dislocs.size(), cid, lastdist, lastdistsq, posx, posy);
}

void nearestDislIndex(const std::vector<disl>& dislocs, size_t beginIndex, size_t endIndex, size_t& cid, double& lastdist, double& lastdistsq, double posx, double posy)
{
    for (size_t i = beginIndex; i < endIndex; ++i) // find the nearest dislocation
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

std::vector<std::vector<double>> c2r_map(const std::vector<std::vector<std::complex<double>>>& in)
{
    int res = static_cast<int>(in.size());

    std::vector<std::vector<double>> ret(res, std::vector<double>(2 * res));
    for (int kx = 0; kx < res; ++kx)
        for (int ky = 0; ky < res; ++ky)
        {
            ret[kx][ky * 2] = in[kx][ky].real();
            ret[kx][ky * 2 + 1] = in[kx][ky].imag();
        }
    return ret;
}

# pragma region Fourier analysis and correlation
// self-written basic implementation
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

// self-written basic implementation
void fouriertr(const std::vector<double>& data_in, std::vector<fftw_complex>& data_out)
{
    int n = static_cast<int>(data_in.size());
    for (int k = 0; k < n; ++k)
    {
        for (int j = 0; j < n; ++j)
        {
            data_out[k][0] += data_in[j] * cos(2 * (M_PI * j * k) / n);
            data_out[k][1] -= data_in[j] * sin(2 * (M_PI * j * k) / n);
        }
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

// effective fftw implementation
std::vector<double> autoCorrelation1D(const std::vector<double>& data)
{
    static int lsize = (int)data.size();
    static int psize = lsize / 2 + 1;
    std::vector<double> ret(lsize, 0);
    fftw_complex* tmp;
    tmp = new fftw_complex[psize];

    fftw_plan a = fftw_plan_dft_r2c_1d(lsize, const_cast<double*>(&data.front()), tmp, FFTW_ESTIMATE);
    fftw_execute(a);
    abs_val2(tmp, psize);

    fftw_plan b = fftw_plan_dft_c2r_1d(lsize, tmp, &ret[0], FFTW_ESTIMATE);
    fftw_execute(b);

    fftw_destroy_plan(b);
    fftw_destroy_plan(a);
    for (int i = 0; i < lsize; ++i)
        ret[i] /= lsize * lsize;

    delete[] tmp;
    tmp = nullptr;
    return ret;
}

// effective fftw implementation
std::vector<double> FourierAbsVal1D(const std::vector<double>& data_in)
{
    const int lsize = (int)data_in.size(); // logical size
    const int psize = lsize / 2 + 1; // physical size
    std::vector<double> ret(psize, 0);
    fftw_complex* tmp;
    tmp = new fftw_complex[psize];

    fftw_plan a = fftw_plan_dft_r2c_1d(lsize, const_cast<double*>(&data_in.front()), tmp, FFTW_ESTIMATE);
    fftw_execute(a);

    fftw_destroy_plan(a);
    for (int i = 0; i < psize; ++i)
        ret[i] = tmp[i][0] * tmp[i][0] + tmp[i][1] * tmp[i][1];

    delete[] tmp;
    tmp = nullptr;
    return ret;
}

// effective fftw implementation
std::vector<fftw_complex> FourierTransform1D(const std::vector<double>& data_in)
{
    const int lsize = (int)data_in.size(); // logical size
    const int psize = lsize / 2 + 1; // physical size
    std::vector<fftw_complex> ret(psize);

    fftw_plan a = fftw_plan_dft_r2c_1d(lsize, const_cast<double*>(&data_in.front()), const_cast<fftw_complex*>(&ret.front()), FFTW_ESTIMATE);
    fftw_execute(a);

    fftw_destroy_plan(a);
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

void addFourierAbsValSq1D(std::vector<double>& F_absValSq, const std::vector<double>& linedensity)
{
    const int lsize = (int)linedensity.size(); // logical size
    const int psize = lsize / 2 + 1; // physical size, size of k
    fftw_complex* lineTransform;
    lineTransform = new fftw_complex[psize];

    fftw_plan a = fftw_plan_dft_r2c_1d(lsize, const_cast<double*>(&linedensity.front()), lineTransform, FFTW_ESTIMATE);
    fftw_execute(a);

    fftw_destroy_plan(a);
    for (int k = 0; k < psize; ++k)
        F_absValSq[k] += (lineTransform[k][0] * lineTransform[k][0] + lineTransform[k][1] * lineTransform[k][1]);
    delete[] lineTransform;
}

#pragma endregion

#pragma region sandbox
void randomfill(std::vector<double>& in)
{
    std::mt19937 generator(1000);
    std::uniform_real_distribution<double> distr(0, 2);
    int n = static_cast<int>(in.size());
    for (int i = 0; i < n; ++i)
        in[i] = distr(generator) / 10. + sin(distr(generator) / 10. + i / 50. * (2 * 3.14159));
}

void test_fourier_and_corr()
{
    const int res = 1024;
    std::vector<double> in(res);    // example data of uncorrelated values

    randomfill(in);                 // fill in the vector data with random values for testing
    std::vector<double> autocorrelation(res, 0);       // to store the autocorrelation calculated by hand
    autocorr(in, autocorrelation);                     // autocorrelation is calculated by hand

    std::vector<fftw_complex> fouriertransform(res);   // to store the fourier components calculated by hand
    fouriertr(in, fouriertransform);                   // fourier components are now calculated by hand

    auto FourierComp_fftw = FourierAbsVal1D(in);   // the abs value of the fourier components of the function
    auto correlation_fftw = autoCorrelation1D(in); // calculates the autocorrelation with fftw
    auto Fourier_fftw = FourierTransform1D(in);    // the fourier components of the function


    std::ofstream ofile("fourier_test.txt");
    if (!ofile)
    {
        std::cerr << "Cannot write file.";
        exit(0);
    }
    else
        ofile << "# linedensity[i]\tfouriertransform[0]\tfouriertransform[1]\tFourierTransform1D[0]\tFourierTransform1D[1]\tautocorrelation[i]/" << res << "\tcorrelation_fftw[i]\tFourierComp_fftw[i]\n" << std::setprecision(16);

    for (size_t i = 0; i < res; ++i)
        ofile << in[i] << "\t"
        << fouriertransform[i][0] << "\t"
        << fouriertransform[i][1] << "\t"
        << Fourier_fftw[i < res - i ? i : res - i][0] << "\t"
        << Fourier_fftw[i < res - i ? i : res - i][1] << "\t"
        << autocorrelation[i] / res << "\t"
        << correlation_fftw[i] << "\t"
        << FourierComp_fftw[i < res - i ? i : res - i] << "\n";
}

void test_dirac(int k)
{
    const int res = 1024;
    std::vector<double> in_bc(res);    // example data of patterned values; bc: like a box counting method
    for (int i = 0; i < res; ++i)
        in_bc[i] = sin((i + 0.5) / res * 2 * M_PI * k) + 1;

    std::vector<double> pos;
    for (int i = 0; i < res; ++i)
    {
        int n = static_cast<int>(in_bc[i] * 100);
        for (int j = 0; j < n; ++j)
            pos.push_back(i + static_cast<double>(j) / n);
    }

    std::vector<fftw_complex> fouriertransform_bc(res);     // to store the fourier components calculated by hand
    fouriertr(in_bc, fouriertransform_bc);                   // fourier components are now calculated by hand

    std::vector<std::complex<double>> fouriertransform_dirac(res);
    for (int k = 0; k < res; ++k)
    {
        for (int j = 0; j < static_cast<int>(pos.size()); ++j)
            fouriertransform_dirac[k] += exp(-2 * (M_PI * pos[j] * k) / res * M_i);
    }

    std::ofstream ofile_c("dirac_test_cont.txt");
    std::ofstream ofile_d("dirac_test_disc.txt");
    std::ofstream ofile_f("dirac_test_fourier.txt");
    if (!ofile_c || !ofile_d || !ofile_f)
    {
        std::cerr << "Cannot write file.";
        exit(0);
    }
    else
    {
        ofile_c << "# density[i], the (dislocation, whatever) density\n";
        ofile_d << "# position[i], the coordinate of a particle (dislocation)\n";
        ofile_f << "# Fourier components (real, tabulator, imag, tabulator, abssq) of the step-like, bc function and of the diract delta sum distribution\n";
    }
    for (auto intensity : in_bc)
        ofile_c << intensity << "\n";

    for (auto xcoord : pos)
        ofile_d << xcoord << "\n";

    for (int i = 0; i < res; ++i)
        ofile_f << fouriertransform_bc[i][0] << "\t"
        << fouriertransform_bc[i][1] << "\t"
        << fouriertransform_bc[i][0] * fouriertransform_bc[i][0] + fouriertransform_bc[i][1] * fouriertransform_bc[i][1] << "\t"
        << fouriertransform_dirac[i].real() << "\t"
        << fouriertransform_dirac[i].imag() << "\t"
        << std::norm(fouriertransform_dirac[i]) << "\n";

}
#pragma endregion