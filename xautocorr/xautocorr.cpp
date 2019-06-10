//
// init_config_gen.cpp : This file contains the 'main' function. Program execution begins and ends there.


#pragma region header

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

#pragma endregion

namespace bpo = boost::program_options;

enum method { na, ws, bc, gc }; // Wigner-Seitz, Box-counting and Gauss-convolving; na: default value

void randomfill(std::vector<double>&); // fill up the whole vector with some arbitrary content
void autocorr(const std::vector<double>&, std::vector<double>&); // calculate the correlation by "hand"
void normalize(std::vector<double>&);  // makes the sum of the vector values 1
void abs_val2(fftw_complex*, int);    //absolute value square for complex numbers
std::vector<double> autoCorrelation1D(const std::vector<double>&); // calculates the autocorrelation using fftw
std::vector<double> FourierAbsVal1D(const std::vector<double>&);   // calculates the fourier values with fftw, and then calculates the absolute values of the elements
void test_fourier_and_corr(std::vector<double>&); // sandbox playing area
using disl = std::tuple<double, double, int>; // a dislocation is a (double, double, int) tuple for (posx,posy,type)
std::istream& operator >> (std::istream&, disl&);
std::ostream& operator << (std::ostream&, const disl&);

#pragma endregion

int main(int argc, char** argv)
{
#pragma region read in variables
    bpo::options_description requiredOptions("Required options"); // must be set from command line or file
    bpo::options_description optionalOptions("Optional options"); // must be set from command line or file

    requiredOptions.add_options()
        ("input-path,i", bpo::value<std::string>(), "The file name of the dislocation configuration file ending with .dconf. If it ends with .txt, a file containing the file names is expected: 1 filename per line.");

    optionalOptions.add_options()
        ("resolution,r", bpo::value<int>()->default_value(1024), "The dislocation density map will be evaluated in this many points and the number of the channels in autocorrelation will be the same. Please use sizes which conforms the suggestions of FFTW. (E.G. power of 2.)")
        ("method,m", bpo::value<std::string>()->default_value("bc"), "The method to calculate the autocorrelation. ws: Wigner-Seitz, bc: box-counting, gc: Gauss-convolving")
        ("sampling-multiplier,s", bpo::value<int>()->default_value(16), "This is a parameter for method ws. It tells how many times should be the mesh denser, on which the the area will be evaluated. (E.G.power of 2.)")
        ("half-width,w", bpo::value<double>()->default_value(0.125), "This is a parameter for method gc. It tells how wide the Gauss-distribution should with which the dirac-delta densities should be concolved.")
        ("output-foldername,o", bpo::value<std::string>()->default_value("xautocorr"), "In which folder should the initial conditions be stored. Symbol . means here.")
        ;

    bpo::options_description options; // the superior container of the options

    options.add(requiredOptions).add(optionalOptions).add_options()
        ("help", "show this help")
        ("hide-copyright,c", "hides the copyright notice from the standard output");

    bpo::variables_map vm; // the storage for the variables

    try
    {
        bpo::store(bpo::parse_command_line(argc, argv, options), vm, false);
    }
    catch (bpo::error& e)
    {
        std::cerr << e.what() << std::endl;
        exit(-1);
    }

    if (!vm.count("hide-copyright"))
    {
        std::cout << "init_config_gen from the 2D4 - a 2D discrete dislocation dynamics simulation program toolset.\n"
            "Copyright (C) Dániel Tüzes <tuzes@metal.elte.hu>\n";
    }
#pragma endregion

#pragma region process variables
    if (vm.count("help")) // if the user is curious 
    {
        std::cout << options << std::endl;
        exit(0);
    }
    else //
    {
        // input-path
        if (0 == vm.count("input-path")) // to test if N is present
        {
            std::cerr << "input-path is missing! Program terminates.\n";
            exit(-1);
        }


    }

    unsigned int res = vm["resolution"].as<int>(); // the resolution of the autocorrelation map
    unsigned int sampl = res * vm["sampling-multiplier"].as<int>(); // sampling rate for method ws
    method um = na; // the used method for the calculation
    double sigma = vm["half-width"].as<double>();
    if (vm["method"].as<std::string>() == "ws")
    {
        std::cout << "Sampling rate for Wigner-Seitz: " << sampl << "\n";
        um = ws;
    }
    else if (vm["method"].as<std::string>() == "bc")
        um = bc;
    else if (vm["method"].as<std::string>() == "gc")
    {
        std::cout << "Helf-width if the Gauss-distribution: " << sigma << "\n";
        um = gc;
    }

#pragma endregion

    std::vector<double> linedensity(res, 0); // a sûrûségre, késõbb majd ez fogja tárolni a korrelációt is
    std::string ifname = vm["input-path"].as<std::string>();
    std::ifstream ifile(ifname);
    if (!ifile)
    {
        std::cerr << "The program couldn't open " << ifname << " for reading. Program terminates." << std::endl;
        return 0;
    }
    
    std::vector<disl> dislocs; // container of the N number of dislocations

    for (disl tmp; ifile >> tmp; dislocs.push_back(tmp));

    if (um == bc) // box counting case
    {
        std::vector < std::vector<double>> map(res, std::vector<double>(res, 0));
        for (auto const& disl : dislocs)
        {

        }
    }

    return 0;
}

#pragma region function definition
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
    double posx = 0, posy=0;
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

#pragma region sandbox
void test_fourier_and_corr(std::vector<double>& in)
{
    const int res = 1024;
    randomfill(in);                 // fill in the vector data with random values for testing
    std::vector<double> correlation(res, 0);
    auto correlation_fftw = autoCorrelation1D(in); // calculate correlation with fftw
    auto FourierComp_fftw = FourierAbsVal1D(in);   // the fourier component of the function
    autocorr(in, correlation);                     // autocorrelation calculated by hand

    std::ofstream ofile("test.txt");
    if (!ofile)
    {
        std::cerr << "Cannot write file.";
        exit(0);
    }
    else
        ofile << std::setprecision(16) << "# linedensity[i]\tcorrelation[i]/" << res << "\tcorrelation_fftw[i]\tFourierComp_fftw[i]\n";

    for (unsigned int i = 0; i < res; ++i)
        ofile << in[i] << "\t" << correlation[i] / res << "\t" << correlation_fftw[i] << "\t" << FourierComp_fftw[i < res - i ? i : res - i] << "\n";

}
#pragma endregion

#pragma endregion