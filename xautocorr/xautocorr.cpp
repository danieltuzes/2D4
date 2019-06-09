//
// init_config_gen.cpp : This file contains the 'main' function. Program execution begins and ends there.



#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments


#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <fftw3.h>
#include <numeric>

namespace bpo = boost::program_options;

enum method { na, ws, bc, gc }; // Wigner-Seitz, Box-counting and Gauss-convolving; na: default value

void randomfill(std::vector<double>&);
void autocorr(const std::vector<double>&, std::vector<double>&);
void normalize(std::vector<double>&);
void abs_val2(fftw_complex*, int); //absolute value square for complex numbers
std::vector<double> autoCorrelation1D(const std::vector<double>&);

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

    try {
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

    int resp = ((res / 2 + (__int64)1) * 2); // 1-gyel nagyobb, mint a res, ha az páros, és 2-vel egyébként

    std::vector<double> linedensity(resp, 0); // a sûrûségre, késõbb majd ez fogja tárolni a korrelációt is
    randomfill(linedensity); // fill in the vector data with random values for testing
    std::vector<double> correlation(resp, 0);
    auto correlation_fftw = autoCorrelation1D(linedensity);
    autocorr(linedensity, correlation);
    
    std::ofstream ofile("results.txt");
    if (!ofile)
    {
        std::cerr << "Cannot write file.";
        return 0;
    }
    else
        ofile << "# linedensity[i]\tcorrelation[i]/1026\tcorrelation_fftw[i]\n";

    for (int i = 0; i < resp; ++i)
    {
        ofile << linedensity[i] << "\t" << correlation[i]/1026 << "\t" << correlation_fftw[i] << "\n";
    }

    if (um == bc)
    {

    }

    return 0;
}

void randomfill(std::vector<double>& in)
{
    std::mt19937 generator(1000);
    std::uniform_real_distribution<double> distr(0, 2);
    int n = static_cast<int>(in.size()) - 1;
    for (int i = 0; i < n; ++i)
        in[i] = distr(generator);
}

void autocorr(const std::vector<double>& data_in, std::vector<double>& data_out)
{
    int n = static_cast<int>(data_in.size()) - 1;
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
    tmp = new fftw_complex[size];

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

void abs_val2(fftw_complex* c, int size) //absolute value square for complex numbers
{
    for (int j = 0; j < size; ++j)
    {
        c[j][0] = c[j][0] * c[j][0] + c[j][1] * c[j][1];
        c[j][1] = 0;
    }
}