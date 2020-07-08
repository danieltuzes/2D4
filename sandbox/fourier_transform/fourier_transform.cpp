// fourier_transform.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define VERSION_fourier_transform 1

/*changelog
# 1
Added functionality: if more than 1 filename is give, program calculates the fourier transformation with dir=1, calcualtes the abs value of the map and add them together.

# 0
First working version.
*/

#pragma region includes
#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#define M_i std::complex<double>(0,1)
#include <fftw3.h>

#include <iomanip>
#include <fstream>
#include <vector>

#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments
#pragma endregion

namespace bpo = boost::program_options;

// Fourier-transforms the 2D array in, returns only the non-trivial part
std::vector<std::vector<std::complex<double>>> fftw_r2c_2d(std::vector<std::vector<double>> in);

// Fourier-transforms the 2D array in, expands it to be hermitian
std::vector<std::vector<double>> fftw_c2r_2d(std::vector<std::vector<std::complex<double>>> in);

// prints out the values in x, y, z order, not in matrix-shape
void print_xyz(std::ostream& o, std::vector<std::vector<std::complex<double>>> in);

// prints out the values in x, y, z order, not in matrix-shape
void print_xyz(std::ostream& o, std::vector<std::vector<double>> in);

// prints out the values z in matrix format
void print_z(std::ostream& o, std::vector<std::vector<std::complex<double>>> in);

// prints out the values z in matrix format
void print_z(std::ostream& o, std::vector<std::vector<double>> in);

// prints out the map to the proper file or stream
template <typename T> void print_results(bpo::variables_map vm, T oVals);

int main(int argc, char** argv)
{
#pragma region read in variables
    bpo::options_description requiredOptions("Required options"); // must be set from command line or file
    bpo::options_description optionalOptions("Optional options"); // must be set from command line or file

    requiredOptions.add_options()
        ("input-fname,i", bpo::value<std::vector<std::string>>()->multitoken(), "The file name(s) of the data to be Fourier-transformed. If more than 1 are given, their Fourier norm are added up.");

    optionalOptions.add_options()
        ("direction,d", bpo::value<int>()->default_value(1), "1: r2c; -1: c2r") // c2r extends the values so that output is real; r2r is such c2r where only r values for c are read
        ("output-fname,o", bpo::value<std::string>(), "The output filename, stdout will be used if not set.")
        ;

    bpo::options_description options; // the superior container of the options

    options.add(requiredOptions).add(optionalOptions).add_options()
        ("help", "shows this help")
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
        std::cout << "fourier_transform (version " << VERSION_fourier_transform << ") from 2D4 - a 2D discrete dislocation dynamics simulation program toolset.\n"
            "Copyright (C) Dániel Tüzes <tuzes@metal.elte.hu>\n";
        std::cout.unsetf(std::ios_base::floatfield);
    }
#pragma endregion

#pragma region read in and process data
    std::vector<std::string> ifnames = vm["input-fname"].as<std::vector<std::string>>();

    std::vector<std::vector<double>> oRvals;                // the onput real values
    std::vector<std::vector<std::complex<double>>> oCvals;  // the onput complex values
    int dir = vm["direction"].as<int>();                    // direction of transformation: 1=r2c; -1=c2r

    for (auto ifname : ifnames)
    {
        std::ifstream ifile(ifnames[0]);
        if (!ifile)
        {
            std::cerr << "cannot open " << ifname << ". Program terminates." << std::endl;
            exit(1);
        }
        std::vector<std::vector<double>> iRvals;                // the input real values
        std::vector<std::vector<std::complex<double>>> iCvals;  // the input complex values

        // read in values
        {
            for (std::string line; std::getline(ifile, line); )
            {
                std::istringstream iss(line);
                if (dir == 1)
                {
                    std::vector<double> linevalR;                // the values stored in one line of the file
                    for (double val; iss >> val; linevalR.push_back(val));
                    iRvals.push_back(linevalR);
                }
                else
                {
                    std::vector<std::complex<double>> linevalC; // the values stored in one line of the file
                    for (double valR, valC; iss >> valR >> valC; linevalC.emplace_back(valR, valC));
                    iCvals.push_back(linevalC);
                }
            }
        }

        for (size_t linenumber = 1; (dir == 1 && linenumber < iRvals.size()) || (dir == -1 && linenumber < iCvals.size()); ++linenumber)
            if ((dir == 1 && iRvals[linenumber].size() < iRvals[0].size()) || (dir == -1 && iCvals[linenumber].size() < iCvals[0].size()))
            {
                std::cerr << "Error: line lengths are not equal: line number 0 and " << linenumber << ". Program terminates." << std::endl;
                exit(1);
            }

        if (ifnames.size() == 1)
        {
            if (dir == 1)
                oCvals = fftw_r2c_2d(iRvals);
            else
                oRvals = fftw_c2r_2d(iCvals);
        }
        else
        {
            auto fourier_map = fftw_r2c_2d(iRvals);
            for (auto& line : fourier_map)
                for (auto& val : line)
                    val = std::norm(val);

            if (oCvals.size() == 0)
                oCvals = fourier_map;
            else
                for (size_t y = 0; y < fourier_map.size(); ++y)
                    for (size_t x = 0; x < fourier_map[0].size(); ++x)
                        oCvals[y][x] += fourier_map[y][x];
        }
    }
    

#pragma endregion

    if (dir == 1)
        print_results(vm, oCvals);
    else
        print_results(vm, oRvals);

    std::cout << "Done.\n";
}

std::vector<std::vector<std::complex<double>>> fftw_r2c_2d(std::vector<std::vector<double>> in)
{
    size_t n0 = in.size();          // the 1st dimension; memory is not censecutive in this direction
    size_t n1 = in[0].size();       // the 2nd dimension; memory __is__ censecutive in this direction

    double* in_fftw;                // to store the input values in consecutive order
    in_fftw = (double*) fftw_malloc (sizeof(double) * n0 * n1);
    for (size_t i = 0; i < n0; ++i)
        for (size_t j = 0; j < n1; ++j)
            in_fftw[i * n1 + j] = in[i][j];

    fftw_complex* out_fftw;         // to store the output values in consecutive order
    out_fftw = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n0 * (n1 / 2 + 1));

    fftw_plan r2c_plan = fftw_plan_dft_r2c_2d(static_cast<int>(n0), static_cast<int>(n1), in_fftw, out_fftw, FFTW_ESTIMATE);
    fftw_execute(r2c_plan);

    std::vector<std::vector<std::complex<double>>> ret(n0, std::vector<std::complex<double>>(n1 / 2 + 1, 0));

    for (size_t i = 0; i < n0; ++i)
        for (size_t j = 0; j < n1 / 2 + 1; ++j)
            ret[i][j] = out_fftw[i * (n1 / 2 + 1) + j][0] + out_fftw[i * (n1 / 2 + 1) + j][1] * M_i;

    fftw_destroy_plan(r2c_plan);
    fftw_free(in_fftw);
    fftw_free(out_fftw);
    return ret;
}

std::vector<std::vector<double>> fftw_c2r_2d(std::vector<std::vector<std::complex<double>>> in)
{
    size_t n0 = in.size();              // the 1st dimension; memory is not censecutive in this direction
    size_t n1 = (in[0].size() - 1) * 2; // the 2nd dimension; memory __is__ censecutive in this direction

    fftw_complex* in_fftw;              // to store the input values in consecutive order
    in_fftw = new fftw_complex[n0 * in[0].size()];
    for (size_t i = 0; i < n0; ++i)
        for (size_t j = 0; j < in[0].size(); ++j)
        {
            in_fftw[i * in[0].size() + j][0] = in[i][j].real();
            in_fftw[i * in[0].size() + j][1] = in[i][j].imag();
        }

    double* out_fftw;                   // to store the output values in consecutive order
    out_fftw = new double[n0 * n1];

    fftw_plan c2r_plan = fftw_plan_dft_c2r_2d(static_cast<int>(n0), static_cast<int>(n1), in_fftw, out_fftw, FFTW_ESTIMATE);
    fftw_execute(c2r_plan);

    std::vector<std::vector<double>> ret(n0, std::vector<double>(n1, 0));

    for (size_t i = 0; i < n0; ++i)
        for (size_t j = 0; j < n1; ++j)
            ret[i][j] = out_fftw[i * n1 + j];

    fftw_destroy_plan(c2r_plan);
    delete[] in_fftw;
    delete[] out_fftw;
    return ret;
}

void print_xyz(std::ostream& o, std::vector<std::vector<std::complex<double>>> in)
{
    for (size_t i = 0; i < in.size(); ++i)
    {
        for (size_t j = 0; j < in[0].size(); ++j)
            o << j << "\t" << i << "\t" << in[i][j].real() << "\t" << in[i][j].imag() << "\n";

        o << "\n";   // each line must end with 0
    }
}

void print_xyz(std::ostream& o, std::vector<std::vector<double>> in)
{
    for (size_t i = 0; i < in.size(); ++i)
    {
        for (size_t j = 0; j < in[0].size(); ++j)
            o << j << "\t" << i << "\t" << in[i][j] << "\n";

        o << "\n";
    }
}

void print_z(std::ostream& o, std::vector<std::vector<std::complex<double>>> in)
{
    for (size_t i = 0; i < in.size(); ++i)
    {
        for (size_t j = 0; j < in[0].size(); ++j)
            o << in[i][j].real() << "\t" << in[i][j].imag() << "\t";

        o << "\n";
    }
}

void print_z(std::ostream& o, std::vector<std::vector<double>> in)
{
    for (size_t i = 0; i < in.size(); ++i)
    {
        for (size_t j = 0; j < in[0].size(); ++j)
            o << in[i][j] << "\t";

        o << "\n";
    }
}

template <typename T> void print_results(bpo::variables_map vm, T oVals)
{
    std::ofstream o_xyz, o_z;               // results in both formats will be printed
    std::string ofname;

    bool use_of = vm.count("output-fname"); // use output file
    if (use_of)
    {
        ofname = vm["output-fname"].as<std::string>();

        auto ofname_xyz = "xyz_" + ofname;
        o_xyz.open(ofname_xyz);
        if (!o_xyz)
        {
            std::cerr << "Cannot create " << ofname_xyz << ". Program terminates." << std::endl;
            exit(1);
        }

        auto ofname_z = "z_" + ofname;
        o_z.open(ofname_z);
        if (!o_z)
        {
            std::cerr << "Cannot create " << ofname_z << ". Program terminates." << std::endl;
            exit(1);
        }
    }

    if (use_of)
    {
        print_xyz(o_xyz, oVals);
        print_z(o_z, oVals);
    }
    else
        print_xyz(std::cout, oVals);
}