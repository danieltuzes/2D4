//
// init_config_gen.cpp : This file contains the 'main' function. Program execution begins and ends there.

#define VERSION_init_config_gen "1.0"

#pragma region header

#define _USE_MATH_DEFINES
#include <boost/program_options.hpp> // to read in program call arguments
#include <boost/program_options/options_description.hpp> // to add descriptions of the program call arguments
#include <boost/random/mersenne_twister.hpp> // random number generator
#include <boost/random/uniform_real_distribution.hpp> // uniform distribution generator

#include <iostream>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <utility>


namespace bpo = boost::program_options;
namespace br = boost::random;

using pair = std::pair<double, double>;
double interpolate_y(pair, pair, double);

#pragma endregion

int main(int argc, char** argv)
{
#pragma region reading in variables
    bpo::options_description requiredOptions("Required options"); // must be set from command line or file
    bpo::options_description optionalOptions("Optional options"); // must be set from command line or file

    requiredOptions.add_options()
        ("N,N", bpo::value<int>(), "The number of dislocations to generate. Must be even, because N/2 number of positive and negative dislocations will be created.");

    optionalOptions.add_options()
        ("seed-start,S", bpo::value<int>()->default_value(1000), "An integer used as an initial seed value for the random number generator.")
        ("seed-end,E", bpo::value<int>(), "An integer used as the last seed value for the random number generator, seed-end > seed_start must hold. If set, seed-end - seed-start number of initial configurations will be created.")
        ("pattern-strength,A", bpo::value<double>()->default_value(0), "Instead of a uniform distribution, the probability density function will be\n1 + A * sin(x * n * 2 pi) with A = pattern-strength for rho_+ and\n1 - A * sin(x * w * 2 pi) for rho_-. A must be in [-1:1].")
        ("linear-wavenumber,n", bpo::value<int>()->default_value(3), "The number of waves in the simulation area [-0.5:0.5].")
        ("pattern-type,T", bpo::value<std::string>()->default_value("s"), "The pattern type in the density distribution.")
        ("unsorted,U", "If set, dislocations will not printed out in order starting with positive Burger's vector and highest value in y, but with alternating Burger's vector and uncorrelated x and y coordinates.")
        ("output-fnameprefix,o", bpo::value<std::string>()->default_value(""), "In which folder should the initial conditions be stored. Symbol ./ means here.")
        ("bare,B", "If set, filenames will not contain the value of the parameter N.")
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
        std::cout << "init_config_gen (version " << VERSION_init_config_gen << ") from 2D4 - a 2D discrete dislocation dynamics simulation program toolset.\n"
            "Copyright (C) Dániel Tüzes <tuzes@metal.elte.hu>\n";
    }
#pragma endregion

#pragma region process variables

    if (vm.count("help")) // if the user is curious 
    {
        std::cout << options << std::endl;
        exit(0);
    }

#pragma region N
    size_t N; // the total number of dislocations
    if (0 == vm.count("N")) // to test if N is present
    {
        std::cerr << "N is missing! Program terminates.\n";
        exit(-1);
    }
    else
    {
        N = vm["N"].as<int>();
        std::cout << "N =\t" << N << std::endl;
    }

    if (N % 2 != 0)  // to test if N is even
    {
        std::cerr << "N must be even. Program terminates." << std::endl;
        exit(-1);
    }
#pragma endregion

#pragma region seed-end, seed-start
    
    int seed_start = seed_start = vm["seed-start"].as<int>(); // the first seed value for the random number generator engine
    int seed_end; // the last seed value for the random number generator engine
    
    std::cout << "seed-start =\t" << seed_start << std::endl;
    if (vm.count("seed-end") != 0)
        seed_end = vm["seed-end"].as<int>();
    else
        seed_end = seed_start + 1;
    std::cout << "seed-end =\t" << seed_end << std::endl;

    if (seed_end <= seed_start)
    {
        std::cerr << "seed-end >= seed-start must hold. Program terminates." << std::endl;
        exit(-1);
    }

#pragma endregion

# pragma region pattern
    double A = vm["pattern-strength"].as<double>(); // pattern strength
    if (abs(A) > 1)
    {
        std::cerr << "Too strong pattern-strength paramter A = " << A << ". A must be in [-1:1]. Program terminates." << std::endl;
        exit(-1);
    }
    else
        std::cout << "A =\t" << A << std::endl;

    std::vector<std::pair<double, double>> cdfr_p, cdfr_n; // x,y values of the cumulative distribution function of rho_p and rho_n
    if (A != 0)
    {
        int res = (int)sqrt(N) * 10 + 1000; // resolution for the inverse function lookup; heuristic guess
        double n = vm["linear-wavenumber"].as<int>(); // linear-wavenumber
        for (int i = 0; i <= res; ++i)
        {
            double x = (double)i / res - 0.5; // the x position in dimensionless units, x∈[-0.5:0.5]
            double yrp = 0.5 + x + A * pow(sin(x * n * 2 * M_PI), 2) / (n * 2 * M_PI); // the y values for the rho_p distribution
            double yrn; // the y values for the rho_n distribution
            if (vm["pattern-type"].as<std::string>() == "s") // signed, i.e. κ has patterns
                yrn = 0.5 + x - A * pow(sin(x * n * 2 * M_PI), 2) / (n * 2 * M_PI);
            else // total, i.e. ρ_t has patterns
                yrn = yrp;
            cdfr_p.push_back(std::pair<double, double>(x, yrp));
            cdfr_n.push_back(std::pair<double, double>(x, yrn));
        }
    }
#pragma endregion

#pragma region output path, filenames and sorting

    std::string of(vm["output-fnameprefix"].as<std::string>()); // the output filename prefix (potentially inculding foldername)
    std::cout << "output-fnameprefix =\t" << of << std::endl;
    
    bool unsorted = false; // leave the output unsorted?
    if (vm.count("unsorted") != 0)
        unsorted = true;
    std::cout << "unsorted =\t" << unsorted << std::endl;
    
    bool bare = false; // shall the output name contain the number of dislocations?
    if (vm.count("bare") != 0)
        bare = true;
    std::cout << "bare=\t" << bare << std::endl;

#pragma endregion

#pragma endregion

#pragma region generate and write out configuration
    for (int seed_val = seed_start; seed_val < seed_end; ++seed_val) // generate configurations with seeds in the range of [seed-start; seed-end]; seed-end has been set at label: seed-end set
    {
        std::string ofname = of + std::to_string(seed_val);
        if (!bare)
            ofname += "_" + std::to_string(N);
        ofname += ".dconf"; // output filename; the file is inside a folder

        std::ofstream ofile(ofname); // the filestream
        if (!ofile) // evaluates to false if file cannot be opened
        {
            std::cerr << "Cannot write " << ofname << ". Program terminates." << std::endl;
            exit(-1);
        }
        std::cout << "Generating configuration with seed value " << seed_val << "." << std::endl;
        ofile.precision(std::numeric_limits<double>::max_digits10); // print out every digits

        br::mt19937 engine(seed_val); // Mersenne twister is a good and expensive random number generator
        br::uniform_real_distribution<> distr(-0.5, 0.5); // uniform distribution on [-0.5; 0.5)
        using disl = std::tuple<double, double, int>; // a dislocation is a (double, double, int) tuple for (posx,posy,type)
        std::vector<disl> dislocs; // container of the N number of dislocations

        //std::ofstream deb_conf_ofile("debug_conf.txt");
        for (int n = 0; n < N; ++n) // generate the N number of dislocations
        {
            double x = distr(engine);
            if (A != 0)
            {
                x += 0.5;
                auto cmp = pair(0, x);
                std::vector<pair>::iterator lb;
                if (n % 2 == 0)
                    lb = std::lower_bound(cdfr_n.begin(), cdfr_n.end(), cmp, [](pair lhs, pair rhs) {return (lhs.second < rhs.second); });
                else
                    lb = std::lower_bound(cdfr_p.begin(), cdfr_p.end(), cmp, [](pair lhs, pair rhs) {return (lhs.second < rhs.second); });

                auto ub = lb - 1;
                //deb_conf_ofile << x << "\t" << lb->first << "\t" << lb->second << "\t" << ub->first << "\t" << ub->second << "\t" << interpolate_y(*lb, *ub, x) << "\n";
                x = interpolate_y(*lb, *ub, x);
            }
            dislocs.push_back(disl(x, distr(engine), (n % 2) * 2 - 1)); // x, y coordinates, and the Burger's vector
        }

        if (!unsorted) // sorting if not told otherwise; first the negative, in increasing y value than positive with increasing y value
            std::sort(dislocs.begin(), dislocs.end(), [](const disl& a, const disl& b) {return (std::get<1>(a) + std::get<2>(a)) > (std::get<1>(b) + std::get<2>(b)); });

        for_each(dislocs.begin(), dislocs.end(), [&ofile](const disl& a) {ofile << std::get<0>(a) << "\t" << std::get<1>(a) << "\t" << std::get<2>(a) << "\n"; }); // print out to ofile
    }
#pragma endregion

    std::cout << "Done.\n";
    return 0;
}

double interpolate_y(pair a, pair b, double y) // returns x for y between a and b
{
    return a.first + (y - a.second) * (b.first - a.first) / (b.second - a.second);
}