//
// init_config_gen.cpp : This file contains the 'main' function. Program execution begins and ends there.

#define VERSION_init_config_gen 2.5
/*changelog
# 2.5
* issue#10: the real pattern had 2*n waves instead of n
* wall-height doesn't have default value anymore

# 2.4
* if shortest-criteria doesn't hold and new y value is assigned, the new value does not depend on the last largest seed number anymore
* wall-height is introduced and implemented

# 2.3
parameter bare and gcc warning are eliminated

# 2.2
* random number engine supports 64 bit randomness eliminating too many too close dislocations
* too close dislocations can be eliminated

# 2.1
Static code analysis suggested to change n to int to avoid size_t (8byte) > int (4byte) concersion

# 2.0
* sorting input information is taken now by the input named sorted and not unsorted. It can be x, y or u. The latter means no sorting. Default is y, so the behaviour is not changed.

# 1.1
* potential bug fix
* version number is not string but number

# 1.0
* first version working with the initial idea
* refactored code
*/

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
#include <iomanip>

namespace bpo = boost::program_options;
namespace br = boost::random;

using pair = std::pair<double, double>;
double interpolate_y(pair, pair, double);

using disl = std::tuple<double, double, int>;       // a dislocation is a (double, double, int) tuple for (posx,posy,type)
size_t closest_yID(const disl& cmp, const std::vector<disl>& list)  // the dislocation that is the closest to the selected cmp one
{
    double smallest = INFINITY;
    size_t ID = list.size() + 1;
    for (size_t i = 0; i < list.size(); ++i)
    {
        double tmp = std::get<1>(list[i]) - std::get<1>(cmp);

        while (tmp >= 0.5)
            tmp -= 1;
        while (tmp < -0.5)
            tmp += 1;

        if (fabs(tmp) < smallest)
        {
            smallest = fabs(tmp);
            ID = i;
        }
    }
    return ID;
}

// transfers the coordinate into the range [-0.5; 0.5)
void normalize(double& n)
{
    while (n < -0.5) // bad predictions can lead to values like n=-10'000 or so and using hyperbolic functions on them is problematic
        n += 1;

    while (n >= 0.5) // in rare cases, if is not enough, and in most cases, this is faster than fprem1: https://stackoverflow.com/questions/58803438/best-way-calculating-remainder-on-floating-points
        n -= 1;
}
#pragma endregion

int main(int argc, char** argv)
{
#pragma region reading in variables

    bpo::options_description requiredOptions("Required options"); // must be set from command line or file
    bpo::options_description optionalOptions("Optional options"); // must be set from command line or file

    requiredOptions.add_options()
        ("N,N", bpo::value<size_t>(), "The number of dislocations to generate. Must be even, because N/2 number of positive and negative dislocations will be created.");

    optionalOptions.add_options()
        ("seed-start,S", bpo::value<int>()->default_value(1000), "An integer used as an initial seed value for the random number generator.")
        ("seed-end,E", bpo::value<int>(), "An integer used as the last seed value for the random number generator, seed-end > seed_start must hold. If set, seed-end - seed-start number of initial configurations will be created.")
        ("shortest-distance,l", bpo::value<double>()->default_value(0), "The allowed smallest distance in y. If the distance is smaller, either a new position is generated or the configuration will be marked. If used, suggested value is 8e-9.")
        ("mark,m", bpo::value<std::string>(), "If set, configuration where shortest-distance criterion doesn't meet will be marked with the given argument in the file name. Otherwise, the renitent newer dislocation's y value will be regenerated.")
        ("pattern-strength,A", bpo::value<double>()->default_value(0), "Instead of a uniform distribution, the probability density function will be\n1 + A * sin(x * n * 2 pi) with A = pattern-strength for rho_+ and\n1 - A * sin(x * w * 2 pi) for rho_-. A must be in [-1:1].")
        ("linear-wavenumber,n", bpo::value<int>()->default_value(3), "The number of waves in the simulation area [-0.5:0.5).")
        ("pattern-type,T", bpo::value<std::string>()->default_value("s"), "The pattern type in the density distribution.")
        ("wall-height,w", bpo::value<size_t>(), "The number of dislocations forming 1 wall. N/2/w must be integer. shortest-distance criteria is checked for the bottom dislocation.")
        ("sorted,s", bpo::value<std::string>()->default_value("y"), "Defines the order of the dislocations in the output file.\n   * x: decreasing Burger's vector, and  decreasing x coordinate\n   * y:  decreasing Burger's vector, and  decreasing y coordinate\n   * u: sign is alternating and x and y coordinates are uncorrelated.")
        ("output-fnameprefix,o", bpo::value<std::string>()->default_value(""), "In which folder should the initial conditions be stored. Symbol ./ means here.")
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
            << "Copyright (C) Dániel Tüzes <tuzes@metal.elte.hu>\n";
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
        N = vm["N"].as<size_t>();
        std::cout << "N =\t" << N << std::endl;
    }

    if (N % 2 != 0)  // to test if N is even
    {
        std::cerr << "N must be even. Program terminates." << std::endl;
        exit(-1);
    }
#pragma endregion

#pragma region seed-end, seed-start

    int seed_start = vm["seed-start"].as<int>(); // the first seed value for the random number generator engine
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

    bool wall = vm.count("wall-height");
    size_t wall_height;                             // the number of dislocations positioned on top of each other
    if (wall)
    {
        wall_height = vm["wall-height"].as<size_t>();
        if (N % (wall_height * 2) != 0)
        {
            std::cerr << "N/2/wall-height must be integer, but it is not. Program terminates." << std::endl;
            exit(-1);
        }
        
    }

#pragma endregion

#pragma region output path, filenames, sorting and restriction

    std::string of(vm["output-fnameprefix"].as<std::string>()); // the output filename prefix (potentially inculding foldername)
    std::cout << "output-fnameprefix =\t" << of << std::endl;

    char sorted; // leave the output unsorted?
    {
        std::string sortedstr = vm["sorted"].as<std::string>();
        if (sortedstr.compare("x") == 0)
            sorted = 'x';
        else if (sortedstr.compare("y") == 0)
            sorted = 'y';
        else if (sortedstr.compare("u") == 0)
            sorted = 'u';
        else
        {
            std::cerr << "Error: unsupported ordering: " << sortedstr << ". Program terminates." << std::endl;
            exit(-1);
        }
    }

    std::cout << "sorted =\t" << sorted << std::endl;

    bool mark = vm.count("mark");   // to mark the output file or to regenerate the problematic dislocation

    double shortestDist = vm["shortest-distance"].as<double>();// The allowed smallest distance in y. If the distance is smaller, either a new position is generated or the configuration will be marked.

#pragma endregion

#pragma endregion

#pragma region cumulative distribution

    std::vector<std::pair<double, double>> cdfr_p, cdfr_n;  // x,y values of the cumulative distribution function of rho_p and rho_n
    if (A != 0)
    {
        int res = (int)sqrt(N) * 10 + 1000;                     // resolution for the inverse function lookup; heuristic guess
        double n = vm["linear-wavenumber"].as<int>();           // linear-wavenumber
        for (int i = 0; i <= res; ++i)
        {
            double x = (double)i / res - 0.5;                   // the x position in dimensionless units, x∈[-0.5:0.5]
            double yrp = 0.5 + x - A * pow(cos(x * n * M_PI), 2) / (n * M_PI); // the y values for the rho_p cdf
            double yrn;                                         // the y values for the rho_n cdf
            if (vm["pattern-type"].as<std::string>() == "s")    // signed, i.e. κ has patterns
                   yrn = 0.5 + x + A * pow(cos(x * n * M_PI), 2) / (n * M_PI);
            else                                                // total, i.e. ρ_t has patterns
                yrn = yrp;
            cdfr_p.push_back(std::pair<double, double>(x, yrp));
            cdfr_n.push_back(std::pair<double, double>(x, yrn));
        }
    }

#pragma endregion

#pragma region generate and write out configuration
    for (int seed_val = seed_start; seed_val < seed_end; ++seed_val) // generate configurations with seeds in the range of [seed-start; seed-end]; seed-end has been set at label: seed-end set
    {
        std::cout << "Generating configuration with seed value " << seed_val << "." << std::endl;

        br::mt19937_64 engine(seed_val);                    // Mersenne twister is a good and expensive random number generator
        br::uniform_real_distribution<> distr(-0.5, 0.5);   // uniform distribution on [-0.5; 0.5)
        std::vector<disl> dislocs;                          // container of the N number of dislocations
        bool marked = false;                                // tells if there is (was) a dislocation with smaller y distance than shortest-distance

        //std::ofstream deb_conf_ofile("debug_conf.txt");
        for (size_t n = 0; n < N; ++n) // generate the N number of dislocations, n = dislocs.size()
        {
            double x = distr(engine);
            if (A != 0) // if no periodic pattern is requested
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
            double y = distr(engine);       // the suggested y value for the dislocation
            int burgersv = (n % 2) * 2 - 1; // the Burgers vector of the dislocation
            if (wall)
                burgersv = ((n / wall_height) % 2) * 2 - 1; // the Burgers vector of the dislocation

            // generate a proper y value
            for (int seed2 = seed_val + 1; shortestDist != 0 && n > 0; seed2++)
            {
                if (seed2 != seed_val + 1)  // if this is not the first time here
                {
                    br::mt19937_64 engine2(seed2);
                    y = distr(engine2);     // y value must be regenerated
                }
                disl tmp(x, y, burgersv);    // the nominated new dislocation
                size_t closestYID = closest_yID(tmp, dislocs);
                disl closestY = dislocs[closestYID];

                double ydist = std::get<1>(closestY) - y;
                normalize(ydist);
                if (fabs(ydist) < shortestDist)
                {
                    std::cout << std::setprecision(25)
                        << "Dislocations were closer than accepted for ID = " << closestYID << " and " << n << ". The positions and types are:\n"
                        << "\t(" << std::get<0>(closestY) << "), (" << std::get<1>(closestY) << ")\t" << std::get<2>(closestY) << "\n"
                        << "\t(" << x << "), (" << y << ")" << "\t" << burgersv << "\n"
                        << "dist = " << ydist << "\n"
                        << "Limit was " << shortestDist << "\n";
                    marked = true;
                    if (!mark)
                    {
                        std::cout << "This y value will be regenerated from a temporary random engine with seed value larger by 1.\n";
                        continue;
                    }
                    else
                        std::cout << "This simulation will be marked with a " << vm["mark"].as<std::string>() << " in the output filename.\n";
                }
                break;
            }
            dislocs.push_back(disl(x, y, burgersv)); // add the dislocation to the list
            for (size_t wc = 1; wall && wc < wall_height; ++wc)
            {
                double ynew = y + wc / sqrt(N); // the y coordinate of the next dislocation above the previous one
                normalize(ynew);
                dislocs.push_back(disl(x, ynew, burgersv)); // add the dislocation to the list
                ++n;
            }
        }

        if (sorted == 'y') // first the negative, in increasing y value than positive with increasing y value
            std::sort(dislocs.begin(), dislocs.end(), [](const disl& a, const disl& b) {return (std::get<1>(a) + std::get<2>(a)) > (std::get<1>(b) + std::get<2>(b)); });
        if (sorted == 'x') // first the negative, in increasing x value than positive with increasing x value
            std::sort(dislocs.begin(), dislocs.end(), [](const disl& a, const disl& b) {return (std::get<0>(a) + std::get<2>(a)) > (std::get<0>(b) + std::get<2>(b)); });

        std::string ofname = of + std::to_string(seed_val);
        if (mark && marked)
            ofname += vm["mark"].as<std::string>();
        ofname += ".dconf"; // output filename; the file is inside a folder

        std::ofstream ofile(ofname); // the filestream
        if (!ofile) // evaluates to false if file cannot be opened
        {
            std::cerr << "Cannot write " << ofname << ". Program terminates." << std::endl;
            exit(-1);
        }
        ofile.precision(std::numeric_limits<double>::max_digits10); // print out every digits

        for_each(dislocs.begin(), dislocs.end(), [&ofile](const disl& a) {ofile << std::get<0>(a) << "\t" << std::get<1>(a) << "\t" << std::get<2>(a) << "\n"; }); // print out to ofile
    }
#pragma endregion

    std::cout << "Done.\n";
    return 0;
}

// returns x for y between a and b
double interpolate_y(pair a, pair b, double y)
{
    return a.first + (y - a.second) * (b.first - a.first) / (b.second - a.second);
}