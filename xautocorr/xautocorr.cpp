//
// xautocorr.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include "xautocorr_utils.h"

int main(int argc, char** argv)
{
#pragma region read in variables
    bpo::options_description requiredOptions("Required options"); // must be set from command line or file
    bpo::options_description optionalOptions("Optional options"); // must be set from command line or file

    requiredOptions.add_options()
        ("input-path,i", bpo::value<std::string>(), "The file name of the dislocation configuration file ending with .dconf. If it ends with .ini, a file containing the file names is expected: 1 filename per line.");

    optionalOptions.add_options()
        ("resolution,r", bpo::value<int>()->default_value(1024), "The dislocation density map will be evaluated in this many points and the number of the channels in autocorrelation will be the same. Please use sizes which conforms the suggestions of FFTW. (E.G. power of 2.)")
        ("method,m", bpo::value<std::string>()->default_value("bc"), "The method to calculate the autocorrelation. wspn: Wigner-Seitz positive and negative, wsts: Wigner-Seitz total and signed, bc: box-counting, gs: Gauss-smoothing")
        ("sub-sampling,s", bpo::value<int>()->default_value(16), "This is a parameter for method ws and gc. It tells how many times should be the mesh denser, on which the density will be evaluated. (E.G.power of 2.)")
        ("half-width,w", bpo::value<double>()->default_value(0.125), "This is a parameter for method gc. It tells how wide the Gauss-distribution should with which the dirac-delta densities should be concolved.")
        ("create-maps", "If set, the program will create the 2D density maps for rho_t and kappa.")
        ("output-fnameprefix,o", bpo::value<std::string>()->default_value("xautocorr/"), "With what output filename prefix should the initial conditions be stored. Symbol ./ means here.")
        ("debug-level,d", bpo::value<int>()->default_value(0), "Debugging information to show.")
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
        std::cout << "xautocorr from the 2D4 - a 2D discrete dislocation dynamics simulation program toolset.\n"
            "Copyright (C) Dániel Tüzes <tuzes@metal.elte.hu>\n";
    }
#pragma endregion

#pragma region process variables

    if (vm.count("help")) // if the user is curious 
    {
        std::cout << options << std::endl;
        exit(0);
    }

    // input-path
    if (0 == vm.count("input-path")) // to test if N is present
    {
        std::cerr << "input-path is missing! Program terminates.\n";
        exit(-1);
    }
    std::string ifname = vm["input-path"].as<std::string>();
    std::ifstream ifile(ifname);
    if (!ifile)
    {
        std::cerr << "The program couldn't open " << ifname << " for reading. Program terminates." << std::endl;
        return 0;
    }

    // output-path
    std::string of(vm["output-fnameprefix"].as<std::string>()); // the output foldername
    std::cout << "Files will be created in " << of << std::endl;

    // methods
    size_t subs = vm["sub-sampling"].as<int>(); // subs-sampling rate
    size_t res = vm["resolution"].as<int>(); // the resolution of the autocorrelation map
    size_t samp = res * subs; // sampling rate for method ws
    method um = na; // the used method for the calculation
    double sigma = vm["half-width"].as<double>();
    std::string methodname;
    if (vm["method"].as<std::string>() == "wspn" || vm["method"].as<std::string>() == "wsts")
    {
        std::cout << "Sampling rate for ";
        if (vm["method"].as<std::string>() == "wspn")
        {
            um = wspn;
            methodname = "Wigner-Seitz positive and negative";
        }
        else
        {
            um = wsts;
            methodname = "Wigner-Seitz total and signed";
        }
        std::cout << methodname << ": " << samp << "\n";
    }
    else if (vm["method"].as<std::string>() == "bc")
    {
        um = bc;
        methodname = "box counting";
    }
    else if (vm["method"].as<std::string>() == "gs")
    {
        std::cout << "Half-width of the Gauss-distribution: " << sigma << "\n";
        if (sigma < (double)1 / samp)
            std::cerr << "Warning: half-width is smaller than the sampling distance. Unwanted side-effect may arise. Consider using larger sampling or larger half-width. Program continues." << std::endl;
        um = gs;
        methodname = "Gauss-smoothing";
    }

    // extra infos
    bool create_maps = false;
    if (vm.count("create-maps"))
    {
        create_maps = true;
        std::cout << "Density maps will be created." << std::endl;
    }

    int debug_level = vm["debug-level"].as<int>();
    std::cout << "Debug level is " << debug_level << ". ";
    if (debug_level == 0)
        std::cout << "No debug information will be shown." << std::endl;
    else
        std::cout << "Debug information will be shown. The name of these files start with deb." << std::endl;

#pragma endregion

    std::vector<double> linedensity(res, 0); // a sûrûségre, késõbb majd ez fogja tárolni a korrelációt is

    std::vector<disl> dislocs; // container of the N number of dislocations
    for (disl tmp; ifile >> tmp; dislocs.push_back(tmp)); // read in dislocations
    std::sort(dislocs.begin(), dislocs.end(),
        [](const disl& a, const disl& b) {return (std::get<1>(a) + std::get<2>(a)) > (std::get<1>(b) + std::get<2>(b)); }); // dislocaions are sorted based on their type and y value
    size_t size = dislocs.size();

#pragma region create maps if asked for
    std::vector<std::vector<std::vector<double>>> maps; // the 2D maps for rho_t, kappa, rho_p, rho_n
    std::vector<std::string> names{ "rho_t", "kappa", "rho_p", "rho_n" };

    std::string ofname_extra = "_" + vm["method"].as<std::string>() + "_r" + std::to_string(res);

    std::vector<std::string> o_maps_fn;
    for (const auto& name : names)
        o_maps_fn.push_back(of + ifname + ofname_extra + "_" + name + ".txt");

    std::vector<std::ofstream> o_maps(4); // output files for the maps;  rho_t, kappa, rho_p, rho_n

    if (create_maps)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            maps.push_back(std::vector<std::vector<double>>(res, std::vector<double>(res, 0)));
            o_maps[i].open(o_maps_fn[i]);
            if (!o_maps[i])
            {
                std::cerr << "Cannot create " << o_maps_fn[i] << ". Program terminates." << std::endl;
                exit(-1);
            }
            o_maps[i] << "# This file contains the density of " << names[i] << " for the file " << ifname << " using the " << methodname << " method.\n";
        }
    }

#pragma endregion

    if (um == bc) // box counting case
    {
        for (size_t i = 0; i < dislocs.size(); ++i)
        {
            int xbin = (int)((std::get<0>(dislocs[i]) + 0.5) * res);
            int ybin = (int)((std::get<1>(dislocs[i]) + 0.5) * res);
            size_t cellnum = res * res;

            maps[0][ybin][xbin] += cellnum;
            if (i < size / 2) // positive dislocation
            {
                maps[2][ybin][xbin] += cellnum;
                maps[1][ybin][xbin] += cellnum;
            }
            else // negative dislocation
            {
                maps[3][ybin][xbin] += cellnum;
                maps[1][ybin][xbin] -= cellnum;
            }
        }
    }

    if (um == gs) // Gauss-smoothing case
    {
        for (size_t y = 0; y < samp; ++y) // samp =  res * subs, it defines a more fine mesh
            for (size_t x = 0; x < samp; ++x)
            {
                pair centerpoint((x + 0.5) / samp - 0.5, (y + 0.5) / samp - 0.5); // centerpoints of the fine mesh
                for (size_t i = 0; i < dislocs.size(); ++i)
                {
                    pair disl_pos_xy(std::get<0>(dislocs[i]), std::get<1>(dislocs[i]));
                    double exponent = exp(-distsq(centerpoint, disl_pos_xy) / (2 * sigma * sigma));

                    if (i < size / 2) // positive dislocation
                        maps[2][y / subs][x / subs] += exponent;
                    else // negative dislocation
                        maps[3][y / subs][x / subs] += exponent;
                }
            }

        normalize(maps[2], res * res * (double)size / 2);
        normalize(maps[3], res * res * (double)size / 2);

        for (size_t y = 0; y < res; ++y) // rho_t and kappa cannot be calculated without the normalizaton of rho_p and rho_n
            for (size_t x = 0; x < res; ++x) // using stl::transform twice on rho_t and kappa would have an unnecessary overhead
            {
                maps[0][y][x] = maps[2][y][x] + maps[3][y][x];
                maps[1][y][x] = maps[2][y][x] - maps[3][y][x];
            }
    }

    if (um == wspn) // Wigner-Seitz positive and negative case
    {
        for (auto& val : dislocs) // zero out the 3rd element, it will be used to measure the area
            std::get<2>(val) = 0;

        measure_area(dislocs, 0, size/2, samp);
        measure_area(dislocs, size/2, size, samp);

        if (create_maps)
        {
            measure_density<true>(dislocs, samp, maps[2], maps[3]);
            
            for (size_t i = 0; i < res; ++i)
                for (size_t j = 0; j < res; ++j)
                {
                    maps[0][i][j] = maps[2][i][j] + maps[3][i][j];
                    maps[1][i][j] = maps[2][i][j] - maps[3][i][j];
                }
        }
    }

    if (um == wsts)
    {
        for (auto& val : dislocs) // zero out the 3rd element, it will be used to measure the area
            std::get<2>(val) = 0;

        measure_area(dislocs, samp);

        if (create_maps)
        {
            measure_density<false>(dislocs, samp, maps[0], maps[1]);
            for (size_t i = 0; i < res; ++i)
                for (size_t j = 0; j < res; ++j)
                {
                    maps[2][i][j] = (maps[0][i][j] + maps[1][i][j]) / 2;
                    maps[3][i][j] = (maps[0][i][j] - maps[1][i][j]) / 2;
                }
        }
    }

    if (debug_level) // creates dislocation area file deb_disl.txt and levels for gnuplot contour
    {
        std::ofstream debo_disl_areafile(of + "deb_disl.txt");
        debo_disl_areafile << "# The positive and negative dislocations with their corresponding area.\n"
            << "# The first half of the dislocations are positive, the rest are negative.\n"
            << "# The area can be calculated separately for different types (in case wspn) or their joined set (in case wsts).\n";

        for (const auto& d : dislocs)
            debo_disl_areafile << d << "\n"; // print out to ofile

        for (size_t i = 0; i < 4; ++i)
            gnuplotlevels(maps[i], of + "deb_levels_" + names[i] + ".txt");
    }

    if (create_maps)
        for (size_t i = 0; i < 4; ++i)
            o_maps[i] << maps[i];

    std::cout << "Done. Program terminates." << std::endl;
    return 0;
}