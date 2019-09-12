//
// xautocorr.cpp : This file contains the 'main' function. Program execution begins and ends there.

#define VERSION_xautocorr 1.0

#include "xautocorr_utils.h"

int main(int argc, char** argv)
{
    // test_dirac(23);
    // return 0;

#pragma region read in variables
    bpo::options_description requiredOptions("Required options"); // must be set from command line or file
    bpo::options_description optionalOptions("Optional options"); // must be set from command line or file

    requiredOptions.add_options()
        ("input-fname,i", bpo::value<std::string>(), "The file name of the dislocation configuration file ending with .dconf. If it ends with .ini, a file containing the file names is expected: 1 filename per line.");

    optionalOptions.add_options()
        ("resolution,r", bpo::value<int>()->default_value(1024), "The dislocation density map will be evaluated in this many points and the number of the Fourier components will be the same (for all methods except df, where it only means the number of Fourier components). Please use sizes which conforms the suggestions of FFTW. (E.G. power of 2.)")
        ("method,m", bpo::value<std::string>()->default_value("bc"), "The method to investigate the pattern.\n - wspn: Wigner-Seitz positive and negative\n - wsts: Wigner-Seitz total and signed\n - bc: box-counting\n - gs: Gauss-smoothing\n - df: direct Fourier")
        ("sub-sampling,s", bpo::value<int>()->default_value(1), "This is a parameter for method gs, wspn and wsts. It tells how many times should be the mesh denser, on which the density will be evaluated. (E.G.power of 2.)")
        ("half-width,w", bpo::value<double>(), "This is a parameter for method gc. It gives the width of the Gauss-distribution with which the Dirac-delta densities are convolved.")
        ("create-maps", "If set, the program will create the 2D density maps for rho_t and kappa.")
        ("output-fnameprefix,o", bpo::value<std::string>()->default_value(""), "With what output filename prefix should the initial conditions be stored. Symbol ./ or empty string means here.")
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
        std::cout << "xautocorr (version " << std::setprecision(1) << std::fixed << (VERSION_xautocorr + VERSION_xautocorr_utils) << ") from 2D4 - a 2D discrete dislocation dynamics simulation program toolset.\n"
            "Copyright (C) Dániel Tüzes <tuzes@metal.elte.hu>\n";
    }
#pragma endregion

#pragma region process variables

    if (vm.count("help")) // if the user is curious 
    {
        std::cout << options << std::endl;
        exit(0);
    }

    // input-fname
    if (0 == vm.count("input-fname")) // to test if N is present
    {
        std::cerr << "input-fname is missing! Program terminates.\n";
        exit(-1);
    }
    std::string ifname = vm["input-fname"].as<std::string>();
    std::vector<std::string> ifnames;

    if (ifname.size() >= 4 && ifname.compare(ifname.size() - 4, 4, ".ini") == 0)
    {
        std::ifstream ifile(ifname);
        if (!ifile)
        {
            std::cerr << "The program couldn't open " << ifname << " for reading. Program terminates." << std::endl;
            return 0;
        }

        for (std::string tmp; ifile >> tmp;)
        {
            if (tmp.size() >= 6 && tmp.compare(tmp.size() - 6, 6, ".dconf") == 0)
                ifnames.push_back(tmp);
            else
            {
                std::cerr << ifname << " contains invalid dislocation configuration filename " << tmp << ". Program termiantes." << std::endl;
                exit(2);
            }
        }
        std::cout << ifname << " as input path contained " << ifnames.size() << " number of files to read." << std::endl;
        if (ifnames.size() == 0)
        {
            std::cerr << "Not enough input files. Program terminates." << std::endl;
            exit(1);
        }
    }
    else if (ifname.size() >= 6 && ifname.compare(ifname.size() - 6, 6, ".dconf") == 0)
    {
        std::cout << ifname << " will be used to read dislocation configuration" << std::endl;
        ifnames.push_back(ifname); // ifnames contains exactly 1 file that can be opened or
    }
    else
    {
        std::cerr << ifname << " is not a valid input filename. Program terminates." << std::endl;
        exit(1);
    }

    // output-path
    std::string of(vm["output-fnameprefix"].as<std::string>()); // the output filename prefix (potentially inculding foldername)
    std::cout << "output-fnameprefix = \t" << of << std::endl;

    // methods
    size_t subs = vm["sub-sampling"].as<int>(); // subs-sampling rate
    size_t res = vm["resolution"].as<int>(); // the resolution of the autocorrelation map
    size_t samp = res * subs; // sampling rate for method ws
    method um = na; // the used method for the calculation
    double sigma;
    size_t half_boxwidth;

    std::string methodname;
    std::string methodnameabbrev = vm["method"].as<std::string>();
    if (methodnameabbrev == "wspn" || methodnameabbrev == "wsts")
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
    else if (methodnameabbrev == "bc")
    {
        um = bc;
        methodname = "box counting";
    }
    else if (methodnameabbrev == "gs")
    {
        if (vm.count("half-width") == 0)
        {
            std::cerr << "Warning: no half-width were set by the user, it will be set automatically so that half-width * sub-sampling * resolution = 4 will hold." << std::endl;
            sigma = (double)4 / samp;
        }
        else
        {
            sigma = vm["half-width"].as<double>();
            if (sigma < (double)1 / samp)
                std::cerr << "Warning: half-width is smaller than the sampling distance. Unwanted side-effect may arise. Consider using larger sampling or larger half-width. Program continues." << std::endl;
        }
        half_boxwidth = std::min(size_t(6 * sigma * samp) + 1, samp / 2);

        std::cout << "Half-width of the Gauss-distribution: " << sigma << ", a dislocation affects 2*" << half_boxwidth << " * 2*" << half_boxwidth << " number of subcells.\n";
        um = gs;
        methodname = "Gauss-smoothing";
    }
    else if (methodnameabbrev == "df")
    {
        um = df;
        methodname = "direct Fourier";
    }

    // extra infos
    bool create_maps = false;
    if (vm.count("create-maps"))
    {
        create_maps = true;
        std::cout << "Maps will be created for direct or k space density." << std::endl;
    }

    int debug_level = vm["debug-level"].as<int>();
    std::cout << "Debug level is " << debug_level << ". ";
    if (debug_level == 0)
        std::cout << "No debug information will be shown." << std::endl;
    else
        std::cout << "Debug information will be shown. The name of these files start with deb." << std::endl;

#pragma endregion

#pragma region define vectors F_absValSqs
    std::vector<std::vector<double>> F_absValSqs; // the abs square of the different Fourier-transformed data; rho_t, kappa, rho_p, rho_n
    std::vector<std::string> names{ "rho_t", "kappa", "rho_p", "rho_n" };  // rho_t, kappa, rho_p, rho_n
    std::string ofname_extra = methodnameabbrev + "_r" + std::to_string(res); // method name and resolution

    std::vector<std::string> o_k_fn; // output filenames for the Fourier components of rho_t, kappa, rho_p and rho_n
    for (const auto& name : names)
        o_k_fn.push_back(of + ofname_extra + "_" + name + "_k.txt");

    std::vector<std::ofstream> o_kf(4); // output files for the k values; rho_t, kappa, rho_p, rho_n

    for (size_t i = 0; i < 4; ++i)
    {
        F_absValSqs.push_back(std::vector<double>(res / 2 + 1));
        o_kf[i].open(o_k_fn[i]);
        if (!o_kf[i])
        {
            std::cerr << "Cannot create " << o_k_fn[i] << ". Program terminates." << std::endl;
            exit(-1);
        }
        o_kf[i] << "# This file contains the Fourier components of " << names[i] << " for the file " << ifname << " using the " << methodname << " method.\n";
    }

#pragma endregion

    for (const auto& ifname : ifnames)
    {
#pragma region create maps
        std::vector<std::vector<std::vector<double>>> maps; // the 2D maps for rho_t, kappa, rho_p, rho_n
        std::vector<std::vector<std::vector<double>>> imaps; // the imaginary 2D maps for rho_t, kappa, rho_p, rho_n for the case df

        std::vector<std::string> o_maps_fn;
        for (const auto& name : names)
            o_maps_fn.push_back(of + ifname + ofname_extra + "_" + name + ".txt");

        std::vector<std::ofstream> o_maps(4); // output files for the maps; rho_t, kappa, rho_p, rho_n

        if (create_maps || um == bc || um == gs || um == df) // bc and gs relies on the maps, but they won't be printed; df needs whole k to be calculated
        {
            for (size_t i = 0; i < 4; ++i)
            {
                if (um != df)
                    maps.push_back(std::vector<std::vector<double>>(res, std::vector<double>(res)));
                else
                    imaps.push_back(std::vector<std::vector<double>>(res, std::vector<double>(res)));
            }
        }

#pragma endregion

#pragma region reading in dislocations
        std::vector<disl> dislocs; // container of the N number of dislocations
        std::ifstream ifile(ifname);
        if (!ifile)
        {
            std::cerr << "Cannot open " << ifname << "for reading. Program terminates." << std::endl;
            exit(1);
        }
        std::cout << "Reading in " << ifname << std::endl;
        for (disl tmp; ifile >> tmp; dislocs.push_back(tmp)); // read in dislocations
        std::sort(dislocs.begin(), dislocs.end(),
            [](const disl& a, const disl& b) {return (std::get<1>(a) + std::get<2>(a)) > (std::get<1>(b) + std::get<2>(b)); }); // dislocaions are sorted based on their type and y value
        size_t size = dislocs.size();
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

            for (size_t i = 0; i < 4; ++i)
                for (size_t linenum = 0; linenum < res; ++linenum)
                    addFourierAbsValSq1D(F_absValSqs[i], maps[i][linenum]);
        }

        if (um == gs) // Gauss-smoothing case
        {
            for (size_t i = 0; i < dislocs.size(); ++i)
            {
                double xcoord = std::get<0>(dislocs[i]);
                double ycoord = std::get<1>(dislocs[i]);
                size_t xid_min = ((size_t)((xcoord + 0.5) * samp - half_boxwidth) + samp) % samp;
                size_t yid_min = ((size_t)((ycoord + 0.5) * samp - half_boxwidth) + samp) % samp;
                pair disl_pos_xy(xcoord, ycoord);
                for (size_t y = yid_min; y < yid_min + 2 * half_boxwidth; ++y)
                {
                    size_t yid = y % samp;
                    for (size_t x = xid_min; x < xid_min + 2 * half_boxwidth; ++x)
                    {
                        size_t xid = x % samp;
                        pair centerpoint((xid + 0.5) / samp - 0.5, (yid + 0.5) / samp - 0.5); // centerpoints of the fine mesh
                        double exponent = exp(-distsq(centerpoint, disl_pos_xy) / (2 * sigma * sigma));

                        if (i < size / 2) // positive dislocation
                            maps[2][yid / subs][xid / subs] += exponent;
                        else // negative dislocation
                            maps[3][yid / subs][xid / subs] += exponent;
                    }
                }
            }

            normalize(maps[2], res * res * (double)size / 2);
            normalize(maps[3], res * res * (double)size / 2);

            for (size_t y = 0; y < res; ++y) // rho_t and kappa cannot be calculated without the normalizaton of rho_p and rho_n
            {
                for (size_t x = 0; x < res; ++x) // using stl::transform twice on rho_t and kappa would have an unnecessary overhead
                {
                    maps[0][y][x] = maps[2][y][x] + maps[3][y][x];
                    maps[1][y][x] = maps[2][y][x] - maps[3][y][x];
                }
                for (size_t i = 0; i < 4; ++i)
                    addFourierAbsValSq1D(F_absValSqs[i], maps[i][y]);
            }
        }

        if (um == wspn) // Wigner-Seitz positive and negative case
        {
            for (auto& val : dislocs) // zero out the 3rd element, it will be used to measure the area
                std::get<2>(val) = 0;

            measure_area(dislocs, 0, size / 2, samp); // measure the area of the positive dislocations
            measure_area(dislocs, size / 2, size, samp); // measure the area of the negative dislocations

            for (size_t y = 0; y < res; ++y)
            {
                std::vector<std::vector<double>> linedensities(4, std::vector<double>(res));
                measure_density<true>(dislocs, y, samp, linedensities[2], linedensities[3]);
                for (size_t x = 0; x < res; ++x)
                {
                    linedensities[0][x] = linedensities[2][x] + linedensities[3][x];
                    linedensities[1][x] = linedensities[2][x] - linedensities[3][x];
                }
                for (size_t i = 0; i < 4; ++i)
                    addFourierAbsValSq1D(F_absValSqs[i], linedensities[i]);

                if (create_maps) // copy linedensity to the map
                {
                    for (size_t i = 0; i < 4; ++i)
                        for (size_t x = 0; x < res; ++x)
                            maps[i][y][x] = linedensities[i][x];
                }
            }
        }

        if (um == wsts)
        {
            for (auto& val : dislocs) // zero out the 3rd element, it will be used to measure the area
                std::get<2>(val) = 0;

            measure_area(dislocs, samp);

            for (size_t y = 0; y < res; ++y)
            {
                std::vector<std::vector<double>> linedensities(4, std::vector<double>(res));
                measure_density<false>(dislocs, y, samp, linedensities[0], linedensities[1]);
                for (size_t x = 0; x < res; ++x)
                {
                    linedensities[2][x] = (linedensities[0][x] + linedensities[1][x]) / 2;
                    linedensities[3][x] = (linedensities[0][x] - linedensities[1][x]) / 2;
                }
                for (size_t i = 0; i < 4; ++i)
                    addFourierAbsValSq1D(F_absValSqs[i], linedensities[i]);

                if (create_maps) // copy linedensity to the map
                {
                    for (size_t i = 0; i < 4; ++i)
                        for (size_t x = 0; x < res; ++x)
                            maps[i][y][x] = linedensities[i][x];
                }
            }
        }

        if (um == df)
        {
            for (int i = 0; i < static_cast<int>(dislocs.size()); ++i)
                for (int kx = 0; kx < res; ++kx)
                    for (int ky = 0; ky < res; ++ky)
                    {
                        maps[0][kx][ky] += cos(2 * (M_PI * ((std::get<0>(dislocs[i]) * kx) + (std::get<1>(dislocs[i]) * ky))) / res);
                        imaps[0][kx][ky] -= sin(2 * (M_PI * ((std::get<0>(dislocs[i]) * kx) + (std::get<1>(dislocs[i]) * ky))) / res);


                        maps[1][kx][ky] += std::get<2>(dislocs[i]) * cos(2 * (M_PI * ((std::get<0>(dislocs[i]) * kx) + (std::get<1>(dislocs[i]) * ky))) / res);
                        imaps[1][kx][ky] -= std::get<2>(dislocs[i]) * sin(2 * (M_PI * ((std::get<0>(dislocs[i]) * kx) + (std::get<1>(dislocs[i]) * ky))) / res);

                        if (std::get<2>(dislocs[i]) == 1) // positive dislocation; if i<dislocs.size()/2, it is positive, and negative otherwise, but no problem, this is still fast
                        {
                            maps[2][kx][ky] += cos(2 * (M_PI * ((std::get<0>(dislocs[i]) * kx) + (std::get<1>(dislocs[i]) * ky))) / res);
                            imaps[2][kx][ky] -= sin(2 * (M_PI * ((std::get<0>(dislocs[i]) * kx) + (std::get<1>(dislocs[i]) * ky))) / res);
                        }
                        else  // negative dislocation
                        {
                            maps[3][kx][ky] += cos(2 * (M_PI * ((std::get<0>(dislocs[i]) * kx) + (std::get<1>(dislocs[i]) * ky))) / res);
                            imaps[3][kx][ky] -= sin(2 * (M_PI * ((std::get<0>(dislocs[i]) * kx) + (std::get<1>(dislocs[i]) * ky))) / res);
                        }
                    }
        }

        if (debug_level) // creates dislocation area file deb_disl.txt and levels for gnuplot contour
        {
            std::ofstream debo_disl_areafile(of + ifname + ofname_extra + "deb_disl.txt");
            debo_disl_areafile << "# The positive and negative dislocations with their corresponding area.\n"
                << "# The first half of the dislocations are positive, the rest are negative.\n"
                << "# The area can be calculated separately for different types (in case wspn) or their joined set (in case wsts).\n";

            if (create_maps && (um == wsts || um == wspn))
            {
                for (const auto& d : dislocs)
                    debo_disl_areafile << d << "\n"; // print out to ofile
                for (size_t i = 0; i < 4; ++i)
                    gnuplotlevels(maps[i], of + ifname + ofname_extra + "_" + names[i] + "levels.txt");
            }
        }

        if (create_maps)
            for (size_t i = 0; i < 4; ++i)
            {
                o_maps[i].open(o_maps_fn[i]);
                if (!o_maps[i])
                {
                    std::cerr << "Cannot create " << o_maps_fn[i] << ". Program terminates." << std::endl;
                    exit(-1);
                }
                if (um != df)
                    o_maps[i] << "# This file contains the density of " << names[i] << " for the file " << ifname << " using the " << methodname << " method.\n";
                else
                    o_maps[i] << "# This file contains the Fourier components (real, tab, imag) of " << names[i] << " for the file " << ifname << " using the " << methodname << " method.\n";

                o_maps[i] << maps[i];
            }

    }

    for (size_t i = 0; i < 4; ++i)
        for (const auto& val : F_absValSqs[i])
            o_kf[i] << val << "\n";

    std::cout << "Done. Program terminates." << std::endl;
    return 0;
}