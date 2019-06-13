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
        ("method,m", bpo::value<std::string>()->default_value("bc"), "The method to calculate the autocorrelation. ws: Wigner-Seitz, bc: box-counting, gc: Gauss-convolving")
        ("sub-sampling,s", bpo::value<int>()->default_value(16), "This is a parameter for method ws and gc. It tells how many times should be the mesh denser, on which the density will be evaluated. (E.G.power of 2.)")
        ("half-width,w", bpo::value<double>()->default_value(0.125), "This is a parameter for method gc. It tells how wide the Gauss-distribution should with which the dirac-delta densities should be concolved.")
        ("create-maps", "If set, the program will create the 2D density maps for rho_t and kappa.")
        ("output-foldername,o", bpo::value<std::string>()->default_value("xautocorr"), "In which folder should the initial conditions be stored. Symbol . means here.")
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
    else //
    {
        // input-path
        if (0 == vm.count("input-path")) // to test if N is present
        {
            std::cerr << "input-path is missing! Program terminates.\n";
            exit(-1);
        }
    }

    unsigned int subs = vm["sub-sampling"].as<int>(); // subs-sampling rate
    unsigned int res = vm["resolution"].as<int>(); // the resolution of the autocorrelation map
    unsigned int samp = res * subs; // sampling rate for method ws
    method um = na; // the used method for the calculation
    double sigma = vm["half-width"].as<double>();
    std::string methodname;
    if (vm["method"].as<std::string>() == "ws")
    {
        std::cout << "Sampling rate for Wigner-Seitz: " << samp << "\n";
        um = ws;
        methodname = "Wigner-Seitz";
    }
    else if (vm["method"].as<std::string>() == "bc")
    {
        um = bc;
        methodname = "box counting";
    }
    else if (vm["method"].as<std::string>() == "gc")
    {
        std::cout << "Half-width of the Gauss-distribution: " << sigma << "\n";
        if (sigma < (double)1 / samp)
            std::cerr << "Warning: half-width is smaller than the sampling distance. Unwanted side-effect may arise. Consider using larger sampling or larger half-width. Program continues." << std::endl;
        um = gc;
        methodname = "Gauss convolution";
    }

    bool create_maps = false;
    if (vm.count("create-maps"))
    {
        create_maps = true;
        std::cout << "Density maps will be created." << std::endl;
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
    for (disl tmp; ifile >> tmp; dislocs.push_back(tmp)); // read in dislocations

    std::string o_rho_t_fn = ifname + "_rho_t.txt";
    std::string o_kappa_fn = ifname + "_kappa.txt";
    std::ofstream o_rho_t;
    std::ofstream o_kappa;

    if (create_maps)
    {
        o_rho_t.open(o_rho_t_fn);
        o_kappa.open(o_kappa_fn);
        if (!o_rho_t)
        {
            std::cerr << "Cannot create " << o_rho_t_fn << ". Program terminates." << std::endl;
            exit(-1);
        }
        if (!o_kappa)
        {
            std::cerr << "Cannot create " << o_kappa_fn << ". Program terminates." << std::endl;
            exit(-1);
        }
        o_rho_t << "# This file contains the frequency of rho_t for the file " << ifname << " using the " << methodname << " method.\n";
        o_kappa << "# This file contains the frequency of kappa for the file " << ifname << " using the " << methodname << " method.\n";
    }

    if (um == bc) // box counting case
    {
        std::vector < std::vector<int>> rho_t(res, std::vector<int>(res, 0)); // the 2D map for rho_total
        std::vector < std::vector<int>> kappa(res, std::vector<int>(res, 0)); // the 2D map for kappa = rho_signed
        for (auto const& disl : dislocs)
        {
            int xbin = (int)((std::get<0>(disl) + 0.5) * res);
            int ybin = (int)((std::get<1>(disl) + 0.5) * res);
            int type = std::get<2>(disl);

            rho_t[ybin][xbin] += 1;
            kappa[ybin][xbin] += type;
        }

        if (create_maps)
        {
            o_rho_t << rho_t;
            o_kappa << kappa;
        }
    }

    if (um == gc) // Gauss convolving case
    {
        std::vector<std::vector<double>> rho_p(res, std::vector<double>(res, 0)); // the 2D map for rho_total
        std::vector<std::vector<double>> rho_n(res, std::vector<double>(res, 0)); // the 2D map for kappa = rho_signed

        for (unsigned int y = 0; y < samp; ++y) // samp =  res * subs, it defines a more fine mesh
            for (unsigned int x = 0; x < samp; ++x)
            {
                pair centerpoint((x + 0.5) / samp - 0.5, (y + 0.5) / samp - 0.5); // centerpoints of the fine mesh
                for (auto const& disl : dislocs)
                {
                    pair disl_pos_xy(std::get<0>(disl), std::get<1>(disl));
                    double exponent = exp(-distsq(centerpoint, disl_pos_xy) / (2 * sigma * sigma));

                    if (std::get<2>(disl) == -1) // negative dislocation
                        rho_n[y / subs][x / subs] += exponent;
                    else // positive dislocation
                        rho_p[y / subs][x / subs] += exponent;
                }
            }

        normalize(rho_p, (double)dislocs.size() / 2);
        normalize(rho_n, (double)dislocs.size() / 2);

        std::vector<std::vector<double>> rho_t(res, std::vector<double>(res, 0)); // the 2D map for rho_total
        std::vector<std::vector<double>> kappa(res, std::vector<double>(res, 0)); // the 2D map for kappa = rho_signed
        for (unsigned int y = 0; y < res; ++y)
            for (unsigned int x = 0; x < res; ++x)
            {
                rho_t[y][x] = rho_p[y][x] + rho_n[y][x];
                kappa[y][x] = rho_p[y][x] - rho_n[y][x];
            }

        if (create_maps)
        {
            o_rho_t << rho_t;
            o_kappa << kappa;
        }
    }

    if (um == ws) // Wigner-Seitz case
    {
        std::sort(dislocs.begin(), dislocs.end(), [](const disl& a, const disl& b) {return (std::get<1>(a) > std::get<1>(b)); }); // dislocaions are sorted based on their y value
        std::vector<int> points(dislocs.size(), 0);
        for (unsigned int y = 0; y < samp; ++y) // samp =  res * subs, it defines a more fine mesh
            for (unsigned int x = 0; x < samp; ++x)
            {
                pair centerpoint((x + 0.5) / samp - 0.5, (y + 0.5) / samp - 0.5); // centerpoints of the fine mesh

            }
    }

    std::cout << "Done." << std::endl;
    return 0;
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

    for (unsigned int i = 0; i < res; ++i)
        ofile << in[i] << "\t" << correlation[i] / res << "\t" << correlation_fftw[i] << "\t" << FourierComp_fftw[i < res - i ? i : res - i] << "\n";

}
#pragma endregion