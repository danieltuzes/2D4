//
// xautocorr.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include "xautocorr_utils.h"

void measure_area(std::vector<disl>&, size_t);
void measure_density(const std::vector<disl>& dislocs, std::vector<double>& line, size_t, size_t);
void measure_density(const std::vector<disl>&, std::vector<std::vector<double>>&, size_t);
void nearestDislIndex(const std::vector<disl>&, size_t&, double&, double&, double, double);


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
    else //
    {
        // input-path
        if (0 == vm.count("input-path")) // to test if N is present
        {
            std::cerr << "input-path is missing! Program terminates.\n";
            exit(-1);
        }
    }

    size_t subs = vm["sub-sampling"].as<int>(); // subs-sampling rate
    size_t res = vm["resolution"].as<int>(); // the resolution of the autocorrelation map
    size_t samp = res * subs; // sampling rate for method ws
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

    int debug_level = vm["method"].as<int>();
    std::cout << "Debug level is " << debug_level << ". ";
    if (debug_level == 0)
        std::cout << "No debug information will be shown." << std::endl;
    else
        std::cout << "Debug information will be shown. The name of these files start with deb." << std::endl;


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
    size_t size = dislocs.size();

#pragma region create maps if asked for
    std::string o_rho_p_fn = ifname + "_rho_p.txt";
    std::string o_rho_n_fn = ifname + "_rho_n.txt";
    std::string o_rho_t_fn = ifname + "_rho_t.txt";
    std::string o_kappa_fn = ifname + "_kappa.txt";
    std::ofstream o_rho_p;
    std::ofstream o_rho_n;
    std::ofstream o_rho_t;
    std::ofstream o_kappa;

    if (create_maps)
    {
        o_rho_p.open(o_rho_p_fn);
        o_rho_n.open(o_rho_n_fn);
        o_rho_t.open(o_rho_t_fn);
        o_kappa.open(o_kappa_fn);
        if (!o_rho_p)
        {
            std::cerr << "Cannot create " << o_rho_p_fn << ". Program terminates." << std::endl;
            exit(-1);
        }
        if (!o_rho_n)
        {
            std::cerr << "Cannot create " << o_rho_n_fn << ". Program terminates." << std::endl;
            exit(-1);
        }
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
        o_rho_p << "# This file contains the frequency of rho_p for the file " << ifname << " using the " << methodname << " method.\n";
        o_rho_n << "# This file contains the frequency of rho_n for the file " << ifname << " using the " << methodname << " method.\n";
        o_rho_t << "# This file contains the frequency of rho_t for the file " << ifname << " using the " << methodname << " method.\n";
        o_kappa << "# This file contains the frequency of kappa for the file " << ifname << " using the " << methodname << " method.\n";
    }

#pragma endregion

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

        for (size_t y = 0; y < samp; ++y) // samp =  res * subs, it defines a more fine mesh
            for (size_t x = 0; x < samp; ++x)
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

        normalize(rho_p, (double)size / 2);
        normalize(rho_n, (double)size / 2);

        std::vector<std::vector<double>> rho_t(res, std::vector<double>(res, 0)); // the 2D map for rho_total
        std::vector<std::vector<double>> kappa(res, std::vector<double>(res, 0)); // the 2D map for kappa = rho_signed
        for (size_t y = 0; y < res; ++y)
            for (size_t x = 0; x < res; ++x)
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
        std::sort(dislocs.begin(), dislocs.end(), [](const disl& a, const disl& b) {return (std::get<1>(a) + std::get<2>(a)) > (std::get<1>(b) + std::get<2>(b)); }); // dislocaions are sorted based on their type and y value
        std::vector<disl> dislocs_p(size / 2); // only the positive dislocations; the 3rd element (int) will be used to measure the area
        std::vector<disl> dislocs_n(size / 2); // only the negative dislocations; the 3rd element (int) will be used to measure the area

        std::copy_n(dislocs.begin(), size / 2, dislocs_p.begin());
        std::copy_n(dislocs.begin() + size / 2, size / 2, dislocs_n.begin());

        // zero out the area of the dislocations
        for (auto& disloc : dislocs_n)
            std::get<2>(disloc) = 0;

        for (auto& disloc : dislocs_p)
            std::get<2>(disloc) = 0;

        measure_area(dislocs_p, samp);
        measure_area(dislocs_n, samp);

        if (debug_level)
        {
            std::ofstream debo_dislpfile("deb_disl_p.txt");
            debo_dislpfile << "# The positive dislocations with their corresponding area.\n";
            std::ofstream debo_dislnfile("deb_disl_n.txt");
            debo_dislnfile << "# The negative dislocations with their corresponding area.\n";

            for_each(dislocs_p.begin(), dislocs_p.end(),
                [&debo_dislpfile](const disl& a) {debo_dislpfile << std::get<0>(a) << "\t" << std::get<1>(a) << "\t" << std::get<2>(a) << "\n"; }); // print out to ofile
            for_each(dislocs_n.begin(), dislocs_n.end(),
                [&debo_dislnfile](const disl& a) {debo_dislnfile << std::get<0>(a) << "\t" << std::get<1>(a) << "\t" << std::get<2>(a) << "\n"; }); // print out to ofile
        }
        
        if (create_maps)
        {
            std::vector<std::vector<double>> rho_p(res, std::vector<double>(res, 0)); // the 2D map for rho_total
            std::vector<std::vector<double>> rho_n(res, std::vector<double>(res, 0)); // the 2D map for kappa = rho_signed
            std::vector<std::vector<double>> rho_t(res, std::vector<double>(res, 0)); // the 2D map for kappa = rho_signed
            std::vector<std::vector<double>> kappa(res, std::vector<double>(res, 0)); // the 2D map for kappa = rho_signed

            measure_density(dislocs_p, rho_p, samp);
            measure_density(dislocs_n, rho_n, samp);

            for (size_t i = 0; i < res; ++i)
                for (size_t j = 0; j < res; ++j)
                {
                    rho_t[i][j] = rho_p[i][j] + rho_n[i][j];
                    kappa[i][j] = rho_p[i][j] - rho_n[i][j];
                }

            o_rho_p << rho_p;
            o_rho_n << rho_n;
            o_rho_t << rho_t;
            o_kappa << kappa;
        }
    }

    std::cout << "Done." << std::endl;
    return 0;
}

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

void measure_density(const std::vector<disl>& dislocs, std::vector<double>& linedensity, size_t yid, size_t samp)
{
    size_t size = linedensity.size(); // number of density points
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
            linedensity[x / subs] += (double)1 / std::get<2>(dislocs[cid]);
            lastdist += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
            lastdistsq = lastdist * lastdist;
        }
        lastdist += (double)1 / samp; // if one moves in x, it increases the distance, this is an overestimation
        lastdistsq = lastdist * lastdist;
    }
}

void measure_density(const std::vector<disl>& dislocs, std::vector<std::vector<double>>& map, size_t samp)
{
    for (size_t yid = 0; yid < map.size(); ++yid)
        measure_density(dislocs, map[yid], yid, samp);
}

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