//
// project_parser.h : contains the function declarations for project_parser.cpp

/*
# 1.1
date and time always exist now, sometimes without the required values, and is printed out
Stops and gives error if (sD->subConfigTimesType == 'b' && sD->subConfigTimes < 1)

# 1.0
* sub-config-times-type is added
* weight-function is added
* initStepSize takes the value of initial-stepsize

# 0.9
* heavisideCutoff is read from heaviside-cutoff
* if (true) type of checks are eliminated

# 0.8
* Restructured program options
* cyclic-external-stress had r short switch overlapping with fixed-rate-external-stress', now it has i
* default values' silly representation is fixed with implicit textual representation

# 0.7
* checks for dipole-precision
* position-precision's default value is 1e-5 again

# 0.6
* prints out umfpack version info
* always prints out version info

# 0.5
subConfigTimes is added

# 0.4
position-precision is decreased by a factor of 20

# 0.3
program call is printed out and an endl is included, therefore the stdout will show the message during the simulation

# 0.2
cutoff multiplier is abbreviated as u, not c

# 0.1
Ther first version tracked file
*/

#ifndef SDDDST_CORE_PROJECT_PARSER_H
#define SDDDST_CORE_PROJECT_PARSER_H

#define VERSION_project_parser 1.1

// define macros for version tracking
#define STR(x) #x
#define XSTR(x) STR(x)

#ifdef __GNUC__
#define COMPILER_VERSION GNUC __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__
#define MACHINE_INFO GCC_MACHINE_INFO
#else
#ifdef _MSC_VER
#define COMPILER_VERSION MSVSC++  _MSC_VER
#define MACHINE_INFO A windows machine
#define USR_COMP_OPTIONS MSVSC++ options
#else
#define COMPILER_VERSION unknown compiler
#define MACHINE_INFO unknown machine
#define USR_COMP_OPTIONS unknown options
#endif
#endif

#if !defined(__DATE__)
#define __DATE__ "__DATE__ isn't supported"
#endif
#if !defined(__TIME__)
#define __TIME__ "__TIME__ isn't supported"
#endif

#include "simulation_data.h"

// must include other header files to provide version numbers during project_parsing
#include "precision_handler.h"
#include "simulation.h"

#include <memory>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

namespace sdddstCore
{
    class ProjectParser
    {
    public:
        ProjectParser(int argc, char** argv);

        std::shared_ptr<SimulationData> getSimulationData();

        void printLicense(int argc, char** argv);

    private:
        void processInput(boost::program_options::variables_map& vm);
        std::shared_ptr<SimulationData> sD;
    };
}

#endif
