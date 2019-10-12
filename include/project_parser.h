
// project_parser.h : contains the function declarations for project_parser.cpp

#ifndef SDDDST_CORE_PROJECT_PARSER_H
#define SDDDST_CORE_PROJECT_PARSER_H

#define VERSION_project_parser 0.1

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
#define USR_COMP_OPTIONS default MSVSC++ options
#else
#define COMPILER_VERSION unknown compiler
#define MACHINE_INFO unknown machine
#define USR_COMP_OPTIONS unknown options
#endif
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

        void printLicense();

    private:
        void processInput(boost::program_options::variables_map& vm);
        std::shared_ptr<SimulationData> sD;
    };
}

#endif
