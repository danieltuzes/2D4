/*
 * SDDDST Simple Discrete Dislocation Dynamics Toolkit
 * Copyright (C) 2015-2019 Gábor Péterffy <peterffy95@gmail.com>, Dániel Tüzes <tuzes@metal.elte.hu> and their friends.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

#include "constants.h"
#include "project_parser.h"
#include "Fields/AnalyticField.h"
#include "Fields/PeriodicShearStressELTE.h"
#include "stress_protocol.h"

#include <iostream>

sdddstCore::ProjectParser::ProjectParser(int argc, char** argv) :
    sD(nullptr)
{
    boost::program_options::options_description requiredOptions("Required Options"); // must be set from command line or file
    boost::program_options::options_description optionalOptions("Optional Options"); // either has value or program runs without them
    boost::program_options::options_description fieldOptions("Available dislocation fields"); // how the interaction force should be calculated: lookup table or calculation
    boost::program_options::options_description externalStressProtocolOptions("External stress protocols (optional)"); // how the external stress should increase

    requiredOptions.add_options()
        ("dislocation-configuration,I", boost::program_options::value<std::string>(), "plain text file path containing dislocation data in (x y b) triplets")
        ("result-dislocation-configuration,O", boost::program_options::value<std::string>(), "path where the result configuration will be stored at the end of the simulation")
        ;

    optionalOptions.add_options()
        ("point-defect-configuration", boost::program_options::value<std::string>(), "plain text file path containing point defect data in (x y) pairs")
        ("logfile-path,L", boost::program_options::value<std::string>(), "path for the plain text log file (it will be overwritten if it already exists)")
        ("time-limit,t", boost::program_options::value<double>(), "the simulation stops if simulation time reached this limit")
        ("step-count-limit,C", boost::program_options::value<unsigned int>(), "the simulation will stop after this many steps")
        ("strain-increase-limit", boost::program_options::value<double>(), "the simulation stops if strain increase reaches this value")
        ("avalanche-detection-limit", boost::program_options::value<unsigned int>(), "the simulation will stop after the threshold was reached with the given number of events from above")
        ("avalanche-speed-threshold", boost::program_options::value<double>()->default_value(1e-3), "speed threshold for counting avalanches")
        ("initial-stepsize", boost::program_options::value<double>()->default_value(1e-6), "first tried step size for the simulation")
        ("cutoff-multiplier,c", boost::program_options::value<double>()->default_value(1e20), "multiplier of the 1/sqrt(N) cutoff parameter")
        ("max-stepsize,M", boost::program_options::value<double>(), "the stepsize can not exceed this value")
        ("calculate-strain,S", "turns on strain calculation for the simulation")
        ("calculate-order-parameter,l", "turns on order parameter calculation during the simulation")
        ("position-precision,P", boost::program_options::value<double>()->default_value(1e-5), "minimum precision for the positions for the adaptive step size protocol")
        ("save-sub-configurations,o", boost::program_options::value<std::string>(), "saves the current configuration after every N successful step to the given destination")
        ("sub-configuration-delay,N", boost::program_options::value<unsigned int>()->default_value(5), "number of successful steps between the sub configurations written out")
        ("sub-configuration-delay-during-avalanche,n", boost::program_options::value<unsigned int>()->default_value(1), "number of successful steps between the sub configurations written out during avalanche if avalanche detection is on")
        ;

    fieldOptions.add_options()
        ("periodic-stress-field-elte,p", boost::program_options::value<std::string>(), "periodic stress field based on ELTE library, the argument should be the folder where the compressed numerical data can be found (default)")
        ("periodic-stress-field-analytic,m", ("analytic periodic stress field (default), number of images in each direction: " + std::to_string(ANALYTIC_FIELD_N)).c_str())
        ;

    externalStressProtocolOptions.add_options()
        ("const-external-stress,s", boost::program_options::value<double>()->default_value(0), "the constant in the external stress during the simulation")
        ("fixed-rate-external-stress,r", boost::program_options::value<double>(), "the slope of external stress - time function (disabled by default)")
        ("cyclic-external-stress,r", boost::program_options::value<double>(), "the time period of the cyclic load")
        ("spring-constant", boost::program_options::value<double>(), "simple model of an experiment where a spring is used, the arg should be the spring constant (it is valid only with the previous option together)")
        ;

    boost::program_options::options_description options;

    options.add(requiredOptions).add(optionalOptions).add(fieldOptions).add(externalStressProtocolOptions).add_options()
        ("help", "shows this help")
        ("hide-copyright,c", "hides the copyright notice from the standard output");

    boost::program_options::variables_map vm;

    try {
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), vm, true);
    }
    catch (boost::program_options::error & e)
    {
        std::cerr << e.what() << std::endl;
        exit(-1);
    }

    if (!vm.count("hide-copyright")) // if the user doesn't hide the bloaty copyright text
        printLicense();

    if (vm.count("help")) // if the user only interested in the help, there is no need to check the variables
    {
        std::cout << "Detailed version and compiler info:\n"
            << "VERSION_analytic_field:\t" << VERSION_analytic_field << "\n"
            << "VERSION_constants:\t" << VERSION_constants << "\n"
            << "VERSION_dislocations:\t" << VERSION_dislocations << "\n"
            << "VERSION_field:\t" << VERSION_field << "\n"
            << "VERSION_periodic_shear_stress_elte:\t" << VERSION_periodic_shear_stress_elte << "\n"
            << "VERSION_point_defect:\t" << VERSION_point_defect << "\n"
            << "VERSION_precision_handler:\t" << VERSION_precision_handler << "\n"
            << "VERSION_project_parser:\t" << VERSION_project_parser << "\n"
            << "VERSION_simulation:\t" << VERSION_simulation << "\n"
            << "VERSION_simulation_data :\t" << VERSION_simulation_data << "\n"
            << "VERSION_stress_protocol:\t" << VERSION_stress_protocol << "\n"
            << "COMPILER_VERSION:\t" << XSTR(COMPILER_VERSION) << "\n"
            << "MACHINE_INFO:\t" << XSTR(MACHINE_INFO) << "\n"
            << "COMPILER_VERSION\t" << XSTR(COMPILER_VERSION) << "\n"
            << "USR_COMP_OPTIONS\t" << XSTR(USR_COMP_OPTIONS) << "\n"
            << options << std::endl;
        exit(0);
    }
    else // check the variables if they are set up properly
        processInput(vm);
}


std::shared_ptr<sdddstCore::SimulationData> sdddstCore::ProjectParser::getSimulationData()
{
    return sD;
}

void sdddstCore::ProjectParser::printLicense() // if used multiple times, should move outside this cpp
{
    double totalVersion = VERSION_analytic_field +
        VERSION_constants +
        VERSION_dislocations +
        VERSION_field +
        VERSION_periodic_shear_stress_elte +
        VERSION_point_defect +
        VERSION_precision_handler +
        VERSION_project_parser +
        VERSION_simulation +
        VERSION_simulation_data +
        VERSION_stress_protocol;
    std::cout << "This is 2D4_sim (version " << totalVersion << "), a 2D discrete dislocation dynamics simulation program toolset based on sdddst. See README.md for copyright.\n";
}

void sdddstCore::ProjectParser::processInput(boost::program_options::variables_map & vm)
{
    // Check for required options
    if (!vm.count("dislocation-configuration"))
    {
        std::cerr << "dislocation-configuration is missing!\n";
        exit(-1);
    }

    if (!vm.count("result-dislocation-configuration"))
    {
        std::cerr << "result-dislocation-configuration is missing!\n";
        exit(-1);
    }

    if (vm.count("point-defect-configuration"))
        sD = std::shared_ptr<SimulationData>(new SimulationData(vm["dislocation-configuration"].as<std::string>(), vm["point-defect-configuration"].as<std::string>()));
    else
        sD = std::shared_ptr<SimulationData>(new SimulationData(vm["dislocation-configuration"].as<std::string>(), ""));


    if (vm.count("time-limit"))
    {
        sD->isTimeLimit = true;
        sD->timeLimit = vm["time-limit"].as<double>();
    }

    if (vm.count("initial-stepsize"))
        sD->stepSize = vm["initial-stepsize"].as<double>();

    if (vm.count("position-precision"))
        sD->prec = vm["position-precision"].as<double>();

    if (vm.count("cutoff-multiplier"))
    {
        sD->cutOffMultiplier = vm["cutoff-multiplier"].as<double>();
        sD->updateCutOff();
    }

    if (vm.count("strain-increase-limit"))
    {
        sD->isStrainIncreaseLimit = true;
        sD->totalAccumulatedStrainIncreaseLimit = vm["strain-increase-limit"].as<double>();
    }

    if (vm.count("max-stepsize"))
    {
        sD->maxStepSizeLimit = vm["max-stepsize"].as<double>();
        sD->isMaxStepSizeLimit = true;
    }

    if (vm.count("point-defect-configuration"))
    {
        sD->A *= 1. / sqrt(sD->dc);
        sD->KASQR *= double(sD->dc);
    }

    if (vm.count("calculate-strain"))
        sD->calculateStrainDuringSimulation = true;

    if (vm.count("calculate-order-parameter"))
        sD->orderParameterCalculationIsOn = true;

    if (vm.count("logfile-path"))
    {
        std::string oLogFname = vm["logfile-path"].as<std::string>();
        sD->standardOutputLog.open(oLogFname);
        if (!sD->standardOutputLog)
        {
            std::cerr << "Error: cannot create logfile for writing with name " << oLogFname << ". Program terminates.\n";
            exit(-1);
        }
    }

    if (vm.count("step-count-limit"))
    {
        sD->isStepCountLimit = true;
        sD->stepCountLimit = vm["step-count-limit"].as<unsigned int>();
    }

    if (vm.count("avalanche-detection-limit"))
    {
        sD->countAvalanches = true;
        sD->avalancheTriggerLimit = vm["avalanche-detection-limit"].as<unsigned int>();
        sD->avalancheSpeedThreshold = vm["avalanche-speed-threshold"].as<double>();
    }

    if (vm.count("save-sub-configurations"))
    {
        sD->isSaveSubConfigs = true;
        sD->subConfigPath = vm["save-sub-configurations"].as<std::string>();
        sD->subConfigDelay = vm["sub-configuration-delay"].as<unsigned int>();
        sD->subConfigDelayDuringAvalanche = vm["sub-configuration-delay-during-avalanche"].as<unsigned int>();
    }

    sD->endDislocationConfigurationPath = vm["result-dislocation-configuration"].as<std::string>();
    {
        if (!std::ofstream(sD->endDislocationConfigurationPath)) // it is better to check now than being surprised after days
        {
            std::cerr << "Cannot create result dislocation configuration file at " << sD->endDislocationConfigurationPath << ". Program terminates.\n";
            exit(-1);
        }
    }

    if (!vm.count("periodic-stress-field-elte"))
        sD->tau = std::unique_ptr<Field>(new AnalyticField());
    else
    {
        std::unique_ptr<sdddstCoreELTE::PeriodicShearStressELTE> tmp(new sdddstCoreELTE::PeriodicShearStressELTE());
        tmp->loadStress(vm["periodic-stress-field-elte"].as<std::string>(), "xy", 1024);
        tmp->loadStress(vm["periodic-stress-field-elte"].as<std::string>(), "xy_diff_x", 1024);
        sD->tau = std::move(tmp);
    }

    double ext_stress = vm["const-external-stress"].as<double>();
    if (vm.count("cyclic-external-stress"))
    {
        if (vm.count("fixed-rate-external-stress"))
            sD->externalStressProtocol = std::unique_ptr<StressProtocol>(new CyclicLoadProtocol(ext_stress, vm["fixed-rate-external-stress"].as<double>(), vm["cyclic-external-stress"].as<double>()));
        else
        {
            std::cerr << "Error: if cyclic-external-stress is set then fixed-rate-external-stress must be also set. Program terminates.\n";
            exit(-1);
        }
    }
    else if (vm.count("fixed-rate-external-stress"))
        sD->externalStressProtocol = std::unique_ptr<StressProtocol>(new FixedRateProtocol(ext_stress, vm["fixed-rate-external-stress"].as<double>()));
    else
        sD->externalStressProtocol = std::unique_ptr<StressProtocol>(new StressProtocol(ext_stress));

    if (vm.count("change-cutoff-to-inf-under-threshold"))
    {
        sD->speedThresholdForCutoffChange = vm["change-cutoff-to-inf-under-threshold"].as<double>();
        sD->isSpeedThresholdForCutoffChange = true;
    }
}

