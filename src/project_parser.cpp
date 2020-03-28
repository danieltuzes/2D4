//
// project_parser.cpp : contains the function definitions for project_parser.h

#include "constants.h"
#include "project_parser.h"
#include "stress_protocol.h"
#include "suitesparse/SuiteSparse_config.h"

#include <iostream>

namespace bpo = boost::program_options;


sdddstCore::ProjectParser::ProjectParser(int argc, char** argv) :
    sD(nullptr)
{
    bpo::options_description requiredOpt("Required options");           // must be set from command line
    bpo::options_description optionalOpt("Group of optional options");  // either has value or program runs without them; can be set from command line

    // long list of grouped optional options
    {
        bpo::options_description IOOpt("Input-output options");         // input - output paths
        bpo::options_description POpt("Precision options");             // precision, also affects program speed
        bpo::options_description LOpt("Simulation limit options");      // limiting simulation by stopping it if a condition is met
        bpo::options_description AOpt("Analysis options");              // in simulation calculation for further investigation
        bpo::options_description SOpt("Ext. stress protocol options");  // how the external stress should increase

        requiredOpt.add_options()
            ("dislocation-configuration,I", bpo::value<std::string>(), "plain text file path containing dislocation data in (x y b) triplets")
            ("result-dislocation-configuration,O", bpo::value<std::string>(), "path where the result configuration will be stored at the end of the simulation")
            ;

        IOOpt.add_options()
            ("logfile-path,L", bpo::value<std::string>(), "path for the plain text log file (it will be overwritten if it already exists)")
            ("save-sub-configurations,o", bpo::value<std::string>(), "saves the current configuration after every N successful step to the given destination")
            ("sub-configuration-delay,N", bpo::value<unsigned int>()->default_value(5), "number of successful steps between the sub configurations written out")
            ("sub-config-times,T", bpo::value<double>()->default_value(0), "subconfigs must be written out at simulation times for all n pos integer (input value T) at\na) T*n \nb) initial-stepsize * T^n")
            ("sub-config-times-type,b", bpo::value<char>()->default_value('b'), "subconfigs must be written out at simulation times for all n pos integer (input value a or b) at\na) T*n \nb) initial-stepsize * T^n")
            ("sub-configuration-delay-during-avalanche,n", bpo::value<unsigned int>()->default_value(1), "number of successful steps between the sub configurations written out during avalanche if avalanche detection is on")
            ("point-defect-configuration", bpo::value<std::string>(), "plain text file path containing point defect data in (x y) pairs");

        POpt.add_options()
            ("position-precision,P", bpo::value<double>()->default_value(1e-5, "1e-5"), "minimum precision for the positions for the adaptive step size protocol")
            ("cutoff-multiplier,u", bpo::value<double>()->default_value(1e20), "multiplier of the 1/sqrt(N) cutoff parameter")
            ("heaviside-cutoff,h", "The weight in the Jacobian is a heavyside step function of the distance with charasteristic value of cutoff-multiplier")
            ("initial-stepsize", bpo::value<double>()->default_value(1e-6, "1e-6"), "first tried step size for the simulation")
            ("max-stepsize,M", bpo::value<double>(), "the stepsize can not exceed this value")
            ("weight-function,w", bpo::value<char>()->default_value('c'), "Weights as the function of the A_ii\nc) coded first, 1/(1+1/s)^2\np) as in paper, 1/(1+1/s)\nm) mathematical hint: (1 - s + 1/(1+s) - 2*exp(-s) ) / (1 - s - 1/(1+s) )")
            ("dipole-precision,p", bpo::value<double>()->default_value(0.05, "0.05"), "minimum precision with respect to the nearest dislocation; use 0 to disable this feature");

        LOpt.add_options()
            ("time-limit,t", bpo::value<double>(), "the simulation stops if simulation time reached this limit")
            ("step-count-limit,C", bpo::value<unsigned int>(), "the simulation will stop after this many steps")
            ("strain-increase-limit", bpo::value<double>(), "the simulation stops if strain increase reaches this value")
            ("avalanche-detection-limit", bpo::value<unsigned int>(), "the simulation will stop after the threshold was reached with the given number of events from above");

        AOpt.add_options()
            ("calculate-strain,S", "turns on strain calculation for the simulation")
            ("avalanche-speed-threshold", bpo::value<double>()->default_value(1e-3), "speed threshold for counting avalanches")
            ("calculate-order-parameter,l", "turns on order parameter calculation during the simulation");

        SOpt.add_options()
            ("const-external-stress,s", bpo::value<double>()->default_value(0), "the constant in the external stress during the simulation")
            ("fixed-rate-external-stress,r", bpo::value<double>(), "the slope of external stress - time function (disabled by default)")
            ("cyclic-external-stress,i", bpo::value<double>(), "the time period of the cyclic load")
            ("spring-constant", bpo::value<double>(), "simple model of an experiment where a spring is used, the arg should be the spring constant (it is valid only with a fixed-rate-external-stress)");

        optionalOpt.add(IOOpt).add(POpt).add(LOpt).add(AOpt).add(SOpt);
    }

    bpo::options_description options("Possible program call options");

    options.add(requiredOpt).add(optionalOpt).add_options()
        ("help", "shows this help")
        ("hide-copyright,c", "hides the copyright notice from the standard output");

    bpo::variables_map vm;

    try {
        bpo::store(bpo::parse_command_line(argc, argv, options), vm, true);
    }
    catch (bpo::error & e)
    {
        std::cerr << e.what() << std::endl;
        exit(-1);
    }

    if (!vm.count("hide-copyright")) // if the user doesn't hide the bloaty copyright text
    {
        printLicense(argc, argv);
    }

    if (vm.count("help")) // if the user only interested in the help, there is no need to check the variables
    {
        std::cout
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

void sdddstCore::ProjectParser::printLicense(int argc, char** argv) // if used multiple times, should move outside this cpp
{
    double totalVersion =
        VERSION_constants +
        VERSION_dislocations +
        VERSION_field +
        VERSION_point_defect +
        VERSION_precision_handler +
        VERSION_project_parser +
        VERSION_simulation +
        VERSION_simulation_data +
        VERSION_stress_protocol;
    std::cout << "This is 2D4_sim (version " << totalVersion;
#ifdef DEBUG_VERSION
    std::cout << ", with DEBUG functions version " << DEBUG_VERSION;
#endif
    std::cout << "),\n"
        << "a 2D discrete dislocation dynamics simulation program toolset based on sdddst.\n"
        << "See README.md for copyright.\n"
        << "Program call: \n";
    for (int i = 0; i < argc; ++i)
        std::cout << argv[i] << " ";
    std::cout << std::endl;

    std::cout
        << "Detailed version and compiler info:\n"
        << "VERSION_constants:                  " << VERSION_constants << "\n"
        << "VERSION_dislocations:               " << VERSION_dislocations << "\n"
        << "VERSION_field:                      " << VERSION_field << "\n"
        << "VERSION_point_defect:               " << VERSION_point_defect << "\n"
        << "VERSION_precision_handler:          " << VERSION_precision_handler << "\n"
        << "VERSION_project_parser:             " << VERSION_project_parser << "\n"
        << "VERSION_simulation:                 " << VERSION_simulation << "\n"
        << "VERSION_simulation_data:            " << VERSION_simulation_data << "\n"
        << "VERSION_stress_protocol:            " << VERSION_stress_protocol << "\n"
        << "USE_IEEE_HYPERBOLIC                 " << USE_IEEE_HYPERBOLIC << "\n"
        << "SUITESPARSE_VERSION:                " << XSTR(SUITESPARSE_MAIN_VERSION) << "." << XSTR(SUITESPARSE_SUB_VERSION) << "." << XSTR(SUITESPARSE_SUBSUB_VERSION) << "\n"
        << "COMPILER_VERSION:                   " << XSTR(COMPILER_VERSION) << "\n"
        << "MACHINE_INFO:                       " << XSTR(MACHINE_INFO) << "\n"
        << "USR_COMP_OPTIONS:                   " << XSTR(USR_COMP_OPTIONS) << std::endl;
}

void sdddstCore::ProjectParser::processInput(bpo::variables_map& vm)
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

    sD->stepSize = vm["initial-stepsize"].as<double>();     // will be modified later
    sD->initStepSize = vm["initial-stepsize"].as<double>(); // constant during the simulation
    sD->prec = vm["position-precision"].as<double>();
    sD->dipole_prec = vm["dipole-precision"].as<double>();
    sD->weightFunc = vm["weight-function"].as<char>();
    if (sD->weightFunc != 'c' && sD->weightFunc != 'p' && sD->weightFunc != 'm')
    {
        std::cerr << "Unsupported weight-function " << sD->weightFunc << ", supported values are c, p and m. Program terminates." << std::endl;
        exit(-1);
    }
    sD->cutOffMultiplier = vm["cutoff-multiplier"].as<double>();
    sD->heavisideCutoff = vm.count("heaviside-cutoff");
    sD->updateCutOff();

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

    sD->calculateStrainDuringSimulation = vm.count("calculate-strain");

    sD->orderParameterCalculationIsOn = vm.count("calculate-order-parameter");

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
        sD->subConfigTimes = vm["sub-config-times"].as<double>();
        sD->subConfigTimesType = vm["sub-config-times-type"].as<char>();
        sD->subConfigDelayDuringAvalanche = vm["sub-configuration-delay-during-avalanche"].as<unsigned int>();

        sD->nextWriteOutTime = sD->getNextWriteOutTime();
    }

    sD->endDislocationConfigurationPath = vm["result-dislocation-configuration"].as<std::string>();
    {
        if (!std::ofstream(sD->endDislocationConfigurationPath)) // it is better to check now than being surprised after days
        {
            std::cerr << "Cannot create result dislocation configuration file at " << sD->endDislocationConfigurationPath << ". Program terminates.\n";
            exit(-1);
        }
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

