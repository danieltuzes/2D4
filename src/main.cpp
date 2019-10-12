
// main.cpp : This file contains the 'main' function for 2D4_sim. Program execution begins and ends there.

#include "project_parser.h"
#include "simulation.h"

int main(int argc, char** argv)
{
    // Parsing the configuration
    sdddstCore::ProjectParser parser(argc, argv);

    // Init the simulation
    sdddstCore::Simulation simulation(parser.getSimulationData());

    // Run the simulation
    simulation.run();

    return 0;
}
