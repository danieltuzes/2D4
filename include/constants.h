
// constants.h : contains some constant expressions from different places of the code. Included in simulation_data.cpp, project_parser.cpp, precision_handler.cpp and AnalyticField.h

#ifndef SDDDST_CORE_CONSTANTS_H
#define SDDDST_CORE_CONSTANTS_H

#define VERSION_constants 0.1

#define SCALE_FACTOR_AALTO 200.0 // 200 b sized system
#define EPS 1e-12
#define ANALYTIC_FIELD_N 4
#define DEFAULT_CUTOFF_MULTIPLIER 1.0
#define DEFAULT_CUTOFF 1.0
#define DEFAULT_PRECISION 1e-8
#define DEFAULT_ITERATION_COUNT 2
#define DEFAULT_TIME_LIMIT 0.0
#define DEFAULT_STEP_SIZE 1.0
#define DEFAULT_SIM_TIME 0.0
#define DEFAULT_KASQR (1.65 * 1.65 * 1e6 / 256)
#define DEFAULT_A (1e-4 * 16)
#define DEFAULT_EXTERNAL_FIELD 0.0

#endif
