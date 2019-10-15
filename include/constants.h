
// constants.h : contains some constant expressions from different places of the code. Included in simulation_data.cpp, project_parser.cpp, precision_handler.cpp and AnalyticField.h

/*changelog
# 0.2 
USE_IEEE_HYPERBOLIC added so AnalyticField.cpp can check against this to decide which method to use. If 1, uses original, IEEE compliant method, and uses a faster but less precise version otherwise

# 0.1
first version with VERSION_constants
*/

#ifndef SDDDST_CORE_CONSTANTS_H
#define SDDDST_CORE_CONSTANTS_H

#define VERSION_constants 0.2

#define ANALYTIC_FIELD_N 4
#define EPS 1e-12
#define DEFAULT_CUTOFF_MULTIPLIER 1.0
#define DEFAULT_CUTOFF 1.0
#define DEFAULT_PRECISION 1e-8
#define DEFAULT_ITERATION_COUNT 2
#define DEFAULT_TIME_LIMIT 0.0
#define DEFAULT_STEP_SIZE 1.0
#define DEFAULT_SIM_TIME 0.0
#define DEFAULT_KASQR (1.65 * 1.65 * 1e6 / 256)
#define DEFAULT_A (1e-4 * 16)

#define USE_IEEE_HYPERBOLIC 1
#if (ANALYTIC_FIELD_N == 4)
#undef USE_IEEE_HYPERBOLIC // I'd like to give a value again
#define USE_IEEE_HYPERBOLIC 0 // if false, program calculates hyperbolic functions with identities: faster but less precise
#if !USE_IEEE_HYPERBOLIC
#define cosh2pi 267.746761483748222245931879901    // cosh(2 * pi)
#define cosh4pi 143375.656570070321051450310133    // cosh(4 * pi)
#define cosh6pi 7.67764676977233502175191858009e7  // cosh(6 * pi)
#define cosh8pi 4.11131577927974976374895337541e10 // cosh(8 * pi)
#define coshTpi 2.20157529303160145057002722946e13 // cosh(10* pi)

#define sinh2pi 267.744894041016514257117449688    // sinh(2 * pi)
#define sinh4pi 143375.656566582978695241314640    // sinh(4 * pi)
#define sinh6pi 7.67764676977233437051070497210e7  // sinh(6 * pi), basically cosh6pi
#define sinh8pi 4.11131577927974976374773721974e10 // sinh(8 * pi), basically cosh8pi
#define sinhTpi 2.20157529303160145057002722719e13 // sinh(10* pi), basically coshTpi

#define cosh__2pi_xp1 (cosh2pix * cosh2pi + sinh2pix * sinh2pi)  // == cosh( 2*pi * (x+1) )
#define cosh__2pi_xm1 (cosh2pix * cosh2pi - sinh2pix * sinh2pi)  // == cosh( 2*pi * (x-1) )
#define cosh__2pi_xp2 (cosh2pix * cosh4pi + sinh2pix * sinh4pi)  // == cosh( 2*pi * (x+2) )
#define cosh__2pi_xm2 (cosh2pix * cosh4pi - sinh2pix * sinh4pi)  // == cosh( 2*pi * (x-2) )
#define cosh__2pi_xp3 (cosh2pix * cosh6pi + sinh2pix * sinh6pi)  // == cosh( 2*pi * (x+3) )
#define cosh__2pi_xm3 (cosh2pix * cosh6pi - sinh2pix * sinh6pi)  // == cosh( 2*pi * (x-3) )
#define cosh__2pi_xp4 (cosh2pix * cosh8pi + sinh2pix * sinh8pi)  // == cosh( 2*pi * (x+4) )
#define cosh__2pi_xm4 (cosh2pix * cosh8pi - sinh2pix * sinh8pi)  // == cosh( 2*pi * (x-4) )
#define cosh__2pi_xp5 (cosh2pix * coshTpi + sinh2pix * sinhTpi)  // == cosh( 2*pi * (x+5) )
#define cosh__2pi_xm5 (cosh2pix * coshTpi - sinh2pix * sinhTpi)  // == cosh( 2*pi * (x-5) )

#define sinh__2pi_xp1 (sinh2pix * cosh2pi + cosh2pix * sinh2pi)  // == sinh( 2*pi * (x+1) )
#define sinh__2pi_xm1 (sinh2pix * cosh2pi - cosh2pix * sinh2pi)  // == sinh( 2*pi * (x-1) )
#define sinh__2pi_xp2 (sinh2pix * cosh4pi + cosh2pix * sinh4pi)  // == sinh( 2*pi * (x+2) )
#define sinh__2pi_xm2 (sinh2pix * cosh4pi - cosh2pix * sinh4pi)  // == sinh( 2*pi * (x-2) )
#define sinh__2pi_xp3 (sinh2pix * cosh6pi + cosh2pix * sinh6pi)  // == sinh( 2*pi * (x+3) )
#define sinh__2pi_xm3 (sinh2pix * cosh6pi - cosh2pix * sinh6pi)  // == sinh( 2*pi * (x-3) )
#define sinh__2pi_xp4 (sinh2pix * cosh8pi + cosh2pix * sinh8pi)  // == sinh( 2*pi * (x+4) )
#define sinh__2pi_xm4 (sinh2pix * cosh8pi - cosh2pix * sinh8pi)  // == sinh( 2*pi * (x-4) )
#define sinh__2pi_xp5 (sinh2pix * coshTpi + cosh2pix * sinhTpi)  // == sinh( 2*pi * (x+5) )
#define sinh__2pi_xm5 (sinh2pix * coshTpi - cosh2pix * sinhTpi)  // == sinh( 2*pi * (x-5) )
#endif
#endif

#endif
