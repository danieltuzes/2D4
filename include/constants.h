
// constants.h : contains some constant expressions from different places of the code. Included in simulation_data.cpp, project_parser.cpp, precision_handler.cpp and AnalyticField.h

/*changelog
# 0.6
DEBUG_CODEPARTS is added to include debugging codeparts

# 0.5
DEFAULT_CUTOFF_MULTIPLIER and DEFAULT_CUTOFF are eliminated, bc they are read from input or has default boost::program_option value

# 0.4
potential bug fix around USE_IEEE_HYPERBOLIC and ANALYTIC_FIELD_N

# 0.3
In case USE_IEEE_HYPERBOLIC is false, program calculates cosh(x) and sinh(x) from the exponential function, sinh = cosh is considered for arguments larger than 6*pi.

# 0.2
USE_IEEE_HYPERBOLIC added so AnalyticField.cpp can check against this to decide which method to use. If 1, uses original, IEEE compliant method, and uses a faster but less precise version otherwise

# 0.1
first version with VERSION_constants
*/

#ifndef SDDDST_CORE_CONSTANTS_H
#define SDDDST_CORE_CONSTANTS_H

#define VERSION_constants 0.5

// uncomment if you want to include and call debugging functions
// modifications in debug functions change the affected file's version only, calls on those functions changes only DEBUG_VERSION number only
#define DEBUG_VERSION 0.2

#define ANALYTIC_FIELD_N 4
#define EPS 1e-12
#define DEFAULT_KASQR (1.65 * 1.65 * 1e6 / 256)
#define DEFAULT_A (1e-4 * 16)

#define USE_IEEE_HYPERBOLIC 0
#if (ANALYTIC_FIELD_N != 4)     // USE_IEEE_HYPERBOLIC must be true, bc the non IEEE version is implemented only for 4 images
#undef USE_IEEE_HYPERBOLIC      // I'd like to give a value again
#define USE_IEEE_HYPERBOLIC 1   // if false, program calculates hyperbolic functions with identities: faster but less precise
#endif

#if !USE_IEEE_HYPERBOLIC
#define cosh2pi 267.746761483748222245931879901    // cosh(2 * pi)
#define cosh4pi 143375.656570070321051450310133    // cosh(4 * pi)
#define cosh6pi 7.67764676977233502175191858009e7  // cosh(6 * pi)
#define cosh8pi 4.11131577927974976374895337541e10 // cosh(8 * pi)
#define coshTpi 2.20157529303160145057002722946e13 // cosh(10* pi)

#define sinh2pi 267.744894041016514257117449688    // sinh(2 * pi)
#define sinh4pi 143375.656566582978695241314640    // sinh(4 * pi)

#define cosh__2pi_x__ ((exp2pix + 1/exp2pix) / 2)
#define cosh__2pi_xp1 (((exp2pix + 1/exp2pix) * cosh2pi + (exp2pix - 1/exp2pix) * sinh2pi) / 2)
#define cosh__2pi_xm1 (((exp2pix + 1/exp2pix) * cosh2pi - (exp2pix - 1/exp2pix) * sinh2pi) / 2)
#define cosh__2pi_xp2 (((exp2pix + 1/exp2pix) * cosh4pi + (exp2pix - 1/exp2pix) * sinh4pi) / 2)
#define cosh__2pi_xm2 (((exp2pix + 1/exp2pix) * cosh4pi - (exp2pix - 1/exp2pix) * sinh4pi) / 2)
#define cosh__2pi_xp3 (cosh6pi * exp2pix)
#define cosh__2pi_xm3 (cosh6pi / exp2pix)
#define cosh__2pi_xp4 (cosh8pi * exp2pix)
#define cosh__2pi_xm4 (cosh8pi / exp2pix)
#define cosh__2pi_xp5 (coshTpi * exp2pix)
#define cosh__2pi_xm5 (coshTpi / exp2pix)

#define sinh__2pi_x__ ((exp2pix - 1/exp2pix) / 2)
#define sinh__2pi_xp1 (((exp2pix - 1/exp2pix) * cosh2pi + (exp2pix + 1/exp2pix) * sinh2pi) / 2)
#define sinh__2pi_xm1 (((exp2pix - 1/exp2pix) * cosh2pi - (exp2pix + 1/exp2pix) * sinh2pi) / 2)
#define sinh__2pi_xp2 (((exp2pix - 1/exp2pix) * cosh4pi + (exp2pix + 1/exp2pix) * sinh4pi) / 2)
#define sinh__2pi_xm2 (((exp2pix - 1/exp2pix) * cosh4pi - (exp2pix + 1/exp2pix) * sinh4pi) / 2)

#endif
#endif

