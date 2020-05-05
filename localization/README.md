# 2D4::localization
The purpose of this program is to tell the level of localization of a dilsocation configuration. The width of the dislocation-ensemble-band in a system seen from a plane parallel to y at x_0 can be defined as

\[d\left( {x_0} \right) = \sum\limits_i {d\left( {{x_i},{x_0}} \right)} /N\]
where $d\left( {{x_i},{x_0}} \right)$ is the distance between the *i*th dislocation and the plane at $x_i$ considering the periodic boundary conditions, and the summation in *i* goes through the dislocation indices from 1 to *N*.

The lowest value of the width for all planes is a unique value and assigned to the dislocation configuration, and denoted by $d = min_{x_0 \in [0; 1)} d(x_0)$. For uniform particle distributions, the *d* is one fourth of the system size, $d = 1/4$, as well as for any odd or even periodic function (with integer wave number larger than 1). The localization is then a measure defined by $η = 1 - 4d$.

## Dependencies
 boost/program_options

## Parameters
* The program takes arguments as filenames for which the value $x_0$ and η will be printed out to stdout.

## Calculating methog

To make a closed formula for the width one can represent the *x* value of each dislocation on the unit circle, calculate the center of mass and then project it onto the circe at $φ_0$. The width multiplied by 2π will given by the formula mentioned above using d to measure angles.