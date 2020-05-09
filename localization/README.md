# 2D4::localization
The purpose of this program is to tell the level of localization of a dilsocation configuration.

## Measured quantities
There are two intuitive ways to define the measure of localization on the range from 0 (no localization) to 1 (full localization). 

### The width of the ensemble
The width of the dislocation-ensemble-band in a system seen from a plane parallel to $y$ at $x_0$ can be defined as

\[d\left( {x_0} \right) = \sum\limits_i {d\left( {{x_i},{x_0}} \right)} /N\]
where $d\left( {{x_i},{x_0}} \right)$ is the distance between the *i*th dislocation and the plane at $x_i$ considering the periodic boundary conditions, and the summation in *i* goes through the dislocation indices from 1 to *N*.

The lowest value of the width for all planes is a unique value and assigned to the dislocation configuration, and denoted by $d = min_{x_0 \in [0; 1)} d(x_0)$. For uniform particle distributions, the *d* is one fourth of the system size, $d = 1/4$. The localization is then a measure defined by $η = 1 - 4d$.

### Center of mass
The center of mass (CoM) with periodic boundary conditions can be defined after putting all points in the range [-0.5; 0.5) onto the circumference of the unit circle. This CoM won't be on the circumference if not all points are at the same place, therefore, needed to be transformed to represent its value on the dislocation systems's plane.
* The distance between the CoM and the center of the circle is the radius *r*, and the measure of localization is the $1-r$.
* The CoM projected onto the circumference is CoM of the dislocation ensemble in their system.

## Dependencies
 boost/program_options

## Parameters
* The program takes arguments as filenames for which the CoM, *r*, position of the plane for which *d* is minimal, and *η* are printed out.

## Calculating method

For the width of the ensemble, the program seeks through 1000 + $N$ number of points equally distributed along the *x* axis. Around the minimum, further 1000 points are investigated between the two neighboring points of the minimum, equally distributed.

The CoM and the radius can be expressed with a closed formula.