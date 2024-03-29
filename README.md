# 2D4
A 2D discrete dislocation dynamics simulation program toolset.

## Short description
This toolkit contains tools to simulate 2D single slip edge dislocation systems under periodic boundary conditions. The integrator is based on an implicit numerical scheme which makes it possible to keep the O(*N*<sup>2</sup>) complexity which arise from the pair interactions while no dislocation annihilation is required and the runtime is greatly decreased.

To decrease the runtime, a matrix inversion on a matrix containing up to *N*<sup>2</sup> number of elements needs to be performed leading to a O(*N*<sup>3</sup>) complexity. By neglecting the gradient of the distant dislocation-pairs, the number of elemnts can be decreased to be proportional with *N* so that the overall complexity of the simulation remains O(*N*<sup>2</sup>). In case of not applying such simplification, the contant factor for the O(*N*<sup>3</sup>) can be still kept low so that the main leading term in runtime for typical system sizes (below 16'000 dislocations) is characterized by the O(*N*<sup>3</sup>) term.

The details of the numerical scheme and the implementation with the achived results [are published](https://iopscience.iop.org/article/10.1088/1361-651X/ab76b2/meta) and a [free arxiv version is also available](https://arxiv.org/abs/1909.05706).

### Tools
The toolkit contains the following tools
1. The [**simulation program** *2D4_sim*](https://github.com/danieltuzes/2D4/tree/master/src) is the main component of this toolset. It evols the dislocation configuration under the prescribed external stress.
2. [**Dislocation system generator** *init_config_gen*](https://github.com/danieltuzes/2D4/init_config_gen) to create the initial configuration of uncorrelated dislocations.
3. **Evaluation programs**

   1. [*xpattern*](xpattern) performs analysis on the simulations obatined looking for patterns
   2. [*conf_compare*](conf_compare) compares simulations and tell if they are the same or not, where the largest deviation is and what the average deviation is.
4. Sandbox programs

    These programs help to identify or measure the effect of one specific change in the code. The following projects are available:

   * [*ieee_hyperbolic*](https://github.com/danieltuzes/2D4/tree/master/sandbox/ieee_hyperbolic): to decide whether the calculation of hyperbolic function using identities are faster
   * [*merge_data_func*](https://github.com/danieltuzes/2D4/tree/master/sandbox/ieee_hyperbolic): merges and averages function evaluated at different values.

The rest of this file belongs to the description of the simulation program 2D4_sim.

## Build & run 2D4_sim
This solution uses several external libraries, such as **umfpack**, **boost**, **FFTW** and furthermore, to keep the code simple, some additional include libraries must be set up for your compiler. A convenient Linux-gcc-cmake built procedure (scenario A) is provided along with a Windows-VS-vcpkg built procedure (scenario B).

### Scenario A
This is the Linux-gcc-cmake case not maintained in this repository. To be able to build the simulator, you will need the following dependencies on your computer:

* g++ or clang
* cmake (at least version 3.8)
* make
* pyhton3 with python3-dev
  * this is only needed if python interface is required in `CMakeLists.txt`
  * this can be found at line 23: `option(BUILD_PYTHON_BINDINGS "Build python interface package" OFF)`
  * a system-wide installation is required by default, otherwise, cmake will not find the libs
  * using a local (user) version of python is also possible: define paths for cmake with the switch `-DPYTHON_LIBRARY=~/.local/lib/libpython` and `-DPYTHON_INCLUDE_DIR=~/.local/include/python3.7m` if you installed python3-dev into your `$HOME/.local` directory
* umfpack from suitesparse
* boost (the program options, and python libraries, if python is required)
* FFTW

The build is fairly simple from the root directory of the source:

```bash
mkdir build
cd build
cmake ..
make
```

The resulting binary which can be used to run the simulations:

```bash
build/src/sdddst
```
You can safely **delete** the files corresponding to scenario B:
* 2D4_sim.vcxproj
* 2D4_sim.vcxproj.filters
* 2D4_sim.vcxproj.user
* 2D4.sln

### Scenario B
This is the Windows-VS-vcpkg case. The project can be compiled for x64 compatible machines only [due to umfpacks's openblas dependency](https://github.com/microsoft/vcpkg/issues/2095), which is available only for x64. This case comes with no python interface yet. Follow these instructions to be able build this project.
1. Install [vcpkg](https://github.com/microsoft/vcpkg), a C++ package manager for Windows.
   1. Open a PowerShell (abbreviated as PS) and first create the directory where you would like to place it. Let's say, `C:\local`, so in PS, execute `mkdir C:\local`. Move there by executing `cd C:\local`.
   2. Download the vcpkg package with git or by downloading it with a browser from [its website](https://github.com/Microsoft/vcpkg).
	3. Extract the files with their parent folder vcpkg-master from the compressed file vcpkg-master.zip. Move it to `C:\local`. Navigate your PS there by `cd vcpkg-master`.
   3. Install the program by executing `.\bootstrap-vcpkg.bat`. The installar may ask for admin privileges.
   4. To use the installed packages automatically, execute `.\vcpkg integrate install`.
2. Install the required dependencies with *vcpkg* by executing
      
      1. `.\vcpkg install boost-program-options:x64-windows` for boost program options,
      2. `.\vcpkg install boost-random:x64-windows` for boost random number generation,
      3. `.\vcpkg install fftw3:x64-windows` for FFTW, and
      4. `.\vcpkg install suitesparse:x64-windows`, this installs the umfpack.
      
	These steps may take at least around 15 minutes.
	
You can safely **delete** building files and folders for scenario A:
* cmake/
* CMakeLists.txt
* sdddstCMakeConfig.h.in
### Configuration files
For the available parameters, check out the help!

```bash
./sdddst --help
```

To be able to run a simulation, data has to be provided in plain text format. The slip planes of the dislocations are parallel with the x axis and the simulation area is between [-0.5, 0.5] in both directions. Based on that an example configuration file which contains dislocation data:

```
-0.1 0.3 1
0.4 0.2 1
-0.2 -0.3 -1
-0.3 0.3 1
...
```

Each line represent a dislocation: in the first column the x coordinates are present, while in the second one the y coordinates can be found. The last column can be either 1 or -1 based on in which direction the dislocation's Burgers vector point.
Point defects represented in files are the same in structure, but without the last column. Point defects are fixed in place during the simulations.

### Field of a dislocation
To be able to simulate dislocation interactions, a field need to be defined. These should be periodic and should reflect the size of the simulation cell. The current default one uses a binary datablob which contains precalculated data. The binary (periodic_stress_xy_1024x1024_bin.dat) need to be in the present working directory, or the path has to be defined with the corresponding option. Do not include the name of the binary at the end of the path!

### External stress
The default protocol for external stress gives 0 for every timestep, but other protocols can be choosen as well, see the help.

### Structure of the log file
If a log file is requested, a file will be continously updated during the simulation, but be aware that it is only written into the file when the buffer needs to be emptied. Each line represent a successful step in the simulation and in order the meaning of the columns:

* simulation time
* number of successful steps
* number of failed steps
* worst error ratio squared
* average speed of the dislocations
* cutoff (used in the implicit scheme)
* order parameter
* value of the external stress
* computation time between the last two successful steps
* accumulated strain
* average *v*<sup>2</sup>
* energy of the system

### Cutoff multiplier
A cutoff parameter is needed for this implicit method. The meaning of the parameter is that if it is infinite the calculation goes like an implicit method was used, but if it is zero, it is like an explicit method. The multiplier multiplied with one on square root N (where N is the number of the dislocations) results in the actual cutoff.

## Python interface
A minimalistic Python interface is also available which makes it possible to run simulations directly from python. A simple description about how to run a simulation can be found below:
### Compile the python module
Just like before, cmake is responsible for the compilation process. The dependencies:

* Just like above
* Boost.Python (already covered by boost)
* Python3 and the developer libraries as well

After the dependencies are in place, the following should be issued in the previously created `build` folder:

```bash
cmake .. -DBUILD_PYTHON_BINDINGS On
make
```

As a result `PySdddstCore.so` named shared library going to be created in the build folder. This should be placed into your path and after that it can be imported into python by issueing the following command (a python package can be created with it as well):

```python
import PySdddstCore as psc
```

### Example
First we need to import the library into python:

```python
import PySdddstCore as psc
```

The next step is to prepare the simulation data. For that we need to use a `SimulationDataObject`. Let's say we are going to simulate the motion of some dislocations, and their data is in `dislocation.dconf` in the current working directory:

```python
simulation_data = psc.SimulationDataObject("dislocation.dconf", "")
```

(Caution: If the paths do not point to a valid file, the program will crash. This will be fixed later on.)

The properties of the simulation will be available through this object during the whole simulation. Set up the simulation parameters:

```python
simulation_data().cutoff_multiplier = 0.5
simulation_data().precision = 1e-6
simulation_data().time_limit = 100
simulation_data().time_limited = True
simulation_data().calculate_strain_during_simulation = True
simulation_data().final_dislocation_configuration_path = "result.dconf"
simulation_data().update_cutoff()
```

The actual values can be obtained by simply invoking the parameters:

```python
simulation_data().cutoff_multiplier
simulation_data().precision
simulation_data().time_limit
simulation_data().time_limited
simulation_data().calculate_strain_during_simulation
simulation_data().final_dislocation_configuration_path
```

Two important thing is still missing. The external stress field (even if it is zero) and the dislocation field model has to be specified. First we specify a zero external stress field:

```python
external_stress = psc.StressProtocolObject()
```

This is not a valid object right now. We should initialize it first:

```python
external_stress.init()
```

Validity can be checked by:

```python
external_stress.valid()
```

We need to provide this object to the simulation as well by:

```python
simulation_data().external_stress = external_stress
```

After that the `StressProtocolObject` becomes invalid again.

The same should be applied in case of the dislocation model:

```python
field = psc.AnalyticFieldObject()
field.init()
simulation_data().tau = field
```

After we are ready, the simulation can be created and run by:

```python
simulation = psc.Simulation.create(simulation_data)
simulation.run()
```

## Impressum, credits, authors and copyright
This is a private fork of the work pgabor/sdddst. The source files imported from that repo may have different copyright, and most of those files start with a copyright section. All other modifications are under my copyright if applicable.