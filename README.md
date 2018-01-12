# LINEFLOW: Library for constructing linear approximations of line flow constraints

## Overview 
LINEFLOW is a library for constructing linear approximations of line flow constraints. The user can choose between inner and outer 
approximations that have a desired quality or number of linear constraints. This repository includes the following files:

1. Source code for the library as well as files required to build it on Windows or Unix systems.
2. Source code for Matlab wrappers, which make the library easily callable from Matlab.
3. Pre-built files for library and Matlab wrappers (in 'lib' directory).
4. Example scripts of calling the library in C and Matlab. 
5. A Matlab script that plots the boundary of the feasible set of the original constraint and its constructed approximation.

- - - -

## Building library in C/C++
LINEFLOW can be built on Windows, Linux or Mac. There are no external dependencies; the only requirement is a C compiler. Below is 
the instruction on how to build the library on different OS.
### Windows (with Visual Studio)
You can use pre-built files from this repository, which were compiled on a 64 bit machine. Alternatively, you can do the following: 

1. Create a dynamic link library (see example [Walkthrough]( https://msdn.microsoft.com/en-us/library/ms235636.aspx)).
2. Add header file (`line_flow.h`), source file (`line_flow.c`) and definition file (`line_flow.def`) to the project.
3. Build the project.

### Unix (Linux/MacOS)
Use Makefile available in this repository.

- - - -

## Installation in Matlab
The library can be called from Matlab using the following interface functions:

* `LF_linearize_line`, which constructs the approximation of the line flow constraint of one line.
* `LF_linearize_system`, which constructs the approximation of all line flow constraints in the system.

These are so-called mex functions, which are written in C for efficiency but can be called from Matlab as regular Matlab functions. 
Compiled mex functions for Windows and MacOS are available from this repository. Alternatively, you can compile them yourself. This 
requires pre-built library files (.o in Unix and .lib in Windows), which you can get from this repository or by building the library 
yourself using the steps described above. Note that compiling a mex function requires that Matlab have a C compiler. The list of 
supported compilers is available [here]( https://ch.mathworks.com/support/compilers.html). To compile function `LF_linearize_XXX`, 
where `XXX` stands for `line` or `system`, follow these steps:

1. Put the following files in a folder on your Matlab path:
    - `line_flow.h`
    - `line_flow.lib` (for Windows) or `line_flow.o` (for Linux/MacOS)
    - `LF_linearize_XXX.c`

2. Make this folder your current folder in Matlab.
3. Type the following command in Matlab command line:
    - for Windows:      `mex LF_linearize_XXX.c line_flow.lib`
    - for Linux/MacOS:  `mex LF_linearize_XXX.c line_flow.o`
4. That's it! Now you can call the compiled function as a regular Matlab function.

- - - -

## Approximation algorithm
The algorithm for constructing the approximation is based on the properties of the feasible set of a line flow constraint. The algorithm 
iteratively increases the number of constructed linear constraints until either the desired approximation accuracy has been achieved or the 
maximum allowed number of linear constraints have been constructed. At each iteration, an attempt is made to equalize the approximation 
errors associated with individual linear constraints in order to increase the approximation accuracy. More information about the algorithm 
is available here. The approximation of a single line flow constraint is constructed 
using Matlab mex function `LF_linearize_line` or C function `LF_construct`. The description of algorithm's inputs and outputs is given below.

### Inputs
The inputs include the parameters of the branch, the side of the line at which the line flow constraint needs to 
be approximated, and algorithm's options.

#### Branch parameters (required)
A branch is modeled by the pi-model. The following values are required by the algorithm (all in p.u.):

- conductance,
- susceptance,
- shunt susceptance,
- transformer's tap ratio (if not transformer, set to 0),
- transformer's phase shift (if not transformer, set to 0),
- value of the thermal limit (if there is no limit, set to 0),
- lower bound on voltage magnitude at the beginning of the line,
- upper bound on voltage magnitude at the beginning of the line,
- lower bound on voltage magnitude at the end of the line,
- upper bound on voltage magnitude at the end of the line.

#### Flow side (required)
Tells the algorithm, the constraint at which side of the branch should be approximated:

- If 1, the line flow constraint at the beginning of the line will be approximated,
- If 2, the line flow constraint at the end of the line will be approximated,
- If 3 (recommended), the line flow constraints at both line ends will be approximated. Note that only the linear constraints relevant 
to the intersection of two feasible sets will be retained, which means that generally the number of retained constraints will be roughly 
the same as in options 1 or 2.

#### Algorithm's options (optional)
Tells the algorithm what the desired parameters of the approximation are and provides control parameters of some low-level functions. 
The following options are available (default values are in the brackets):

- Approximation type: inner or outer (inner),
- Maximum approximation error in the current magnitude in percent (5.0),
- Maximum number of constructed linear constraints for approximating one side of the boundary surface (15),
- Computation mode - 1 or 2 (2):
    - If 1, the algorithm constructs the approximation with the maximum given number of linear constraints, regardless of the resulting 
      approximation error,
    - If 2, the algorithm iteratively increases the number of constructed linear constraints until the desired approximation accuracy has 
      been achieved, the maximum allowed number of linear constraints have been constructed, or the change of the maximum error at two 
      consecutive iterations is smaller than a given value.
- Threshold value of the change of the maximum error at two consecutive iterations. If the actual value is smaller than the threshold, the 
  algorithm stops (0.1),
- Maximum number of adjustments that are carried out to equalize the approximation errors associated with individual linear constraints (4),
- Threshold value of the error ratio, which is defined as the minimum error associated with an individual linear constraint divided by the 
  maximum error associated with individual constraint. If actual ratio is higher than the threshold, the adjustments stop (0.9),
- Transformer model type - 0 or 1 (0):
    - If 0, the given shunt susceptance is equally split between two ends of the line,
    - If 1, the given shunt susceptance is only put at the beginning of the line. This helps model the reactive power loss associated with 
      the magnetizing current. 
- Threshold value of phase angle difference in degrees. The approximation is only constructed for the points on the boundary surface, at 
  which the angle value does not exceed the threshold (85.0),
- Maximum number of iterations in the bisection algorithm and Newton-Raphson algorithm, which are used by low-level functions (25),
- Threshold value of the change of the step in the bisection algorithm and Newton-Raphson algorithm. If the actual step is smaller than 
  the threshold, the algorithms stop (0.0001).

It is recommended to keep the last three options at their default values.

### Outputs

#### Matrix of constraints normals
This corresponds to A in Ax<b.

#### Vector of constraints 'offsets'
This corresponds to b in Ax<b.

#### Number of constructed constraints

#### Extimate of the maximum approximation error

#### Output flag
Shows the result of the approximation algorithm:

- 0 - The nonlinear constraint cannot become binding,
- 1 - The nonlinear constraint is infeasible,
- 2 - The approximation was successfully constructed,
- 3 - There was an error in the input branch parameters,
- 4 - There was an error in the input algorithm's options,
- 5 - There is no thermal limit for the given branch,
- 6 - Other (currently not used).

The approximation is only constructed if output flag = 2.

- - - -

## Example usage
This repository contains a number of examples of how to use the library in Matlab and C/C++:

1. Matlab:
    - `approximate_single_line.m` shows how to approximate a single line flow constraint. It contains a description of inputs and outputs of 
      mex function `LF_linearize_line`.
    - `approximate_system.m` shows how to construct the approximation of all line flow constraints in the system and put it into a format 
      that can be easily fed into a solver. It contains a description of inputs and outputs of mex function `LF_linearize_system`.
    - `construct_and_plot.m` plots the boundary surface of the actual nonlinear line flow constraint and its approximation. The actual plotting 
      is done using the provided Matlab function `plot_constraints.m`.
    - `use_in_OPF` shows how to use the constructed approximation for a given system in the AC OPF using 
      [MATPOWER](http://www.pserc.cornell.edu/matpower/). Note that the best results were observed when the AC OPF with linearized line flow
      constraints was solved with IPOPT.
2. C/C++:
    - `approximate_single_line.c` shows how to approximate a single line flow constraint. It contains a description of inputs and outputs of
      the main interface function `LF_construct` of the library. A brief description of other functions in the library is given in `line_flow.h`.

- - - -

## External packages using LINEFLOW
LINEFLOW has been integrated into a modeling library [PFNET](https://github.com/ttinoco/PFNET).

- - - -

## License
BSD 3-clause license.

- - - -

## Authors
Dmitry Shchetinin

- - - -

## Thanks
[Tomas Tinoco De Rubira](https://ttinoco.github.io/) and Ivo Caduff for their help with developing the library.

