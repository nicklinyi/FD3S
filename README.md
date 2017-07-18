This directory contains the source code for high-order finite difference
3D elastic wave propagation on staggered grids for spherical geometry.

The serial algorithm has been written to run on the
Pacific Sierra F90 to F77 translator.


### COMPILATION

To compile the code, type
```
$ make main.exe
```

to run the `makefile`.

### STRUCTURE AND COMPONENTS OF THE ALGORITHM

The makefile connects the different subfiles of the code.

The declaration of global parameters and variables is accomplished by

`common_global.h`

The main program is called `fd3s_main.f90` including all other subroutines
by call statements.

At first, parameters are read in from the file `Par` in the
routine `fd3s_input.f90`.

`fd3s_init.f90` initializes all further values such as model space,
source and receiver parameters, absorbing boundaries, and FD operator
specifications.

The subroutine `fd3s_model.f90` computes the background model as chosen from:
- Homogeneous half space
- Two-layered model (discontinuity at 240km depth)
- isotropic part of PREM (Dziewonski and Anderson, 1981)
- Slab model read in from files for each node as perturbation values for PREM

`fd3s_check.f90` writes all important setup features interactively and to
log-files, including a file designed to be loaded by Matlab.

`fd3s_evolution.f90` solves the velocity-stress formalism of the elastodynamic
equations in spherical geometry by finite difference derivatives of
variable operator length and second-order time extrapolation.
Two different surface treatments are involved as well as
switches for elastic interpolation, and computation of curl and divergence
for the velocity field.  
Interpolation and partial derivatives are explicitely calculated in
subroutine `fd3s_oper.f90`.

`fd3s_output` saves seismograms of the velocity components and/or curl and
divergence to file. Additionally, snapshots of the velocity field, stress
tensor trace, curl and divergenve of velocities can be written out
on the surface and two vertical planes, usually through the hypocenter.

It is recommended to insert delta functions for seismogram simulations
to be able to alter the frequencies afterwards. In this case, snapshots
are automatically excluded.


### Parameter control

Parameter settings are undertaken within the file `Par`, together with
`params.h` for amount of grid points.

The amount of grid points should be derived together with the range of
model space and number of nodes used by applying Equation (2.26) of
the Diplom Thesis.

After defining a background model, the time increment can be calculated from
the numerical stability criterion for spherical geometry (Equ. (2.25)) which
is also written out while running the program. The highest value for the
stability factor (usually from the phi-axis of the lowermost node)
should then act as a constraint on the time increment.
From that value, the total amount of time steps can be chosen according
to desired total length of wave propagation.
Values for the approximate frequency resolution when using a delta pulse
are written out as well.

It is therefore recommended to run several simulations of one time step
to account for the best fit between time increment and frequency resolution.

Changing parameters of `params.h` requires a new compilation, while
changing the `Par` file can be done upon the same executable.
Especially due to the times spent on compiling on the Hitachi, the values
of `params.h` should be chosen at first.


### Running the Program

The compilation produces the executable `main.exe`

To run this file type

```
$ make test
```

### File storage and additional folders

Object files of the compilation process are stored in a directory
called `OBJECTS` as seen in the makefile.

In the case of the specific Par file shown here, the directory containing the
output data is called `DATA` as seen in the definition of the parameter
`seisfile`.

The log files are written into the folder `output` as seen from the
parameter `outfile` in `fd3s_main.f90`.

All those folders must obviously be created if these settings are used.

The output files of seismograms are `filename0_x`, `filename0_y`, `filename0_z`.



This algorithm has been developed during the Diploma Thesis of
Tarje Nissen-Meyer on the basis of a code written by Prof. Heiner Igel,
both at the Institute of Geophysics of the Ludwig-Maximilians-University
Munich, Germany until August 2001.

Copyright by Tarje Nissen-Meyer, 2002

Further information: tarje@geophysik.uni-muenchen.de

Last change: March 9th 2002

### What I have done
- rewrite the makefiles in root directory and `MODELS` directory.
- edit the `model_parameters` in `MODELS` directory.
- modify the `README.md` to suit changes made by me.

### Notification
**<font color=red>You must run this code with >= 13 CPU cores, or you will get a run-time error</font>**.   
