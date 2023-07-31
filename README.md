<p align="center">
  <img src="./Logo_Unicamp__0.jpg" />
</p>

# LAMC - Linear alkanes Monte Carlo
This repository centralizes the code implemented for research during my Master's of Science in Chemical Engineering at Unicamp.

## Structure
This code is divided into folders:
1. `src`: source code files, main components for running the program, implementation of Monte Carlo algorithms as well as properties calculations and input/output routines;
2. `scripts`: `python` files used for post-processing data, running multiple simulations on parallel and other routines;
3. `data`: files used for simulating the data presented on the dissertation;
4. `notebooks`: additional `jupyter-notebook` files used for analysing, post-processing the data and creating input files.

## Code compilation and utilization
The main program is written in `C` and requires an adequate compiler (such as `gcc`). It can be compiled with the `MakeFile` available inside `src` folder:
```
cd src
make
```
This will generate several `.o` files and an executable `LAMC.exe`. For running the program it is necessary to parse the input file with extension `.inp`. 

An input file is a text file specifying the simulation parameters, with the following structure:
```
SIMULATION_ENSEMBLE NVT
EQUILIBRATION_ENSEMBLE NVT
NUBMER_OF_MOLECULES 1000
CHAIN_SIZE 1
TEMPERATURE 600.0
BOX_X_SIZE 133.24475446068766
BOX_Y_SIZE 133.24475446068766
BOX_Z_SIZE 30.0
NUMBER_TOTAL_CYCLES 60000
NUMBER_EQUILIBRATION_CYCLES 30000
ACCEPTANCE_RATIO 0.4
CYCLES_TO_CALCULATE_PROPERTIES 3
CYCLES_TO_RECORD_CONFIGURAITON 3000
CYCLES_TO_GET_PROFILES 3000
CLOSED_BOX 1
DISPLACEMENT_PROBABILITY 0.33
INSERTION_PROBABILITY 0.33
DELETION_PROBABILITY 0.33
FUGACITY 1000000.0
REFERENCE_POTENTIAL LENNARD_JONES
PERTURBED_POTENTIAL LENNARD_JONES
```

Each field is described below:
- `SIMULATION_ENSEMBLE` - ensemble for running the production steps, must be one of `NVT`, `muVT` or `NPT`;
- `EQUILIBRATION_ENSEMBLE` - ensemble for running the equilibration steps, must be one of `NVT`, `muVT` or `NPT`. More specifically, if `EQUILIBRATION_ENSEMBLE` and `SIMULATION_ENSEMBLE` different than the first half the equilibration steps are run in the equilibration ensemble and the latter half on the simulation ensemble;
- `NUBMER_OF_MOLECULES` - number of molecules for the simulation, must be an integer. Fixed number for `NVT` or `NPT` simulations or initial number for `muVT`;
- `CHAIN_SIZE` - integer specifying the size of the alkane chain;
- `TEMPERATURE` - simulation temperature in Kelvin;
- `BOX_X_SIZE` - float with the size of the `X` axis of the box;
- `BOX_Y_SIZE` - float with the size of the `Y` axis of the box;
- `BOX_Z_SIZE` - float with the size of the `Z` axis of the box;  
- `NUMBER_TOTAL_CYCLES` - integer with the number of cycles for the total simulation, during each cycle each molecule is submitted to one movement attempt (on average);
- `NUMBER_EQUILIBRATION_CYCLES` - integer with the number of cycles for the equilibration part, same definition of cycle as before;
- `ACCEPTANCE_RATIO` - target acceptance ratio for displacement attempts;
- `CYCLES_TO_CALCULATE_PROPERTIES` - integer with the number of cycles to update properties calculation;
- `CYCLES_TO_RECORD_CONFIGURAITON` - integer with the number of cycles to record configuration;
- `CYCLES_TO_GET_PROFILES` - integer with the number of cycles to record raw profile data, only used for closed box simulations;
- `CLOSED_BOX`: either 1 for closed box or 0 for open;
- `DISPLACEMENT_PROBABILITY` - probability for performing a displacement attempt on `muVT` simulation;
- `INSERTION_PROBABILITY` - probability for performing a molecule insertion attempt on `muVT` simulation;
- `DELETION_PROBABILITY` - probability for performing a molecule deletion attempt on `muVT` simulation;
- `FUGACITY` - fugacity for `muVT` simulation;
- `REFERENCE_POTENTIAL` - reference potential, either one or two values: `LENNARD_JONES` or `HARD_SPHERE`;
- `PERTURBED_POTENTIAL` - perturbation potential, either one or two values: `LENNARD_JONES` or `HARD_SPHERE`;



## References
The development of this code has been made possible by several research materials mainly, but not restricted, to:
1. Frenkel, Daan, and Berend Smit. Understanding molecular simulation: from algorithms to applications. Vol. 1. Elsevier, 2001.
2. Allen, Michael P., and Dominic J. Tildesley. Computer simulation of liquids. Oxford university press, 2017.

Additional references are enumerated on the dissertation associated with this work.
