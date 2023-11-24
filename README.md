# qtip4pf_water
Showcasing some of the code written to run classical and quantum mechanical water simulations during my DPhil in C++.
The main purpose of this code was to provide a general framework for setting up and running these water simulations to be used by other members of my research group to be used to generate and test "fast" centroid potentials of mean force for CMD calculations and to test them.
The actual code provided here does not calculate much in the way of simulation data, only the average kinetic and potential energies in order to check the simulation is running reasonably.
This is because the radial distribution functions used to calculate these potentials of mean force were written more generally by the other research group members and then interfaced with this code.

The code was not necessary designed to be read by people who are not me or the members of the research group who I was collaborating with, and as such may not be particularly easy to follow - especially for someone without prior experience of path integral simulations.
A brief overview of the code is as follows:

The main function is contained within water.cpp, and deals with reading in data from the input file (here config.txt) and then determining how the simulation needs to be run from this.

The simulation is then run inside simulation.cpp, which uses update.cpp to propagate the trajectories.

Force calculations can be found inside forces.cpp

Thermostats for NVT simulations are contained in thermostat.cpp

Other files are used such as initialise.cpp to generate the starting positions, neighbour.cpp to calculate verlet neighbour lists, path_integrals.cpp to calculate useful quantities from path integral simulations and fourier.cpp to perform the fast fourier transforms used in path integral simulations

Parameters used throughout the code can be found inside qtip4pfparams.h, and helper functions used to intuitively index arrays inside indexing.h
