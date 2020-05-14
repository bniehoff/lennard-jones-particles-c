# molecular-dynamics
Molecular Dynamics simulation of a microcanonical ensemble of ~10000 particles subject to pairwise Lennard-Jones interactions, within a cubical box with periodic boundary conditions.  See the Jupyter notebook Molecular-Dynamics.ipynb for analysis the simulation and the output it produces.

The main program is written in C under md-simulate.c, md-simulate.h.  The usage (with defaults) is

md-simulate --cellcount 5 --density 0.8 --temperature 0.9 --output-directory "data" --prefix "rho_0.8/T_0.9"

The numbers given are just some reasonable values that will run quickly.  Change them to whatever you like.  Density and temperature are dimensionless, which means physical constants (like Boltzmann's) have already been absorbed.  "Cellcount" refers to the number of face-centered-cubic crystal cells, along one side of the cubical box, which are used to set up the initial positions of the atoms.  Thus the total number of atoms is

Number of atoms = 4 * (cellcount)^3

since there are 4 atoms in a single unit cell of an FCC lattice.  The output-directory must already exist.  The prefix is a sub-directory of output-directory where the data files of this run will be stored, and it must *also* already exist!  The C program does no checking at all of these conditions.

There is a Python wrapper which runs md-simulate, ensures that the appropriate subdirectories exist, and facilitates sweeping over the parameter range in density and temperature.  Its usage, with defaults, is

python3 batch-simulate.py --cellcount 5 --density-linspace "[0.8, 0.8, 1]" --temperature-linspace "[0.1, 0.9, 5]" --output-directory "data"

There is an optional flag --remove-data which will delete all files in output-directory.  The --density-linspace and --temperature-linspace parameters are in the format of arguments which would be passed to numpy.linspace().  For example "[0.1, 0.9, 5]" means to take 5 samples evenly spaced from 0.1 to 0.9.  The Python script automatically runs md-simulate for each value of density and temperature in the combined sample space, storing the results of individual runs under subdirectories of output-directory.

Files produced:

In output-directory: the file thermo_measurements.csv contains the energy, heat capacity (at constant volume), and pressure, as a function of density and temperature.  The data in this file should be accumulated over many runs.

In output-directory/prefix:  Three files:

final_state.csv:  The positions, velocities, and speeds of all particles at the end of simulation.

summary_info.csv:  Various global parameters of the simulation, such as the size of the box and the number of atoms.

time_series.csv:  The temperature, potential energy, total energy, and mean squared displacement of the atoms as a function of time, over the entire simulation.

By using this output data, one can infer the existence of a phase transition in the system from solid (at low temperatures) to liquid (at higher temperatures.  This is exhibited by a jump in the graph of energy-vs-temperature, as well as by a qualitative difference in the arrangement of the particles at the end of the simulation (a regular crystalline array, vs a random distribution).  See the included Mathematica file, which has also been printed to pdf.  (Will create a Jupyter notebook soon that does the same thing, but with more visualizations.)

Some brief discussion of the algorithms:

We integrate the motion of the particles with Velocity Verlet, which preserves energy.  During the initial stage of the simulation, we "kick" the particles to adjust their total energy, until the system reaches equilibrium.  Reaching equilibrium is the longest part of the simulation.  After equilibrium is reached, we measure thermodynamic quantities of interest.

The Lennard-Jones potential induces pairwise forces between all the particles, which is in principle very expensive to compute.  We save time by dividing our box into "cells".  We approximate the L-J potential by cutting it off at some appropriate maximum range, thus limiting the distance at which particles can influence other particles.  The cell size is chosen to be twice the maximum range of the modified L-J potential.  Therefore, the particles within a given cell can only interact with that cell, and the 26 adjacent cells.  This saves a great deal of computational time, while still clearly showing the qualitative physics of the phase transition from solid to liquid.
