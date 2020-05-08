# molecular-dynamics
Molecular Dynamics simulation (written in C99) of a microcanonical ensemble of ~10000 particles subject to Lennard-Jones interactions, within a cubical box with periodic boundary conditions.  It outputs numerous files containing the results of the simulation, such as the energy as a function of temperature, pressure as a function of temperature, velocity distribution of the particles, and final position of the particles.

By using this output data, one can infer the existence of a phase transition in the system from solid (at low temperatures) to liquid (at higher temperatures.  This is exhibited by a jump in the graph of energy-vs-temperature, as well as by a qualitative difference in the arrangement of the particles at the end of the simulation (a regular crystalline array, vs a random distribution).  See the included Mathematica file, which has also been printed to pdf.

One can adjust the (dimensionless) density and temperature, as well as the number of particles (see below), by editing project1.c and re-compiling.  In the future this may be fixed to allow command-line input instead.

One does not choose the number of particles directly, but instead the number of unit cells along one dimension of the box.  The number of particles is then 4 times the unit cell number cubed.  The larger the simulation, the longer the run time.

Some brief discussion of the algorithms:

We integrate the motion of the particles with Velocity Verlet, which preserves energy.  During the initial stage of the simulation, we "kick" the particles to adjust their total energy, until the system reaches equilibrium.  Reaching equilibrium is the longest part of the simulation.  After equilibrium is reached, we measure thermodynamic quantities of interest.

The Lennard-Jones potential induces pairwise forces between all the particles, which is in principle very expensive to compute.  We save time by dividing our box into "cells".  We approximate the L-J potential by cutting it off at some appropriate maximum range, thus limiting the distance at which particles can influence other particles.  The cell size is chosen to be twice the maximum range of the modified L-J potential.  Therefore, the particles within a given cell can only interact with that cell, and the 26 adjacent cells.  This saves a great deal of computational time, while still clearly showing the qualitative physics of the phase transition from solid to liquid.
