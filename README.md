# md
The program mdeval.py / main_mdtraj.py reads and analyzes lammps-data from molecular dynamics simulations of single-particle impingements of low-temperature plasma particles onto a metal surface.

With this program the particle trajectories are read, then the collision events with the surface are counted and analyzed.
From there, particle states are defined and their respective populations as well as transition fluxes are calculated.
Then the transition rates are computed and used to construct the analytical solution of the rate equation model for particle popuations. This enables us to make long-time predictions for the adsorbed particle number to gain a better understanding of the behavior of plasmas near a surface.
Also, the energy distributions for adatoms are calculated and the average residence time can be obtained from the analytical solution.

The theoretical foundation is based on the work by A. Filinov, M. Bonitz, and D. Loffhagen (https://doi.org/10.1088/1361-6595/aac61e, https://doi.org/10.1088/1361-6595/aac620).
