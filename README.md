## PyrticleSim :8ball: ##

Alexandros Papagiannakis, Stanford University, 2025

This repository includes two classes that can be used to simulate freely diffusing particles, or particle trajectories within rod-shaped cells, with and without nucleoid exclusion.
A Jupyter notebook includes an example where the two classes are implemented to simulate particle trajectories.

"Simulations_of_free_particles.py" script includes the "unbound_particle_fast_simulation" class which takes the frame interval and the number of particle positions as initial input.
Once the class is initialized, the use can use the "run_simulation" attribute to perform simulations for: 
1. a specified number of particles
2. sampling from a specified diffusion coefficient array
3. with a specific pseudo-random seed
4. and the option to show or hide the simulated traces. 

"Simulations_of_confined_particles.py" script includes the "bound_particle_simulation" class which takes the particle radius and the number of particle positions as initial input.
Once the class is initialized, the use can use the "run_simulation" function to perform simulations for:
1. a specified number of particles
2. a specified cell length
3. a specified cell width
4.  specified nuleoid width
5.  a specified nuleoid length
6.  a specified number of time steps
7.  sampling from a specified diffision coefficient array
8.  for a specified time lag
9.  for a specific particle radius
10.  initializing the simulations from a speficic pseudo-random seed
11.  selecting a binary input for nucleoid confinement
12.  showing the cell-bound particle trajectories for a specific interval (e.g., every 100 traces)
13.  with the possibility to save images of the particle traces at a specified path.


![Particle simulation example](https://github.com/alexSysBio/PyrticleSim/blob/main/Example_maps_movie_nucleoid_confinement_short.gif)
