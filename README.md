## PyrticleSim :8ball: ##

Alexandros Papagiannakis, Stanford University, 2025

This repository includes two classes that can be used to simulate freely diffusing particles, or particle trajectories within rod-shaped cells, with and without nucleoid exclusion.
A Jupyter notebook includes an example where the two classes are implemented to simulate particle trajectories.

'Simulations_of_free_particles.py' script includes the "unbound_particle_fast_simulation" class which takes the frame interval and the number of particle positions as initial input.
Once the class is initialized, the use can use the "run_simulation" attribute to perform simulations for a specified number of particles, sampling from a specified diffusion_coefficient_array, with a specific pseudo-random seed and the option to show or hide the simulated traces. 


![Particle simulation example](https://github.com/alexSysBio/PyrticleSim/blob/main/Example_maps_movie_nucleoid_confinement_short.gif)
