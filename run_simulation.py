import numpy as np

from simulation_tools.planet_to_sph import convert_planetary_model_to_SPH

from analysis_tools.plot_system import systemplotter
from analysis_tools.animation import Animator

from amuse.community.fi.interface import Fi
from amuse.lab import nbody_system

import pickle

# create planet
number_of_sph_particles = 100
target_core_mass = 1 | units.MEarth
csv_file = 'earth_profile/PREM.csv' # file with earth profiles

# for now create the model with planet_to_sph,
# can only run after the profiles are complete in earth_profile/PREM.csv

model = convert_planetary_model_to_SPH(
        None,
        number_of_sph_particles,
        with_core_particle = True
        seed=12345,
        csv_file=csv_file,
        #        base_grid_options = dict(type = "glass", target_rms = 0.01),
        target_core_mass = target_core_mass
    )

core, gas_without_core, core_radius = \
        model.core_particle, model.gas_particles, model.core_radius

# make a plot with the planet and hydroparticles
systemplotter(gas_without_core + core)

# add the hydroparticles to a hydrocode and run for some time to check if it's stable

converter = nbody_system.nbody_to_si(10|units.MSun, core_radius)

hydro_code = Fi(converter, mode='openmp', redirection='none')
hydro_code.parameters.epsilon_squared = core_radius**2
hydro_code.parameters.n_smooth_tol = 0.01
hydro_code.gas_particles.add_particles(gas_without_core)
hydro_code.dm_particles.add_particle(core)

path = 'simulation_results/test_planet/'
index = 0
while (hydro_code.model_time < 1 | units.day):
    hydro_code.evolve_model(1 | units.hour)

    # save the x and y positions for plotting
    x, y = hydro_code.gas_particles.x, hydro_code.gas_particles.y
    x.append(hydro_code.dm_particles.x)
    y.append(hydro_code.dm_particles.y)
    
    with open(path+f'x_positions_{index}', 'w') as file:
        pickle.dump(x, file)

    with open(path+f'y_positions_{index}', 'w') as file:
        pickle.dump(y, file)

    index += 1


# make an animation of the system over time

animator = Animator(path, xlabel='x', ylabel='y')
animator.make_animation(show=True)


