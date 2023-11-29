import numpy as np

from analysis_tools.plot_system import systemplotter
from analysis_tools.animation import Animator
import matplotlib.pyplot as plt

from amuse.ext.star_to_sph import (pickle_stellar_model, convert_stellar_model_to_SPH,)
from amuse.community.fi.interface import Fi
from amuse.lab import nbody_system
from amuse.units import units
from amuse.io import write_set_to_file
from amuse.datamodel import Particles

import pickle

# create planet
number_of_sph_particles = 1000
target_core_mass = 5 | units.MSun # is this a possible unit?
pickle_file = 'profiles/super_giant_stellar_structure.pkl' # insert planet profile name from mesa

# for now create the model with planet_to_sph,
# can only run after the profiles are complete in earth_profile/PREM.csv

model = convert_stellar_model_to_SPH(
        None,
        number_of_sph_particles,
        with_core_particle = True,
        seed=12345,
        pickle_file=pickle_file,
        #        base_grid_options = dict(type = "glass", target_rms = 0.01),
        target_core_mass = target_core_mass
    )

core, gas_without_core, core_radius = \
        model.core_particle, model.gas_particles, model.core_radius

# make a plot with the planet and hydroparticles

plt.scatter(gas_without_core.x.value_in(units.AU), gas_without_core.y.value_in(units.AU),c="blue", marker='o')
plt.scatter(core.x.value_in(units.AU), core.y.value_in(units.AU),c="r", marker='o')
plt.savefig('simulation_results/planetbefore.png')

#systemplotter(merged_dataset).plot(save=True)

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
    hydro_code.evolve_model(hydro_code.model_time + (1 | units.hour))
    
    # save the gas and dm particles for plotting
    # x, y = hydro_code.gas_particles.x, hydro_code.gas_particles.y
    # x.append(hydro_code.dm_particles.x)
    # y.append(hydro_code.dm_particles.y)
    write_set_to_file(hydro_code.gas_particles, path+f'gas_particles_{index}.hdf5', overwrite_file=True)
    write_set_to_file(hydro_code.dm_particles, path+f'dm_particles_{index}.hdf5', overwrite_file=True)

    index += 1


systemplotter(hydro_code.gas_particles, xlabel='x', ylabel='y').plot(save=True, c="blue",s=2)
# plt.scatter(hydro_code.gas_particles.x, hydro_code.gas_particles.y,c="blue", marker='o')
# plt.scatter(hydro_code.dm_particles.x, hydro_code.dm_particles.y,c="r", marker='o')
# plt.savefig('simulation_results/planetafter.png')
# make an animation of the system over time
hydro_code.stop()
animator = Animator(path, xlabel='x', ylabel='y')
animator.make_animation(save_path='simulation_results/animation.mp4')

