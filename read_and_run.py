from analysis_tools.plot_system import systemplotter
from analysis_tools.animation import Animator

from simulation_tools import inject_energy

from amuse.ext.star_to_sph import (pickle_stellar_model, convert_stellar_model_to_SPH,)
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree
from amuse.couple.bridge import Bridge
from amuse.lab import nbody_system
from amuse.units import units
from amuse.units import constants
from amuse.io import write_set_to_file, read_set_from_file
from amuse.datamodel import Particles

import os
import pickle
import matplotlib.pyplot as plt
import numpy as np

#-----------------------INPUT PARAMETERS-----------------------
# characteristics of the planet
number_of_sph_particles = 30000
target_core_mass = 12 | units.MEarth # 0.5 | units.MJupiter
pickle_file = 'simulation_tools/profiles/jupiter_like_planet_structure.pkl' # insert planet profile name from mesa

# characteristics of the moon
# # moon_mass = 0.008 | units.MEarth # mass Europa, previously: 0.012*planet_mass
# moon_mass = 5 | units.MEarth # 10 x mass Io
# # distance_planet_moon = 671000 | units.km # distance Europa from jupiter
# distance_planet_moon = 2 * 621600 | units.km # distance Io and Jupiter
# eccentricity_moon = 0 # 0.009 for Europa, 0.004 for Io, but close enough to 0.

# # explosion
# outer_fraction = 0.3
# explosion_energy = 2.0e+42|units.erg

# model evolution
timestep = 6 | units.hour
simulation_duration = 8 | units.day

# -------- PREVIOUS MODEL --------------------------------------

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

print(core)
# -----------------READ IN PARTICLES----------------------

counter = 0
path_results = 'simulation_results/jupiterlike_planet/'
for path in os.listdir(path_results):
    if os.path.isfile(os.path.join(path_results, path)):
        counter += 1

counter = counter//3 # 3 different kind of files in directory
gas_particles = read_set_from_file(path_results+f'gas_particles_{counter-1}.hdf5')
dm_particle = read_set_from_file(path_results+f'dm_particles_{counter-1}.hdf5')
gravity_particle = read_set_from_file(path_results+f'gravity_particles_{counter-1}.hdf5')


# -----------------PHYSICS SETUP----------------------

# add the hydroparticles to a hydrocode and run for some time to check if it's stable

converter = nbody_system.nbody_to_si(1|units.MSun, core_radius)

#setup hydro code
hydro_code = Fi(converter, mode='openmp', redirection='none')
hydro_code.parameters.epsilon_squared = core_radius**2
hydro_code.parameters.n_smooth_tol = 0.01
hydro_code.gas_particles.add_particles(gas_particles)
hydro_code.dm_particles.add_particles(dm_particle)

#setup gravity code
gravity_code = BHTree(converter)
gravity_code.parameters.epsilon_squared = core_radius**2
gravity_code.particles.add_particles(gravity_particle)

#setup the bridge
bridge = Bridge(use_threading=False)
bridge.add_system(gravity_code, (hydro_code,)) #makes sure planet/atmosphere particles are affected by gravity
bridge.add_system(hydro_code, (gravity_code,)) #makes sure moon is also affected by hydrodynamics and stays in orbit

# Set the timestep for the bridge
bridge.timestep = timestep

path_results = 'simulation_results/jupiterlike_planet/'
if not os.path.exists(path_results):
    os.mkdir(path_results)


#----------------------BEFORE PLOT-----------------------
# systemplotter(hydro_code.gas_particles, xlabel='x', ylabel='y').plot(save="simulation_results/planetbefore.png", c="b",s=2, close=False)
# systemplotter(gravity_code.particles, xlabel='x', ylabel='y').plot(save="simulation_results/planetbefore.png", c="r",s=40, close=False)
# systemplotter(hydro_code.dm_particles, xlabel='x', ylabel='y').plot(save="simulation_results/testplanetBefore.png", c="r",s=20)


#----------------------RUN THE SIMULATION----------------------
initial_energy = gravity_code.get_total_energy() + hydro_code.get_total_energy()
while (hydro_code.model_time < simulation_duration):
    #hydro_code.evolve_model(hydro_code.model_time + bridge.timestep)
    bridge.evolve_model(hydro_code.model_time + bridge.timestep)
    
    dE = initial_energy/(gravity_code.get_total_energy() + hydro_code.get_total_energy())
    print('Time=', hydro_code.model_time.in_(units.hour))
    print('dE=', dE)
    
    write_set_to_file(hydro_code.gas_particles, path_results+f'gas_particles_{counter}.hdf5', overwrite_file=True)
    write_set_to_file(hydro_code.dm_particles, path_results+f'dm_particles_{counter}.hdf5', overwrite_file=True)
    write_set_to_file(gravity_code.particles, path_results+f'gravity_particles_{counter}.hdf5', overwrite_file=True)

    counter += 1



#----------------------AFTER PLOT-----------------------
# systemplotter(hydro_code.gas_particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="b",s=2, close=False)
# systemplotter(gravity_code.particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="r",s=40, close=False)
# systemplotter(hydro_code.dm_particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="r",s=20)


#----------------------STOP THE SIMULATION-----------------------
hydro_code.stop()
gravity_code.stop()


path = 'simulation_results/jupiterlike_planet/'
animator = Animator(path, xlabel='x', ylabel='y', xlim=0.03, ylim=0.03)
animator.make_animation(save_path='simulation_results/animation_jup.mp4')



