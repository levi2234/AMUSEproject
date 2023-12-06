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
from amuse.io import write_set_to_file
from amuse.datamodel import Particles

import os
import pickle
import matplotlib.pyplot as plt
import numpy as np

#-----------------------INPUT PARAMETERS-----------------------
# characteristics of the planet
number_of_sph_particles = 1000
target_core_mass = 12 | units.MEarth # 0.5 | units.MJupiter
pickle_file = 'simulation_tools/profiles/jupiter_like_planet_structure.pkl' # insert planet profile name from mesa

# characteristics of the moon
# moon_mass = 0.008 | units.MEarth # mass Europa, previously: 0.012*planet_mass
moon_mass = 0.015 | units.MEarth # mass Io
# distance_planet_moon = 671000 | units.km # distance Europa from jupiter
distance_planet_moon = 421600 | units.km # distance Io and Jupiter
eccentricity_moon = 0 # 0.009 for Europa, 0.004 for Io, but close enough to 0.

# explosion
outer_fraction = 0.1
explosion_energy = 1.0e+46|units.erg

# model evolution
timestep = 6 | units.hour
simulation_duration = 5 | units.day

# ----------------------CREATE THE PLANET----------------------

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



#---------------------Create a Moon --------------------
#calcualte the planets total mass
planet_mass = core.mass.sum() + gas_without_core.mass.sum()

#create binary system with the planet and the moon from orbital elements

system = new_binary_from_orbital_elements(planet_mass, moon_mass, distance_planet_moon, eccentricity_moon, G = constants.G)

#move most massive object to 0,0,0
system.position = system.position -system[0].position

#add moon to own particle set (could be done better but this only works for 2 objects)
moon = Particles(1)
moon.mass = moon_mass
moon.position = system[1].position
moon.velocity = system[1].velocity
print(moon.velocity.in_(units.kms))

#----------------------BEFORE PLOT-----------------------
# plt.scatter(gas_without_core.x.value_in(units.AU), gas_without_core.y.value_in(units.AU),c="blue", marker='o')
# plt.scatter(core.x.value_in(units.AU), core.y.value_in(units.AU),c="r", marker='o')
# plt.savefig('simulation_results/planetbefore.png')

#----------------------PHYSICS SETUP----------------------

# add the hydroparticles to a hydrocode and run for some time to check if it's stable

converter = nbody_system.nbody_to_si(1|units.MSun, core_radius)

#setup hydro code
hydro_code = Fi(converter, mode='openmp', redirection='none')
hydro_code.parameters.epsilon_squared = core_radius**2
hydro_code.parameters.n_smooth_tol = 0.01
hydro_code.gas_particles.add_particles(gas_without_core)
hydro_code.dm_particles.add_particle(core)

#setup gravity code
gravity_code = BHTree(converter)
gravity_code.parameters.epsilon_squared = core_radius**2
gravity_code.particles.add_particles(moon)

#setup the bridge
bridge = Bridge(use_threading=False)
bridge.add_system(gravity_code, (hydro_code,)) #makes sure planet/atmosphere particles are affected by gravity
bridge.add_system(hydro_code, (gravity_code,)) #makes sure moon is also affected by hydrodynamics and stays in orbit

# Set the timestep for the bridge
bridge.timestep = timestep



path_results = 'simulation_results/jupiterlike_planet/'
if not os.path.exists(path_results):
    os.mkdir(path_results)

index = 0
triggered_injection = False

#----------------------BEFORE PLOT-----------------------
systemplotter(hydro_code.gas_particles, xlabel='x', ylabel='y').plot(save="simulation_results/planetbefore.png", c="b",s=2, close=False)
systemplotter(gravity_code.particles, xlabel='x', ylabel='y').plot(save="simulation_results/planetbefore.png", c="r",s=40, close=False)
systemplotter(hydro_code.dm_particles, xlabel='x', ylabel='y').plot(save="simulation_results/testplanetBefore.png", c="r",s=20)


#----------------------RUN THE SIMULATION----------------------
while (hydro_code.model_time < simulation_duration):
    #hydro_code.evolve_model(hydro_code.model_time + bridge.timestep)
    bridge.evolve_model(hydro_code.model_time + bridge.timestep)

    if (hydro_code.model_time.value_in(units.hour) >= 24) & (triggered_injection == False): #do something at the 6th timestep( 6 hours in)   
        #inject energy
        triggered_injection = True
        print("injecting energy")
        planet_radius = max([np.sqrt(pos.length_squared()) for pos in gas_without_core.position])
        
        explosion_radius = planet_radius - (planet_radius * outer_fraction)
        print(explosion_radius.in_(units.RJupiter))

        inject_energy.inject_explosion_energy(hydro_code.gas_particles,explosion_energy=explosion_energy, exploding_region=explosion_radius)
        
    
    write_set_to_file(hydro_code.gas_particles, path_results+f'gas_particles_{index}.hdf5', overwrite_file=True)
    write_set_to_file(hydro_code.dm_particles, path_results+f'dm_particles_{index}.hdf5', overwrite_file=True)
    write_set_to_file(gravity_code.particles, path_results+f'gravity_particles_{index}.hdf5', overwrite_file=True)

    index += 1



#----------------------AFTER PLOT-----------------------
# systemplotter(hydro_code.gas_particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="b",s=2, close=False)
# systemplotter(gravity_code.particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="r",s=40, close=False)
# systemplotter(hydro_code.dm_particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="r",s=20)


#----------------------STOP THE SIMULATION-----------------------
hydro_code.stop()
gravity_code.stop()


path = 'simulation_results/jupiterlike_planet/'
animator = Animator(path, xlabel='x', ylabel='y', xlim=0.01, ylim=0.01)
animator.make_animation(save_path='simulation_results/animation_jup.mp4')



