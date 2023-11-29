import numpy as np

from analysis_tools.plot_system import systemplotter
from analysis_tools.animation import Animator
import matplotlib.pyplot as plt
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


import pickle

# ----------------------CREATE THE PLANET----------------------
number_of_sph_particles = 50000
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


#---------------------CREATE MOON--------------------
#System Parameters
planet_mass = core.mass.sum() + gas_without_core.mass.sum()
moon_mass = 0.012*planet_mass
semimajor_axis = 5.2 | units.AU
eccentricity = 0

#Create Binary
system = new_binary_from_orbital_elements(planet_mass,moon_mass,semimajor_axis, eccentricity, G = constants.G)


#Translate system to geocentric frame
system.position = system.position - system[0].position

#add moon to own particle set (could be done better but this only works for 2 objects)
moon = Particles(1)
moon.mass = moon_mass
moon.position = system[1].position
moon.velocity = system[1].velocity



#----------------------PHYSICS SETUP----------------------

#setup converter
converter = nbody_system.nbody_to_si(10|units.MSun, core_radius)

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
bridge.timestep = 1 | units.hour



path = 'simulation_results/test_planet/'
index = 0
triggered_injection = False

#----------------------BEFORE PLOT-----------------------
systemplotter(hydro_code.gas_particles, xlabel='x', ylabel='y').plot(save="simulation_results/planetbefore.png", c="b",s=2, close=False)
systemplotter(gravity_code.particles, xlabel='x', ylabel='y').plot(save="simulation_results/planetbefore.png", c="r",s=40, close=False)
systemplotter(hydro_code.dm_particles, xlabel='x', ylabel='y').plot(save="simulation_results/testplanetBefore.png", c="r",s=20)


#----------------------RUN THE SIMULATION----------------------
while (hydro_code.model_time < 10 | units.day):
    
    bridge.evolve_model(hydro_code.model_time + bridge.timestep)

    if (hydro_code.model_time.value_in(units.hour) >= 9) & (triggered_injection == False): #do something at the 6th timestep( 6 hours in)   
        #inject energy
        triggered_injection = True
        print("injecting energy")
        inject_energy.inject_explosion_energy(hydro_code.gas_particles,explosion_energy=2e50|units.erg,exploding_region=10|units.RSun)
        
    
    write_set_to_file(hydro_code.gas_particles, path+f'gas_particles_{index}.hdf5', overwrite_file=True)
    write_set_to_file(hydro_code.dm_particles, path+f'dm_particles_{index}.hdf5', overwrite_file=True)
    write_set_to_file(gravity_code.particles, path+f'gravity_particles_{index}.hdf5', overwrite_file=True)

    index += 1



#----------------------AFTER PLOT-----------------------
#this code needs to run before the simulation is stopped
systemplotter(hydro_code.gas_particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="b",s=2, close=False)
systemplotter(gravity_code.particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="r",s=40, close=False)
systemplotter(hydro_code.dm_particles, xlabel='x', ylabel='y').plot(save="simulation_results/testplanetAfter.png", c="r",s=20)


#----------------------STOP THE SIMULATION-----------------------
hydro_code.stop()
gravity_code.stop()

# #----------------------ANIMATION-----------------------
# animator = Animator(path, xlabel='x', ylabel='y')
# animator.make_animation(save_path='simulation_results/testplanetanim.mp4')



