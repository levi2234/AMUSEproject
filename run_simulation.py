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
number_of_sph_particles = 1500
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


#---------------------Create a Moon --------------------
#calcualte the planets total mass
planet_mass = core.mass.sum() + gas_without_core.mass.sum()
moon_mass = 0.012*planet_mass

#create binary system with the planet and the moon from orbital elements
system = new_binary_from_orbital_elements(planet_mass,moon_mass,4 | units.AU, 0, G = constants.G)
#move most massive object to 0,0,0
system.position = system.position -system[0].position

#add moon to own particle set (could be done better but this only works for 2 objects)
moon = Particles(1)
moon.mass = moon_mass
moon.position = system[1].position
moon.velocity = system[1].velocity

#----------------------BEFORE PLOT-----------------------
plt.scatter(gas_without_core.x.value_in(units.AU), gas_without_core.y.value_in(units.AU),c="blue", marker='o')
plt.scatter(core.x.value_in(units.AU), core.y.value_in(units.AU),c="r", marker='o')
plt.savefig('simulation_results/planetbefore.png')

#----------------------PHYSICS SETUP----------------------

# add the hydroparticles to a hydrocode and run for some time to check if it's stable

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
# Create the bridge
bridge = Bridge(use_threading=False)
bridge.add_system(gravity_code, (hydro_code,))
bridge.add_system(hydro_code, (gravity_code,))

# Set the timestep for the bridge
bridge.timestep = 1 | units.hour



path = 'simulation_results/test_planet/'
index = 0
triggered_injection = False


#----------------------RUN THE SIMULATION----------------------
while (hydro_code.model_time < 4 | units.day):
    #hydro_code.evolve_model(hydro_code.model_time + bridge.timestep)
    bridge.evolve_model(hydro_code.model_time + bridge.timestep)

    if (hydro_code.model_time.value_in(units.hour) >= 6) & (triggered_injection == False): #do something at the 6th timestep( 6 hours in)   
        #inject energy
        triggered_injection = True
        print("injecting energy")
        inject_energy.inject_explosion_energy(hydro_code.gas_particles,explosion_energy=1.0e+50|units.erg,exploding_region=10|units.RSun)
        
    
    write_set_to_file(hydro_code.gas_particles, path+f'gas_particles_{index}.hdf5', overwrite_file=True)
    write_set_to_file(hydro_code.dm_particles, path+f'dm_particles_{index}.hdf5', overwrite_file=True)
    write_set_to_file(gravity_code.particles, path+f'gravity_particles_{index}.hdf5', overwrite_file=True)

    index += 1



#----------------------AFTER PLOT-----------------------
systemplotter(hydro_code.gas_particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="b",s=2, close=False)
systemplotter(gravity_code.particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="r",s=40, close=False)
systemplotter(hydro_code.dm_particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="r",s=20)

hydro_code.stop()
gravity_code.stop()
animator = Animator(path, xlabel='x', ylabel='y',xlim=6,ylim=6)
animator.make_animation(save_path='simulation_results/animation.mp4')



