from analysis_tools.plot_system import systemplotter
from analysis_tools.animation import Animator

from simulation_tools import create_planet_and_atmosphere as CPA
from simulation_tools.profiles import atmospheric_profiles as AP
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
import sys
import pickle
import matplotlib.pyplot as plt
import numpy as np

# model evolution
timestep = 1 | units.hour
simulation_duration = 10 | units.day

# ----------------------CREATE THE PLANET----------------------

# for now create the model with planet_to_sph,
# can only run after the profiles are complete in earth_profile/PREM.csv

planet_mass = 1 | units.MEarth
planet_radius = 1 | units.REarth
atmosphere_profile = AP.true_profile #wether this works is dependant on your system and the initialization of the particles. If it doesn't work, try the other profile such as AP.barometric_profile
atmosphere_height = 100|units.km
atmosphere_density = 1.225 | units.kg * units.m** (-3)
num_points = 50000
gamma = 1.4

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 2:
    print("Usage: python my_script.py integer_argument")
    sys.exit(1)

# Extract the integer argument from the command line
try:
    arg1 = int(sys.argv[1])
except ValueError:
    print("Error: Argument must be an integer.")
    sys.exit(1)
atmosphere_speed = int(sys.argv[1])
print("Speed:", atmosphere_speed)

planet, atmosphere = CPA.create_planet_and_atmosphere(planet_mass=planet_mass , 
                                                planet_radius= planet_radius, 
                                                atmosphere_profile=atmosphere_profile, 
                                                atmosphere_height= atmosphere_height, 
                                                atmosphere_density=atmosphere_density, 
                                                num_points=num_points,
                                                gamma=gamma)


#adding radial velocities to atmosphere
#get x,y,z positions and normalize
x = atmosphere.x.value_in(units.m)
y = atmosphere.y.value_in(units.m)
z = atmosphere.z.value_in(units.m)
r = np.sqrt(x**2 + y**2 + z**2)
x = x/r
y = y/r
z = z/r

#set velocities

vx = x * atmosphere_speed + np.random.uniform(-10,10) | units.m / units.s
vy = y * atmosphere_speed + np.random.uniform(-10,10)| units.m / units.s
vz = z * atmosphere_speed + np.random.uniform(-10,10)| units.m / units.s
atmosphere.vx = vx 
atmosphere.vy = vy
atmosphere.vz = vz


#---------------------Create a Moon --------------------
#calcualte the planets total mass
planet_mass = planet.mass.sum()+ atmosphere.mass.sum()
print('planet mass:', planet_mass.in_(units.MEarth))

# moon_mass = 0.008 | units.MEarth # mass Europa, previously: 0.012*planet_mass
moon_mass = 0.012 | units.MEarth #
# distance_planet_moon = 671000 | units.km # distance Europa from jupiter
distance_planet_moon = 621600 | units.km # distance Io and Jupiter
eccentricity_moon = 0 # 0.009 for Europa, 0.004 for Io, but close enough to 0.
#create binary system with the planet and the moon from orbital elements

system = new_binary_from_orbital_elements(planet_mass, moon_mass, distance_planet_moon, eccentricity_moon, G = constants.G)

#move most massive object to 0,0,0
system.position = system.position -system[0].position

#add moon to own particle set (could be done better but this only works for 2 objects)
moon = Particles(1)
moon.mass = moon_mass
moon.position = system[1].position
moon.velocity = system[1].velocity


#add sun from geocentric frame
system = new_binary_from_orbital_elements(planet_mass, 1 | units.MSun, 1| units.AU, 0, G = constants.G)
#center around earth
system.position = system.position -system[0].position
#add sun to own particle set
sun = Particles(1)
sun.mass = 1 | units.MSun
sun.position = system[1].position
sun.velocity = system[1].velocity

# #----------------------BEFORE PLOT-----------------------
fig, ax = plt.subplots()
planet_radius_temp = planet.radius.value_in(units.m)
circle = plt.Circle((0, 0), planet_radius_temp[0], color='green',alpha= 0.2, fill=True)
ax.scatter(atmosphere.x.value_in(units.m), atmosphere.y.value_in(units.m), s=0.01)
ax.add_artist(circle)
ax.set_aspect('equal')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('Planet and Atmosphere before simulation \n (Earth Atmosphere Distribution)')
plt.savefig('simulation_results/planetbefore.png')

#----------------------PHYSICS SETUP----------------------

# add the hydroparticles to a hydrocode and run for some time to check if it's stable

converter = nbody_system.nbody_to_si(1e24 | units.kg, 6.371e6 | units.m)

#setup hydro code
hydro_code = Fi(converter, mode='openmp', redirection='none')
hydro_code.parameters.epsilon_squared = planet_radius**2
hydro_code.parameters.n_smooth_tol = 0.01
hydro_code.gas_particles.add_particles(atmosphere)
hydro_code.dm_particles.add_particle(planet)

#setup gravity code
gravity_code = BHTree(converter)
gravity_code.parameters.epsilon_squared = planet_radius**2
gravity_code.particles.add_particles(moon)
gravity_code.particles.add_particles(sun)

#setup the bridge
bridge = Bridge(use_threading=False)
bridge.add_system(gravity_code, (hydro_code,)) #makes sure planet/atmosphere particles are affected by gravity
bridge.add_system(hydro_code, (gravity_code,)) #makes sure moon is also affected by hydrodynamics and stays in orbit

# Set the timestep for the bridge
bridge.timestep = timestep



path_results = f'simulation_results/Earth_planet_{atmosphere_speed}ms_long/'
if not os.path.exists(path_results):
    os.mkdir(path_results)

index = 0
triggered_injection = False

# #----------------------BEFORE PLOT-----------------------
# systemplotter(hydro_code.gas_particles, xlabel='x', ylabel='y').plot(c="b",s=0.01, close=False)
# systemplotter(gravity_code.particles, xlabel='x', ylabel='y').plot( c="r",s=40, close=False)
# systemplotter(hydro_code.dm_particles, xlabel='x', ylabel='y').plot(save="simulation_results/testplanetBefore.png", c="r",s=20)

try:
    previous_time = None
    #----------------------RUN THE SIMULATION----------------------
    while (hydro_code.model_time < simulation_duration):
        if previous_time == None:
            previous_time = hydro_code.model_time.value_in(units.hour)
        #hydro_code.evolve_model(hydro_code.model_time + bridge.timestep)
        bridge.evolve_model(hydro_code.model_time + bridge.timestep)
        if previous_time != hydro_code.model_time.value_in(units.hour):
            previous_time = hydro_code.model_time.value_in(units.hour)
        else:
            print("Didn't evolve")
            break
        print('Time=', hydro_code.model_time.in_(units.hour))

        write_set_to_file(hydro_code.gas_particles, path_results+f'gas_particles_{index}.hdf5', overwrite_file=True)
        write_set_to_file(hydro_code.dm_particles, path_results+f'dm_particles_{index}.hdf5', overwrite_file=True)
        write_set_to_file(gravity_code.particles, path_results+f'gravity_particles_{index}.hdf5', overwrite_file=True)
    
        index += 1
except Exception as e:
    print("Speed failed:", atmosphere_speed)
    pass


# #----------------------AFTER PLOT-----------------------
# systemplotter(hydro_code.gas_particles, xlabel='x', ylabel='y').plot(save=f"simulation_results/endplot_{atmosphere_speed}.png", c="b",s=0.01, close=False)
# systemplotter(gravity_code.particles, xlabel='x', ylabel='y').plot(save="simulation_results/endplot.png", c="r",s=40, close=False)
# systemplotter(hydro_code.dm_particles, xlabel='x', ylabel='y').plot(save=f"simulation_results/endplot_{atmosphere_speed}.png", c="r",s=20)


#----------------------STOP THE SIMULATION-----------------------
hydro_code.stop()
gravity_code.stop()


# path = f'simulation_results/Earth_planet_{atmosphere_speed}ms/'
# animator = Animator(path, xlabel='x', ylabel='y', xlim=0.005, ylim=0.005, s_atmosphere=0.01)
# animator.make_animation(save_path=f'simulation_results/Earth_planet_{atmosphere_speed}.mp4')



