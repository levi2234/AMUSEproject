"""This file contains the functions to create a planet.
The funciton returns a planet object with the following attributes:
The planet at its core and the atmosphere around it.

The planet is created as an AMUSE hydrodynamics particle set.

The function returns an AMUSE Particle Object with the planet and the atmosphere
""" 

import numpy as np
from amuse.datamodel import Particles
from simulation_tools import create_atmosphere as CA
from simulation_tools import create_bare_planet as CBP
from simulation_tools.profiles import atmospheric_profiles as AP
#import module from parent directory
from amuse.units import units
#example of combined units

def create_planet_and_atmosphere(planet_mass, planet_radius, atmosphere_profile, atmosphere_height, atmosphere_density, num_points=10000 , **kwargs):
    
    #create the atmosphere particle
    outer_radius = planet_radius + atmosphere_height
    atmosphere = CA.create_atmosphere_from_profile(profile=atmosphere_profile, num_points=num_points, inner_radius=planet_radius, outer_radius=planet_radius + atmosphere_height, rho0= atmosphere_density, **kwargs)
    
    planet = CBP.create_planet(planet_mass, planet_radius, **kwargs)
    return planet , atmosphere

