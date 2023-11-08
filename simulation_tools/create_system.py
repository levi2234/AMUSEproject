from amuse.ext.orbital_elements import new_binary_from_orbital_elements
import amuse.units as u
import numpy as np
import create_planet as cp

    
#all variables for our system. This can later be oofloaded to a config file or something
moon_mass = 0.0123 | u.MEarth


#first we create the planet
cp.create_planet() #this function returns an AMUSE Particle Object with the planet and the atmosphere

#then we create a binary system with the planet and moon

