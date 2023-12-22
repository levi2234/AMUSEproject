import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl

from amuse.units import units, constants
from amuse.datamodel import Particles

#load h5 file
from amuse.io import read_set_from_file
from amuse.io import write_set_to_file

def speed(particle_set, unit=units.m/units.s):
    x_term = particle_set.vx.value_in(unit)**2
    y_term = particle_set.vy.value_in(unit)**2
    z_term = particle_set.vz.value_in(unit)**2
    return np.sqrt(x_term + y_term + z_term)

def relative_speed(object1, object2, unit=units.m/units.s):
    x_term = (object1.vx.value_in(unit) - object2.vx.value_in(unit))**2
    y_term = (object1.vy.value_in(unit) - object2.vy.value_in(unit))**2
    z_term = (object1.vz.value_in(unit) - object2.vz.value_in(unit))**2
    return np.sqrt(x_term + y_term + z_term)

def v_esc(particle1, particle2, unit=units.m/units.s):
    mass = (particle2.mass.value_in(units.kg)) # + particle1.mass.value_in(units.kg)
    distance_ = distance(particle1, particle2) 
    G = constants.G.value_in(units.m**3 / units.kg / units.s**2)
    v_esc = np.sqrt(2*G*mass/distance_)
    return v_esc

def distance(object1, object2, unit=units.m):
    x_term = (object1.x.value_in(unit) - object2.x.value_in(unit))**2
    y_term = (object1.y.value_in(unit) - object2.y.value_in(unit))**2
    z_term = (object1.z.value_in(unit) - object2.z.value_in(unit))**2
    return np.sqrt(x_term + y_term + z_term)

def hill_radius(m, M, a, e): 
    return (a*(1-e))*((m/(3*M))**(1/3))

def add_hill_sphere_attributes(gas_particles, rh_planet, rh_moon, dm_part, gravity_part):
    
    # calculate what particles are in the hill sphere of the planet
    distance_to_planet = distance(gas_particles, dm_part[0])
    gas_particles.in_hill_planet = np.zeros(len(gas_particles))
    mask = distance_to_planet < rh_planet
    gas_particles.in_hill_planet[mask] = 1
        
    # calculate what particles are in the hill sphere of the moon
    distance_to_moon = distance(gas_particles, gravity_part[0])
    gas_particles.in_hill_moon = np.zeros(len(gas_particles))
    mask = distance_to_moon < rh_moon
    gas_particles.in_hill_moon[mask] = 1
    # escape_velocities_planet = v_esc(gravity_part[0], dm_part[0])
    radii_moon = gravity_part[0].radius 
    if np.sum(mask)!=0:
        hill_sphere_moon = gas_particles.select(
            lambda in_hill_sphere: in_hill_sphere==1, ['in_hill_moon'])
        escape_velocities_moon = v_esc(hill_sphere_moon, gravity_part[0])
        speeds_particles = speed(hill_sphere_moon)
        v_rel_particles = relative_speed(hill_sphere_moon, gravity_part[0])
        escape_velocities_planet = v_esc(hill_sphere_moon, dm_part[0])

        if np.sum(v_rel_particles<escape_velocities_moon)!=0:
            print(np.sum(v_rel_particles<escape_velocities_moon))
            mask = v_rel_particles < escape_velocities_moon
            if np.sum(speeds_particles[mask]>escape_velocities_planet[mask])!=0:
                print('particle should be bound now')
        crossing_times_moon = radii_moon.value_in(units.m) / speeds_particles
        crossing_times_moon = (crossing_times_moon | units.s).value_in(units.hour)
        # print('crossing_time_moon:', crossing_times_moon)
        return v_rel_particles, escape_velocities_moon, speeds_particles, escape_velocities_planet
    return None, None, None, None