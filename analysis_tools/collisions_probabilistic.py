from pathlib import Path
from amuse.units import units as u
from amuse.units import constants as c
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from scipy.stats import normaltest
from amuse.datamodel import Particles
import os
from amuse.io import read_set_from_file
from amuse.io import write_set_to_file

def distance(object1, object2, unit=u.m):
    '''
    A function that computes the distance between two objects
    
    Args:
        object1: gas particle
        object2: moon
    
    Returns:
        distance in meters
    '''
    x_term = (object1.x.value_in(unit) - object2.x.value_in(unit))**2
    y_term = (object1.y.value_in(unit) - object2.y.value_in(unit))**2
    z_term = (object1.z.value_in(unit) - object2.z.value_in(unit))**2
    return np.sqrt(x_term + y_term + z_term)

def velocity_modulus(object, unit=u.m*u.s**-1):
    '''
    A function that computes the modulus of the velocity of a particle
    
    Args:
        object: particle
    
    Returns:
        velocity: the velocity of the particle in m/s
    '''
    vx=object.vx.value_in(unit)
    vy=object.vy.value_in(unit)
    vz=object.vz.value_in(unit)
    velocity=np.sqrt(vx**2+vy**2+vz**2)
    return velocity


def check_collisions(object1,object2,radius_object_2):
    '''
    Calculates the number of particles within a region 
    and their velocities
    
    Args: 
        object1: gas particles from the planet atmosphere
        object2: moon
        radius_object_2: moon radius / radius of its Hill sphere
    Returns:
        number_of_collisions: number of collisions found
        velocities: velocity of each particle that collided
    '''
    distance_moon_particle=distance(object1,object2)
    number_of_collisions=len(distance_moon_particle[distance_moon_particle<=radius_object_2])
    indexes_particles=np.where(distance_moon_particle[distance_moon_particle<=radius_object_2])
    velocities=velocity_modulus(object1[indexes_particles])
    return number_of_collisions,velocities

def velocity_cdf(velocity_list):
    '''
    A function that calculates the probability distribution of the velocities
    distribution, a Gaussian fit of that probability and the probability of having 
    0 velocity

    Args:
        velocity_list: list of velocities
    
    Returns:
        sorted_velocities: velocities from smaller to biggest
        empirical_CDF: the true probability distribution
        normal_pdf: a Gaussian fit of the distribution
        p_value: the p-value of a normality test
        probability_zero_val: the probability of having 0 velocity.
    '''
    velocities_array=np.array(velocity_list)
    p_value=np.round(normaltest(velocities_array)[1],2)
    mu,std=norm.fit(velocities_array)
    N=len(velocities_array)
    sorted_velocities = np.sort(velocities_array)
    empirical_CDF = np.array(range(N))/float(N)
    normal_CDF = norm.cdf(sorted_velocities,mu,std)
    probability_zero_vel=(norm.cdf(0,mu,std))
    return sorted_velocities, empirical_CDF, normal_CDF,p_value,probability_zero_vel

