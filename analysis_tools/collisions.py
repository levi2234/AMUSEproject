from pathlib import Path

from amuse.units import units as u
from amuse.units import constants as c
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from scipy.stats import normaltest
from amuse.datamodel import Particles
import os
#load h5 file
from amuse.io import read_set_from_file
from amuse.io import write_set_to_file

def distance(object1, object2, unit=u.m):
    x_term = (object1.x.value_in(unit) - object2.x.value_in(unit))**2
    y_term = (object1.y.value_in(unit) - object2.y.value_in(unit))**2
    z_term = (object1.z.value_in(unit) - object2.z.value_in(unit))**2
    return np.sqrt(x_term + y_term + z_term)

def velocity_modulus(object, unit=u.m*u.s**-1):
    vx=object.vx.value_in(unit)
    vy=object.vy.value_in(unit)
    vz=object.vz.value_in(unit)
    velocity=np.sqrt(vx**2+vy**2+vz**2)
    return velocity


def check_collisions(object1,object2,radius_object_2):
    distance_moon_particle=distance(object1,object2)
    number_of_collisions=len(distance_moon_particle[distance_moon_particle<=radius_object_2])
    indexes_particles=np.where(distance_moon_particle[distance_moon_particle<=radius_object_2])
    velocities=velocity_modulus(object1[indexes_particles])
    return number_of_collisions,velocities

def velocity_cdf(velocity_list):
    velocities_array=np.array(velocities_list)
    p_value=np.round(normaltest(velocities_array)[1],2)
    mu,std=norm.fit(velocities_array)
    N=len(velocities_array)
    sorted_velocities = np.sort(velocities_array)
    empirical_CDF = np.array(range(N))/float(N)
    normal_CDF = norm.cdf(sorted_velocities,mu,std)
    probability_zero_vel=(norm.cdf(0,mu,std))
    return sorted_velocities, empirical_CDF, normal_CDF,p_value,probability_zero_vel

if __name__ == "__main__": 
    
    #navigate to the directory where the h5 files is located relative to this file
    os.chdir("../simulation_results/jupiterlike_planettest_fast") #change the folder to the folder where the h5 file is located
    files = os.listdir() #list all files in the directory


    #get the range of numbers in the file names
    file_numbers = []
    for file in files:
        file_numbers.append(int(file.split("_")[2].split(".")[0]))
    file_numbers = np.sort(file_numbers)
    total_collisions=0
    velocities_list=[]
    for j in np.unique(file_numbers):
        try:
            dm_part = read_set_from_file(f"dm_particles_{j}.hdf5", "hdf5")
            gas_part = read_set_from_file(f"gas_particles_{j}.hdf5", "hdf5")
            gravity_part = read_set_from_file(f"gravity_particles_{j}.hdf5", "hdf5")
        except:
            print(f"file {j} not found for either dm, gas or gravity particles")
            continue
        collisions,velocities=check_collisions(gas_part,gravity_part[0],1.18*10**8)
        total_collisions+=collisions
        for particle_velocity in velocities:
            velocities_list.append(particle_velocity)


        #print(collisions)
    print(total_collisions)
    #path_results = '/data2/presa/AMUSEproject/simulation_results/probability_distributions/'
    #path_results = Path.cwd().parent/'probability_distributions'
    #if not os.path.exists(path_results):
    #    os.mkdir(path_results)
    velocities,empirical_cdf,normal_cdf,pvalue,probability_v0=velocity_cdf(velocities_list)
    plt.plot(velocities, empirical_cdf,label='Simulation results')
    plt.plot(velocities, normal_cdf,label='Gaussian fit')
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Cumulative distribution function')
    plt.title('Velocity of particles that reach the Hill sphere of the moon')
    plt.text(0.6, 0.5, f'p-value$={pvalue}$', transform=plt.gca().transAxes)
    plt.legend()
    plt.show()
 

