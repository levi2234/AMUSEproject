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


        print(collisions)
    from scipy.stats import norm
    Z=np.array(velocities_list)
    mu,std=norm.fit(Z)
    print(total_collisions)
    plt.hist(velocities_list,bins=10,density=True)
    xmin, xmax = plt.xlim() 
    x = np.linspace(xmin, xmax, 100) 
    p = norm.pdf(x, mu, std)
    plt.plot(x,p,'k-') 
    plt.show()
    N=len(Z)

    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    
    p_value=np.round(normaltest(Z)[1],2)
    plt.plot(X2, F2,label='Simulation results')
    plt.plot(X2, norm.cdf(X2,mu,std),label='Gaussian fit')
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Cumulative distribution function')
    plt.title('Velocity of particles that reach the Hill sphere of the moon')
    plt.text(0.6, 0.5, f'p-value$={p_value}$', transform=plt.gca().transAxes)
    plt.legend()
    plt.show()
    print(norm.cdf(0,mu,std))

