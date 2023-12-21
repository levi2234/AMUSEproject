from amuse.units import units as u
from amuse.units import constants as c
import numpy as np
from amuse.datamodel import Particles
import os
#load h5 file
from amuse.io import read_set_from_file
from amuse.io import write_set_to_file

def potential_energy(m1,m2,r):
    G =  c.G.value_in(u.m**3 / (u.kg * u.s**2))
    return (-G*m1*m2/r)

def kinetic_energy(m,v):
    return (0.5*m*v**2)

def binding_energy(m1,m2,r, v): #here m1 is the mass of the planet, m2 is the mass of the moon, r is the distance between the planet and the moon, v is the velocity of the moon
    U =potential_energy(m1,m2,r)
    
    T = kinetic_energy(m2,v)
    binding_energy = U+T
   
    return binding_energy

def speed(particle_set, unit=u.m/u.s):
    x_term = particle_set.vx.value_in(unit)**2
    y_term = particle_set.vy.value_in(unit)**2
    z_term = particle_set.vz.value_in(unit)**2
    return np.sqrt(x_term + y_term + z_term)

def distance(object1, object2, unit=u.m):
    x_term = (object1.x.value_in(unit) - object2.x.value_in(unit))**2
    y_term = (object1.y.value_in(unit) - object2.y.value_in(unit))**2
    z_term = (object1.z.value_in(unit) - object2.z.value_in(unit))**2
    return np.sqrt(x_term + y_term + z_term)   

def add_binding_energy_attributes(gas_particles):
    
    gas_speed = speed(gas_particles)
    
    distance_to_planet = distance(gas_particles, dm_part[0])
    E_binding_to_planet = binding_energy(m1=dm_part[0].mass.value_in(u.kg), \
                                         m2=gas_particles.mass.value_in(u.kg),\
                                         r=distance_to_planet,\
                                         v=gas_speed)
    
    distance_to_moon = distance(gas_particles, gravity_part[0])
    E_binding_to_moon = binding_energy(m1=gravity_part[0].mass.value_in(u.kg),
                                       m2=gas_particles.mass.value_in(u.kg),
                                       r=distance_to_moon,
                                       v=gas_speed)
    
    gas_particles.bound_to_planet = np.zeros(len(gas_particles))
    mask = (E_binding_to_planet < 0) & (E_binding_to_planet < E_binding_to_moon)
    gas_particles.bound_to_planet[mask] = 1

    gas_particles.bound_to_moon = np.zeros(len(gas_particles))
    mask = (E_binding_to_moon < 0) & (E_binding_to_planet > E_binding_to_moon)
    gas_particles.bound_to_moon[mask] = 1

if __name__ == "__main__": 
    
    #navigate to the directory where the h5 files is located relative to this file
    os.chdir("simulation_results/rocky_results/Earth_planet_12550/") #change the folder to the folder where the h5 file is located
    files = os.listdir() #list all files in the directory


    #get the range of numbers in the file names
    file_numbers = []
    for file in files:
        file_numbers.append(int(file.split("_")[2].split(".")[0]))
    file_numbers = np.sort(file_numbers)

    for j in file_numbers:
        try:
            dm_part = read_set_from_file(f"dm_particles_{j}.hdf5", "hdf5")
            gas_part = read_set_from_file(f"gas_particles_{j}.hdf5", "hdf5")
            gravity_part = read_set_from_file(f"gravity_particles_{j}.hdf5", "hdf5")
        except:
            print(f"file {j} not found for either dm, gas or gravity particles")
            continue
        
        add_binding_energy_attributes(gas_part)

       
        write_set_to_file(gas_part, f"gas_particles_{j}.hdf5", "hdf5", overwrite_file=True)
                
    # snapshot = read_set_from_file(f"gas_particles_0.hdf5", "hdf5")
    # print(snapshot.bound_to_planet