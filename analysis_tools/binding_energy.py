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
    
    
if __name__ == "__main__": 
    
    #navigate to the directory where the h5 files is located relative to this file
    os.chdir("../simulation_results/test_planet") #change the folder to the folder where the h5 file is located
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
        
    
        for i in gas_part: 
            
            gas_speed = np.sqrt(i.vx**2 + i.vy**2 + i.vz**2).value_in(u.m/u.s)
            
            #calculate the binding energy of the gas to the planet
            r_to_planet = np.sqrt((i.x-dm_part[0].x)**2 + (i.y-dm_part[0].y)**2 + (i.z-dm_part[0].z)**2)
            r_to_planet = r_to_planet.value_in(u.m)
            E_binding_to_planet = binding_energy(dm_part[0].mass.value_in(u.kg),i.mass.value_in(u.kg),r_to_planet, gas_speed)
            
                        
            #calculate the binding energy of the gas to the moon
            r_to_moon = np.sqrt((i.x-gravity_part[0].x)**2 + (i.y-gravity_part[0].y)**2 + (i.z-gravity_part[0].z)**2)
            r_to_moon = r_to_moon.value_in(u.m)
            E_binding_to_moon = binding_energy(gravity_part[0].mass.value_in(u.kg),i.mass.value_in(u.kg),r_to_moon, gas_speed)
            
            if (E_binding_to_planet < 0) & (E_binding_to_planet < E_binding_to_moon):
                i.bound_to_planet = 1
            else:
                i.bound_to_planet = 0
            

            if (E_binding_to_planet < 0) & (E_binding_to_planet > E_binding_to_moon):
                i.bound_to_moon = 1
            else:
                i.bound_to_moon = 0
                
        write_set_to_file(gas_part, f"gas_particles_{j}.hdf5", "hdf5", overwrite_file=True)
                
    # snapshot = read_set_from_file(f"gas_particles_0.hdf5", "hdf5")
    # print(snapshot.bound_to_planet