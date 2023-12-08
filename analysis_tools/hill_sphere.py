from amuse.units import units as u
import numpy as np
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

def hill_radius(m,M,a,e): 
    return (a*(1-e))*((m/(3*M))**(1/3))

def add_hill_sphere_attributes(gas_particles, rh_planet, rh_moon):
 
    distance_to_planet = distance(gas_particles, dm_part[0])
    gas_particles.in_hill_planet = np.zeros(len(gas_particles))
    mask = distance_to_planet < rh_planet
    gas_particles.in_hill_planet[mask] = 1
        
    distance_to_moon = distance(gas_particles, gravity_part[0])
    gas_particles.in_hill_moon = np.zeros(len(gas_particles))
    mask = distance_to_moon < rh_moon
    gas_particles.in_hill_moon[mask] = 1
    
if __name__ == "__main__": 
    
    #navigate to the directory where the h5 files is located relative to this file
    os.chdir("../simulation_results/jupiterlike_planet") #change the folder to the folder where the h5 file is located
    files = os.listdir() #list all files in the directory


    #get the range of numbers in the file names
    file_numbers = []
    for file in files:
        file_numbers.append(int(file.split("_")[2].split(".")[0]))
    file_numbers = np.sort(file_numbers)

    for j in file_numbers:
            
        dm_part = read_set_from_file(f"dm_particles_{j}.hdf5", "hdf5")
        gas_part = read_set_from_file(f"gas_particles_{j}.hdf5", "hdf5")
        gravity_part = read_set_from_file(f"gravity_particles_{j}.hdf5", "hdf5")
        

        #calculate hill radius parameters and hill radii
        a = np.sqrt((dm_part[0].x-gravity_part[0].x)**2 + (dm_part[0].y-gravity_part[0].y)**2 + (dm_part[0].z-gravity_part[0].z)**2)
        a = a.value_in(u.m)
        
        M = dm_part[0].mass.value_in(u.kg) + gravity_part[0].mass.value_in(u.kg)
        eccentricity = 0
        rh_planet = hill_radius(dm_part[0].mass.value_in(u.kg), M, a, eccentricity) #(beware to set the eccentricity to an appropriate value)
        rh_moon = hill_radius(gravity_part[0].mass.value_in(u.kg), M, a, eccentricity) #(beware to set the eccentricity to an appropriate value)
        
        add_hill_sphere_attributes(gas_part, rh_planet, rh_moon)
        
        write_set_to_file(gas_part, f"gas_particles_{j}.hdf5", "hdf5", overwrite_file=True)



    
    
    
    # #plot hill spheres
    # import matplotlib.pyplot as plt
    

    # #plot circle for hill sphere
    # print(f"Planet Hill radius: {rh_planet} m")
    # circle1 = plt.Circle((dm_part[0].x.value_in(u.m),dm_part[0].y.value_in(u.m)),radius=rh_planet,color="b",fill=False)
    # circle2 = plt.Circle((gravity_part[0].x.value_in(u.m),gravity_part[0].y.value_in(u.m)),radius=rh_moon,color="b",fill=False)
    # fig = plt.gcf()
    # ax = fig.gca()
    # ax.add_artist(circle1)
    # ax.add_artist(circle2)
    
    # #plot particles
    # # ax.scatter(gas_part.x.value_in(u.m),gas_part.y.value_in(u.m),s=0.05,c="b")
    # ax.scatter(dm_part.x.value_in(u.m),dm_part.y.value_in(u.m),s=2,c="r")
    # ax.scatter(gravity_part.x.value_in(u.m),gravity_part.y.value_in(u.m),s=2,c="r")
    # ax.axis('equal')
    
    # plt.savefig("../test.png")
    
    