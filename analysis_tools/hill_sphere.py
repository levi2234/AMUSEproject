from amuse.units import units as u
import numpy as np
from amuse.datamodel import Particles
import os
#load h5 file
from amuse.io import read_set_from_file
from amuse.io import write_set_to_file


def hill_radius(m,M,a,e): 
    return (a*(1-e))*((m/(3*M))**(1/3))

    
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
            
        dm_part = read_set_from_file(f"dm_particles_{j}.hdf5", "hdf5")
        gas_part = read_set_from_file(f"gas_particles_{j}.hdf5", "hdf5")
        gravity_part = read_set_from_file(f"gravity_particles_{j}.hdf5", "hdf5")
        

        #calculate hill radius parameters and hill radii
        a = np.sqrt((dm_part[0].x-gravity_part[0].x)**2 + (dm_part[0].y-gravity_part[0].y)**2 + (dm_part[0].z-gravity_part[0].z)**2)
        a = a.value_in(u.m)
        
        M = dm_part[0].mass.value_in(u.kg) + gravity_part[0].mass.value_in(u.kg)
        rh_planet = hill_radius(dm_part[0].mass.value_in(u.kg),M,a,0) #(beware to set the eccentricity to an appropriate value)
        rh_moon = hill_radius(gravity_part[0].mass.value_in(u.kg),M,a,0) #(beware to set the eccentricity to an appropriate value)
        
        for i in gas_part:
            
            r_to_planet = np.sqrt((i.x-dm_part[0].x)**2 + (i.y-dm_part[0].y)**2 + (i.z-dm_part[0].z)**2)
            r_to_planet = r_to_planet.value_in(u.m)
            if r_to_planet < rh_planet:
                i.in_hill_planet = 1
            else:
                i.in_hill_planet = 0
            
            r_to_moon = np.sqrt((i.x-gravity_part[0].x)**2 + (i.y-gravity_part[0].y)**2 + (i.z-gravity_part[0].z)**2)
            r_to_moon = r_to_moon.value_in(u.m)
            if r_to_moon < rh_moon:
                i.in_hill_moon = 1
            else:
                i.in_hill_moon = 0
        
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
    
    