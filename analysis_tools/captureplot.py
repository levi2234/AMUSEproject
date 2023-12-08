from amuse.units import units as u
import numpy as np
from amuse.datamodel import Particles
import os
#load h5 file
from amuse.io import read_set_from_file
from amuse.io import write_set_to_file

if __name__ == "__main__": 
    
    #navigate to the directory where the h5 files is located relative to this file
    os.chdir("../simulation_results/test_planet") #change the folder to the folder where the h5 file is located
    files = os.listdir() #list all files in the directory


    #get the range of numbers in the file names
    file_numbers = []
    for file in files:
        file_numbers.append(int(file.split("_")[2].split(".")[0]))
    file_numbers = np.sort(file_numbers)

    planet_hill_contained_percentage = []
    moon_hill_contained_percentage = []
    
    planet_bound_percentage = []
    moon_bound_percentage = []
    
    #loop over all files
    for j in file_numbers:
        try:
            gas_part = read_set_from_file(f"gas_particles_{j}.hdf5", "hdf5")
        except:
            print(f"file {j} not be found for gas particles")
            #remove the file number from the list of file numbers
            file_numbers = np.delete(file_numbers, np.where(file_numbers == j))
            continue
        
        planet_hill_contained= np.sum(gas_part.mass[np.array(gas_part.in_hill_planet).astype(bool)])/np.sum(gas_part.mass) #percentage
        moon_hill_contained = np.sum(gas_part.mass[np.array(gas_part.in_hill_moon).astype(bool)])/np.sum(gas_part.mass) #percentage
        
        planet_hill_contained_percentage.append(planet_hill_contained)
        moon_hill_contained_percentage.append(moon_hill_contained)
        
        planet_bound = np.sum(gas_part.mass[np.array(gas_part.bound_to_planet).astype(bool)])/np.sum(gas_part.mass) #percentage
        moon_bound = np.sum(gas_part.mass[np.array(gas_part.bound_to_moon).astype(bool)])/np.sum(gas_part.mass)#percentage
        
        planet_bound_percentage.append(planet_bound)
        moon_bound_percentage.append(moon_bound)
        
    
    #plot the results
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10,10))
    plt.plot(file_numbers, planet_hill_contained_percentage, label="planet hill contained")
    plt.plot(file_numbers, moon_hill_contained_percentage, label="moon hill contained")
    plt.plot(file_numbers, planet_bound_percentage, label="planet bound")
    plt.plot(file_numbers, moon_bound_percentage, label="moon bound")
    plt.legend()
    plt.xlabel("file number")
    plt.ylabel("percentage")
    plt.title("percentage of gas particles bound to the planet and the moon")
    plt.savefig("../binding_percentage.png")
        
        
    
        
        