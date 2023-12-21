from analysis_tools.collisions_probabilistic import *
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

explosion_energies=[11800,12000,12100,12100,12350,12375,
                    12400,12450,12475,12500,12525,12550,12575,12600,12700,12900,13200]
N_particles=50000 #number of particles in the hydro code
R_hill=6.15*10**7 #hill radius of the moon in the gas giant system

current_dir=os.getcwd()
accreted_fraction={}
if __name__ == "__main__": 
    for impact_energy in explosion_energies:
        #navigate to the directory where the h5 files is located relative to this file
        os.chdir("simulation_results/rocky_results/Earth_planet_{}ms".format(impact_energy,'g'))  
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
            # calculation of the number of particles that reach the Hill sphere and
            # their velocities
            collisions,velocities=check_collisions(gas_part,gravity_part[0],R_hill)
            total_collisions+=collisions
            for particle_velocity in velocities:
                velocities_list.append(particle_velocity)


        os.chdir(current_dir)
        
        # Calculation of the probability distributions and Gaussian fits
        if total_collisions != 0:
            velocities,empirical_cdf,normal_cdf,pvalue,probability_v0=\
                velocity_cdf(velocities_list)

            if pvalue > 0.2:
                accreted_fraction[format(impact_energy,'g')]=\
                    probability_v0*total_collisions/N_particles
                plt.plot(velocities, empirical_cdf,label='Simulation results')
                plt.plot(velocities, normal_cdf,label='Gaussian fit')
                plt.xlabel('Velocity [m/s]')
                plt.ylabel('Cumulative distribution function')
                plt.text(0.6, 0.5, f'p-value$={pvalue}$', transform=plt.gca().transAxes)
                plt.legend()
                plt.savefig('simulation_results/earth_0.01ratio/CDF{}.png'.format(impact_energy,'g'))
                plt.clf()
                plt.hist(velocities,density=True)
                plt.xlabel('Velocity [m/s]')
                plt.ylabel('Density')
                plt.title('{}'.format(impact_energy,'g'))
                plt.savefig('simulation_results/earth_0.01ratio/PDF{}.png'.format(impact_energy,'g'))
                plt.clf()           
    energies = list(map(float, accreted_fraction.keys()))
    fractions = list(accreted_fraction.values())

    plt.semilogy(energies, fractions, marker='o', linestyle='-',color='g')
    plt.axvspan(0, 12449, alpha=0.3, color='grey')  
    plt.axvspan(13201, 100000, alpha=0.3, color='grey')
    plt.axvspan(12449, 13200, alpha=0.3, color='green')
    plt.text(0.48, 0.5, 'Accretion' ,transform=plt.gca().transAxes,weight='bold',color='g')
    plt.text(0.48, 0.45, 'window' ,transform=plt.gca().transAxes,weight='bold',color='g')
    plt.text(0.02, 0.8, 'Low velocity region' ,transform=plt.gca().transAxes,weight='bold')
    plt.text(0.03, 0.6, 'Collisionless regime' ,transform=plt.gca().transAxes,weight='normal')
    plt.text(0.73, 0.6, 'Escape regime' ,transform=plt.gca().transAxes,weight='normal')
    plt.text(0.675, 0.8, 'High velocity region' , transform=plt.gca().transAxes,weight='bold')
    plt.xlabel('Velocity of the ejected particles [m/s]')
    plt.ylabel('Captured mass fraction')
    plt.title('Fraction of the ejected atmosphere captured by the moon') 
    plt.xlim(11500,14100)
    plt.savefig('simulation_results/earth_0.01ratio/fraction_earth.png')
    