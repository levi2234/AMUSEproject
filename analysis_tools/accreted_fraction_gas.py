from collisions_probabilistic import *
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
from hill_sphere import hill_radius

explosion_energies=np.arange(1,10.1,0.3)*1e42|u.erg
N_particles=2000 #number of particles in the hydro code
m_moon=20*5.972*10**24
m_planet=300*5.972*10**24
semi_major_axis=421600000
R_hill=hill_radius(m_moon,m_planet,semi_major_axis,0)
current_dir = os.getcwd()
accreted_fraction = {}
if __name__ == "__main__": 
    for impact_energy in explosion_energies:
        #navigate to the directory where the h5 files is located relative to this file
        os.chdir("simulation_results/gas_results/{:.1e}_erg".format(impact_energy.value_in(u.erg),'g')) 
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
        
        # check if directories exist that we want to save our plots in
        if not os.path.exists('simulation_results/earth_0.01ratio'):
            os.mkdir('simulation_results/earth_0.01ratio')
        if not os.path.exists('simulation_results/jupiter_0.05ratio'):
            os.mkdir('simulation_results/jupiter_0.05ratio')
        
        # calculation of the probability distributions and Gaussian fits
        if total_collisions != 0:
            velocities,empirical_cdf,normal_cdf,pvalue,probability_v0=\
                velocity_cdf(velocities_list)

            if pvalue > 0.2:
                # Calculation of the accreted fraction
                accreted_fraction[format(impact_energy.value_in(u.erg),'g')]=\
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
                plt.savefig('simulation_results/jupiter_0.05ratio/PDF{}.png'.format(impact_energy,'g'))
                plt.clf()           
    energies = list(map(float, accreted_fraction.keys()))
    fractions = list(accreted_fraction.values())

    plt.rcParams['figure.dpi'] = 150    
    plt.loglog(energies, fractions, marker='o', linestyle='-',color='g')
    plt.axvspan(0, 2.3 * 10**42, alpha=0.3, color='grey')  
    plt.axvspan(6.5 * 10**42, 3*10**43, alpha=0.3, color='grey')
    plt.axvspan(2.3 * 10**42, 6.5*10**42, alpha=0.3, color='green')
    plt.text(0.5, 0.5, 'Accretion' ,transform=plt.gca().transAxes,weight='bold',color='g')
    plt.text(0.5, 0.45, 'window' ,transform=plt.gca().transAxes,weight='bold',color='g')
    plt.text(0.02, 0.8, 'Low velocity region' ,transform=plt.gca().transAxes,weight='bold')
    plt.text(0.03, 0.6, 'Collisionless regime' ,transform=plt.gca().transAxes,weight='normal')
    plt.text(0.73, 0.6, 'Escape regime' ,transform=plt.gca().transAxes,weight='normal')
    plt.text(0.675, 0.8, 'High velocity region' , transform=plt.gca().transAxes,weight='bold')
    plt.xlabel('Energy of the explosion [erg]')
    plt.ylabel('Captured mass fraction')
    plt.title('Fraction of the ejected atmosphere captured by the moon') 
    plt.xlim(8e41,1.85e43)
    custom_ticks = [1e42, 5e42, 1e43]
    plt.xticks(custom_ticks, [f'{tick:.0e}' for tick in custom_ticks])
    plt.savefig('simulation_results/jupiter_0.05ratio/fraction_jupiter.png')

    