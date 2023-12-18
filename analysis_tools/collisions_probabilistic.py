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

#explosion_energies=np.arange(1,10.1,0.3)*1e42|u.erg
explosion_energies=[11800,12000,12100,12100,12350,12375,
                    12400,12450,12475,12500,12525,12550,12575,12600,12700,12900,13200]
current_dir=os.getcwd()
accreted_fraction={}
if __name__ == "__main__": 
    for impact_energy in explosion_energies:
        #navigate to the directory where the h5 files is located relative to this file
        #os.chdir("../simulation_results/energies_results/{}".format(impact_energy,'g')) 
        os.chdir("../simulation_results/earth_planet/Earth_planet_{}ms".format(impact_energy,'g'))  
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
            collisions,velocities=check_collisions(gas_part,gravity_part[0],6.15*10**7)
            total_collisions+=collisions
            for particle_velocity in velocities:
                velocities_list.append(particle_velocity)


            #print(collisions)
        os.chdir(current_dir)
        print(total_collisions)
        if total_collisions != 0:
            velocities,empirical_cdf,normal_cdf,pvalue,probability_v0=velocity_cdf(velocities_list)
            plt.plot(velocities, empirical_cdf,label='Simulation results')
            plt.plot(velocities, normal_cdf,label='Gaussian fit')
            plt.xlabel('Velocity [m/s]')
            plt.ylabel('Cumulative distribution function')
            #plt.title('{}'.format(impact_energy,'g'))
            plt.text(0.6, 0.5, f'p-value$={pvalue}$', transform=plt.gca().transAxes)
            plt.legend()
            plt.show()
            #plt.savefig('../simulation_results/jupiter_0.05ratio/CDF{}.png'.format(impact_energy,'g'))
            plt.clf()
            # uncomment to show histograms as well
            #plt.hist(velocities,density=True)
            #plt.xlabel('Velocity [m/s]')
            #plt.ylabel('Density')
            #plt.title('{}'.format(impact_energy,'g'))
            #plt.savefig('../simulation_results/jupiter_0.05ratio/PDF{}.png'.format(impact_energy,'g'))
            #plt.clf()
            if pvalue > 0.2:
                #accreted_fraction[format(impact_energy.value_in(u.erg),'g')]=probability_v0*total_collisions/2000
                accreted_fraction[format(impact_energy,'g')]=probability_v0*total_collisions/50000
            print(accreted_fraction)
    energies = list(map(float, accreted_fraction.keys()))
    fractions = list(accreted_fraction.values())

    #plt.rcParams['figure.dpi'] = 150    
    #plt.loglog(energies, fractions, marker='o', linestyle='-',color='g')
    #plt.axvspan(0, 2.3 * 10**42, alpha=0.3, color='grey')  
    #plt.axvspan(6.5 * 10**42, 3*10**43, alpha=0.3, color='grey')
    #plt.axvspan(2.3 * 10**42, 6.5*10**42, alpha=0.3, color='green')
    #plt.text(0.5, 0.5, 'Accretion' ,transform=plt.gca().transAxes,weight='bold',color='g')
    #plt.text(0.5, 0.45, 'window' ,transform=plt.gca().transAxes,weight='bold',color='g')
    #plt.text(0.02, 0.8, 'Low velocity region' ,transform=plt.gca().transAxes,weight='bold')
    #plt.text(0.03, 0.6, 'Collisionless regime' ,transform=plt.gca().transAxes,weight='normal')
    #plt.text(0.73, 0.6, 'Escape regime' ,transform=plt.gca().transAxes,weight='normal')
    #plt.text(0.675, 0.8, 'High velocity region' , transform=plt.gca().transAxes,weight='bold')
    #plt.xlabel('Energy of the explosion [erg]')
    #plt.ylabel('Captured mass fraction')
    #plt.title('Fraction of the ejected atmosphere captured by the moon') 
    #plt.xlim(8e41,1.85e43)
    #custom_ticks = [1e42, 5e42, 1e43]
    #plt.xticks(custom_ticks, [f'{tick:.0e}' for tick in custom_ticks])
    #plt.show()
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
    plt.savefig('../simulation_results/earth_0.01ratio/fraction_earth.png')