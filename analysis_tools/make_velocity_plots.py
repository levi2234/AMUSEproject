import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl

from amuse.units import units, constants
from amuse.datamodel import Particles

#load h5 file
from amuse.io import read_set_from_file
from amuse.io import write_set_to_file
from hill_sphere import *

def make_velocities_plot(relative_velocities, escape_velocities, path, option='moon'):
    total_relative_velocities = [v_rel[i] for v_rel in relative_velocities for i in range(len(v_rel))]   
    total_escape_velocities = [v_esc[i] for v_esc in escape_velocities for i in range(len(v_esc))]

    fig, ax = plt.subplots()
    fig.suptitle(f'Velocities with respect to the {option}')
    ax.scatter(total_relative_velocities, total_escape_velocities)
    ax.plot(np.arange(0, 40000, 1), np.arange(0, 40000, 1), color='grey', linestyle='--')
    # ax.set_xlim(0, 3000)
    # ax.set_ylim(0, 800)
    ax.set_xlabel('Relative velocity (m/s)')
    ax.set_ylabel(f'Escape velocity {option} (m/s)')
    fig.savefig(path)

def make_plot_all_simulations_rocky(all_relative_velocities, all_escape_velocities, save_results_path, start_velocities, option='moon'):
    cm = plt.cm.get_cmap('plasma', len(start_velocities))
    fig, ax = plt.subplots()
    fig.suptitle(f'Velocities with respect to the {option}', fontsize=15)
    
    for j in range(len(start_velocities)):
        relative_velocities = all_relative_velocities[j]
        escape_velocities = all_escape_velocities[j]
        total_relative_velocities = [v_rel[i] for v_rel in relative_velocities for i in range(len(v_rel))]   
        total_escape_velocities = [v_esc[i] for v_esc in escape_velocities for i in range(len(v_esc))]
        ax.scatter(np.array(total_relative_velocities), np.array(total_escape_velocities), s=4, label=r'$v_\mathrm{eject}$='+str(start_velocities[j])+' m/s', color=cm(j))
    
    ax.plot(np.arange(0, 1600, 1), np.arange(0, 1600, 1), color='grey', linestyle='--')
    ax.set_xlabel('Relative velocity (km/s)', size=15)
    ax.set_ylabel(f'Escape velocity {option} (km/s)', size=15)

    norm = mpl.colors.Normalize(vmin=min(start_velocities), vmax=max(start_velocities))
    cb = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    cb = fig.colorbar(cb, ax=ax)
    cb.set_label(r'$v_\mathrm{eject}$ (m/s)', size=15)

    fig.tight_layout()
    fig.savefig(save_results_path)

def make_plot_all_simulations_gas(all_relative_velocities, all_escape_velocities, save_results_path, explosion_energies, option='moon'):
    cm = plt.cm.get_cmap('plasma', len(explosion_energies))
    fig, ax = plt.subplots()
    fig.suptitle(f'Velocities with respect to the {option}', fontsize=15)
    
    for j in range(len(explosion_energies)):
        relative_velocities = all_relative_velocities[j]
        escape_velocities = all_escape_velocities[j]
        total_relative_velocities = [v_rel[i] for v_rel in relative_velocities for i in range(len(v_rel))]   
        total_escape_velocities = [v_esc[i] for v_esc in escape_velocities for i in range(len(v_esc))]
        ax.scatter(np.array(total_relative_velocities)/1000, np.array(total_escape_velocities)/1000, s=4, label=r'E ='+str(explosion_energies[j]), color=cm(j))

    ax.plot(np.arange(0, 45, 1), np.arange(0, 45, 1), color='grey', linestyle='--')
    ax.set_xlabel('Relative velocity (km/s)', size=15)
    ax.set_ylabel(f'Escape velocity {option} (km/s)', size=15)

    norm = mpl.colors.Normalize(vmin=min(explosion_energies).value_in(units.erg), vmax=max(explosion_energies).value_in(units.erg))
    cb = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    cb = fig.colorbar(cb, ax=ax)
    cb.set_label('Explosion energy (erg)', size=15)

    fig.tight_layout()
    fig.savefig(save_results_path)

def make_velocity_plots_rocky():
    start_velocities = [12475,12500,12525,12550,12575,12600,12700,12900,13200]
    path = 'simulation_results/rocky_results/Earth_planet_'

    all_relative_velocities = []
    all_escape_velocities = []
    all_velocities = []
    all_escape_vel_planet = []
    for v_start in start_velocities:
        print('Ejecting velocity = {} m/s'.format(v_start))
        path_simulation = path + str(v_start) + 'ms/'
        files = os.listdir(path_simulation)

        # get the range of numbers in the file names
        file_numbers = []
        for file in files:
            if file.endswith('.hdf5'):
                file_numbers.append(int(file.split("_")[2].split(".")[0]))
        file_numbers = np.sort(file_numbers)

        relative_velocities_particles = []
        escape_velocities_moon = []
        velocities_particles = []
        escape_velocities_planet = []
        for j in file_numbers:
                
            dm_part = read_set_from_file(path_simulation+f"dm_particles_{j}.hdf5", "hdf5")
            gas_part = read_set_from_file(path_simulation+f"gas_particles_{j}.hdf5", "hdf5")
            gravity_part = read_set_from_file(path_simulation+f"gravity_particles_{j}.hdf5", "hdf5")
            
            # calculate hill radius parameters and hill radii
            a = distance(dm_part[0], gravity_part[0], unit=units.m)         
            M = dm_part[0].mass.value_in(units.kg) + gravity_part[0].mass.value_in(units.kg)
            eccentricity = 0
            rh_planet = hill_radius(dm_part[0].mass.value_in(units.kg), M, a, eccentricity) 
            rh_moon = hill_radius(gravity_part[0].mass.value_in(units.kg), M, a, eccentricity) 
            
            v_rel_particles, v_esc_moon, v_particles, v_esc_planet = add_hill_sphere_attributes(gas_part, rh_planet, rh_moon, dm_part, gravity_part)
            if v_rel_particles is not None:
                relative_velocities_particles.append(v_rel_particles)
                escape_velocities_moon.append(v_esc_moon)
                velocities_particles.append(v_particles)
                escape_velocities_planet.append(v_esc_planet)
            # write_set_to_file(gas_part, f"gas_particles_{j}.hdf5", "hdf5", overwrite_file=True)
        
        all_relative_velocities.append(relative_velocities_particles)
        all_escape_velocities.append(escape_velocities_moon)
        all_velocities.append(velocities_particles)
        all_escape_vel_planet.append(escape_velocities_planet)
        
        save_result_path = f'simulation_results/rocky_results/velocities_in_hill_sphere_moon_{v_start}.png'
        save_result_path2 = f'simulation_results/rocky_results/velocities_in_hill_sphere_planet_{v_start}.png'
        make_velocities_plot(relative_velocities_particles, escape_velocities_moon, save_result_path)
        make_velocities_plot(velocities_particles, escape_velocities_planet, save_result_path2, option='planet')

    save_results_path = 'velocities_in_hill_sphere_moon_rocky.png'
    save_results_path2 = 'velocities_in_hill_sphere_planet_rocky.png'
    make_plot_all_simulations_rocky(all_velocities, all_escape_vel_planet, save_results_path2, start_velocities, option='planet')
    make_plot_all_simulations_rocky(all_relative_velocities, all_escape_velocities, save_results_path, start_velocities)

def make_velocity_plots_gas():
    explosion_energies = np.arange(1,10.1,0.3)*1e42|units.erg
    path = 'simulation_results/gas_results/'

    all_relative_velocities = []
    all_escape_velocities = []
    all_velocities = []
    all_escape_vel_planet = []
    for explosion_energy in explosion_energies:
        print('explosion energy = {:.1e} erg'.format(explosion_energy.value_in(units.erg)))
        path_simulation = path+'{:.1e}.erg/'.format(explosion_energy.value_in(units.erg), 'g')
        files = os.listdir(path_simulation)

        # get the range of numbers in the file names
        file_numbers = []
        for file in files:
            if file.endswith('.hdf5'):
                file_numbers.append(int(file.split("_")[2].split(".")[0]))
        file_numbers = np.sort(file_numbers)

        relative_velocities_particles = []
        escape_velocities_moon = []
        velocities_particles = []
        escape_velocities_planet = []
        for j in file_numbers:
                
            dm_part = read_set_from_file(path_simulation+f"dm_particles_{j}.hdf5", "hdf5")
            gas_part = read_set_from_file(path_simulation+f"gas_particles_{j}.hdf5", "hdf5")
            gravity_part = read_set_from_file(path_simulation+f"gravity_particles_{j}.hdf5", "hdf5")
            
            # calculate hill radius parameters and hill radii
            a = distance(dm_part[0], gravity_part[0], unit=units.m)           
            M = dm_part[0].mass.value_in(units.kg) + gravity_part[0].mass.value_in(units.kg)
            eccentricity = 0
            rh_planet = hill_radius(dm_part[0].mass.value_in(units.kg), M, a, eccentricity) 
            rh_moon = hill_radius(gravity_part[0].mass.value_in(units.kg), M, a, eccentricity) 
            
            v_rel_particles, v_esc_moon, v_particles, v_esc_planet = add_hill_sphere_attributes(gas_part, rh_planet, rh_moon, dm_part, gravity_part)
            if v_rel_particles is not None:
                relative_velocities_particles.append(v_rel_particles)
                escape_velocities_moon.append(v_esc_moon)
                velocities_particles.append(v_particles)
                escape_velocities_planet.append(v_esc_planet)
            # write_set_to_file(gas_part, f"gas_particles_{j}.hdf5", "hdf5", overwrite_file=True)
        all_relative_velocities.append(relative_velocities_particles)
        all_escape_velocities.append(escape_velocities_moon)
        all_velocities.append(velocities_particles)
        all_escape_vel_planet.append(escape_velocities_planet)
        
        save_result_path = f'simulation_results/gas_results/velocities_in_hill_sphere_moon_{explosion_energy.value_in(units.erg)}_erg.png'
        save_result_path2 = f'simulation_results/gas_results/velocities_in_hill_sphere_planet_{explosion_energy.value_in(units.erg)}_erg.png'
        make_velocities_plot(relative_velocities_particles, escape_velocities_moon, save_result_path)
        make_velocities_plot(velocities_particles, escape_velocities_planet, save_result_path2, option='planet')

    save_results_path = 'velocities_in_hill_sphere_moon_gas.png'
    save_results_path2 = 'velocities_in_hill_sphere_planet_gas.png'
    make_plot_all_simulations_gas(all_relative_velocities, all_escape_velocities, save_results_path, explosion_energies)
    make_plot_all_simulations_gas(all_velocities, all_escape_vel_planet, save_results_path2, explosion_energies, option='planet')
    
if __name__ == "__main__": 
    print('rocky:')
    make_velocity_plots_rocky()
    print('gassy:')
    make_velocity_plots_gas()
    # check_long_run()
    