import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl

from amuse.units import units, constants
from amuse.datamodel import Particles

#load h5 file
from amuse.io import read_set_from_file
from amuse.io import write_set_to_file

def speed(particle_set, unit=units.m/units.s):
    x_term = particle_set.vx.value_in(unit)**2
    y_term = particle_set.vy.value_in(unit)**2
    z_term = particle_set.vz.value_in(unit)**2
    return np.sqrt(x_term + y_term + z_term)

def relative_speed(object1, object2, unit=units.m/units.s):
    x_term = (object1.vx.value_in(unit) - object2.vx.value_in(unit))**2
    y_term = (object1.vy.value_in(unit) - object2.vy.value_in(unit))**2
    z_term = (object1.vz.value_in(unit) - object2.vz.value_in(unit))**2
    return np.sqrt(x_term + y_term + z_term)

def v_esc(particle1, particle2, unit=units.m/units.s):
    mass = (particle2.mass.value_in(units.kg))# + particle1.mass.value_in(units.kg))
    distance_ = distance(particle1, particle2) 
    G = constants.G.value_in(units.m**3 / units.kg / units.s**2)
    v_esc = np.sqrt(2*G*mass/distance_)
    return v_esc

def distance(object1, object2, unit=units.m):
    x_term = (object1.x.value_in(unit) - object2.x.value_in(unit))**2
    y_term = (object1.y.value_in(unit) - object2.y.value_in(unit))**2
    z_term = (object1.z.value_in(unit) - object2.z.value_in(unit))**2
    return np.sqrt(x_term + y_term + z_term)

def hill_radius(m, M, a, e): 
    return (a*(1-e))*((m/(3*M))**(1/3))

def add_hill_sphere_attributes(gas_particles, rh_planet, rh_moon, dm_part, gravity_part):
    
    # calculate what particles are in the hill sphere of the planet
    distance_to_planet = distance(gas_particles, dm_part[0])
    gas_particles.in_hill_planet = np.zeros(len(gas_particles))
    mask = distance_to_planet < rh_planet
    gas_particles.in_hill_planet[mask] = 1
        
    # calculate what particles are in the hill sphere of the moon
    distance_to_moon = distance(gas_particles, gravity_part[0])
    gas_particles.in_hill_moon = np.zeros(len(gas_particles))
    mask = distance_to_moon < rh_moon
    gas_particles.in_hill_moon[mask] = 1
    # escape_velocities_planet = v_esc(gravity_part[0], dm_part[0])
    radii_moon = gravity_part[0].radius 
    if np.sum(mask)!=0:
        hill_sphere_moon = gas_particles.select(
            lambda in_hill_sphere: in_hill_sphere==1, ['in_hill_moon'])
        escape_velocities_moon = v_esc(hill_sphere_moon, gravity_part[0])
        speeds_particles = speed(hill_sphere_moon)
        v_rel_particles = relative_speed(hill_sphere_moon, gravity_part[0])
        escape_velocities_planet = v_esc(hill_sphere_moon, dm_part[0])

        if np.sum(v_rel_particles<escape_velocities_moon)!=0:
            print(np.sum(v_rel_particles<escape_velocities_moon))
            if np.sum(speeds_particles>escape_velocities_planet):
                print('particle should be bound now')
        crossing_times_moon = radii_moon.value_in(units.m) / speeds_particles
        crossing_times_moon = (crossing_times_moon | units.s).value_in(units.hour)

        return v_rel_particles, escape_velocities_moon, speeds_particles, escape_velocities_planet
    return None, None, None, None

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
    cm = plt.cm.get_cmap('hot', len(start_velocities))
    fig, ax = plt.subplots()
    fig.suptitle(f'Velocities with respect to the {option}')
    
    for j in range(len(start_velocities)):
        relative_velocities = all_relative_velocities[j]
        escape_velocities = all_escape_velocities[j]
        total_relative_velocities = [v_rel[i] for v_rel in relative_velocities for i in range(len(v_rel))]   
        total_escape_velocities = [v_esc[i] for v_esc in escape_velocities for i in range(len(v_esc))]
        # ax.scatter(total_relative_velocities, total_escape_velocities, s=4, label=r'$v_\mathrm{start}$='+str(start_velocities[j])+' m/s')
        ax.scatter(np.array(total_relative_velocities), np.array(total_escape_velocities), s=4, label=r'$v_\mathrm{eject}$='+str(start_velocities[j])+' m/s', color=cm(j))
    ax.plot(np.arange(0, 1600, 1), np.arange(0, 1600, 1), color='grey', linestyle='--')
    ax.set_xlabel('Relative velocity (km/s)')
    ax.set_ylabel(f'Escape velocity {option} (km/s)')
    # ax.legend()
    norm = mpl.colors.Normalize(vmin=min(start_velocities), vmax=max(start_velocities))
    cb = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    # cb.set_label('Explosion energy (erg)')
    fig.colorbar(cb, ax=ax, label=r'$v_\mathrm{eject}$ (m/s)')
    fig.savefig(save_results_path)

def make_plot_all_simulations_gas(all_relative_velocities, all_escape_velocities, save_results_path, explosion_energies, option='moon'):
    cm = plt.cm.get_cmap('hot', len(explosion_energies))
    fig, ax = plt.subplots()
    fig.suptitle(f'Velocities with respect to the {option}')
    
    for j in range(len(explosion_energies)):
        relative_velocities = all_relative_velocities[j]
        escape_velocities = all_escape_velocities[j]
        total_relative_velocities = [v_rel[i] for v_rel in relative_velocities for i in range(len(v_rel))]   
        total_escape_velocities = [v_esc[i] for v_esc in escape_velocities for i in range(len(v_esc))]
        # ax.scatter(total_relative_velocities, total_escape_velocities, s=4, label=r'$v_\mathrm{start}$='+str(start_velocities[j])+' m/s')
        ax.scatter(np.array(total_relative_velocities)/1000, np.array(total_escape_velocities)/1000, s=4, label=r'E ='+str(explosion_energies[j]), color=cm(j))
    ax.plot(np.arange(0, 45, 1), np.arange(0, 45, 1), color='grey', linestyle='--')
    ax.set_xlabel('Relative velocity (km/s)')
    ax.set_ylabel(f'Escape velocity {option} (km/s)')
    # ax.legend()
    norm = mpl.colors.Normalize(vmin=min(explosion_energies).value_in(units.erg), vmax=max(explosion_energies).value_in(units.erg))
    cb = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    # cb.set_label('Explosion energy (erg)')
    fig.colorbar(cb, ax=ax, label='Explosion energy (erg)')
    fig.savefig(save_results_path)

def make_velocity_plots_rocky():
    start_velocities = [12475,12500,12525,12550,12575,12600,12700,12900,13200]
    path = 'simulation_results/rocky_results/Earth_planet_'
    path = '/net/vdesk/data2/vanes/AMUSEproj/AMUSEproj/AMUSEproject/simulation_results/Earth_planet_'

    all_relative_velocities = []
    all_escape_velocities = []
    all_velocities = []
    all_escape_vel_planet = []
    for v_start in start_velocities:
        os.chdir(path+str(v_start)+'ms')
        files = os.listdir()

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
                
            dm_part = read_set_from_file(f"dm_particles_{j}.hdf5", "hdf5")
            gas_part = read_set_from_file(f"gas_particles_{j}.hdf5", "hdf5")
            gravity_part = read_set_from_file(f"gravity_particles_{j}.hdf5", "hdf5")
            
            # calculate hill radius parameters and hill radii
            a = np.sqrt((dm_part[0].x-gravity_part[0].x)**2 + (dm_part[0].y-gravity_part[0].y)**2 + (dm_part[0].z-gravity_part[0].z)**2)
            a = a.value_in(units.m)
            
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
        save_result_path = f'/net/vdesk/data2/evdijk/SMA/AMUSEproject-1/simulation_results/velocities_in_hill_sphere_{v_start}.png' 
        save_result_path2 = f'/net/vdesk/data2/evdijk/SMA/AMUSEproject-1/simulation_results/velocities_wrtplanet_in_hill_sphere_{v_start}.png'
        make_velocities_plot(relative_velocities_particles, escape_velocities_moon, save_result_path)
        make_velocities_plot(velocities_particles, escape_velocities_planet, save_result_path2, option='planet')

    save_results_path = 'velocities_in_hill_sphere_moon_rocky.png'
    save_results_path2 = 'velocities_in_hill_sphere_planet_rocky.png'
    save_results_path = '/net/vdesk/data2/evdijk/SMA/AMUSEproject-1/velocities_in_hill_sphere_earthlike.png'
    save_results_path2 = '/net/vdesk/data2/evdijk/SMA/AMUSEproject-1/velocities_wrtplanet_in_hill_sphere_earthlike.png'
    make_plot_all_simulations_rocky(all_velocities, all_escape_vel_planet, save_results_path2, start_velocities, option='planet')

def make_velocity_plots_gas():
    explosion_energies = np.arange(1,10.1,0.3)*1e42|units.erg
    path = 'simulation_results/gas_results/'
    path = '/net/vdesk/data2/presa/AMUSEproject/simulation_results/gas_results/'

    all_relative_velocities = []
    all_escape_velocities = []
    all_velocities = []
    all_escape_vel_planet = []
    for explosion_energy in explosion_energies:
        # os.chdir(path+'{:.1e}.erg/'.format(explosion_energy.value_in(units.erg), 'g))
        os.chdir(path+'{}/'.format(explosion_energy,'g'))
        files = os.listdir()

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
                
            dm_part = read_set_from_file(f"dm_particles_{j}.hdf5", "hdf5")
            gas_part = read_set_from_file(f"gas_particles_{j}.hdf5", "hdf5")
            gravity_part = read_set_from_file(f"gravity_particles_{j}.hdf5", "hdf5")
            
            # calculate hill radius parameters and hill radii
            a = np.sqrt((dm_part[0].x-gravity_part[0].x)**2 + (dm_part[0].y-gravity_part[0].y)**2 + (dm_part[0].z-gravity_part[0].z)**2)
            a = a.value_in(units.m)
            
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
        save_result_path = f'/net/vdesk/data2/evdijk/SMA/AMUSEproject-1/simulation_results/velocities_in_hill_sphere_{explosion_energy.value_in(units.erg)}_erg.png'
        save_result_path2 = f'/net/vdesk/data2/evdijk/SMA/AMUSEproject-1/simulation_results/velocities_wrtplanet_in_hill_sphere_{explosion_energy.value_in(units.erg)}_erg.png'
        make_velocities_plot(relative_velocities_particles, escape_velocities_moon, save_result_path)
        make_velocities_plot(velocities_particles, escape_velocities_planet, save_result_path2, option='planet')

    save_results_path = 'velocities_in_hill_sphere_moon_gas.png'
    save_results_path2 = 'velocities_in_hill_sphere_planet_gas.png'
    save_results_path = '/net/vdesk/data2/evdijk/SMA/AMUSEproject-1/velocities_in_hill_sphere_jupiterlike.png'
    save_results_path2 = '/net/vdesk/data2/evdijk/SMA/AMUSEproject-1/velocities_wrtplanet_in_hill_sphere_jupiterlike.png'
    make_plot_all_simulations_gas(all_relative_velocities, all_escape_velocities, save_results_path, explosion_energies)
    make_plot_all_simulations_gas(all_velocities, all_escape_vel_planet, save_results_path2, explosion_energies, option='planet')
    

if __name__ == "__main__": 

    make_velocity_plots_rocky()
    make_velocity_plots_gas()
    