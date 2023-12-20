import numpy as np
import matplotlib.pyplot as plt
import os

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
    mass = (particle1.mass.value_in(units.kg) + particle2.mass.value_in(units.kg))
    distance_ = distance(particle1, particle2) 
    G = constants.G.value_in(units.m**3 / units.kg / units.s**2)
    v_esc = np.sqrt(2*G*mass/distance_)
    return v_esc

def distance(object1, object2, unit=units.m):
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
    
    radii_moon = gravity_part.radius # 1000 | units.km # [1000, 2000, 3000, 4000] | units.km
    if np.sum(mask)!=0:
        hill_sphere_moon = gas_particles.select(
            lambda in_hill_sphere: in_hill_sphere==1, ['in_hill_moon'])
        escape_velocities_moon = v_esc(hill_sphere_moon, gravity_part)
        speeds_particles = speed(hill_sphere_moon)
        v_rel_particles = relative_speed(hill_sphere_moon, gravity_part[0])
        print(np.sum(escape_velocities_moon>speeds_particles))
        print(v_rel_particles)
        print(escape_velocities_moon)
        crossing_times_moon = radii_moon.value_in(units.m) / speeds_particles
        crossing_times_moon = (crossing_times_moon | units.s).value_in(units.hour)
        print('crossing time moon:', crossing_times_moon)
        # print(speeds_particles)
        # print(escape_velocities_moon)
        return v_rel_particles, escape_velocities_moon
    return None, None

def make_velocities_plot(relative_velocities, escape_velocities, path):
    total_relative_velocities = [v_rel[i] for v_rel in relative_velocities for i in range(len(v_rel))]   
    total_escape_velocities = [v_esc[i] for v_esc in escape_velocities for i in range(len(v_esc))]

    fig, ax = plt.subplots()
    ax.scatter(total_relative_velocities, total_escape_velocities)
    ax.set_xlabel('Relative velocity (m/s)')
    ax.set_ylabel('Escape velocity moon (m/s)')
    fig.savefig(path)

if __name__ == "__main__": 
    
    #navigate to the directory where the h5 files is located relative to this file
    os.chdir("../simulation_results/jupiterlike_planet") #change the folder to the folder where the h5 file is located
    files = os.listdir() #list all files in the directory


    #get the range of numbers in the file names
    file_numbers = []
    for file in files:
        file_numbers.append(int(file.split("_")[2].split(".")[0]))
    file_numbers = np.sort(file_numbers)
    relative_velocities_particles = []
    escape_velocities_moon = []
    for j in file_numbers:
            
        dm_part = read_set_from_file(f"dm_particles_{j}.hdf5", "hdf5")
        gas_part = read_set_from_file(f"gas_particles_{j}.hdf5", "hdf5")
        gravity_part = read_set_from_file(f"gravity_particles_{j}.hdf5", "hdf5")
        

        #calculate hill radius parameters and hill radii
        a = np.sqrt((dm_part[0].x-gravity_part[0].x)**2 + (dm_part[0].y-gravity_part[0].y)**2 + (dm_part[0].z-gravity_part[0].z)**2)
        a = a.value_in(units.m)
        
        M = dm_part[0].mass.value_in(units.kg) + gravity_part[0].mass.value_in(units.kg)
        eccentricity = 0
        rh_planet = hill_radius(dm_part[0].mass.value_in(units.kg), M, a, eccentricity) #(beware to set the eccentricity to an appropriate value)
        rh_moon = hill_radius(gravity_part[0].mass.value_in(units.kg), M, a, eccentricity) #(beware to set the eccentricity to an appropriate value)
        
        v_rel_particles, v_esc_moon = add_hill_sphere_attributes(gas_part, rh_planet, rh_moon)
        if v_rel_particles is not None:
            relative_velocities_particles.append(v_rel_particles)
            escape_velocities_moon.append(v_esc_moon)
        write_set_to_file(gas_part, f"gas_particles_{j}.hdf5", "hdf5", overwrite_file=True)
    
    save_result_path = '../velocities_in_hill_sphere.png'
    make_velocities_plot(relative_velocities_particles, escape_velocities_moon, save_result_path)


    
    
    
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
    
    