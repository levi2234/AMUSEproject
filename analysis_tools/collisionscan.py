import collisiondetector as cd
import numpy as np
from amuse.io import read_set_from_file
from amuse.units import units
import os
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor

def collisionscan(velocity, radius):
    print("Collisionscan for velocity: ", velocity)
    #CHANGE THESE values depending on the simulation
    moon_radius =radius #in meters
    

    folder_location = f'simulation_results/rocky_results/Earth_planet_{velocity}ms'

    number_of_files = len(os.listdir(folder_location))/3

    # Read gas_particles file for each frame

    collisions = []
    for frame in range(int(number_of_files)-1):
        print(f"Velocity: {velocity},  Snapshot:", frame)
        #read atmosphere particles
        atmosphere_before = read_set_from_file(f'{folder_location}/gas_particles_{frame}.hdf5', 'hdf5')
        atmosphere_after = read_set_from_file(f'{folder_location}/gas_particles_{frame+1}.hdf5', 'hdf5')
        
        #read moon particles
        moon_before = read_set_from_file(f'{folder_location}/gravity_particles_{frame}.hdf5', 'hdf5')[0]
        moon_after = read_set_from_file(f'{folder_location}/gravity_particles_{frame+1}.hdf5', 'hdf5')[0]

    #calculate average moon location
        average_moon_location = (moon_before.position.value_in(units.m) + moon_after.position.value_in(units.m)) /2
        
        #calculate average moon velocity
        collision_bool = [cd.collision_detector(atmosphere_before[i].position.value_in(units.m), atmosphere_after[i].position.value_in(units.m), average_moon_location, moon_radius) for i in range(len(atmosphere_before))]
        collisions.append(collision_bool)
        print("Number of collisions: ", np.sum(collision_bool))
        # print("Collisionbool: ", collision_bool)

    #add the velocity and number of collisions for each timestep to collisions.txt in a row as follows: velocity number_of_collisions_timestep1 number_of_collisions_timestep2 etc.
    collisions = np.array(collisions)
    number_of_collisions = np.sum(collisions, axis=1)


    with open(f'simulation_results/{radius}m.csv', 'a') as f:
        f.write(f'{velocity},{np.sum(number_of_collisions)},{np.sum(number_of_collisions)/len(collisions[0])},{number_of_collisions}\n')
    



    
#make a parallelized version of the code above

for radius in np.linspace(1737.1e3, 384567780,10):
    #make a csv file for each radius
    if not os.path.isfile(f'simulation_results/{radius}m.csv'):
        with open(f'simulation_results/{radius}m.csv', 'w') as f:
            f.write(f'Velocity,Number of collisions total,Fraction of collisions,Number of collisions per timestep\n')
    
    velocities = [11800,12000,12100,12350, 12375, 12400, 12450,12475,12500,15250,12550,12575,12600,12700,12900,113200]
    with ProcessPoolExecutor(max_workers=6) as executor:
        executor.map(collisionscan,velocities,radius*np.ones(len(velocities)))

    

    
    
    
# print(collisions)
# print("Total number of collisions: ", np.sum(collisions))
# print("Collsion percentage: ", np.sum(collisions)/len(collisions[0]))
# collisions = np.array(collisions)
# number_of_collisions = np.sum(collisions, axis=1)
# plt.plot(number_of_collisions)
# plt.xlabel('Frame number')
# plt.ylabel('Number of collisions')
# plt.title('Number of collisions per frame')
# plt.savefig('number_of_collisions.png')
# plt.show()