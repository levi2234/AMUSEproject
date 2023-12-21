from amuse.io import read_set_from_file
from amuse.units import units
import numpy as np
import matplotlib.pyplot as plt

from simulation_tools import create_planet_and_atmosphere as CPA
from simulation_tools.profiles import atmospheric_profiles as AP
from amuse.io import read_set_from_file
from amuse.units import units
import os
import os
from matplotlib.animation import FuncAnimation

#import datafile
def update(frame):
    # Read gas_particles file for each frame
    filename = f'simulation_results/jupiterlike_planet/gas_particles_{frame}.hdf5'
    atmosphere = read_set_from_file(filename, 'hdf5')

    planet_radius = 1 | units.REarth
    #calculate density from mass and radius of each shell
    atmosphere.density = atmosphere.mass / (4/3 * np.pi * atmosphere.radius**3)
    atmosphere.radius = np.sqrt(atmosphere.x**2 + atmosphere.y**2 + atmosphere.z**2) 

    # Clear the previous plot
    plt.clf()

    # Plot density vs radius
    plt.subplot(1, 2, 1)
    plt.scatter(atmosphere.density.value_in(units.kg/units.m**3), atmosphere.radius.value_in(units.km)- planet_radius.value_in(units.km), s=0.1)
    plt.xlabel('Density (kg/m^3)')
    plt.ylabel('Height (km)')
    plt.title('Density vs Radius')
    plt.grid(True)  # Add grid lines

    # Customize the plot appearance
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)

    # Plot temperature vs radius
    plt.subplot(1, 2, 2)
    plt.scatter(atmosphere.u.value_in(atmosphere.u.unit), atmosphere.radius.value_in(units.km) - planet_radius.value_in(units.km), s=0.1)
    plt.xlabel('Internal Energy (J) (per unit mass)')
    plt.ylabel('Height (km)')
    plt.title('Internal Energy vs Radius \n (Temperature profile Analog)')
    plt.grid(True)  # Add grid lines

    # Customize the plot appearance
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)

    plt.tight_layout()

# Create the animation
animation = FuncAnimation(plt.gcf(), update, frames=range(20), interval=300)

# Save the animation as a GIF
animation.save('density_temperature_vs_radius_animation_jup.mp4', writer='imagemagick')

