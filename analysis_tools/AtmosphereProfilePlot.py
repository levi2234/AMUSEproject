from amuse.io import read_set_from_file
from amuse.units import units
import numpy as np
import matplotlib.pyplot as plt

from simulation_tools import create_planet_and_atmosphere as CPA
from simulation_tools.profiles import atmospheric_profiles as AP

planet_mass = 1 | units.MEarth
planet_radius = 1 | units.REarth
atmosphere_profile = AP.true_profile
atmosphere_height = 100|units.km
atmosphere_density = 1.225 | units.kg * units.m** (-3)
num_points = 50000
gamma = 1.4


planet, atmosphere = CPA.create_planet_and_atmosphere(planet_mass=planet_mass , 
                                                planet_radius= planet_radius, 
                                                atmosphere_profile=atmosphere_profile, 
                                                atmosphere_height= atmosphere_height, 
                                                atmosphere_density=atmosphere_density, 
                                                num_points=num_points,
                                                gamma=gamma)
atmosphere.radius = np.sqrt(atmosphere.x**2 + atmosphere.y**2 + atmosphere.z**2) 

#calculate density from mass and radius of each shell
atmosphere.density = atmosphere.mass / (4/3 * np.pi * atmosphere.radius**3)

# Plot density vs radius
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.scatter(atmosphere.density.value_in(units.kg/units.m**3), atmosphere.radius.value_in(units.km)- planet_radius.value_in(units.km), s=0.1)
plt.xlabel('Density (kg/m^3)')
plt.ylabel('Height (km)')
plt.title('Density vs Radius')
plt.grid(True)  # Add grid lines

# Customize the plot appearance
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
# plt.xlim(0, atmosphere.density.value_in(units.kg/units.m**3).max() + 10)
# plt.ylim(0, atmosphere.radius.value_in(units.km).max() + 10)

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
# plt.xlim(0, atmosphere.u.value_in(atmosphere.u.unit).max() + 10)
# plt.ylim(0, atmosphere.radius.value_in(units.km).max() + 10)

plt.tight_layout()
plt.savefig('density_temperature_vs_radius.png')
