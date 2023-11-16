import numpy as np
from amuse.datamodel import Particles
from amuse.units import units
import matplotlib.pyplot as plt
from simulation_tools import create_planet_and_atmosphere
from simulation_tools import atmospheric_profiles as AP

planet, atmosphere = create_planet_and_atmosphere(planet_mass= 1e24 | units.kg , planet_radius= 6.371e6 | units.m , atmosphere_profile=AP.true_profile, atmosphere_height= 90|units.km , atmosphere_density=1.225 | units.kg * units.m** (-3) , gamma=1.4)

merged = Particles()
merged.add_particles(planet)
merged.add_particles(atmosphere)

#print(sum(atmosphere.mass.value_in(units.kg))) #should be around 5,15 x 10^18 kg
#scale atmosphere mass to sum of atmosphere mass
atmosphere.mass = atmosphere.mass *  5.15e18/sum(atmosphere.mass.value_in(units.kg))
#systemplotter(atmosphere, "x", "y").plot(savefig="planet_and_atmosphere.png")
print(sum(atmosphere.mass.value_in(units.kg)))
fig, ax = plt.subplots()
planet_radius = planet.radius.value_in(units.m)

circle = plt.Circle((0, 0), planet_radius[0], color='green',alpha= 0.2, fill=True)
ax.scatter(atmosphere.x.value_in(units.m), atmosphere.y.value_in(units.m), s=1)

#plot circle with radius of planet

ax.add_artist(circle)
ax.set_aspect('equal')
ax.set_xlabel('x ')
ax.set_ylabel('y ')

plt.savefig('atmosphere.png')
