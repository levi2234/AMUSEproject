import numpy as np
from matplotlib import pyplot as plt
from amuse.datamodel import Particles
from amuse.units import units
from simulation_tools import create_planet_and_atmosphere
from simulation_tools import atmospheric_profiles as AP
from amuse.lab import nbody_system
from amuse.community.fi.interface import Fi
from amuse.community.hermite0.interface import Hermite
from amuse.couple import bridge


def plot(planet,atmosphere, savefig = 'simulation_results/simulation.png'):
    fig, ax = plt.subplots()
    planet_radius = planet.radius.value_in(units.m)
    circle = plt.Circle((0, 0), planet_radius[0], color='green',alpha= 0.2, fill=True)
    ax.scatter(atmosphere.x.value_in(units.m), atmosphere.y.value_in(units.m), s=1)
    ax.add_artist(circle)
    ax.set_aspect('equal')
    ax.set_xlabel('x ')
    ax.set_ylabel('y ')
    plt.savefig(savefig)

    

print("Creating planet and atmosphere")
planet, atmosphere = create_planet_and_atmosphere(planet_mass= 1e24 | units.kg , 
                                                  planet_radius= 6.371e6 | units.m , 
                                                  atmosphere_profile=AP.true_profile, 
                                                  atmosphere_height= 90|units.km , 
                                                  atmosphere_density=1.225 | units.kg * units.m** (-3) , 
                                                  num_points=5000,
                                                  gamma=1.4)

plot(planet,atmosphere, savefig='simulation_results/simulation_before.png')

print("Constructing hydro simulation")
converter = nbody_system.nbody_to_si(1e24 | units.kg, 6.371e6 | units.m)
hydro = Fi(converter, mode='openmp')
hydro.gas_particles.add_particles(atmosphere)

# Initialize the gravity code
gravity = Hermite(converter)
gravity.particles.add_particles(planet)

# Create a bridge
system_bridge = bridge.Bridge(use_threading=False)
system_bridge.add_system(gravity, (hydro,) )
#system_bridge.add_system(hydro, (gravity,) )

print(hydro.gas_particles)
# Define a function to prevent atmosphere particles from penetrating the planet
def prevent_atmosphere_penetration():
    for p in hydro.gas_particles:
        distance_from_planet_center = (p.position - planet.position).length()
        if distance_from_planet_center < planet.radius:
            # Handle the penetration (e.g., reflect the particle)
            # This is a placeholder; you'll need a physical model for this
            p.velocity = -p.velocity



print("Evolving the system")
# Evolving the system
end_time = 0.5 | units.day
while hydro.model_time < end_time:
    print("Time: ", hydro.model_time)
    dt = 0.05 | units.day  # Time step
    system_bridge.evolve_model(hydro.model_time + dt)
   
    prevent_atmosphere_penetration()
    hydro.model_time = hydro.model_time + dt
print("Done evolving the system")
hydro_particles = hydro.gas_particles.copy()
planet = gravity.particles.copy()
hydro.stop()
gravity.stop()

plot(planet,hydro_particles, savefig='simulation_results/simulation_after.png')
