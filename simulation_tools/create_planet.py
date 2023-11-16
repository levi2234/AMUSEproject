    
from amuse.datamodel import Particles
from amuse.units import units

def create_planet(planet_mass, planet_radius, **kwargs):
    """
    Creates a single planet particle with the given mass and radius at the center of the simulation.
    
    Args:
        planet_mass (float): Mass of the planet particle.
        planet_radius (float): Radius of the planet particle.
    
    Returns:
        planet (Particles): A single particle representing the planet.
    """
    #create the center planet particle
    planet = Particles(1)
    planet.mass = planet_mass
    planet.radius = planet_radius
    planet.x = 0 | planet_radius.unit
    planet.y = 0 | planet_radius.unit
    planet.z = 0 | planet_radius.unit
    planet.vx = 0 | units.m / units.s
    planet.vy = 0 | units.m / units.s
    planet.vz = 0 | units.m / units.s
    planet.name = 'planet'
    
    return planet