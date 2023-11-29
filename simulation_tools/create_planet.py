"""This file contains the functions to create a planet.
The funciton returns a planet object with the following attributes:
The planet at its core and the atmosphere around it.

The planet is created as an AMUSE hydrodynamics particle set.

The function returns an AMUSE Particle Object with the planet and the atmosphere
""" 

from amuse.ext.star_to_sph import (pickle_stellar_model, convert_stellar_model_to_SPH,)
from amuse.test.amusetest import get_path_to_results
from amuse.community.mesa_r2208.interface import MESA
from amuse.datamodel import Particles

import os

def setup_stellar_evolution_model():
    out_pickle_file = os.path.join(get_path_to_results(), 
                                   "super_giant_stellar_structure.pkl")

    stellar_evolution = MESA(redirection="none")
    stars = Particles(1)
    stars.mass = 1.0 | units.MJupiter
    stellar_evolution.particles.add_particles(stars)
    stellar_evolution.commit_particles()

    print(
        "Evolving a MESA star with mass:",
        stellar_evolution.particles[0].mass
    )
    try:
        while stellar_evolution.model_time<0.12|units.Myr:
            stellar_evolution.evolve_model()
            print("star:", stellar_evolution.particles[0].stellar_type, stellar_evolution.model_time.in_(units.Myr))
    except AmuseException as ex:
        print("Evolved star to", stellar_evolution.particles[0].age)
        print("Radius:", stellar_evolution.particles[0].radius)
    
    pickle_stellar_model(stellar_evolution.particles[0], out_pickle_file)
    stellar_evolution.stop()
    return out_pickle_file

pickle_file = setup_stellar_evolution_model()
#pickle_file = './super_giant_stellar_structure.pkl'
print("Star generated.")



    