# copied from running_hydrodynamics.ipynb from the AMUSE tutorial
# changes made by: Esther van Dijk - inner region -> outer region

from amuse.units import units, constants


def inject_explosion_energy(gas_particles, 
                            explosion_energy=1.0e+51|units.erg,
                            exploding_region=10|units.RSun):
    outer = gas_particles.select(
        lambda pos: pos.length_squared() > exploding_region**2,
        ["position"])
    print(len(outer), "outermost particles selected.")
    print("Adding", explosion_energy / outer.total_mass(), "of explosion " \
        "(specific internal) energy to each of the n=", len(outer), "SPH particles.")
    outer.u += explosion_energy / outer.total_mass()
    return 