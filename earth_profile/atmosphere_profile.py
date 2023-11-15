import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units,constants

def get_atmospheric_profile(n_points=50):
    '''
    Returns a dictionary with the Earth's atmospheric profile
        Inputs:
            - n_points: number of points to interpolate the profile.
                        Default 50 

        Returns:
            Dictionary with the following keys:
            - 'altitude': radial profile of atmosphere in units of km
            - 'density': density profile of atmosphere in units of kg*m**-3
            - 'temperature': temperature profile of atmosphere in units of K
            - 'mmw': mean molecular weight in units of gmol-1
    '''
    altitude=np.arange(0,90,5)#|units.km
    temperature=np.array([288.150,255.650,223.150,216.650,216.650,221.650,226.650,237.050,251.050,
                      265.050,270.650,259.450,245.450,231.450,217.450,206.650,196.650,186.946])#|units.K
    density=np.array([1.225,0.736116,0.412707,0.193674,0.088,0.0395,0.018,0.0082,0.0039,0.0019,0.000977,0.000537,0.000288,0.000149,
                  0.0000742,0.0000349,0.0000157,0.00000677])#| units.kg * units.m**(-3)
    
    new_altitude = np.linspace(0, 85, n_points)
    mean_molecular_weight=28.96*np.ones(len(new_altitude))|units.g * units.mol
    new_temp=np.interp(new_altitude,altitude,temperature)|units.K
    new_density=np.interp(new_altitude,altitude,density)|units.kg * units.m**(-3)

    profile={'altitude':new_altitude|units.km,'density':new_density,
             'temperature':new_temp,'mmw':mean_molecular_weight}
    return profile