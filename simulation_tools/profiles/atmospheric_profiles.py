import numpy as np
from amuse.units import units

#------------------- Atmospheric Density Profiles -------------------#
def barometric_profile(r, **kwargs):
    """

    
    A barometric density profile.
    
    Args:
        r (float): The distance from the center of the planet.
        r0 (float): The radius at which the density is rho0.
        rho0 (float): The density at r0.
        gamma (float): The adiabatic index.
    
    Returns:
        rho (float): The density at r.
    """
    gamma = kwargs.get('gamma') #adiabatic index
    r0 = kwargs.get('inner_radius') #inner radius of the planet
    rho0 = kwargs.get('rho0') #density at r0
    
    
    
    if r0 is None:
        raise ValueError('r0 is not given for barometric profile (earth radius is 6.371e6 m))')
    if rho0 is None:
        raise ValueError('rho0 is not given for barometric profile (earth density is 1.225 kg/m^3))')
    if gamma is None:
        raise ValueError('gamma is not given for barometric profile (earth gamma is 1.4))')
    
    density = rho0 * (r0 / r)**gamma
    temperature = 288.15 * (r0 / r)**gamma | units.K
    mmw = 28.96 * (r0 / r)**gamma | units.g / units.mol
    values = {'density': density, 'altitude': r, "temperature":temperature,"mmw":mmw}

    
    return values

def exponential_profile(r,  **kwargs):
    """
    An exponential density profile.
    
    Args:
        r (float): The distance from the center of the planet.
        r0 (float): The radius at which the density is rho0.
        rho0 (float): The density at r0.
        gamma (float): The adiabatic index.
    
    Returns:
        rho (float): The density at r.
    """
    #check if r0 is given and alert if not
    r0 = kwargs.get('inner_radius') #inner radius of the planet
    rho0 = kwargs.get('rho0') #density at r0
    #raise error if r0 is not given
    if r0 is None:
        raise ValueError('r0 is not given for exponential profile (earth radius is 6.371e6 m))')
    #raise error if rho0 is not given
    if rho0 is None:
        raise ValueError('rho0 is not given for exponential profile (earth density is 1.225 kg/m^3))')
    density = rho0 * np.exp(-(r - r0) / r0)
    temperature = 288.15 * np.exp(-(r - r0) / r0)
    mmw = 28.96 * np.exp(-(r - r0) / r0)
    values = {'density': density, 'altitude': r, "temperature":temperature,"mmw":None}

    return values


def true_profile(r, **kwargs):
    """
    Returns a dictionary containing atmospheric profile data for a planet.

    Parameters:
    -----------
    r : astropy.units.quantity.Quantity
        Distance from the center of the planet.
    **kwargs : dict
        Keyword arguments that can be passed to the function. The following
        keyword argument is used:
            - inner_radius : astropy.units.quantity.Quantity
                Inner radius of the planet.

    Returns:
    --------
    profile : dict
        A dictionary containing the following keys:
            - altitude : astropy.units.quantity.Quantity
                Altitude above the planet's surface.
            - density : astropy.units.quantity.Quantity
                Density of the atmosphere at the given altitude.
            - temperature : astropy.units.quantity.Quantity
                Temperature of the atmosphere at the given altitude.
            - mmw : astropy.units.quantity.Quantity
                Mean molecular weight of the atmosphere.
    """
    
    r = r.in_(units.m)
    r0 = kwargs.get('inner_radius') #inner radius of the planet
    
    r = r - r0
    
    r_stripped = r.value_in(units.m)
    altitude=np.arange(0,90000,5000)#|units.m
    temperature=np.array([288.150,255.650,223.150,216.650,216.650,221.650,226.650,237.050,251.050,
                      265.050,270.650,259.450,245.450,231.450,217.450,206.650,196.650,186.946])#|units.K
    density=np.array([1.225,0.736116,0.412707,0.193674,0.088,0.0395,0.018,0.0082,0.0039,0.0019,0.000977,0.000537,0.000288,0.000149,
                  0.0000742,0.0000349,0.0000157,0.00000677])#| units.kg * units.m**(-3)
    
    
    mean_molecular_weight=28.96*np.ones(len(r_stripped))|units.g / units.mol
    new_temp=np.interp(r_stripped,altitude,temperature)|units.K
    new_density=np.interp(r_stripped,altitude,density)|units.kg * units.m**(-3)

    profile={'altitude':r,'density':new_density,
             'temperature':new_temp,'mmw':mean_molecular_weight}
    return profile

def uniform_profile(r,**kwargs):
    """
    A uniform density profile.
    
    Args:
        r (float): The distance from the center of the planet.
        r0 (float): The radius at which the density is rho0.
        rho0 (float): The density at r0.
        gamma (float): The adiabatic index.
    
    Returns:
        rho (float): The density at r.
    """
    density = 1.225| units.kg * units.m**(-3)
    temperature = 288.150 * np.ones(len(r)) | units.K 
    mmw = 28.96* np.ones(len(r)) | units.g / units.mol
    values = {'density': density, 'altitude': r, "temperature":temperature,"mmw":mmw}

    return values