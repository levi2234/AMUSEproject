
import numpy as np
from amuse.datamodel import Particles
from amuse.units import units   
from amuse.units import constants
from simulation_tools import atmospheric_profiles 

#------------------- Sample Density Profiles ------------------------#
def generate_points_in_shell(num_points, inner_radius, outer_radius, units=True,**kwargs):
    """
    Generates random points in a 3D shell using Rejection Sampling
    
    Args:
        num_points (int): The number of points to generate.
        inner_radius : The inner radius of the shell
        outer_radius : The outer radius of the shell
    
    Returns:
        points (array): A 2D array of shape (num_points, 3) containing the generated points. (units are in same units as inner_radius)
    """
    try:
        stripped_outer_radius = outer_radius.value_in(inner_radius.unit)
    except:
        print('outer_radius units could not be converted to inner_radius units')
        
    stripped_inner_radius = inner_radius.value_in(inner_radius.unit)
    
    points = np.random.uniform(-stripped_outer_radius, stripped_outer_radius, size=(num_points, 3)) 
    
    # Calculate the distance of each point from the origin
    distances = np.linalg.norm(points, axis=1)
    
    # Select only the points that are in the shell
    points = points[(distances >= stripped_inner_radius) & (distances <= stripped_outer_radius)]
    
    # If we don't have enough points, recurse
    if len(points) < num_points:
        points = np.vstack((points, generate_points_in_shell(num_points - len(points), inner_radius, outer_radius, units = False)))
    
    # If we have too many points, trim
    if len(points) > num_points:
        points = points[:num_points]

    if units == False:
        return points
    else:
       
        return points | inner_radius.unit

def create_atmosphere_from_profile(profile=None, **kwargs):
    """
    Creates a particle atmosphere with a given density profile.

    Args:
        profile (callable): A function that takes a distance and returns a density. (distance given is from the center of the planet (change if you think it is more usefull to make it from surface))
        num_points (int): The number of points to generate in the atmosphere.
        inner_radius (float): The inner radius of the spherical shell in which to generate points.(default is earth radius(m))
        outer_radius (float): The outer radius of the spherical shell in which to generate points. (default is 1e9 (m))
        **kwargs: Additional keyword arguments to pass to the density profile function.

    Returns:
        atmosphere (Particles): A particle set representing the atmosphere.
        
    Notes: 
        Default values are for Earth.
    """
    
    num_points = kwargs.get('num_points')
    inner_radius = kwargs.get('inner_radius')
    outer_radius = kwargs.get('outer_radius')
    
    #converting units
    if inner_radius is not None:
        inner_radius = inner_radius.in_(units.m)
    if outer_radius is not None:
        outer_radius = outer_radius.in_(units.m)
    
    
    
    
    if num_points is None:
        print('num_points is not given for create_atmosphere_from_profile (default is 10000)')
        num_points = 10000
        
    if inner_radius is None:
        raise ValueError('inner_radius is not given for create_atmosphere_from_profile (earth radius is 6.371e6 m))')
    
    if outer_radius is None:
        print('outer_radius is not given for create_atmosphere_from_profile (default is 20km above earth surface (1e9 m) )')
        outer_radius = inner_radius + (2e4 | units.m).value_in(inner_radius.unit)
        
    # Generate random points in a spherical shell
    
    points = generate_points_in_shell(num_points, inner_radius , outer_radius).value_in(units.m)
    
    atmosphere = Particles(len(points))
    # Calculate the density at each point
    if profile is None:
        print('profile is not given for create_atmosphere_from_profile (default is uniform density)')
        density = np.ones(num_points)
    else:
        point_distances = np.linalg.norm(points, axis=1) | units.m
        density = profile(point_distances, **kwargs)["density"] #this returns the density at each point (array length num_points)
        #atmosphere mass is the density at each point multiplied by the volume of the point
        atmosphere.mass = density * (4/3) * np.pi * (outer_radius**3 - inner_radius**3)
        
        temperature = profile(point_distances, **kwargs)["temperature"]
        mmw = profile(point_distances, **kwargs)["mmw"]
        

        
        #calculate internal energy from temperature and mmw
        if temperature.unit == units.K and mmw.unit == units.g / units.mol:
            v = 1 | units.mol
            internal_energy = (3/2) * ((1/mmw) * (1 | units.mol**(-1))  * constants.kB * temperature)
        else:
            internal_energy = None
        
        
    
    # Create the atmosphere (check units!)
    
    atmosphere.x = points[:, 0] | units.m
    atmosphere.y = points[:, 1] | units.m
    atmosphere.z = points[:, 2] | units.m
    atmosphere.vx = 0. | units.m / units.s
    atmosphere.vy = 0. | units.m / units.s
    atmosphere.vz = 0. | units.m / units.s
    atmosphere.u = internal_energy


    #add name label to the particles
    atmosphere.name = 'atmosphere'
    
    return atmosphere


#uncomment to test the function

#results = create_atmosphere_from_profile(profile=true_profile, num_points=10000, inner_radius=6.371e6 | units.m, outer_radius=6.461e6 | units.m, rho0=1.225 | units.kg * units.m**(-3),gamma=1.4)
#results = create_atmosphere_from_profile(profile= exponential_profile,inner_radius=6.371e6, outer_radius=1e9, rho0=1.225, gamma=1.4)

# print(results)
# #plot the results

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(results.x.value_in(units.m), results.y.value_in(units.m), results.z.value_in(units.m), s=1)
# ax.set_xlabel('x ')
# ax.set_ylabel('y ')
# ax.set_zlabel('z ')
# plt.savefig('atmosphere.png')
# #plot the density profile

# plt.figure()
# positions = results.position.value_in(units.m)
# plt.plot(np.linalg.norm(positions, axis=1), results.mass.value_in(results.mass.unit), '.')
# plt.xlabel('Distance from center (m)')
# #plot temperature profile
# #plt.plot(np.linalg.norm(positions, axis=1), results.U.value_in(results.U.unit), '.')
# plt.savefig('density_profile.png')


