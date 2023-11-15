
import numpy as np
from amuse.datamodel import Particles


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
    gamma = kwargs.get('gamma')
    r0 = kwargs.get('inner_radius') #inner radius of the planet
    rho0 = kwargs.get('rho0') #density at r0
    
    if r0 is None:
        raise ValueError('r0 is not given for barometric profile (earth radius is 6.371e6 m))')
    if rho0 is None:
        raise ValueError('rho0 is not given for barometric profile (earth density is 1.225 kg/m^3))')
    if gamma is None:
        raise ValueError('gamma is not given for barometric profile (earth gamma is 1.4))')
    
    return rho0 * (r0 / r)**gamma

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
    
    return rho0 * np.exp(-(r - r0) / r0)


#------------------- Sample Density Profiles ------------------------#
def generate_points_in_shell(num_points, inner_radius, outer_radius, **kwargs):
    """
    Generates random points in a 3D shell using Rejection Sampling
    
    Args:
        num_points (int): The number of points to generate.
        inner_radius (float): The inner radius of the shell (remember the units!)
        outer_radius (float): The outer radius of the shell. (remember the units!)
    
    Returns:
        points (array): A 2D array of shape (num_points, 3) containing the generated points.
    """
    
    # Generate random points in a cube
    points = np.random.uniform(-outer_radius, outer_radius, size=(num_points, 3))
    
    # Calculate the distance of each point from the origin
    distances = np.linalg.norm(points, axis=1)
    
    # Select only the points that are in the shell
    points = points[(distances >= inner_radius) & (distances <= outer_radius)]
    
    # If we don't have enough points, recurse
    if len(points) < num_points:
        points = np.vstack((points, generate_points_in_shell(num_points - len(points), inner_radius, outer_radius)))
    
    # If we have too many points, trim
    if len(points) > num_points:
        points = points[:num_points]
    
    return points

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
    
    
    if num_points is None:
        print('num_points is not given for create_atmosphere_from_profile (default is 10000)')
        num_points = 10000
        
    if inner_radius is None:
        raise ValueError('inner_radius is not given for create_atmosphere_from_profile (earth radius is 6.371e6 m))')
    
    if outer_radius is None:
        print('outer_radius is not given for create_atmosphere_from_profile (default is 20km above earth surface (1e9 m) )')
        outer_radius = inner_radius + 2e4
        
    # Generate random points in a spherical shell
    points = generate_points_in_shell(num_points, inner_radius, outer_radius)
    
    # Calculate the density at each point
    if profile is None:
        print('profile is not given for create_atmosphere_from_profile (default is uniform density)')
        density = np.ones(num_points)
    else:
        point_distances = np.linalg.norm(points, axis=1)
        density = profile(point_distances, **kwargs) #this returns the density at each point (array length num_points)
    
    # Create the atmosphere (check units!)
    atmosphere = Particles(len(points))
    atmosphere.x = points[:, 0] 
    atmosphere.y = points[:, 1] 
    atmosphere.z = points[:, 2] 
    atmosphere.mass = (density * (outer_radius - inner_radius) / num_points)  #the masses of the particles are the density times the volume of the shell divided by the number of particles (double check)
    
    #add name label to the particles
    atmosphere.name = 'atmosphere'
    
    return atmosphere


#uncomment to test the function

#results = create_atmosphere_from_profile(profile=exponential_profile, num_points=10000, inner_radius=6.371e6, outer_radius=1e9, rho0=1.225,gamma=1.4)
# results = create_atmosphere_from_profile(profile= exponential_profile,inner_radius=6.371e6, outer_radius=1e9, rho0=1.225, gamma=1.4)

# print(results)
# #plot the results
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(results.x, results.y, results.z, s=1)
# ax.set_xlabel('x ')
# ax.set_ylabel('y ')
# ax.set_zlabel('z ')
# plt.savefig("test.jpg")


