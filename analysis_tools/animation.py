# to do:
# add options to plot gas, dm and gravity particles, and read that in more efficiently
# add options to more consistently work with units

from analysis_tools.plot_system import systemplotter

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

from amuse.io import read_set_from_file
from amuse.units import units

import numpy as np

import pickle

import os

class Animator:
    """
    A class that can be used for showing or creating an animation

    Attributes:

    path_simulation_results: path where the x and y data for the different frames are stored
    data should be stored in the format: 'gas_particles_{index}.hdf5' or 'dm_particles_{index}.hdf5'

    
    xlabel: str - None
    ylabel: str - None


    """
    def __init__(self, path_simulation_results, xlabel, ylabel, center='planet', xlim=0.05, ylim=0.05):
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xlim = (-xlim, xlim)
        self.ylim = (-ylim, ylim)
        self.data_path = path_simulation_results
        self.center = center
        self.n_frames = self.count_files()
    
    def count_files(self):
        # could add an extra check for the right filenames

        dir_path = self.data_path
        counter = 0

        for path in os.listdir(dir_path):
            if os.path.isfile(os.path.join(dir_path, path)):
                counter += 1

        return counter
    
    def animate(self, i):
        self.ax.clear()

        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)

        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

        filename_gas = self.data_path + f'gas_particles_{i}.hdf5'
        filename_dm = self.data_path + f'dm_particles_{i}.hdf5'
        filename_gravity = self.data_path + f'gravity_particles_{i}.hdf5'
        
        gas_particles = read_set_from_file(filename_gas)
        dm_particles = read_set_from_file(filename_dm)
        gravity_particles = read_set_from_file(filename_gravity)
        
        if self.center=="moon":
            gas_particles.position -= gravity_particles[0].position
            dm_particles.position -= gravity_particles[0].position
            gravity_particles[1].position -= gravity_particles[0].position
            gravity_particles[0].position = [0, 0, 0] | units.m
        # print(len(x_positions), len(y_positions))
        self.ax.scatter(dm_particles.x.value_in(units.AU), dm_particles.y.value_in(units.AU), color='darkred', label='planet core', s=10)
        self.ax.scatter(gas_particles.x.value_in(units.AU), gas_particles.y.value_in(units.AU), color='indianred', label='planet atmosphere/envelope', s=5)
        self.ax.scatter(gravity_particles.x.value_in(units.AU), gravity_particles.y.value_in(units.AU), color='blue', label='moon')
        self.ax.legend()

    def make_animation(self, show=False, save_path=None):
        self.fig, self.ax = plt.subplots()

        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        
        # self.animate(0)
        frames = np.arange(0, int(self.n_frames/3)-1)

        animation = FuncAnimation(self.fig, self.animate, frames)
        
        
        if show:
            self.fig.show()

        if save_path is not None:
            animation.save(save_path)


