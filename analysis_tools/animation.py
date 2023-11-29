from analysis_tools.plot_system import systemplotter

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

from amuse.io import read_set_from_file

import numpy as np

import pickle

import os

class Animator:
    """
    A class that can be used for showing or creating an animation

    Attributes:

    path_simulation_results: path where the x and y data for the different frames are stored
    data should be stored in the format: 'x_positions_{index}' or 'y_positions_{index}'

    
    xlabel: str - None
    ylabel: str - None


    """
    def __init__(self, path_simulation_results, xlabel, ylabel):
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.data_path = path_simulation_results
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
        filename_gas = self.data_path + f'gas_particles_{i}.hdf5'
        filename_dm = self.data_path + f'dm_particles_{i}.hdf5'
        
        print(filename_gas)
        
        gas_particles = read_set_from_file(filename_gas)
        dm_particles = read_set_from_file(filename_dm)

        x_positions = gas_particles.x.append(dm_particles.x)
        y_positions = gas_particles.y.append(dm_particles.y)

        self.ax.scatter(x_positions, y_positions)

    def make_animation(self, show=False, save_path=None):
        self.fig, self.ax = plt.subplots()

        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

        frames = np.arange(0, int(self.n_frames/2))

        animation = FuncAnimation(self.fig, self.animate, frames)

        if show:
            self.fig.show()

        if save_path is not None:
            animation.save(save_path)


