from plot_system import systemplotter

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

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
        filename_x = self.data_path + f'x_positions_{i}'
        with open(filename_x, 'r') as f:
            x_positions = pickle.load(f)

        filename_y = self.data_path + f'y_positions_{i}'
        with open(filename_y, 'r') as f:
            y_positions = pickle.load(f)

        self.ax.scatter(x_positions, y_positions)

    def plot_animation(self, show=False, save_path=None):
        self.fig, self.ax = plt.subplots()

        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

        frames = np.arange(0, self.n_frames)

        animation = FuncAnimation(self.fig, self.animate, frames)

        if show:
            self.fig.show()

        if save_path is not None:
            animation.save(save_path)


