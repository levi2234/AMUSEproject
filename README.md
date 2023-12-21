# 1. Astronomical Simulation Project


This repository contains the code and documentation for an astronomical simulation project that aims to simulate the ejection of a planet's atmosphere and the subsequent capture of that atmosphere by the planet's moon. The project is being conducted as part of the fulfillment of the AMUSE (Astronomy Multi-purpose Software Environment) course requirements.

## 1.1. Project Team
- Levi van Es - s2115409  
- Andres Presa - s3643751
- Esther van Dijk - s2034174

## 1.2. Table of Contents
- [1. Astronomical Simulation Project](#1-astronomical-simulation-project)
  - [1.1. Project Team](#11-project-team)
  - [1.2. Table of Contents](#12-table-of-contents)
  - [1.3. Introduction](#13-introduction)
  - [1.4. Simulation Details](#14-simulation-details)
  - [1.5. Installation](#15-installation)
  - [1.6. Project Setup](#16-project-setup)
    - [1.6.1. Degrees of Freedom](#161-degrees-of-freedom)
    - [1.6.2. Constraints](#162-constraints)
  - [1.7. Simulation Design](#17-simulation-design)
    - [1.7.1. Hydrodynamics](#171-hydrodynamics)
    - [1.7.2. Gravitational Interaction](#172-gravitational-interaction)
    - [1.7.3. Hydro-Gravity Bridge](#173-hydro-gravity-bridge)
  - [1.8. Implementing the Simulation](#18-implementing-the-simulation)
    - [1.8.1. Simulating the MESA model](#181-simulating-the-mesa-model)
    - [1.8.2. Importing the MESA model into AMUSE](#182-importing-the-mesa-model-into-amuse)
    - [1.8.3. Making the Moon and the Planet](#183-making-the-moon-and-the-planet)
    - [1.8.4. Setting up the Hydrodynamics](#184-setting-up-the-hydrodynamics)
    - [1.8.5. Setting up the Gravity](#185-setting-up-the-gravity)
    - [1.8.6. Coupling the Hydrodynamics and Gravity](#186-coupling-the-hydrodynamics-and-gravity)
    - [1.8.7. Running the Simulation](#187-running-the-simulation)
  - [1.9. Exploring different scenarios](#19-exploring-different-scenarios)
  - [1.10. Results](#110-results)

## 1.3. Introduction

The goal of this project is to develop an astronomical simulation that models the fascinating phenomenon of a planet's atmosphere being ejected and subsequently captured by its moon. Through this simulation, we aim to explore the dynamics and consequences of such an event in a celestial system and the likelihood of such an event occurring for a given set of parameters. We will evaluate a set of 2 systems, one discribinga a gas giant and its moon, and the other a rocky planet and its moon. These two situations will be simulated using the AMUSE framework, a versatile tool for astrophysical simulations. 

## 1.4. Simulation Details



## 1.5. Prerequisites

The following software is required to run the simulation:
- [AMUSE](https://amusecode.github.io/) (version 13.0 or higher)
- [MESA](http://mesa.sourceforge.net/) (version 10398 or higher)
- [Python](https://www.python.org/) (version 3.7 or higher)
- [Matplotlib](https://matplotlib.org/) (version 3.3.4 or higher)
- [Numpy](https://numpy.org/) (version 1.20.1 or higher)

## 1.6. Running the Simulation
After git cloning the repository, the project contains two main files aimed to run the gas_giant simulation and the rocky_planet simulation. The project is segmented into sections each adressing a different step in the complete project. 

Initially the simulations will need to be run. The running of these simulations is devided up into two different files since their approaches vary greatly. Running the gas_giant simulation is done by running the _run_simulation_gas.py_ file. This file incorporates the MESA model and the AMUSE simulation. Running this simulation for different  Running the rocky_planet simulation is done by running the _run_simulation_rocky.py_ file. This file incorporates the AMUSE simulation and the MESA model. The top of the file contains the parameters that can be changed to run the simulation for different scenarios.  The gas giant simulation is run as follows:

```bash
    python3 run_simulation_gas.py
```	

The rocky planet can be run in a simular fashion. However this approach requires the user to set the atmospheric radial velocity. Changing this value is different than the gas giant simulation and requires to set this velocity in m/s as forllows

```bash
    python3 run_simulation_rocky.py 12550
```
Here the 12550 is the radial velocity in m/s.

## 1.7 Post Simulation Analysis
After running the simulation, the results are stored in the _results_ folder. The results are stored in a .hdf5 file. From this location all the following analyses will be done.

This brings us to the _analyis tools_ folder which contains the following files:

-  _animation.py_  - This file creates an animation of the simulation. The animation is stored alongside the simulation folders in _simulation_results_.
- AtmosphereProfilePlot.py - This file creates a atmosphere identical to the atmosphere of the rocky planet simulation. This atmosphere is then plotted in a density vs radius plot. This plot is stored in the _simulation_results_ folder under the velocity its name. This plot aims as a check to see if the atmosphere is created correctly.
- _binding_energy.py_ - This file calculates the binding energy of the atmosphere and the moon. The binding energy is calculated for each timestep and stored in a new column in the .hdf5 file of the specific timestep. In this column a boolean is stored indicating wether a particle is bound to the moon or not, and similarly for the planet and the atmsophere.
- _collisiondetector.py_ - This file contains the collision detector function. This function detects wether the line between two particle locations intersects with a sphere. This is used to calculate wether a particle has collided with the moon between or in a timestep. 
- collisionscan.py - This file puts the _collsiondetector.py_ to use and scans the simulation for collisions. It stores the results in a .csv file in the _simulation_results_ folder under the velocity its name.
- _collisions_probabalistic.py_ - Andres fill this please with a discription
- _hill_sphere.py_ - Esther fill this please with a discription
- ParticleProfileAnim.py - This file creates an animation of the particles in a simulation. Changing the parameters in the file allows for different animations. The animation returns a videofile showing the density distribution and internal energy distribution of the particles. This file is stored in the _simulation_results_ folder under the velocity its name.
- _plot_system.py_ - This file contains the code to plot a AMUSE particle set. It builds on the matplotlib library and is used to plot the system at a given timestep. This is mainly for visualisation during the development of the simulation.

## Instruction typical simulation

1. Run the simulation

```bash
    python3 run_simulation_gas.py
```
or 
```bash
    python3 run_simulation_rocky.py 12550
```

2. Analyse the simulation

Run any of the following files to analyse the simulation

```bash
    python3 analysis_tools/animation.py
    python3 analysis_tools/binding_energy.py
    python3 analysis_tools/collisionscan.py
    python3 analysis_tools/collisions_probabalistic.py
    python3 analysis_tools/hill_sphere.py
    python3 analysis_tools/plot_system.py
```
IMPORTANT! it is important to run this from the main directory and not from the analysis_tools directory. The analysis tools are designed to be run from the main directory because of their relative imports. 



