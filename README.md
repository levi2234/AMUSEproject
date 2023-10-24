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

The goal of this project is to develop an astronomical simulation that models the fascinating phenomenon of a planet's atmosphere being ejected and subsequently captured by its moon. Through this simulation, we aim to explore the dynamics and consequences of such an event in a celestial system and the likelihood of such an event occurring for a given set of parameters. The number 

## 1.4. Simulation Details

Our simulation is based on the AMUSE framework, a versatile tool for astrophysical simulations. It uses numerical methods and principles of astrophysics to model the complex interactions between celestial bodies. Here are some key details of the simulation:

- **Programming Language**: Python
- **Simulation Environment**: AMUSE
- **Simulation Scenarios**: Ejection of the planet's atmosphere and capture by the moon
- **Parameters**: [List important parameters here]

## 1.5. Installation

To setup all dependencies for your venv, follow these steps:

1. Clone this repository to your local machine:

   ```bash
   git clone https://github.com/your-username/your-repo.git

2. Install the required dependencies (recommended to use a virtual environment)):

   ```bash
   pip install -r requirements.txt
   
3. Load the AMUSE module (*working with version --2023.5.1--* ):

   ```bash
   module load AMUSE/2023.5.1

4. startup the AMUSE environment:

## 1.6. Project Setup

Before starting the construction of the simulation, we first conducted some research to understand the phenomenon of a planet's atmosphere being ejected and captured by its moon. We also explored the AMUSE framework and its capabilities to determine the best approach to simulate the phenomenon. 

### 1.6.1. Degrees of Freedom
A good way to frame the problem is to list the free parameters of the simulation. These are the parameters that can be varied to explore the dynamics of the system and the consequences of the phenomenon. Here are the free parameters of our simulation:

1. **Mass of the planet**: The mass of the planet from which the atmosphere is ejected
2. **Mass of the moon**: The mass of the moon that captures the ejected atmosphere
3. **Distance between the planet and the moon**: The distance between the planet and the moon at the time of ejection
4. **Velocity of the ejected atmosphere**: The velocity of the ejected atmosphere at the time of ejection
5. **Mass of the ejected atmosphere**: The mass of the ejected atmosphere
6. **Location of the ejected atmosphere**: The location of the ejected atmosphere at the time of ejection w.r.t the planet and the moon
7. **Atmospheric parameters**: The atmospheric parameters of the ejected atmosphere (e.g. temperature, pressure, etc.)

Looking at these parameters we immediately notice that there is a lot of freedom in the simulation. However, because there are so many parameters, it is not feasible to explore the entire parameter space. Therefore, we will focus on a few key parameters and explore the dynamics of the system for different values of these parameters. To do this we will have to constrain the other parameters.

### 1.6.2. Constraints
To constrain the parameters of the simulation, we will use the following assumptions:

1. Constrain the mass of the planet to be similar to that of Earth
2. Constrain the mass of the moon to be similar to that of Earth's moon
3. Constrain the distance between the planet and the moon to be similar to that of Earth and its moon
4. Constrain the mass of the ejected atmosphere through the physics of the phenomenon
5. Constrain the location of the ejected atmosphere to be uniform across the surface of the planet (for initial tests)


## 1.7. Simulation Design
With the aforementioned constraints, we can now design the simulation. The situation requires two types of physics. 

1. **Hydrodynamics**: The hydrodynamics of the ejected atmosphere and its interaction with the planets surface (optionally also the moon's surface)
2. **Gravitational**: The gravitational interaction between the planet, the moon, and the ejected atmosphere
3. **Hydro-Gravity Bridge**: The coupling of the hydrodynamics and gravitational interaction

These two types of physics can be simulated separately. The hydrodynamics can be simulated using the hydrodynamics module of AMUSE. The gravitational interaction can be simulated using the gravitational module of AMUSE. However, because the hydrodynamics and gravitational interaction are coupled, we will have to combine the two simulations. This can be done by using the bridge module of AMUSE. The bridge module allows us to combine two simulations and exchange information between them. 

This kind of coupling is not trivial and requires some experimentation. Therefore, we will first simulate the hydrodynamics of the ejected atmosphere and its interaction with the planet's surface. Once we have a working simulation, we will add the gravitational interaction with the moon. 

### 1.7.1. Hydrodynamics
The hydrodynamics of the ejected atmosphere and its interaction with the planet's surface can be simulated using the hydrodynamics module of AMUSE. The hydrodynamics module allows us to simulate the hydrodynamics of a fluid and its interaction with a solid surface. The module uses the Smoothed Particle Hydrodynamics (SPH) method to simulate the fluid. To get the hydrodynamics module to work, we will have to define the following:

1. **The fluid**: The fluid that is being simulated (in this case the ejected atmosphere)
2. **The solid**: The solid surface that the fluid interacts with (in this case the planet's surface)
3. **The initial conditions**: The initial conditions of the fluid and the solid
4. **The hydrodynamics parameters**: The hydrodynamics parameters of the fluid and the solid
5. **The hydrodynamics simulation**: The hydrodynamics simulation that simulates the fluid and the solid
6. **The hydrodynamics bridge**: The hydrodynamics bridge that couples the fluid and the solid (check if this is necessary and we can just include the solid in the atmosphere simulation)

Applying the aforementioned steps to our simulation, we get the following:
1. **The fluid**: to construct the fluid atmosphere there are various ways to do this. The first way to do this is to use the SPH module of AMUSE to construct the fluid atmosphere. The SPH module allows us to construct a fluid atmosphere using the SPH method. The SPH method uses particles to represent the fluid. The particles are distributed across the fluid and have various properties such as mass, position, velocity, etc. The SPH module allows us to construct the fluid atmosphere by defining the atmosphere parameters and atmosphere particle locations. Simulating these however is difficult and is best done with the use of a pre-existing functions from AMUSE ``star_to_sph()`` and ``convert_stellar_model_to_sph()``. By using a simulation from MESA we are able to convert a MESA stellar model object to an amuse SPH object which then can be used in our HYDRO simulations.<br>


2. **The Solid**: the before mentioned functions ``star_to_sph()`` and ``convert_stellar_model_to_sph()``, do not only export the atmosphere of a model but also allow us to convert a solid center core  MESA stellar model object to an hydro SPH object representation in AMUSE . This object can then be used to represent the solid surface of the planet. <br>
<!--A crude way to construct the central planet is to use the ``new_sph_particle()`` function from AMUSE to create a central SPH particle. This function allows us to construct a particle with the desired properties. We can then this particle in combination with the aforementioned atmosphere particles to construct the planet with atmosphere.-->

3. **Initial Conditions** The initial conditions as mentioned before can be included in the MESA simulation and from then on are integraded in the final MESA solution <br>

4. **Hydrodynamics Parameters** The hydrodynamics parameters are the parameters that define the hydrodynamics of the fluid and the solid. These parameters include the viscosity, the density, the pressure, etc. These parameters can be obtained from the MESA stellar model object. <br>






### 1.7.2. Gravitational Interaction
The gravitational interaction between the planet, the moon, and the ejected atmosphere can be simulated using the gravitational module of AMUSE. The gravitational module allows us to simulate the gravitational interaction between celestial bodies. The module uses the N-body method to simulate the gravitational interaction. To get the gravitational module to work, we will have to define the following:

1. **The celestial bodies**: The celestial bodies that are being simulated (in this case the planet, the moon, and the ejected atmosphere)
2. **The initial conditions**: The initial conditions of the celestial bodies i.e. Masses, positions, velocities, etc.

Applying the aforementioned steps to our simulation, we get the following:

1. **The celestial bodies**: The celestial bodies in this problem are just the earth, moon and host star. The earth and moon can be represented as SPH particles and the host star can be represented as a point particle. Gravitational between these gravity and SPH particles is solved in a later step through the Hydro-Gravity bridge<br>

2. **The initial conditions**: The initial conditions of the planet can be obtained from the MESA stellar model object. The initial conditions of the moon can be constructed by hand through the use of the ``new_sph_particle()`` function from AMUSE. This particle is initialized to have parameters (position,velocity, mass, radius) similar to that of a/the moon. 

### 1.7.3. Hydro-Gravity Bridge

The Hydro-Gravity Bridge is the magic that combines N-body gravitaitonal effects with the fluid-like behaviour of the atmosphere. The bridge works by coupling the two simulations and exchanging information between them. The bridge allows us to simulate the hydrodynamics of the ejected atmosphere and its interaction with the planet's surface and the gravitational interaction between the planet, the moon, and the ejected atmosphere. To get the bridge to work, we will have to define the following:

1. **The bridge**: The bridge that couples the hydrodynamics and gravitational interaction
2. **The coupling parameters**: The coupling parameters that define the coupling between the hydrodynamics and gravitational interaction

Applying the aforementioned steps to our simulation, we get the following:

1. **The bridge**: The bridge can be constructed using the ``bridge()`` function from AMUSE. This function allows us to construct a bridge between two simulations. The function takes the two simulations as input and returns a bridge object. The bridge object can then be used to couple the two simulations. The direction this bridge works in is very important for the runtime and carefull thought needs to go into the significance of the direction of the bridge. Since the atmosphere likely has a marginal impact on gravitational effects it is likely not needed to couple the hydro to the gravity. Because the gravity likely has a large effect on the behaviour of the hydro code this bridge direction is very significant and should not be omitted <br>

2. **The coupling parameters**: The coupling parameters are the parameters that define the coupling between the hydrodynamics and gravitational interaction. These parameters include the coupling time, the coupling radius, etc. These parameters are less trivial and some thought needs to go into these parameters to ensure the coupling is done correctly and for the code to even run.

<!-- Say something about how we chose these coupling parameters>  -->

## 1.8. Implementing the Simulation

### 1.8.1. Simulating the MESA model
<!--Insert mesa model making-->

### 1.8.2. Importing the MESA model into AMUSE
<!--Insert mesa model importing-->

### 1.8.3. Making the Moon and the Planet
<!--Insert moon and planet making-->

### 1.8.4. Setting up the Hydrodynamics
<!--Insert hydrodynamics setup-->

### 1.8.5. Setting up the Gravity
<!--Insert gravity setup-->

### 1.8.6. Coupling the Hydrodynamics and Gravity
<!--Insert coupling setup-->

### 1.8.7. Running the Simulation
<!--Insert simulation running-->

## 1.9. Exploring different scenarios
<!--Insert different scenarios-->

## 1.10. Results
<!--Insert results-->
