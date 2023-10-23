# Astronomical Simulation Project


This repository contains the code and documentation for an astronomical simulation project that aims to simulate the ejection of a planet's atmosphere and the subsequent capture of that atmosphere by the planet's moon. The project is being conducted as part of the fulfillment of the AMUSE (Astronomy Multi-purpose Software Environment) course requirements.

## Project Team
- Levi van Es - s2115409
- Andres Presa - s3643751
- Esther

## Table of Contents
- [Introduction](#introduction)
- [Simulation Details](#simulation-details)
- [Installation](#installation)
- [Project Setup](#Project-Setup)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)

## Introduction

The goal of this project is to develop an astronomical simulation that models the fascinating phenomenon of a planet's atmosphere being ejected and subsequently captured by its moon. Through this simulation, we aim to explore the dynamics and consequences of such an event in a celestial system and the likelihood of such an event occurring for a given set of parameters. The number 

## Simulation Details

Our simulation is based on the AMUSE framework, a versatile tool for astrophysical simulations. It uses numerical methods and principles of astrophysics to model the complex interactions between celestial bodies. Here are some key details of the simulation:

- **Programming Language**: Python
- **Simulation Environment**: AMUSE
- **Simulation Scenarios**: Ejection of the planet's atmosphere and capture by the moon
- **Parameters**: [List important parameters here]

## Installation

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

## Project Setup

Before starting the construction of the simulation, we first conducted some research to understand the phenomenon of a planet's atmosphere being ejected and captured by its moon. We also explored the AMUSE framework and its capabilities to determine the best approach to simulate the phenomenon. 

#### Degrees of Freedom
A good way to frame the problem is to list the free parameters of the simulation. These are the parameters that can be varied to explore the dynamics of the system and the consequences of the phenomenon. Here are the free parameters of our simulation:

1. **Mass of the planet**: The mass of the planet from which the atmosphere is ejected
2. **Mass of the moon**: The mass of the moon that captures the ejected atmosphere
3. **Distance between the planet and the moon**: The distance between the planet and the moon at the time of ejection
4. **Velocity of the ejected atmosphere**: The velocity of the ejected atmosphere at the time of ejection
5. **Mass of the ejected atmosphere**: The mass of the ejected atmosphere
6. **Location of the ejected atmosphere**: The location of the ejected atmosphere at the time of ejection w.r.t the planet and the moon
7. **Atmospheric parameters**: The atmospheric parameters of the ejected atmosphere (e.g. temperature, pressure, etc.)

Looking at these parameters we immediately notice that there is a lot of freedom in the simulation. However, because there are so many parameters, it is not feasible to explore the entire parameter space. Therefore, we will focus on a few key parameters and explore the dynamics of the system for different values of these parameters. To do this we will have to constrain the other parameters.

#### Constraints
To constrain the parameters of the simulation, we will use the following assumptions:

1.
2.
3.
4.
5.
6.


