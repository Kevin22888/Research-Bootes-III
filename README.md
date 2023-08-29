# Analysis of the dwarf galaxy Bootes III: a system likely undergoing tidal disruption

This is a repository for my work in AST425, the astrophysics research course at 
U of T. My supervisor is Professor Ting Li.

The analysis is divided into the following parts: 
1. Finding member stars of Bootes III. Calculating its systemic parameters.
2. Simulating orbit and tidal stream of Bootes III. Investigate the changes in the
model results when we add the LMC into the Galactic potential model.
3. Finding tidal tail feature of Bootes III. Searching for evidence of connection to
the Styx stellar stream.
4. Compare Typhon stream to Bootes III, in angular momentum (actions) and energy
parameter spaces. 

Relevant papers:\
[Boötes III is a Disrupting Dwarf Galaxy Associated with the Styx Stellar Stream](https://iopscience.iop.org/article/10.3847/1538-4357/aad8c1)\
[Four new stellar debris streams in the galactic halo](https://iopscience.iop.org/article/10.1088/0004-637X/693/2/1118)


## Table of Contents
* [General Description](#general-description)
* [Python packages](#python-packages)
* [Data](#data)

## General Description
This section needs update.

The folder `main-notebooks` contains the Jupyter notebooks that accomplish the analysis tasks. Inside the folder, `Bootes III Part 1.ipynb` is about member identification and calculation of parameters. `Bootes III Part 2.ipynb` does the orbital integration and stream simulation (This includes all the figures about orbit and stream). `Bootes III Part 3.ipynb` searches for tidal tail features using Gaia DR3 catalogue, and uses the stream model from Part 2 as a reference. Additionally, the file `analysis_functions_v2.py` has custom functions for data reduction and analysis (such as proper motion data cutting and CMD filtering), and "plot_functions.py" is just a collection of helper functions that make plotting easier. These are required for the main notebooks to run successfully.

The main members of Bootes III (within 5 times half-light radius) are in the data files `boo3_main_members.dat` and `boo3_main_rrl.dat`. The first one contains 21 stars (not matched to the Gaia DR3 RR Lyrae catalogue) and the second one has 4 stars (matched to Gaia RRL). There are also many data files named `orbit_*.obj`, these are the stream simulation results from notebook Part 2. There are 14 of those files, and they numbered by pairs (because leading and trailing arms of the stream), from 1 to 7. The numbers refer to the different changes in the simulation input parameters. See below for the description. The first pair, numbered by 1, are the fiducial parameters used. See the paper for description of these default parameters (Paper in progress!). Then, 
each of the results from numbers 2 to 7 only modified one single value from the default parameters in number 1. The following list describes the different changes in the leading arm stream data (the files are named `orbit_leading*.obj`. Their trailing arm counterparts also have their changes.
1. `orbit_leading1.obj`: No change, just the fiducial parameters
2. `orbit_leading2.obj`: LMC mass x2
3. `orbit_leading3.obj`: LMC rhm x2
4. `orbit_leading4.obj`: MW halo mass x2
5. `orbit_leading5.obj`: v_phi = 252 km/s (instead of v_phi = 232)
6. `orbit_leading6.obj`: without LMC
7. `orbit_leading7.obj`: without LMC and with vphi = 252

The folder `additional` contains an `orbit_tutorial.ipynb` that gives a brief tutorial on using galpy to do orbit integration (which is what notebook Part 2 is about).

## Python packages
First, the commonly used packages in research: `numpy`, `scipy`, `matplotlib` and `astropy`. 

For orbit integration and stream simulation––[galpy](https://docs.galpy.org/en/v1.8.1/index.html).

## Data
The survey data used in this project are:\
[Southern Stellar Stream Spectroscopic Survey (S5)](https://s5collab.github.io)\
[Gaia DR3](https://www.cosmos.esa.int/web/gaia/data-release-3)\
[The Dark Energy Camera Legacy Survey (DECaLS)](https://www.legacysurvey.org/dr9/description/)
