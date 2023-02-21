# Research Project: Analyzing the tidally disrupting dwarf galaxy Bootes III

This is a repository for my work in AST425, the astrophysics research course at 
U of T. My project is to study the dwarf galaxy Bootes III. The analysis is divided 
into the following parts: 
1. Identify member stars of Bootes III and measure its physical parameters
2. Orbit integration and tidal stream simulation of Bootes III. Particularly, we investigate 
the effects on the model result when we add the LMC into the Galactic potential model.
3. A further search of Bootes III members using the tidal stream model as a guide. For this
part, we want to search for tidal disruption features of Bootes III.
4. Compare with observation data on known streams that are postulated to be associated 
with Bootes III (in particular, Styx stream and Typhon stream). 

My project supervisor is Prof. Ting Li.

Relevant papers:\
[Boötes III is a Disrupting Dwarf Galaxy Associated with the Styx Stellar Stream](https://iopscience.iop.org/article/10.3847/1538-4357/aad8c1)\
[Four new stellar debris streams in the galactic halo](https://iopscience.iop.org/article/10.1088/0004-637X/693/2/1118)

The collection of Python notebooks and data files summarize my work so far. 
Although the data files are not uploaded to this repository, the notebooks
themselves contain code for the whole analysis procedure.

(README is not complete!)

## Table of Contents
* [General Description](#general-description)
* [Python packages](#python-packages)
* [Data](#data)

## General Description
The folder `main-notebooks` contains the Jupyter notebooks that accomplish the analysis tasks. Inside the folder, `Bootes III Part 1.ipynb` is about member identification and calculation of parameters. `Bootes III Part 2.ipynb` does the orbital integration and stream simulation (This includes all the figures about orbit and stream). `Bootes III Part 3.ipynb` searches for tidal tail features using Gaia DR3 catalogue, and uses the stream model from Part 2 as a reference. Additionally, the file `analysis_functions_v2.py` has custom functions for data reduction and analysis (such as proper motion data cutting and CMD filtering), and "plot_functions.py" is just a collection of helper functions that make plotting easier. These are required for the main notebooks to run successfully.

The main members of Bootes III (within 5 times half-light radius) are in the data files `boo3_main_members.dat` and `boo3_main_rrl.dat`. The first one contains 21 stars (not matched to the Gaia DR3 RR Lyrae catalogue) and the second one has 4 stars (matched to Gaia RRL). There are also many data files named `orbit_*.obj`, these are the stream simulation results from notebook Part 2.

The folder `additional` contains an `orbit_tutorial.ipynb` that gives a brief tutorial on using galpy to do orbit integration (which is what notebook Part 2 is about).

## Python packages
First, the commonly used packages in research: `numpy`, `scipy`, `matplotlib` and `astropy`. 

For orbit integration and stream simulation––[galpy](https://docs.galpy.org/en/v1.8.1/index.html).

## Data
The survey data used in this project are:\
[Southern Stellar Stream Spectroscopic Survey (S5)](https://s5collab.github.io)\
[Gaia DR3](https://www.cosmos.esa.int/web/gaia/data-release-3)\
[The Dark Energy Camera Legacy Survey (DECaLS)](https://www.legacysurvey.org/dr9/description/)
