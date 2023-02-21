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
* [Outline](#outline)
* [Python packages](#python-packages)
* [Data](#data)
* [Additional notebooks](#additional)

## Outline
The folder `main-notebooks` contains the Jupyter notebooks that accomplish the analysis tasks. Inside the folder, `Bootes III Part 1.ipynb` is about member identification and calculation of parameters. `Bootes III Part 2.ipynb` does the orbital integration and stream simulation (This includes all the figures about orbit and stream). `Bootes III Part 3.ipynb` searches for tidal tail features using Gaia DR3 catalogue, and uses the stream model from Part 2 as a reference. Additionally, the file `analysis_functions_v2.py` has custom functions for data reduction and analysis (such as proper motion data cutting and CMD filtering), and "plot_functions.py" is just a collection of helper functions that make plotting easier. These are required for the main notebooks to run successfully.

## Python packages
For orbit and stream modelling, we use 
[galpy](https://docs.galpy.org/en/v1.8.1/index.html)

## Data

## Additional notebooks



