# Research Project: Analyzing the tidally disrupting dwarf galaxy Bootes III

This is a repository for my work in AST425, the astrophysics research course at 
U of T. My project is to study the dwarf galaxy Bootes III. Specifically, I 
find member stars and updated physical parameters of this galaxy, and use the
paramaters to run simulations for orbit and stream. The simulations are then 
compared with observation data on known streams that are postulated to be
associated with Bootes III (in particular, Styx stream and Typhon stream). My
project supervisor is Prof. Ting Li.

Relevant papers:\
[Boötes III is a Disrupting Dwarf Galaxy Associated with the Styx Stellar Stream](https://iopscience.iop.org/article/10.3847/1538-4357/aad8c1)\
[Four new stellar debris streams in the galactic halo](https://iopscience.iop.org/article/10.1088/0004-637X/693/2/1118)

The collection of Python notebooks and data files summarize my work so far. 
Although the data files are not uploaded to this repository, the notebooks
themselves contain code for the whole analysis procedure.

(README is not complete!)

## Table of Contents
* [Project outline](#project-outline)
* [Python packages](#python-packages)
* [Data](#data)
* [Additional notebooks](#additional)

## Project outline
The main Jupyter notebooks are separated into many parts, each covering 
a different stage of the analysis of Bootes III (Boo III). We find member 
stars of Boo III and use the members to compute kinematics and metallicity 
parameters of Boo III in Part 1. We then use the parameters to simulate 
the orbit of Boo III in Part 2. We also simulate a stellar stream 
given Boo 
III kinematics. These are compared with known stellar streams that are 
suspected to be associated with Boo III. In Part 3, we search for more 
member stars in a larger area from the centre of Boo III, using our stream 
simulation as a 
guide, with the aim of 
finding tidal tail features.

## Python packages
For orbit and stream modelling, we use 
[galpy](https://docs.galpy.org/en/v1.8.1/index.html)

## Data

## Additional notebooks



