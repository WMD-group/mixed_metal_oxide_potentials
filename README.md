mixed metal oxide interatomic potentials
============================

A python script that takes GULP input file and calculates the effective k2, k4, C6 and BSM values for oxygen in mixed-metal oxides and plots 
Buckingham, Lennard-Jones, Coulomb potentials, as well as the total potential on a single graph.

Requirements
------------
Python

Current Status
------------
- Reads specified gulp input file 
- Defaults to GULP.gin 
- Prints out all effective oxygen values as well as all potentials and spring constants which can be copied into GUP input file
- Plots the Buckingham potential, Lennard-Jones and corresponding Coulombic interaction on same graph along with the total potential for the specified interaction
- Saves plots as .eps files, named after the atoms involved in that potential

Execution 
------------
python metal_oxide_mix.py -f file_name

Short-term goals
------------
- To allow for other metals, not just Zn, Sn and In

Disclaimer
----------
This file is not affiliated with *GULP*. Feel free to use and modify, but do so at your own risk.
