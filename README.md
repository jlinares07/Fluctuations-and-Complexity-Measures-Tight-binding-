# Fluctuations-and-Complexity-Measures-Tight-binding-

This repository contains the reference implementation and numerical codes developed to reproduce the results presented in the manuscript titled *"Fluctuations, Complexity and Statistical Measures for Particles in a Tight-Binding Lattice"*.

You could find here the Python 3.9.25 scripts that we used to generate the Figures 3 and 4 of the paper. 

For **Figure 3**, we present the system in a distributed configuration over the lattice sites satisfying -M ≤ l+m ≤ M (with M = 3). The relevant scrits are:

- 'dloc1.py': non-interacting system (U = 0.0)
- 'dloc2.py': interacting system (U = 100λ = 1.0)

These scripts output are the compressed Numpy archive files named 'Delocalized_v2(...).npz'

Alternatively, for **Figure 4**, the scripts model a two-particle system highly localized where, in the plane (l,m), the particles are located in (-1,1) and (1,-1). The scripts are:

- 'locfinal0.py': non-interacting system (U = 0.0)
- 'locfinal1.py': interacting system (U = 100λ = 1.0)
- 'locfinalminus.py' : interacting system (U = -100λ = -1.0)

Additionally, the repository includes a brief guide "Guide npz files" explaining how to work with the .npz file format. All simulation outputs were subsequently visualized and processed in separate plotting scripts to produce the final versions of Figures 3 and 4 as shown in the publication.

