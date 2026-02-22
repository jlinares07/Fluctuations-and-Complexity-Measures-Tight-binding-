# Fluctuations-and-Complexity-Measures-Tight-binding-

This repository contains the reference implementation and numerical codes developed to reproduce the results presented in the manuscript titled *"Fluctuations, Complexity and Statistical Measures for Particles in a Tight-Binding Lattice"*.

You could find here the Python 3.9.25 scripts that we used to generate the for Figures 3 and 4 of the paper. 

For **Figure 3**, we present the system in a distributed configuration over the lattice sites satisfying -M ≤ l+m ≤ M (with M = 3). The relevant scryts are:

- 'dloc1.py': non-interacting system (U = 0.0)
- 'dloc2.py': interacting system (U = 100λ = 1.0)

These scripts output are the compressed Numpy archive files named 'Delocalized_v2(...).npz'

Alternatively, you could find the codes "locfinal0.py", "locfinal1.py" and "locfinalminus.py" which were used to build Figure 4. We model the two-particle system highly localized follwing that the particles are located in [-1,1] and [1,-1], respectively, in a (l,m) plane. Again, the corresponding output are those "Localized(...).npz".

Additionally, the "Guide npz files" serves as a brief guide to work with this file extension. All the results were plotted later in a new Python code, resulting in the Figures 3 and 4 that we presented in our work.

