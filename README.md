# Fluctuations-and-Complexity-Measures-Tight-binding-

This repository contains the reference implementation and codes built to present our results in "Fluctuations, Complexity and Statistical Measures for Particles in a Tight-Binding Lattice".

You could find here the Python 3.9.25 files that we used to build Figures 3 and 4. 

In Figure 3, we present the system in a distributed way, following that -M ≤ l+m ≤ M, with M=3. The codes used to obtain that are named "dloc1.py" and "dloc2.py", such codes corresponds to without interaction (U=0.0) and with interaction (U=100λ=1.0), respectively. The output are those Delocalized_v2(...).npz files (Numpy files).
Alternatively, you could find the codes "locfinal0.py", "locfinal1.py" and "locfinalminus.py" which were used to build Figure 4. We model the two-particle system highly localized follwing that the particles are located in [-1,1] and [1,-1], respectively, in a (l,m) plane. Again, the corresponding output are those "Localized(...).npz".

Additionally, the "Guide npz files" serves as a brief guide to work with this file extension. All the results were plotted later in a new Python code, resulting in the Figures 3 and 4 that we presented in our work.

