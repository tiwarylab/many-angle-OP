# ManyAngle

This repository demonstrates how to use ManyAngle.cpp with urea simulation produced by Pablo M. Piaggi and Michele Parrinello in [this paper](https://www.pnas.org/content/115/41/10251). The simulation files including inputParameters.mdp, melt450K-B.gro, topol.top, and centers.dat are taken from Pablo Piaggi. For more details, please refer to his github or website, or email him directly.

# How to perform simulation
```
sh createTpr.sh // This will generate md.tpr
gmx_mpi mdrun -deffnm md -plumed plumed_many.dat -cpi md.cpt
```

The plumed_many.dat will compile ManyAngle.cpp and generate ManyAngle.o & ManyAngle.so. At the same time, this plumed file will also generate COLVAR and HILLS files with metadynamics adding bias as function of two angles calculated by the ManyAngle.cpp.
