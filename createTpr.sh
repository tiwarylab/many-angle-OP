# Remember to module load plumed first.

/opt/packages/gromacs-5.1.4-plumed/bin/gmx_mpi grompp -f inputParameters.mdp -c melt450K-B.gro -p topol.top -o md.tpr
