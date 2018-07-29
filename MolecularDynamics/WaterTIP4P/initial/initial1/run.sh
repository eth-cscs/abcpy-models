#!/bin/bash 

module load GROMACS

#export PLUMED_KERNEL="/users/fbrotz/software/plumed2-ves_FM_daint/src/lib/libplumedKernel.so"

#module swap PrgEnv-cray/6.0.4  PrgEnv-gnu/6.0.4
#module load CrayGNU/17.08
#module load GSL/2.4-CrayGNU-17.08
#module load cray-libsci/17.06.1

#source /users/fbrotz/software/gromacs-5.1.4_FM_daint/bin/GMXRC


#Solvation, preparing box
gmx solvate -cs /users/duttar/abcpy/initial/initial2/tip4p.gro  -o topol.pdb -box 2.5 2.5 2.5 -p empty.top
rm -f \#empty.top.1*
mv empty.top topol.top

#Energy minimization
gmx grompp -c topol.pdb -f em.mdp -p topol.top -o em.tpr
gmx mdrun -deffnm em -nt 1
rm -f em.trr mdout.mdp em.tpr em.edr em.log 

#NPT
gmx grompp -f npt.mdp -c em.gro -p topol.top -o npt.tpr -maxwarn 2
rm -f mdout.mdp em.mdp em.gro *itp
gmx mdrun -deffnm npt
rm -f npt.trr npt.cpt npt_prev.cpt npt.log npt.edr npt.gro npt.mdp *pdb *top

#Radial distribution function OO.Output:rdf_OO.xvg
gmx rdf -f npt.xtc -s npt.tpr -n /users/duttar/abcpy/initial/initial2/index.ndx -b 20 -e 100 -o rdf_OO.xvg << EOF
3
3
EOF
#Radial distribution function OH
#take the fit at x interval of 0.1-0.7. output: rdf_OH.xvg
gmx rdf -f npt.xtc -s npt.tpr -n /users/duttar/abcpy/initial/initial2/index.ndx -b 20 -e 100 -o rdf_OH.xvg << EOF
3
4
EOF

#Mean square displacement for calculation of Diffusion constant, 
#and for low value of MSD saving. Output: msd.xvg
echo 3 3 |gmx msd -f npt.xtc -s npt.tpr -b 20 -e 100 -n /users/duttar/abcpy/initial/initial2/index.ndx  

