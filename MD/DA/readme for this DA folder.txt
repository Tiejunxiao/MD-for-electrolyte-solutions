(1)This folder contains input files for MD simulation of diatomic solutes.
Three input files (CONFIG, CONTROL, FIELD ) are used for MD simulation.
All the simulations share the same CONTROL file.

In the cubic MD simulation box, there are  254 Cl-, 254 Na+, a diatomic solute 1 with site charges (z1,z2) 
and a diatomic solute 2 with site charges (z3,z4).
 Solute 1 consists of site SA1 and SA2, while  solute 2 consists of site SA3 and SA4.
The site separation distance for  two sites of solute 1 is R12. 
The site separation distance for two sites of solute 2 R34. 
We keep  R12=R34, z3=-z1 and z4=-z2 in the simulations.

The "FIELD" files and "CONFIG" files depend on R12 and (z1,z2).
The following keys are used for simulations.
key        R12(Angstrom)     (z1,z2)
DA1       1.0                            (-1,0)
DA2        2.0                           (-1,0)
DA3        3.0                          (-1,0)
DA4        4.0                          (-1,0)
DA4b      4.0                          (-1,1)
DA5         5.0                          (-1,0) 
DA5b       5.0                         (-1,1)

When used for MD simulation,  the content of the file "FIELD-key"  is copied to the file "FIELD" and the content of the file "CONFIG-key"  is copied to the file "CONFIG".

(2)The fortran file "solvation-DA.f90" and  "solvation-DAb.f90" are used to  evaluate mean  electric potential of the sites of the diatomic solute using the Ewald summation method.
A  "HISTORY" file is the input of the file "solvation-DA.f90" and "solvation-DAb.f90" .

When the  "HISTORY" file is generated  from MD simulation with key= DA1,DA2,DA3,DA4,DA5, (wherein (z1=-1,z2=0))
the  "HISTORY" file is used as input of  "solvation-DA.f90".

When the  "HISTORY" file is generated  from MD simulation with key= DA4b,DA5b, (wherein (z1=-1,z2=1))
the  "HISTORY" file is used as input of  "solvation-DAb.f90".