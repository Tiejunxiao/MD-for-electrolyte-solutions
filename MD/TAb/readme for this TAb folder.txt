(1)This folder contains input files for MD simulation of nonlinear triatomic solutes.
Three input files (CONFIG, CONTROL, FIELD ) are used for MD simulation.

In the cubic MD simulation box, there are  254 Cl-, 254 Na+, one nonlinear triatomic solute 1 with site charge numbers (z1,z2,z3) 
and a counter ion SA4 with  charge number z4=-(z1+z2+z3).
 Solute 1 consists of three sites which are tagged as SA1, SA2  and SA3.
The site separation distance for  two sites SAi  and SAj of solute 1 is Rij (i,j=1,2,3). 
We keep  R12=R23=2.5 Angstrom and R13= 4.330 Angstrom for all simulations.
So that all the simulations of this work share the same "CONTROL" file and "CONFIG" file.

The "FIELD" files  depend on site charges (z1,z2,z3).
The following keys are used for simulations.
key        (z1,z2,z3)
TA1b        (1,0,0)
TA2b        (0,1,0)
TA3b        (1,0,-1)
TA4b        (1,-1,1)
TA5b        (1,1,-1) 

When used for MD simulation, please copy  the file "FIELD-key"  to the file "FIELD".

(2)The fortran file "solvation-key.f90" is used to  evaluate mean  electric potential of the sites of the solute using the Ewald summation method.
A  "HISTORY" file is the input of the file "solvation-key.f90".

When the  "HISTORY" file is generated  from MD simulation with key= TA1b,TA2b,TA3b,TA4b,TA5b,
the  "HISTORY" file is used as input of  "solvation-key.f90".
