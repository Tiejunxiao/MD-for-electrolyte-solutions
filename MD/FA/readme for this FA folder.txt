(1)This folder contains input files for MD simulation of tetraatomic solutes.
Three input files (CONFIG, CONTROL, FIELD ) are used for MD simulation.

In the cubic MD simulation box, there are  254 Cl-, 254 Na+, one linear tetraatomic solute 1 with site charge numbers (z1,z2,z3,z4) .
 Solute 1 consists of three sites which are tagged as SA1, SA2, SA3  and SA4.
The site separation distance for  two sites SAi  and SAj of solute 1 is Rij (i,j=1,2,3,4). 
We keep  z1+z2+z3+z4=0,  R12=R23=R34=2.5 Angstrom,   R13=R24= 5.0 Angstrom  and R14=7.5 Angstrom for  all simulations.
So that all the simulations of this work share the same "CONTROL" file and "CONFIG" file.

The "FIELD" files  depend on the solute site charges (z1,z2,z3,z4).
The following keys are used for simulations.
key        (z1,z2,z3,z4)
FA1        (-1,0,1,0)
FA2        (-1,0,0,1)
FA3        (-1,-1,1,1)
FA4        (-1,1,-1,1)

When used for MD simulation,  the file "FIELD-key"  is copied to the file "FIELD".

(2)The fortran file "solvation-key.f90" is used to  evaluate mean  electric potential of the sites of the solute using the Ewald summation method.
A  "HISTORY" file is the input of the file "solvation-key.f90".

When the  "HISTORY" file is generated  from MD simulation with key= FA1,FA2,FA3,FA4,
the  "HISTORY" file is used as input of  "solvation-key.f90".
