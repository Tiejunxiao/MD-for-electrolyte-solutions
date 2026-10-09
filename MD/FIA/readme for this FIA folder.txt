(1)This folder contains input files for MD simulation of tetraatomic solutes.
Three input files (CONFIG, CONTROL, FIELD ) are used for MD simulation.

In the cubic MD simulation box, there are  253 Cl-, 253 Na+, one nonlinear tetraatomic solute 1 with site charge numbers (z1,z2,z3,z4,z5) and a counter solute 2 with charge number z6.
 Solute 1 consists of three sites which are tagged as SA1, SA2, SA3, SA4  and SA5.
The site separation distance for  two sites SAi  and SAj of solute 1 is Rij (i,j=1,2,3,4,5). 
We keep  z1+z2+z3+z4 +z5 +z6=0,  R12=R13=R14=R15=3.0 Angstrom,   R23=R24=R25=R34=R35=R45= 4.899 Angstrom for  all simulations.
So that all the simulations of this work share the same "CONTROL" file and "CONFIG" file.

The "FIELD" files  depend on the solute site charges (z1,z2,z3,z4,z5).
The following keys are used for simulations.
key        (z1,z2,z3,z4,z5)
FIA1        (-1,0,0,0,0)
FIA2        (0,-1,0,0,0)
FIA3        (0,-1,1,0,0)
FIA4        (-1,-1,1,1,0)
FIA5        (0,-1,1,-1,1)

When used for MD simulation,  the file "FIELD-key"  is copied to the file "FIELD".

(2)The fortran file "solvation-key.f90" is used to  evaluate mean  electric potential of the sites of the solute using the Ewald summation method.
A  "HISTORY" file is the input of the file "solvation-key.f90".

When the  "HISTORY" file is generated  from MD simulation with key= FIA1,FIA2,FIA3,FIA4,FIA5,
the  "HISTORY" file is used as input of  "solvation-key.f90".
