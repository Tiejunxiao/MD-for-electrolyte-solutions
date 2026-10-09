(1)This SA folder contains input files for MD simulation of spherical solutes.

Three input files (CONFIG, CONTROL, FIELD ) are used for MD simulation.

In the cubic MD simulation box, there are  255 Cl-, 255 Na+, one SA1+  and one SA2+ ion.
The Na+ and Cl- are taken as solvent species, while the SA1- and SA2+ are taken as solutes.

All the simulations share the same CONTROL file and CONFIG file.
The LJ energy parameter is fixed at 1.0KJ/mol. 
The LJ diameter \sigma_{so} for solute-solvent interaction is a tunable parameter.
The FIELD file depends on the  \sigma_{so} parameter and six different FIELD files are needed.
The following keys are used for simulations.
key          \sigma_{so} (Angstrom)     
SA2.5             2.5
SA3                3.0
SA3.5             3.5
SA4                4.0  
SA4.5             4.5             
SA5                5.0   
SA5.5             5.5

For  \sigma_{so}= 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5 Angstrom, the FIELD files are saved in the files
FIELD-SA2.5, FIELD-SA3, FIELD-SA3.5, FIELD-SA4, FIELD-SA4.5, FIELD-SA5 and FIELD-SA5.5, respectively.

When key=SA3.5, the radial distribution function and  the electrostatic energy of the solution reduces to the pure solvent case.

When used for MD simulation, please copy the file "FIELD-key"  to the file "FIELD".

(2)The fortran file "solvation.f90" is used to  evaluate mean  electric potential or equivalently the induced potential of the solute ions using the Ewald summation method.
A  "HISTORY" file is the input of the code "solvation.f90".
When the  "HISTORY" file is generated  from MD simulation with key= SA2.5,SA3,SA3.5,SA4,SA4.5,SA5, SA5.5,
the  "HISTORY" file is used as input of  "solvation.f90". 

