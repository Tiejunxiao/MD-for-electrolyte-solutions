(1)This  folder contains input files for MD simulation of a NaCl solution (pure solvent).
In the cubic MD simulation box, there are  256 Cl- and 256 Na+.
Three input files (CONFIG, CONTROL, FIELD ) are used for MD simulation.

(2) The fortran code "cavity.f90" is used to generate the cavity distribution function C(R).
 The "HISTORY" file from MD simulation of  the NaCl solution  is used as input of code  "cavity.f90". 

C(R)  represents the probability of ﬁnding a cavity of radius R and is further used to evaluate the cavity exclusion probability P(R) via P(R)= 1- \int_0^R C(r)dr.
Then the cavity formation energy mu^{cav}(R)= - kb*T*lnP(R)  can be evaluated as a function of the cavity radius R. 
[refer to: Fernando Bresme and Jose Alejandre,Cavities in ionic liquids , J Chem. Phys. 118,4134,2003；
 Cavities in Molecular Liquids and the Theory of Hydrophobic Solubilities, A. Pohorille and L. R. Pratt, J. Am. Chem. Soc. 112,5066,1990 ] 


