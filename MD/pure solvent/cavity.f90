!CC*************************************************************************************************************************
!CC*   Copyright (c) 2026 [Xueyu Song and Tiejun Xiao]. 
!CC*  This code is distributed under the terms of the Creative Commons Attribution 4.0 International  License (CC BY 4.0).
!CC*  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
!CC*
!CC*  Author: Xueyu Song
!CC*  Email:   [xsong@iastate.edu]
!CC*  Author: Tiejun Xiao
!CC*  Email:   [tjxiao@gznc.edu.cn]
!CC*  Date:    2026-10-04
!CC*************************************************************************************************************************
!!****************************************************************************************************************
!!   This  code is used to  compute the cavity  distribution function  of a  NaCl aqueous solution                      ***
!  wherein water is treated as a dielectric background with  relative  dielectric constant  78                              ***
! ref: Fernando Bresme and Jose Alejandre,Cavities in ionic liquids , J Chem. Phys. 118,4134,2003;                   ***
!!* The "HISTORY" file from MD simulations of pure solvent (using DL_POLY4 software) is the input ;              ***
!!  the cavity radius distribution function  is the output and is saved in the file 'CR.dat'                                        ***
!!****************************************************************************************************************
    program cavity
    implicit double precision (a-h,o-z)  
	  parameter(npmx1=100000)  !! number of CONFIGs
	  parameter(n_mol=512)     !! number of ions In the cell
	  parameter(n_smp=100)    !! number of insertion for each CONFIG
    Parameter(Maxx = 50)    !! bin number used for output   
    Parameter(Delta=0.20d0)  !! bin size
    dimension RR0(3), RR( n_mol,3) , r12(n_mol),NGg(0:100) 
    character*5  char0(12),char1(12) 
    
     dbox= 41.380d0   !! length of the cubic MD simulation box 
     open(10, file='HISTORY',status='old')
     
     read(10,*)  
     read(10,*)         
     
        rmax=dbox/2.0d0  ! the maximum radius is half of the box length
        ikey=0
        NGg(:)=0
    do i=1, npmx1   !!! loop over CONFIGS
      read(10,*)     
      read(10,*)                
      read(10,*)                    
      read(10,*)                
       
	    do j=1, n_mol   !!! loop over ions

	      read(10,*) 
	      read(10,*)   RR(j,1), RR(j,2), RR(j,3)   !! position in x, y, z derections 
           enddo  !! end for j, loop over ions
        
           
         do k= 1, n_smp  
           ikey = ikey+1   
             call RANDOM_SEED( )
             call RANDOM_NUMBER(tx)
             call RANDOM_NUMBER(ty)
             call RANDOM_NUMBER(tz)
             
          RR0(1)= dbox*(tx-0.5d0)
          RR0(2)= dbox*(ty-0.5d0)
          RR0(3)= dbox*(tz-0.5d0)
          
          do jj =1, n_mol
        dRRx = RR0(1) -RR(jj,1) 
        dRRy = RR0(2) -RR(jj,2) 
        dRRz = RR0(3) -RR(jj,3) 

        dRRx = dRRx   - ANINT ( dRRx /dbox)*dbox
        dRRy = dRRy   - ANINT ( dRRy /dbox)*dbox
        dRRz = dRRz   - ANINT ( dRRz /dbox)*dbox
         
         r12(jj)= sqrt(dRRx**2 + dRRy**2 + dRRz**2)
          enddo   !! end for jj
         
          rcav = minval(r12)          
          
          if(rcav.le.rmax) then
             ig=int(rcav/Delta)
             NGg(ig)=NGg(ig)+1 
         endif
                      
         enddo  !! end for k
            
     enddo  !! end for i, loop over CONFIGS
  
!! output the cavity  distribution function C(R) 
!! C(R)  represents the probability of ﬁnding a cavity of radius  R

     open(30, file='CR.dat', status='unknown')
     Ntot = npmx1*n_smp
     do i= 0, Maxx
     Rx= (i*1.0+0.5d0)*Delta
     CR=  NGg(i)*1.0/Ntot/Delta
     write(30,*) Rx , CR
     enddo
     
    close(30) 
    close(10)

  stop

end program cavity        
        
        
