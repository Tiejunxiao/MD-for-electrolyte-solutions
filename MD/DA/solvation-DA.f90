CC*************************************************************************************************************************
CC*   Copyright (c) 2026 [Xueyu Song and Tiejun Xiao]. 
CC*  This code is distributed under the terms of the Creative Commons Attribution 4.0 International  License (CC BY 4.0).
CC*  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
CC*
CC*  Author: Xueyu Song
CC*  Email:   [xsong@iastate.edu]
CC*  Author: Tiejun Xiao
CC*  Email:   [tjxiao@gznc.edu.cn]
CC*  Date:    2026-10-04
CC*************************************************************************************************************************

Program Ewaldsum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    This code is used to calculate the mean electric potential of a site of a diatomic solute  solvated in  a NaCl aqueous solution          !!!!
!!     using the Ewald summation method given  the configurations  from MD simulation                                                                                  !!!!
!!     The input file is a "HISTORY" file with 10^5 configurations                                                                                                                               !!!!
!!     The "HISTORY" file is  generated with the MD simulation software DL_POLY4                                                                                             !!!!
!!******************************************************************************************************************************************
!!     There are six kinds of species in the  electrolyte solution, Cl-, Na+, SA1, SA2, SA3 and SA4                                                                     !!!!
!!      SA1 and SA2 are two sites of a diatomic solute 1,  SA3 and SA3 are two sites of a diatomic solute 2                                                    !!!!
!!     Denote q1 and q2 as the two site charge numbers of solute 1;    q3 and q4 as the two site charge numbers of solute 2                   !!!!
!!    we take q3=-q1, q4=-q2 to keep the MD simulation box neutral.                                                                                                                    !!!!
!*******************************************************************************************************************************************
!!   IMPORTANT: without modification, this code  applies to the case of q1=-1, q2=0                                                                                       !!!!
!!  specifically, the HISTORY file generated for MD with key (DAj, j=1,2,3,4,5) can be used as input of this code.                                        !!!!
!!    When q1,2 take other values, one needs to modify the related site charge numbers  zi(3:6)  in line 59                                                   !!!!
!*******************************************************************************************************************************************
!!   IMPORTANT: This code outputs the mean electric potential of each solute site rather than the electrostatic energy of the solute.  !!!!
!!     Xavej (j=1,2) is the mean electric potential felt by the SAj site                                                                                                                         !!!!
!!     Xvarj (j=1,2) is the standard deviation of the electric potential felt by SAj site                                                                                             !!!!
!!     The mean electric potential and its standard deviation are saved in the output file 'potential.dat'.                                                            !!!!
!!******************************************************************************************************************************************
!!  The mean electric potentials can further be used to evaluate the electrostatic energy of the solute.                                                       !!!! 
!!  For the diatomic solute 1 with site charges  Qk=qk*e_0(k=1,2)  and site separation distance R12,                                                            !!!!    
!!   the electrostatic energy of the solute can  be evaluated   as  u^e= \sum_{k=1,2} Qk*phi_k^ind/2 ,                                                           !!!!
!   with phi_k^ind the induced potential for site k of the solute, and                                                                                                                      !!!!
!!       phi_1^ind= Xave1 -  Q2/(4*PI* Epsilon_0*epsilon_r* R12),                                                                                                                               !!!!
!!       phi_2^ind= Xave2 -  Q1/(4*PI* Epsilon_0*epsilon_r* R12),                                                                                                                              !!!!
!!   with  Epsilon_0  the permittivity of the vacuum,  e_0  the elementary charge .                                                                                               !!!!
!!   and epsilon_r= 78.0  the relative dielectric constant of the aqueous electrolyte solution.                                                                          !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

 Implicit double precision (a-h,o-z)
 parameter(npmx1=100000 )     !!  number of sampling 
 parameter(N=512)        !! number of atoms in MD
 parameter(maxk=1000 )
 dimension RR( N,3), dRR(N,3)    
 dimension zi(6), fkm(maxk) , Zq(N) , Xtot1(npmx1),Xtot2(npmx1) , Xtot3(npmx1),Xtot4(npmx1)
 character*5   char1 
 real*8   Phyr1, Phyk1, Phys1, Phytot1, Xave1, xvar1
 real*8   Phyr2, Phyk2, Phys2, Phytot2, Xave2, xvar2
 real*8   Phyr3, Phyk3, Phys3, Phytot3, Xave3, xvar3
 real*8   Phyr4, Phyk4, Phys4, Phytot4, Xave4, xvar4

 common/fxmat/dbox, xkappa, Np
 
 zi=(/-1d0, 1.d0, -1.d0, 0.d0, 1.d0, 0.d0/)
!!  zi(1:6) are the charge numbers of the Cl-, Na+,SA1, SA2, SA3 and SA4

!! En_trans is an energy  coefficient for an electrostatic system.
!! e_0 is the elementary charge
!! Epsilon_0 is the permittivity of the vacuum.
!! epsilon_r= 78.0 is the relative dielectric constant of the aqueous electrolyte solution. 
!! the electrostatic energy of two elementary charges separated by 1 Angstrom is about 1389.35  kJ/mol.
!!! i.e.,  1 e_0^2/(4 *PI * Epsilon_0* 1 A) =1389.35 kJ/mol
 epsilon_r=78.0
 En_tran= 1389.35d0/epsilon_r  
  Np=N
         
          j1= 509    !! tag for SA1 site of the diatomic molecule 1
          j2= 510    !! tag for SA2  site of  the diatomic molecule 1
          j3= 511    !! tag for SA3 site of the diatomic molecule 2
    	  j4= 512    !! tag for SA4 site of the diatomic molecule 2

       open(20, file='potential.dat',status='unknown')   
 
!!!****************************************************************************************
 open(10, file='HISTORY',status='old')
	 read(10,*)   
         read(10,*)   

!!*****************************************************************************************        
        npmx3= npmx1
	 do i=1, npmx3

              read(10,*)                 !! information of simulation (timestep, number of molecules, etc.)  
	      read(10,*)  dbox1   !! length of  the cubic box   
	      read(10,*)                
	      read(10,*)  
 
              dbox =dbox1
              xkappa=  6.5d0/dbox
            call stk( fkm)
 
	    do j=1,  N
 
 
         read(10,*)   char1 , no2 , xnn 
	 read(10,*)   RR(j,1), RR(j,2), RR(j,3)               !! position in x, y, z derections 
  
 
          if(char1=="Cl-")then 
         ns1=1

      elseif(char1=="Na+")then
         ns1=2

      elseif(char1=="SA1")then
         ns1=3

      elseif(char1=="SA2")then
         ns1=4

      elseif(char1=="SA3")then
         ns1=5

      elseif(char1=="SA4")then
         ns1=6

      endif

        Zq(j)=zi(ns1)
 
	    enddo  !! end for j
          
!!*********************************************************************************
!! use the Ewald summation method to compute the  electric potential .
!! Phytotj (j=1,2,3,4) is the total electric potential felt by the SAj site  

        j0=  j1
        z0=zq(j0)
        do k =1, N
        dRR(k,1) = RR(k,1) -RR(j0,1) 
        dRR(k,2) = RR(k,2) -RR(j0,2) 
        dRR(k,3) = RR(k,3) -RR(j0,3) 

        dRR(k,1) = dRR(k,1)   - ANINT ( dRR(k,1) /dbox)*dbox
        dRR(k,2) = dRR(k,2)   - ANINT ( dRR(k,2) /dbox)*dbox
        dRR(k,3) = dRR(k,3)   - ANINT ( dRR(k,3) /dbox)*dbox

        enddo
 
      call EWSR(Phyr1, dRR, Zq     )        !! r-space contribution
      call EWSK(Phyk1, dRR, Zq, fkm)       !! k-space contribution   
      call EWSS(Phys1, z0)             !! solute image contribution
      
      Phytot1= Phyr1 + Phyk1 - Phys1
!**********************************************************************************
         j0=  j2
        z0=zq(j0)
        do k =1, N
        dRR(k,1) = RR(k,1) -RR(j0,1) 
        dRR(k,2) = RR(k,2) -RR(j0,2) 
        dRR(k,3) = RR(k,3) -RR(j0,3) 

        dRR(k,1) = dRR(k,1)   - ANINT ( dRR(k,1) /dbox)*dbox
        dRR(k,2) = dRR(k,2)   - ANINT ( dRR(k,2) /dbox)*dbox
        dRR(k,3) = dRR(k,3)   - ANINT ( dRR(k,3) /dbox)*dbox

        enddo

 
      call EWSR(Phyr2, dRR, Zq     )        !! r-space contribution
      call EWSK(Phyk2, dRR, Zq, fkm)       !! k-space contribution   
      call EWSS(Phys2, z0)             !! solute image contribution
      
      Phytot2= Phyr2 + Phyk2 - Phys2
!**********************************************************************************
        j0=  j3
        z0=zq(j0)
        do k =1, N
        dRR(k,1) = RR(k,1) -RR(j0,1) 
        dRR(k,2) = RR(k,2) -RR(j0,2) 
        dRR(k,3) = RR(k,3) -RR(j0,3) 

        dRR(k,1) = dRR(k,1)   - ANINT ( dRR(k,1) /dbox)*dbox
        dRR(k,2) = dRR(k,2)   - ANINT ( dRR(k,2) /dbox)*dbox
        dRR(k,3) = dRR(k,3)   - ANINT ( dRR(k,3) /dbox)*dbox

        enddo
 
      call EWSR(Phyr3, dRR, Zq     )        !! r-space contribution
      call EWSK(Phyk3, dRR, Zq, fkm)       !! k-space contribution   
      call EWSS(Phys3, z0)             !! solute image contribution
      
      Phytot3= Phyr3 + Phyk3 - Phys3
!**********************************************************************************
        j0=  j4
        z0=zq(j0)
        do k =1, N
        dRR(k,1) = RR(k,1) -RR(j0,1) 
        dRR(k,2) = RR(k,2) -RR(j0,2) 
        dRR(k,3) = RR(k,3) -RR(j0,3) 

        dRR(k,1) = dRR(k,1)   - ANINT ( dRR(k,1) /dbox)*dbox
        dRR(k,2) = dRR(k,2)   - ANINT ( dRR(k,2) /dbox)*dbox
        dRR(k,3) = dRR(k,3)   - ANINT ( dRR(k,3) /dbox)*dbox

        enddo
 
      call EWSR(Phyr4, dRR, Zq     )        !! r-space contribution
      call EWSK(Phyk4, dRR, Zq, fkm)       !! k-space contribution   
      call EWSS(Phys4, z0)             !! solute image contribution
      
      Phytot4= Phyr4 + Phyk4 - Phys4

!**********************************************************************************

       Xtot1(i)=  Phytot1 * En_tran       
       Xtot2(i)=  Phytot2 * En_tran  
       Xtot3(i)=  Phytot3 * En_tran       
       Xtot4(i)=  Phytot4 * En_tran 

!************************************************************************************
   enddo   !! end for i
        
!************************************************************************************
!! evaluate the mean value (Xave) and standard deviation (Xvar) of the electric potential
!! Xavej (j=1,2,3,4) is the mean electric potential felt by the SAj site
!! Xvarj  (j=1,2,3,4) is the standard deviation of the electric potential felt by the SAj site
!! Xavej and Xvarj (j=1,2,3,4) are saved in the output file 'Energy.dat'

        Xave1=SUM(   Xtot1(1:npmx3)           )/Real(npmx3)
        Var1= SUM( ( Xtot1(1:npmx3) -Xave1)**2 )/Real(npmx3)
        xvar1= dsqrt(Var1)

        Xave2=SUM(   Xtot2(1:npmx3)           )/Real(npmx3)
        Var2= SUM( ( Xtot2(1:npmx3) -Xave2)**2 )/Real(npmx3)
        xvar2= dsqrt(Var2)

        Xave3=SUM(   Xtot3(1:npmx3)           )/Real(npmx3)
        Var3= SUM( ( Xtot3(1:npmx3) -Xave3)**2 )/Real(npmx3)
        xvar3= dsqrt(Var3)

        Xave4=SUM(   Xtot4(1:npmx3)           )/Real(npmx3)
        Var4= SUM( ( Xtot4(1:npmx3) -Xave4)**2 )/Real(npmx3)
        xvar4= dsqrt(Var4)

       write(*,*) Xave1, xvar1
       write(*,*) Xave2, xvar2   
       write(*,*) Xave3, xvar3
       write(*,*) Xave4, xvar4

       write(20,*)xave1, xvar1
       write(20,*)xave2, xvar2 
       write(20,*)xave3, xvar3
       write(20,*)xave4, xvar4 

close(10)
 
stop
end program Ewaldsum

 

!!*****************************************************************************************************
!! r-space contribution for Ewald sum
!! PhyR is the r-space contribution to potential in the center of solute 
subroutine EWSR(Phyr, dRR, Z)
Implicit double precision (a-h,o-z)
integer  i,j, Np
real*8 PI,dbox, xkappa, Phyr,  RIJ, uij
parameter(PI=3.141592654)
real*8  dRR(Np,3) ,  Z(Np)
common/fxmat/dbox, xkappa, Np

           PhyR = 0.0
          
        DO   i = 1, Np   
         
              RIJ   = dSQRT ( dRR(i,1)**2+ dRR(i,2)**2 + dRR(i,3)**2)
              
			  if(RIJ.NE.0.)then			  
			  uij= z(i) / RIJ  * ERFC (xkappa * RIJ )
			  endif

			  if(RIJ==0.)then               !!  RIJ=0 is related to i=j0 term
			  uij =0.d0  
              endif
              
			 Phyr  = Phyr +   uij                         !! r-space contribution
!            if(i<10)then
!            write(*,*) RIJ
!            endif  
            enddo  !!  ii

         
   return
end subroutine EWSR



!!*******************************************************************************************************
!! k-space contribution to Ewald sum
!! Phyk is the k-space contribution to potential in the center of solute
subroutine EWSK(Phyk, dRR, Z, fkm)
Implicit double precision (a-h,o-z)
integer  i,j, km, ktot, Np, kx,ky,kz 
parameter(km=5, ksqmax=27,maxk=1000)
real*8  dRR(Np,3) , fkm(maxk), Z(Np)
real*8  PI, dbox, xkappa, Phyk
COMPLEX*8   EIKX(1:Np,   0:km)
COMPLEX*8   EIKY(1:Np, -km:km)
COMPLEX*8   EIKZ(1:Np, -km:km)
complex*8   EIKR(Np),sum1 
common/fxmat/dbox, xkappa, Np  
parameter(PI=3.141592654)


         DO  i = 1, Np

           EIKX(i, 0) = (1.0, 0.0)
           EIKY(i, 0) = (1.0, 0.0)
           EIKZ(i, 0) = (1.0, 0.0)

           do j=1, km
 
            EIKX(i, j) = CMPLX ( dCOS ( 2*PI*Real(j)/dbox * dRR(i,1) ) , dSIN ( 2*PI*Real(j)/dbox * dRR(i,1) ) )
            EIKY(i, j) = CMPLX ( dCOS ( 2*PI*Real(j)/dbox * dRR(i,2) ) , dSIN ( 2*PI*Real(j)/dbox * dRR(i,2) ) )
            EIKZ(i, j) = CMPLX ( dCOS ( 2*PI*Real(j)/dbox * dRR(i,3) ) , dSIN ( 2*PI*Real(j)/dbox * dRR(i,3) ) )

            EIKY(i, -j) = CONJG ( EIKY(i, j) )
            EIKZ(i, -j) = CONJG ( EIKZ(i, j) )
          enddo  !! j
 		  enddo !! i

      Phyk=0.d0
	  ktot = 0
   do  kx= 0, km
     
          IF ( KX .EQ. 0 ) THEN
              FACTOR = 1.0
           ELSE
              FACTOR = 2.0
           ENDIF

      do ky= -km, km      
	    do kz= -km ,km
            KSQ = KX * KX + KY * KY + KZ * KZ
            IF ( ( KSQ .LT. KSQMAX ) .AND. ( KSQ .NE. 0 ) ) THEN

             ktot = Ktot + 1
             SUM1  = (0.0, 0.0)                  
			  DO  i = 1, Np
               EIKR(i) = EIKX(i, KX) * EIKY(i, KY) * EIKZ(i, KZ)
               SUM1     = SUM1 + Z(i) * EIKR(i)
              enddo ! i
             Phyk =   Phyk + 2.* FACTOR * fkm(ktot) * (CONJG ( SUM1 ) + SUM1 )/2.

            ENDIF
           enddo  ! kz
          enddo  ! ky
         enddo  ! kx
       

  return
end subroutine EWSK

!!***************************************************************************************************
!! image contribution of the solute in Ewald sum
!! Phys is the contribution to potential in the center of solute due to image of the original solute
subroutine EWSS(Phys, Z0)
Implicit double precision (a-h,o-z)
real*8  ESQPI, Phys, Z0
parameter(RSQPI = 0.5641896)
common/fxmat/dbox, xkappa, Np  


Phys= 2.*RSQPI * xkappa *Z0

end subroutine EWSS
 

!!***************************************************************************************************
subroutine  stk( fkm)   !! initiate the k-function before Ewald summation
Implicit double precision (a-h,o-z)
real*8 rvk, dbox, xkappa, PI
integer kx, ky, kz, km, ktot
parameter(PI=3.141592654,ksqmax=27, km=5, maxk=1000)
real*8 fkm(ksqmax), vk(3) 
common/fxmat/dbox, xkappa, Np
 
!    ** LOOP OVER K-VECTORS. NOTE KX IS NON-NEGATIVE ** 
             ktot=0
        DO   KX = 0, km

           vk(1) = 2.0* PI* REAL ( KX )/dbox

           DO   KY = -km, km

              vk(2) = 2.0* PI * REAL ( KY )/dbox

              DO KZ = -km, km

                 vk(3) = 2.0*PI * REAL ( KZ )/dbox

                 rvk = vk(1)**2 + vk(2)**2 +vk(3)**2
                 KSQ = KX * KX + KY * KY + KZ * KZ

                 IF ( ( KSQ .LT. KSQMAX ) .AND. ( rvk .NE. 0 ) ) THEN
                     Ktot= ktot +1
                   fkm(ktot) = 2.0* PI/dbox**3   * dEXP ( -rvk/4./xkappa**2 ) / rvk

                 IF ( KTOT .GT. MAXK ) STOP 'fkm is too small'

                 ENDIF

                enddo  ! kz
			 enddo    ! ky
		 enddo      ! kx

  
 !       WRITE( *, ' ( '' EWALD SUM SETUP COMPLETE ''     ) ' )
 !       WRITE( *, ' ( '' NUMBER OF WAVEVECTORS IS '', I5 ) ' ) ktot
    
    return
end subroutine stk
!****************************************************************************************

!!***************************************************************************************
        REAL*8 FUNCTION ERFC1 ( X )

!    *******************************************************************
!    ** APPROXIMATION TO THE COMPLEMENTARY ERROR FUNCTION             **
!    **                                                               **
!    ** REFERENCE:                                                    **
!    **                                                               **
!    ** ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,    **
!    **    NATIONAL BUREAU OF STANDARDS, FORMULA 7.1.26               **
!    *******************************************************************

        REAL*8        A1, A2, A3, A4, A5, P

        PARAMETER ( A1 = 0.254829592, A2 = -0.284496736 )
        PARAMETER ( A3 = 1.421413741, A4 = -1.453152027 )
        PARAMETER ( A5 = 1.061405429, P  =  0.3275911   )

        REAL*8        T, X, XSQ, TP

!    *******************************************************************

        T  = 1.0 / ( 1.0 + P * X )
        XSQ = X * X

        TP = T * ( A1 + T * ( A2 + T * ( A3 + T * ( A4 + T * A5 ) ) ) )

        ERFC1 = TP * dEXP ( -XSQ )

        RETURN
        END FUNCTION ERFC1 

!******************************************************************************************
