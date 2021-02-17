!  ------------------------------------------------------------------------------------------------------
!                                                               
!    Program HAMS for the diffraction and radiation of waves 
!    by 3D structures.
! 
!  License:
! 
!    This routine is part of HAMS.
!
!    HAMS is a free software framework: you can redistribute it and/or modify it 
!    under the terms of the Apache License, Version 2.0 (the "License"); you may 
!    not use these subroutines except in compliance with the License. The software
!    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
!
!    You should have received a copy of the Apache License, Version 2.0, along with 
!    HAMS. If not, see <http://www.apache.org/licenses/LICENSE-2.0>.
!
!  Code Original Author:
!
!    Yingyi Liu
!
!  ------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!       Calculate the motion response in frequency domain by panel model.
!---------------------------------------------------------------------------------------------
!
      SUBROUTINE SolveMotion(W1,TP,OUFR,BETA,AMP,AMAS,BDMP,VDMP,EXFC,DSPL)
      USE HAMS_mod
      USE Const_mod
      USE Body_mod
      IMPLICIT   NONE

      REAL*8,INTENT(IN)::  W1,TP,OUFR,BETA,AMP
      REAL*8,INTENT(IN):: AMAS(6,6),BDMP(6,6),VDMP(6,6)
      COMPLEX*16,INTENT(IN):: EXFC(6)
      COMPLEX*16,INTENT(OUT):: DSPL(6)
      
      INTEGER INFO,IPV(6),MD,MEXP

      REAL*8 MOD,PHS,REL,IMG,NREL,NIMG,NFAC
      COMPLEX*16 LEFT(6,6),RIGHT(6)
!
! ========================================================

      LEFT=-W1**2*(MATX+AMAS)-CI*W1*(BDMP+VDMP)+CRS+KSTF
      RIGHT=EXFC

      CALL ZGESV( 6, 1, LEFT, 6, IPV, RIGHT, 6, INFO )
      DSPL=RIGHT
!
!   =================================================== 
!    Write WAMIT-style output files
!
      DO MD=1,6
          
       IF (MD.LE.3) THEN
        MEXP=2
       ELSEIF (MD.GE.4) THEN
        MEXP=3
       ENDIF

       NFAC=(RHO*G*AMP)*REFL**MEXP
       
       !print*,'NFAC',(RHO*G*AMP),REFL,MEXP
       !pause
       
       REL=REAL(DSPL(MD))/NFAC
       IMG=IMAG(DSPL(MD))/NFAC
       MOD=SQRT(REL**2+IMG**2) !ABS(EXFC(IMD))/NFAC
       NREL= REL
       NIMG=-IMG
       PHS=ATAN2D(NIMG,NREL)
       
       IF (ABS(TP+1.D0).GT.1.E-6.AND.ABS(TP).GT.1.E-6) THEN
        WRITE(63,1030)  OUFR,BETA*180.0D0/PI,MD,MOD,PHS,NREL,NIMG
       ENDIF
           
      ENDDO
!   =================================================== 
1030 FORMAT(2ES14.6,I6,4ES14.6)
     
      RETURN        
      END SUBROUTINE SolveMotion