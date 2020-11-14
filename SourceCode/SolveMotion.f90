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
      SUBROUTINE SolveMotion(WK,W1,TP,WL,AMP,AMAS,BDMP,VDMP,EXFC,DSPL)
	  USE HAMS_mod
	  USE Const_mod
	  USE Body_mod
      IMPLICIT   NONE

      REAL*8,INTENT(IN)::  WK,W1,TP,WL,AMP
      REAL*8,INTENT(IN):: AMAS(6,6),BDMP(6,6)
      COMPLEX*16,INTENT(IN):: EXFC(6)
      COMPLEX*16,INTENT(OUT):: DSPL(6)
      
	  INTEGER I,J,K,L,NSTP,INFO,IPV(6)

      REAL*8 RERR,PHAS(6),VDMP(6,6)
      COMPLEX*16 EXFC2(6),LEFT(6,6),RIGHT(6),DSPL1(6),DX(6)
!
! ========================================================

      DO I=1,6
      EXFC2(I)=CMPLX(-IMAG(EXFC(I)),-REAL(EXFC(I)))
      ENDDO

      RERR=100.D0
      LEFT=-W1**2*(MATX+AMAS)+CRS+KSTF-CI*W1*BDMP
      RIGHT=EXFC

      CALL ZGESV( 6, 1, LEFT, 6, IPV, RIGHT, 6, INFO )

      DSPL=RIGHT

      DO 100 WHILE (RERR.GT.0.0001D0)

        DO I=1,6
        DO J=1,6
         VDMP(I,J)=0.D0 !BDMP(I,J)+VDMP(I,J)*8.D0/3.D0/PI*W1*ABS(DSPL(J))
        ENDDO
        ENDDO
          
        LEFT=-W1**2*(MATX+AMAS)-CI*W1*(BDMP+VDMP)+CRS+KSTF
        RIGHT=EXFC
        
        CALL ZGESV( 6, 1, LEFT, 6, IPV, RIGHT, 6, INFO )
        DSPL1=RIGHT
        
        DX=DSPL1-DSPL
        RERR=0.D0
        DO K=1,6
        RERR=RERR+ABS(DX(K))/ABS(DSPL1(K))
        ENDDO
        
        DSPL=DSPL1  !+0.75D0*DX

100   CONTINUE

      RETURN        
      END SUBROUTINE SolveMotion