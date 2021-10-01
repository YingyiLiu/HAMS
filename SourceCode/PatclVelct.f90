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
MODULE PatcVelct

   USE HAMS_mod
   USE Const_mod 
   USE WaveDyn_mod
   USE PanelMesh_mod
   
   IMPLICIT NONE

  ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: VINP
   PUBLIC :: DINP

CONTAINS
!------------------------------------------------- 
!           For Incident potential 
!-------------------------------------------------
!
        COMPLEX*16 FUNCTION VINP(X,Y,Z,XW,YW,SITA)
        IMPLICIT  NONE

        REAL*8,INTENT(IN):: X,Y,Z,XW,YW,SITA
        REAL*8:: WKX,DUM,W0,WK0
        
        IF (ABS(V).LT.1.E-8) THEN
         W0=1.E-20
         WK0=1.E-20
        ELSEIF (ABS(V+1.D0).LT.1.E-8) THEN
         W0=1.E+20
         WK0=1.E+20
        ELSE
         W0=W1
         WK0=WK
        ENDIF
 
        IF (Z.GT.0.D0) WRITE(6,*) 'Z>0, IN VINP', 'Z=',Z 
        CI=CMPLX(0.0D0,1.0D0)
 
        WKX=WK0*((X-XW)*COS(SITA)+(Y-YW)*SIN(SITA))

	    IF (H .LT. 0.0D0 .or. H .GT. 500.0D0) THEN
          VINP=-CI*AMP*G/W0*EXP(WK0*Z)*EXP(CI*WKX)
	    ELSE
          DUM=-AMP*W0/(WK0*SINH(WK0*H))
          VINP=CI*DUM*COSH( WK0*(Z+H) )*EXP(CI*WKX)
        ENDIF

        VINP=VINP/CI
        RETURN
        END FUNCTION VINP

 
!    ----------------------------------------------------------------------------------------
!       For derivatives of incident wave potential, i.e. incident wave orbital acceleration
!    ----------------------------------------------------------------------------------------
!
        SUBROUTINE  DINP(X,Y,Z,XW,YW,SITA,DPOX,DPOY,DPOZ)
	    IMPLICIT    NONE
	  
        REAL*8,INTENT(IN):: X,Y,Z,XW,YW,SITA
        COMPLEX*16,INTENT(OUT)::  DPOX,DPOY,DPOZ

        REAL*8:: WKX,THKDH,THKDV,DUM,W0,WK0
        
        IF (ABS(V).LT.1.E-8) THEN
         W0=1.E-20
         WK0=1.E-20
        ELSEIF (ABS(V+1.D0).LT.1.E-8) THEN
         W0=1.E+20
         WK0=1.E+20         
        ELSE
         W0=W1
         WK0=WK     
        ENDIF 
        
        IF (Z.GT.0.D0) WRITE(6,*) 'Z>0, IN DINP', 'Z=',Z

        WKX=WK0*((X-XW)*COS(SITA)+(Y-YW)*SIN(SITA))

	    IF (H .LT. 0.0D0) THEN
          DPOX=    W0*AMP*EXP(WK0*Z)*EXP(CI*WKX)*COS(SITA)
          DPOY=    W0*AMP*EXP(WK0*Z)*EXP(CI*WKX)*SIN(SITA)
          DPOZ=-CI*W0*AMP*EXP(WK0*Z)*EXP(CI*WKX) 
        ELSE
!         write(*,*)  z+H,COSH(WK*(z+H)),SINH(WK*H)           
          THKDH=   COSH(WK0*(Z+H))/SINH(WK0*H)
          THKDV=   SINH(WK0*(Z+H))/SINH(WK0*H)  
          DPOX=    W0*AMP*THKDH*EXP(CI*WKX)*COS(SITA)
          DPOY=    W0*AMP*THKDH*EXP(CI*WKX)*SIN(SITA)
          DPOZ=-CI*W0*AMP*THKDV*EXP(CI*WKX) 
	    ENDIF

          DPOX= DPOX/CI
          DPOY= DPOY/CI
          DPOZ= DPOZ/CI
        RETURN
        END SUBROUTINE  DINP
!-------------------------------------------------------------------------------
END MODULE PatcVelct
!*******************************************************************************





