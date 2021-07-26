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
MODULE HydroStatic

   USE HAMS_mod
   USE Body_mod
   USE Const_mod 
   USE WaveDyn_mod

   IMPLICIT NONE

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Initialisation
   PUBLIC :: ReadHydroStatic

CONTAINS
! ---------------------------------------------------------------------------------------------
!      Initialisation of the matrices, variables, and parameters
! ---------------------------------------------------------------------------------------------
! 
      SUBROUTINE Initialisation
      IMPLICIT   NONE
      
      XW=0.D0
      AMP=1.D0
      SYBO=0
        
      XB=0.D0
      MATX=0.D0
      CRS=0.D0
      VDMP=0.D0
      KSTF=0.D0
      
      EXFC=CMPLX(0.D0,0.D0)
      AMAS=0.D0
      BDMP=0.D0
      DSPL=0.D0
      
      RETURN
      END SUBROUTINE Initialisation
      
      
! ---------------------------------------------------------------------------------------------
!      Read the data of gravity center, mass matrix and restoring matrix.
! ---------------------------------------------------------------------------------------------
! 
      SUBROUTINE ReadHydroStatic
      IMPLICIT   NONE  

      INTEGER I,J

      READ(4,*)
      READ(4,*) XG(1),XG(2),XG(3)
      READ(4,*)
      DO I=1,6
       READ(4,120) (MATX(I,J), J=1, 6)
      ENDDO
      READ(4,*)
      DO I=1,6
       READ(4,120) (VDMP(I,J), J=1, 6)
      ENDDO
      READ(4,*)
      DO I=1,6
       READ(4,120) (CRS(I,J), J=1, 6)
      ENDDO
      READ(4,*)
      DO I=1,6
       READ(4,120) (KSTF(I,J), J=1, 6)
      ENDDO
      
      MATX(1,5)=MATX(1,1)*XG(3)
      MATX(1,6)=-MATX(1,1)*XG(2)
      MATX(2,4)=-MATX(1,1)*XG(3)
      MATX(2,6)=MATX(1,1)*XG(1)
      MATX(3,4)=MATX(1,1)*XG(2)
      MATX(3,5)=-MATX(1,1)*XG(1)
      
      MATX(5,1)=MATX(1,5)
      MATX(6,1)=MATX(1,6)
      MATX(4,2)=MATX(2,4)
      MATX(6,2)=MATX(2,6)
      MATX(4,3)=MATX(3,4)
      MATX(5,3)=MATX(3,5)
      
      DO I=1,6
       DO J=1,6
       WRITE(65,130) I,J,CRS(I,J)/(RHO*G)
       ENDDO
      ENDDO
 
120   FORMAT(6(2x,E12.5))
130   FORMAT(2I6,2X,ES14.6)
      
      RETURN        
      END SUBROUTINE ReadHydroStatic
!-------------------------------------------------------------------------------
END MODULE HydroStatic
!*******************************************************************************