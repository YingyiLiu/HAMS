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
MODULE BodyIntgr

   USE HAMS_mod
   USE Body_mod
   USE Const_mod
   USE WaveDyn_mod
   USE PanelMesh_mod
   USE Inerfs_mod
   USE PatcVelct
   USE Potentials_mod
   USE omp_lib
   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................
   
   PUBLIC :: BODINT_LEFT
   PUBLIC :: RBC_RIGHT
   PUBLIC :: DBC_RIGHT
   
CONTAINS

!   --------------------------------------------------------
!     Integration inside an element, which is not singular 
!   --------------------------------------------------------

      SUBROUTINE BODINT_LEFT(IS,IEL,JEL,TINDP,FLAG)
      IMPLICIT NONE 
      
      INTEGER,INTENT(IN):: IS,IEL,JEL,FLAG
      COMPLEX*16,INTENT(OUT):: TINDP(4)
      
      REAL*8::   RKN(4),ENV(6),EAR
      COMPLEX*16  GRN(4)

      EAR=DS(JEL)
      ENV(:)=DXYZ_P(JEL,:)

      RKN(:)=RKBN(IEL,JEL,IS,:)
      GRN(:)=CGRN(IEL,JEL,IS,:)


      IF (FLAG.EQ.1) THEN

        TINDP(IS)=(RKN(2)+GRN(2)*EAR)*ENV(1)   &
                 +(RKN(3)+GRN(3)*EAR)*ENV(2)   &
                 +(RKN(4)+GRN(4)*EAR)*ENV(3)

      ELSE

        TINDP(IS)=(GRN(2)*ENV(1)+GRN(3)*ENV(2)+GRN(4)*ENV(3))*EAR
        
      ENDIF 
       
      RETURN
      END SUBROUTINE BODINT_LEFT
 
!   --------------------------------------------------------
!     Integration inside an element, which is not singular 
!   --------------------------------------------------------

      SUBROUTINE RBC_RIGHT(IS,IEL,JEL,TINRD,FLAG)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN):: IS,IEL,JEL,FLAG
      COMPLEX*16,INTENT(OUT):: TINRD(4,6,4)
      
      INTEGER  IP
      REAL*8  XQ,YQ,ZQ
      COMPLEX*16  DUM,DPOX,DPOY,DPOZ

!  Above----
!  Field point changes, source element keeps the same as in the first quadrant
       
      IF (FLAG.EQ.1) THEN
        DUM=RKBN(IEL,JEL,IS,1)+CGRN(IEL,JEL,IS,1)*DS(JEL)
      ELSE
        DUM=CGRN(IEL,JEL,IS,1)*DS(JEL)
      ENDIF

     DO 200 IP=1, NSYS
       IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
        XQ=SY(IS,IP)*XYZ_P(JEL,1)
        YQ=SX(IS,IP)*XYZ_P(JEL,2)
        ZQ=           XYZ_P(JEL,3)
        TINRD(IS,1,IP)=TINRD(IS,1,IP)+DUM*SY(IS,IP)*DXYZ_P(JEL,1)
        TINRD(IS,2,IP)=TINRD(IS,2,IP)+DUM*SX(IS,IP)*DXYZ_P(JEL,2)
        TINRD(IS,3,IP)=TINRD(IS,3,IP)+DUM*DXYZ_P(JEL,3)          
        TINRD(IS,4,IP)=TINRD(IS,4,IP)+DUM*SX(IS,IP)*DXYZ_P(JEL,4)
        TINRD(IS,5,IP)=TINRD(IS,5,IP)+DUM*SY(IS,IP)*DXYZ_P(JEL,5)
        TINRD(IS,6,IP)=TINRD(IS,6,IP)+DUM*SY(IS,IP)*SX(IS,IP)*DXYZ_P(JEL,6)
       ELSE
        XQ=SX(IS,IP)*XYZ_P(JEL,1)
        YQ=SY(IS,IP)*XYZ_P(JEL,2)
        ZQ=          XYZ_P(JEL,3)
        TINRD(IS,1,IP)=TINRD(IS,1,IP)+DUM*SX(IS,IP)*DXYZ_P(JEL,1)
        TINRD(IS,2,IP)=TINRD(IS,2,IP)+DUM*SY(IS,IP)*DXYZ_P(JEL,2)
        TINRD(IS,3,IP)=TINRD(IS,3,IP)+DUM*DXYZ_P(JEL,3)          
        TINRD(IS,4,IP)=TINRD(IS,4,IP)+DUM*SY(IS,IP)*DXYZ_P(JEL,4)
        TINRD(IS,5,IP)=TINRD(IS,5,IP)+DUM*SX(IS,IP)*DXYZ_P(JEL,5)
        TINRD(IS,6,IP)=TINRD(IS,6,IP)+DUM*SX(IS,IP)*SY(IS,IP)*DXYZ_P(JEL,6)
       ENDIF
200  CONTINUE
       
      RETURN
      END SUBROUTINE RBC_RIGHT
      
      
!   --------------------------------------------------------
!     Integration inside an element, which is not singular 
!   --------------------------------------------------------

      SUBROUTINE DBC_RIGHT(IS,IEL,JEL,TINRD,FLAG)
      IMPLICIT NONE 
      
      INTEGER,INTENT(IN):: IS,IEL,JEL,FLAG
      COMPLEX*16,INTENT(OUT):: TINRD(4,4)
      
      INTEGER  IP
      REAL*8 XQ,YQ,ZQ
      COMPLEX*16  DUM(2),DPOX,DPOY,DPOZ

!  Above----
!  Field point changes, source element keeps the same as in the first quadrant

      IF (FLAG.EQ.1) THEN
        DUM(1)=RKBN(IEL,JEL,IS,1)+CGRN(IEL,JEL,IS,1)*DS(JEL)
      ELSE
        DUM(1)=CGRN(IEL,JEL,IS,1)*DS(JEL)
      ENDIF

     DO 200 IP=1, NSYS
       IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
        XQ=SY(IS,IP)*XYZ_P(JEL,1)
        YQ=SX(IS,IP)*XYZ_P(JEL,2)
        ZQ=           XYZ_P(JEL,3)
        CALL DINP(XQ,YQ,ZQ,XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)
        DUM(2)=SY(IS,IP)*DPOX*DXYZ_P(JEL,1)+SX(IS,IP)*DPOY*DXYZ_P(JEL,2)+DPOZ*DXYZ_P(JEL,3)        
        TINRD(IS,IP)=TINRD(IS,IP)-DUM(1)*DUM(2)
       ELSE
        XQ=SX(IS,IP)*XYZ_P(JEL,1)
        YQ=SY(IS,IP)*XYZ_P(JEL,2)
        ZQ=          XYZ_P(JEL,3)
        CALL DINP(XQ,YQ,ZQ,XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)
        DUM(2)=SX(IS,IP)*DPOX*DXYZ_P(JEL,1)+SY(IS,IP)*DPOY*DXYZ_P(JEL,2)+DPOZ*DXYZ_P(JEL,3)        
        TINRD(IS,IP)=TINRD(IS,IP)-DUM(1)*DUM(2)
       ENDIF

200  CONTINUE
       
       RETURN
       END SUBROUTINE DBC_RIGHT
!-------------------------------------------------------------------------------
END MODULE BodyIntgr
!*******************************************************************************