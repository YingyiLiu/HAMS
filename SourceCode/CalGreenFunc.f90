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
MODULE CalGreenFunc

   USE HAMS_mod
   USE Body_mod
   USE Const_mod
   USE WaveDyn_mod
   USE PanelMesh_mod
   USE Inerfs_mod

   USE SingularIntgr
   USE INFG3D_Open
   USE FinGreen3D_Open
   USE Potentials_mod
   USE omp_lib
   IMPLICIT NONE

  ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: CALGREEN

   
CONTAINS
!   ----------------------------------------------------------------------------
!      Calculate WAVE TERM's value of Green function for arbitray two points
!   ----------------------------------------------------------------------------

      SUBROUTINE CALGREEN
      IMPLICIT NONE 
      
      INTEGER  IEL,JEL,IS,IP,IRR,FLAG
      REAL*8   XQ,YQ,ZQ,XP,YP,ZP,SIJ,DIJ(3),DIST
      COMPLEX*16  GRN(4),DPOX,DPOY,DPOZ

      IRR=1

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,IS,XQ,XP,YQ,YP,ZQ,ZP,SIJ,DIJ,GRN,FLAG)
         
      DO IEL=1, NELEM

        DO JEL=1, NELEM
            
         XQ=XYZ_P(JEL,1)     ! JEL: source point,  IEL: field point
         YQ=XYZ_P(JEL,2)
         ZQ=XYZ_P(JEL,3)

         DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF
         
        DO IS=1, NSYS

         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          XP=SY(IS,1)*XYZ_P(IEL,1)
          YP=SX(IS,1)*XYZ_P(IEL,2)
          ZP=         XYZ_P(IEL,3)
         ELSE
          XP=SX(IS,1)*XYZ_P(IEL,1)
          YP=SY(IS,1)*XYZ_P(IEL,2)
          ZP=         XYZ_P(IEL,3)
         ENDIF

         IF (NCN(JEL).EQ.3) THEN

          CALL SGLINTBD_TRI(IS,IEL,JEL,SIJ,DIJ,1)

         ELSEIF (NCN(JEL).EQ.4) THEN

          CALL SGLINTBD_QUAD(IS,IEL,JEL,SIJ,DIJ,1)

         ENDIF

         IF (H.LT.0.D0) THEN
           CALL INFGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,GRN,FLAG)
         ELSE
          CALL FINGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,WVN,NK,H,GRN,FLAG)
         ENDIF

         RKBN(IEL,JEL,IS,1)=SIJ
         RKBN(IEL,JEL,IS,2)=DIJ(1)
         RKBN(IEL,JEL,IS,3)=DIJ(2)
         RKBN(IEL,JEL,IS,4)=DIJ(3)
         CGRN(IEL,JEL,IS,:)=GRN(:)
      
        ENDDO
        ENDDO
      ENDDO
        
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      RETURN
      END SUBROUTINE CALGREEN 
      
!   ----------------------------------------------------------------------------
!      Calculate WAVE TERM's value of Green function for arbitray two points
!   ----------------------------------------------------------------------------

      SUBROUTINE CALGREEN_IRR
      IMPLICIT NONE
      
      INTEGER  IEL,JEL,IS,IP,IRR,FLAG
      REAL*8   XQ,YQ,ZQ,XP,YP,ZP,SIJ,DIJ(3),DIST
      COMPLEX*16  GRN(4),DPOX,DPOY,DPOZ

      IRR=1
      
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,IS,XQ,XP,YQ,YP,ZQ,ZP,SIJ,DIJ,GRN,FLAG)
         
      DO IEL=1, NELEM

        DO JEL=1, NELEM
            
         XQ=XYZ_P(JEL,1)
         YQ=XYZ_P(JEL,2)
         ZQ=XYZ_P(JEL,3)

         DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF
         
        DO IS=1, NSYS

         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          XP=SY(IS,1)*XYZ_P(IEL,1)
          YP=SX(IS,1)*XYZ_P(IEL,2)
          ZP=          XYZ_P(IEL,3)
         ELSE
          XP=SX(IS,1)*XYZ_P(IEL,1)
          YP=SY(IS,1)*XYZ_P(IEL,2)
          ZP=         XYZ_P(IEL,3)
         ENDIF

         IF (NCN(JEL).EQ.3) THEN

          CALL SGLINTBD_TRI(IS,IEL,JEL,SIJ,DIJ,IRR)

         ELSEIF (NCN(JEL).EQ.4) THEN

          CALL SGLINTBD_QUAD(IS,IEL,JEL,SIJ,DIJ,IRR)

         ENDIF 

         IF (H.LT.0.D0) THEN
           CALL INFGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,GRN,FLAG)
         ELSE
          CALL FINGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,WVN,NK,H,GRN,FLAG)
         ENDIF

         RKBN(IEL,JEL,IS,1)=SIJ
         RKBN(IEL,JEL,IS,2)=DIJ(1)
         RKBN(IEL,JEL,IS,3)=DIJ(2)
         RKBN(IEL,JEL,IS,4)=DIJ(3)
         CGRN(IEL,JEL,IS,:)=GRN(:)

        ENDDO
        ENDDO
      ENDDO
        
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

      IRR=3

!$OMP PARALLEL DO PRIVATE(IEL,JEL,IS,XQ,XP,YQ,YP,ZQ,ZP,SIJ,DIJ,GRN,FLAG)

      DO IEL=1, iNELEM
          
        DO IS=1, NSYS

        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XP=SY(IS,1)*iXYZ_P(IEL,1)
         YP=SX(IS,1)*iXYZ_P(IEL,2)
         ZP=          iXYZ_P(IEL,3)
        ELSE
         XP=SX(IS,1)*iXYZ_P(IEL,1)
         YP=SY(IS,1)*iXYZ_P(IEL,2)
         ZP=         iXYZ_P(IEL,3)
        ENDIF

        DO JEL=1, NELEM

         XQ=XYZ_P(JEL,1)
         YQ=XYZ_P(JEL,2)
         ZQ=XYZ_P(JEL,3)
         
         DIST=SQRT((XP-XQ)**2+(YP-YQ)**2+(ZP-ZQ)**2)
         IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF
         
         IF (NCN(JEL).EQ.3) THEN
             
          CALL SGLINTBD_TRI(IS,IEL,JEL,SIJ,DIJ,IRR)
          
         ELSEIF (NCN(JEL).EQ.4) THEN
             
          CALL SGLINTBD_QUAD(IS,IEL,JEL,SIJ,DIJ,IRR)

         ENDIF 

         IF (H.LT.0.D0) THEN
           CALL INFGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,GRN,FLAG)
         ELSE
          CALL FINGREEN3D(XQ,XP,YQ,YP,ZQ,ZP,V,WVN,NK,H,GRN,FLAG)
         ENDIF

         PKBN(IEL,JEL,IS,1)=SIJ
         PKBN(IEL,JEL,IS,2)=DIJ(1)
         PKBN(IEL,JEL,IS,3)=DIJ(2)
         PKBN(IEL,JEL,IS,4)=DIJ(3)
         DGRN(IEL,JEL,IS,:)=GRN(:)
         
        ENDDO
        ENDDO
      ENDDO

!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE CALGREEN_IRR
!-------------------------------------------------------------------------------
END MODULE CalGreenFunc
!*******************************************************************************