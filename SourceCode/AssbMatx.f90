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
MODULE AssbMatx

   USE HAMS_mod
   USE Body_mod
   USE Const_mod 
   USE WaveDyn_mod
   USE PanelMesh_mod

   USE BodyIntgr
   USE PatcVelct
   IMPLICIT NONE

   PRIVATE
   
  ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ASSB_LEFT
   PUBLIC :: ASSB_RBC
   
   PUBLIC :: ASSB_DBC
   PUBLIC :: RADIATION_SOLVER
   PUBLIC :: DIFFRACTION_SOLVER
   
CONTAINS
! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the left-hand
!    side coefficient matrix [A]
! ------------------------------------------------------------------- 
      SUBROUTINE ASSB_LEFT(AMAT,IPIV,NELEM,NSYS)
      IMPLICIT  NONE
!     
      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER,INTENT(OUT):: IPIV(NELEM,NSYS)
      COMPLEX*16,INTENT(OUT):: AMAT(NELEM,NELEM,NSYS)
      
      INTEGER IEL,JEL,IS,IP,IRR,MD,FLAG,INFO
      REAL*8 DIST
      COMPLEX*16 TINDP(4)

      IRR=1
      AMAT=CMPLX(0.0D0,0.0D0)

!$OMP PARALLEL NUM_THREADS(NTHREAD)          
!$OMP DO PRIVATE(IEL,JEL,IP,IS,FLAG,DIST,TINDP) !$OMP REDUCTION(+:AMAT)
      
      DO  1000 IEL=1,  NELEM
        
       DO 100  IP=1,  NSYS

        AMAT(IEL,IEL,IP)=2.D0*PI
        
100    CONTINUE

       DO 200 JEL=1,  NELEM
            
        DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
          
        TINDP=(0.0D0, 0.0D0)
        
        DO  200   IS=1,  NSYS
        
         CALL BODINT_LEFT(IS,IEL,JEL,TINDP,FLAG)

         DO IP=1, NSYS
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           AMAT(IEL,JEL,IP)=AMAT(IEL,JEL,IP)+RXY(IS,IP)*TINDP(IS)
          ELSE
           AMAT(IEL,JEL,IP)=AMAT(IEL,JEL,IP)+RXY(IS,IP)*TINDP(IS)
          ENDIF
         ENDDO
        
200    CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$omp parallel do private(NTHREAD)
       DO IP=1, NSYS
          CALL ZGETRF( NELEM, NELEM, AMAT(:,:,IP), NELEM, IPIV(:,IP), INFO )
       ENDDO
!$omp end parallel do

       RETURN
      END SUBROUTINE ASSB_LEFT
      
      
! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the right-hand
!    side vector[B] and solve the equations
! ------------------------------------------------------------------- 
      SUBROUTINE ASSB_RBC(BRMAT,NELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,NSYS
      COMPLEX*16,INTENT(OUT):: BRMAT(NELEM,6,NSYS)
      
      INTEGER IEL,JEL,IS,IP,IRR,MD,FLAG
      REAL*8 DIST
      COMPLEX*16 TINRD(4,6,4),BTMP(6,4)

      IRR=1
      BRMAT=CMPLX(0.0D0,0.0D0)

!$OMP PARALLEL NUM_THREADS(NTHREAD)           
!$OMP DO PRIVATE(IEL,JEL,MD,IP,IS,FLAG,DIST,TINRD,BTMP) !$OMP REDUCTION(+:BRMAT)
      
      DO  1000 IEL=1,  NELEM
            
       BTMP=CMPLX(0.0D0,0.0D0)

        DO 200 JEL=1,  NELEM
            
         DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF

         TINRD=CMPLX(0.0D0, 0.0D0)
        
         DO  200   IS=1,  NSYS
        
          CALL RBC_RIGHT(IS,IEL,JEL,TINRD,FLAG)
        
          DO MD=1,  6
          DO IP=1, NSYS
            BTMP(MD,IP)=BTMP(MD,IP)+TINRD(IS,MD,IP)
          ENDDO
          ENDDO
        
200     CONTINUE
        
         DO  300  MD=1,  6
         DO  300  IP=1, NSYS

         DO  300  IS=1, NSYS
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY(IP,IS)*BTMP(MD,IS)
          ELSE
           BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY(IP,IS)*BTMP(MD,IS)
          ENDIF
300      CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL

       RETURN
      END SUBROUTINE ASSB_RBC
      
! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the right-hand
!    side vector[B] and solve the equations
! ------------------------------------------------------------------- 
      SUBROUTINE ASSB_DBC(BDMAT,NELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,NSYS
      COMPLEX*16,INTENT(OUT):: BDMAT(NELEM,NSYS)
      
      INTEGER IEL,JEL,IS,IP,IRR,MD,FLAG
      REAL*8 XP,YP,ZP,DIST
      COMPLEX*16 TINRD(4,4),BTMP(4)

      IRR=1
      MD=7
      BDMAT=CMPLX(0.0D0,0.0D0)

!$OMP PARALLEL NUM_THREADS(NTHREAD)           
!$OMP DO PRIVATE(XP,YP,ZP,IEL,JEL,IP,IS,FLAG,DIST,TINRD,BTMP) !$OMP REDUCTION(+:BDMAT)
      
      DO  1000 IEL=1,  NELEM
            
       BTMP=CMPLX(0.0D0,0.0D0)
       
       IF (ISOL.EQ.2) THEN
            
         DO 100  IP=1,  NSYS
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           XP=RX(IP,1)*XYZ_P(IEL,1)
           YP=RX(IP,2)*XYZ_P(IEL,2)
           ZP=         XYZ_P(IEL,3)
          ELSE
           XP=RY(IP,1)*XYZ_P(IEL,1)
           YP=RY(IP,2)*XYZ_P(IEL,2)
           ZP=         XYZ_P(IEL,3)
          ENDIF
          BTMP(IP)=4.D0*PI*VINP(XP,YP,ZP,XW(1),XW(2),BETA)
100      CONTINUE
         
       ELSEIF (ISOL.EQ.1) THEN       

        DO 200 JEL=1,  NELEM
            
         DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF
         TINRD=CMPLX(0.0D0, 0.0D0)
        
         DO  200   IS=1,  NSYS
        
           CALL DBC_RIGHT(IS,IEL,JEL,TINRD,FLAG)
        
           DO IP=1, NSYS
             BTMP(IP)=BTMP(IP)+TINRD(IS,IP)
           ENDDO
        
200     CONTINUE
        
       ELSE
           
        PRINT*,"  Error: The input for ISOL should be either 1 or 2."
        STOP
        
       ENDIF
       
        
         DO  300  IP=1, NSYS

         DO  300  IS=1, NSYS
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           BDMAT(IEL,IP)=BDMAT(IEL,IP)+RXY(IP,IS)*BTMP(IS) 
          ELSE
           BDMAT(IEL,IP)=BDMAT(IEL,IP)+RXY(IP,IS)*BTMP(IS)  
          ENDIF       
300      CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      RETURN
      END SUBROUTINE ASSB_DBC
      
! -----------------------------------------------------------------------------
!     Solve the radiation potentials for not removing irregular frequencies
! -----------------------------------------------------------------------------

      SUBROUTINE RADIATION_SOLVER(AMAT,BRMAT,IPIV,MXPOT,NELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER,INTENT(IN):: IPIV(NELEM,NSYS)
      COMPLEX*16,INTENT(IN):: AMAT(NELEM,NELEM,NSYS),BRMAT(NELEM,6,NSYS)
      COMPLEX*16,INTENT(OUT):: MXPOT(NELEM,7,NSYS)
      
      INTEGER IEL,IS,IP,MD,INFO
      COMPLEX*16,ALLOCATABLE:: ATMAT(:,:,:),BRTMAT(:,:,:)
      
      ALLOCATE(ATMAT(NELEM,NELEM,NSYS),BRTMAT(NELEM,6,NSYS))
      
      ATMAT=AMAT
      BRTMAT=BRMAT
      
!$omp parallel do private(NTHREAD)
      DO IP=1, NSYS
           CALL ZGETRS( 'No transpose', NELEM, 6, ATMAT(:,:,IP), NELEM, IPIV(:,IP), BRTMAT(:,:,IP), NELEM, INFO )
      ENDDO
!$omp end parallel do
      
      DO 1400 MD=1, 6
      DO 1400 IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,IS) !$OMP REDUCTION(+:MXPOT)
      DO 1450 IEL=1, NELEM
         MXPOT(IEL,MD,IP)=CMPLX(0.0D0, 0.0D0)
      DO 1460 IS=1, NSYS
         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+BRTMAT(IEL,MD,IS)*RXY(IP,IS)
         ELSE
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+BRTMAT(IEL,MD,IS)*RXY(IP,IS)
         ENDIF
1460  CONTINUE
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)/NSYS
1450  CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
1400  CONTINUE
      
      DEALLOCATE(ATMAT,BRTMAT)
      
      RETURN
      END SUBROUTINE RADIATION_SOLVER
      
! -----------------------------------------------------------------------------
!     Solve the diffraction potentials for not removing irregular frequencies
! -----------------------------------------------------------------------------

      SUBROUTINE DIFFRACTION_SOLVER(AMAT,BDMAT,IPIV,MXPOT,NELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER,INTENT(IN):: IPIV(NELEM,NSYS)
      COMPLEX*16,INTENT(IN):: AMAT(NELEM,NELEM,NSYS),BDMAT(NELEM,NSYS)
      COMPLEX*16,INTENT(OUT):: MXPOT(NELEM,7,NSYS)
      
      INTEGER IEL,IS,IP,MD,INFO
      COMPLEX*16,ALLOCATABLE:: ATMAT(:,:,:),BDTMAT(:,:)
      
      ALLOCATE(ATMAT(NELEM,NELEM,NSYS),BDTMAT(NELEM,NSYS))
      
      ATMAT=AMAT
      BDTMAT=BDMAT

      MD=7
      
!$omp parallel do private(NTHREAD)
      DO IP=1, NSYS
           CALL ZGETRS( 'No transpose', NELEM, 1, ATMAT(:,:,IP), NELEM, IPIV(:,IP), BDTMAT(:,IP), NELEM, INFO )
      ENDDO
!$omp end parallel do
      
       DO 400 IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,IS) !$OMP REDUCTION(+:MXPOT)
       DO 450 IEL=1, NELEM
         MXPOT(IEL,MD,IP)=CMPLX(0.0D0, 0.0D0)
       DO 460 IS=1, NSYS
         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+BDTMAT(IEL,IS)*RXY(IP,IS)
         ELSE
          MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+BDTMAT(IEL,IS)*RXY(IP,IS)
         ENDIF
460    CONTINUE
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)/NSYS
450    CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
400    CONTINUE

      DEALLOCATE(ATMAT,BDTMAT)
      
      RETURN
      END SUBROUTINE DIFFRACTION_SOLVER
!-------------------------------------------------------------------------------
END MODULE AssbMatx
!*******************************************************************************