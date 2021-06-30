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

MODULE PressureElevation
    
   USE HAMS_mod
   USE Body_mod
   USE Const_mod
   USE WaveDyn_mod
   USE PanelMesh_mod
   USE Inerfs_mod
   USE Potentials_mod

   USE INFG3D_Open
   USE FinGreen3D_Open
   USE ImplementSubs
   USE SingularIntgr
   USE PatcVelct
   USE FieldOutput_mod

   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: CalPotential
   PUBLIC :: CalPressure
   PUBLIC :: CalElevation
   
CONTAINS
! -----------------------------------------------------------------------------------------------------
!       Compute the velocity potential due to the separate radiation or diffraction modes
! -----------------------------------------------------------------------------------------------------
      SUBROUTINE CalPotential(XET,RDFLG,MD,POT)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD
      REAL*8,INTENT(IN):: XET(3)
      COMPLEX*16,INTENT(OUT):: POT
      CHARACTER(*),INTENT(IN)::RDFLG

      INTEGER  JEL,IS,IRR,FLAG
      REAL*8   XQ(3),XP(3),XT(3),DIST,EAR,RKN(4),ENV(3),ENT(3),SLD
      COMPLEX*16  F0,DPOX,DPOY,DPOZ,DINCP,TERM1,TERM2,GRN(4),DUM(2)
      COMPLEX*16, ALLOCATABLE:: XPOT(:)
      
      ALLOCATE(XPOT(NELEM))
        
      IRR=1
      XPOT=DCMPLX(0.0D0, 0.0D0)

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(JEL,XQ,ENV,EAR,IS,XP,DIST,FLAG,RKN,GRN,DUM,XT,ENT,DPOX,DPOY,DPOZ,DINCP,TERM1,TERM2) !$OMP REDUCTION(+:XPOT)
      
      DO JEL=1, NELEM

       XQ(1)=XYZ_P(JEL,1)     ! XQ: source point,  XP: field point
       XQ(2)=XYZ_P(JEL,2)
       XQ(3)=XYZ_P(JEL,3)
       ENV(1)=DXYZ_P(JEL,1)
       ENV(2)=DXYZ_P(JEL,2)
       ENV(3)=DXYZ_P(JEL,3)
       EAR=DS(JEL)

       DO IS=1, NSYS

        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XP(1)=SY(IS,1)*XET(1)
         XP(2)=SX(IS,1)*XET(2)
         XP(3)=         XET(3)
        ELSE
         XP(1)=SX(IS,1)*XET(1)
         XP(2)=SY(IS,1)*XET(2)
         XP(3)=         XET(3)
        ENDIF

        DIST=SQRT((XP(1)-XQ(1))**2+(XP(2)-XQ(2))**2+(XP(3)-XQ(3))**2)
        IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF

         IF (NCN(JEL).EQ.3) THEN
          CALL SGLINTBD_TRI2(XP,JEL,RKN,IRR)
         ELSEIF (NCN(JEL).EQ.4) THEN
          CALL SGLINTBD_QUAD2(XP,JEL,RKN,IRR)
         ENDIF

         IF (H.LT.0.D0) THEN
          CALL INFGREEN3D(XQ(1),XP(1),XQ(2),XP(2),XQ(3),XP(3),V,GRN,FLAG)
         ELSE
          CALL FINGREEN3D(XQ(1),XP(1),XQ(2),XP(2),XQ(3),XP(3),V,WVN,NK,H,GRN,FLAG)
         ENDIF

         IF (FLAG.EQ.1) THEN
          DUM(1)= RKN(1)+GRN(1)*EAR
          DUM(2)=(RKN(2)+GRN(2)*EAR)*ENV(1)   &
                +(RKN(3)+GRN(3)*EAR)*ENV(2)   &
              +(RKN(4)+GRN(4)*EAR)*ENV(3)
         ELSE
          DUM(1)= GRN(1)*EAR
          DUM(2)=(GRN(2)*ENV(1)+GRN(3)*ENV(2)+GRN(4)*ENV(3))*EAR
         ENDIF

         IF (MD.EQ.7) THEN
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN   
           XT(1)=SY(1,IS)*XQ(1)
           XT(2)=SX(1,IS)*XQ(2)
           XT(3)=         XQ(3)
           ENT(1)=SY(1,IS)*ENV(1)
           ENT(2)=SX(1,IS)*ENV(2)
           ENT(3)=         ENV(3)
          ELSE
           XT(1)=SX(1,IS)*XQ(1)
           XT(2)=SY(1,IS)*XQ(2)
           XT(3)=         XQ(3)
           ENT(1)=SX(1,IS)*ENV(1)
           ENT(2)=SY(1,IS)*ENV(2)
           ENT(3)=         ENV(3)
          ENDIF
          CALL DINP(XT(1),XT(2),XT(3),XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)   ! Calculate the diffraction potential
          DINCP= DPOX*ENT(1)+DPOY*ENT(2)+DPOZ*ENT(3)
          TERM1=-DUM(1)*DINCP
          TERM2= DUM(2)*MXPOT(JEL,MD,IS)
         ELSEIF (MD.EQ.1.or.MD.EQ.5) THEN
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           TERM1=DUM(1)*SY(1,IS)*DXYZ_P(JEL,MD)
          ELSE
           TERM1=DUM(1)*SX(1,IS)*DXYZ_P(JEL,MD)
          ENDIF
          TERM2=DUM(2)*MXPOT(JEL,MD,IS)
         ELSEIF (MD.EQ.2.or.MD.EQ.4) THEN
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           TERM1=DUM(1)*SX(1,IS)*DXYZ_P(JEL,MD)
          ELSE
           TERM1=DUM(1)*SY(1,IS)*DXYZ_P(JEL,MD)
          ENDIF
          TERM2=DUM(2)*MXPOT(JEL,MD,IS)
         ELSEIF (MD.EQ.3) THEN
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           TERM1=DUM(1)         *DXYZ_P(JEL,MD)
          ELSE
           TERM1=DUM(1)         *DXYZ_P(JEL,MD)
          ENDIF
          TERM2=DUM(2)*MXPOT(JEL,MD,IS)
         ELSEIF (MD.EQ.6) THEN
          IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
           TERM1=DUM(1)*SX(1,IS)*SY(1,IS)*DXYZ_P(JEL,MD)
          ELSE
           TERM1=DUM(1)*SX(1,IS)*SY(1,IS)*DXYZ_P(JEL,MD)
          ENDIF
          TERM2=DUM(2)*MXPOT(JEL,MD,IS)
         ENDIF

        IF (ISOL.EQ.1) THEN
         XPOT(JEL)=XPOT(JEL)+TERM1-TERM2
        ELSEIF (ISOL.EQ.2) THEN
         XPOT(JEL)=XPOT(JEL)-TERM2 
        ENDIF

       ENDDO
      ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL
      
      POT=0.D0
      DO JEL=1, NELEM
        POT=POT+XPOT(JEL)
      ENDDO

      SLD=4.D0*PI
      POT=POT/SLD
      
      IF (MD.EQ.7) THEN
        IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
         F0=0.D0
        ELSE
         F0=VINP(XET(1),XET(2),XET(3),XW(1),XW(2),BETA)
        ENDIF
        POT=POT+F0
      ENDIF

      DEALLOCATE(XPOT)
      
      RETURN
      END SUBROUTINE CalPotential

! -----------------------------------------------------------------------------------------------------
!       Compute the pressure at a point due to the separate radiation or diffraction modes
! ----------------------------------------------------------------------------------------------------- 
      SUBROUTINE CalPressure(XP,RDFLG,MD,PRS)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD
      REAL*8,INTENT(IN):: XP(3)
      COMPLEX*16,INTENT(OUT):: PRS
      CHARACTER(*),INTENT(IN)::RDFLG
      
      COMPLEX*16  XPOT
      
      CALL CalPotential(XP,RDFLG,MD,XPOT)
      
      IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
       PRS=RHO*XPOT
      ELSE
       PRS=CI*W1*RHO*XPOT
      ENDIF

      RETURN
      END SUBROUTINE CalPressure

! -----------------------------------------------------------------------------------------------------
!    Compute the free-surface elevation at a point due to the separate radiation or diffraction modes
! ----------------------------------------------------------------------------------------------------- 
      SUBROUTINE CalElevation(XP,RDFLG,MD,ELV)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD
      REAL*8,INTENT(IN):: XP(3)
      COMPLEX*16,INTENT(OUT):: ELV
      CHARACTER(*),INTENT(IN)::RDFLG
      
      COMPLEX*16  XPOT

      CALL CalPotential(XP,RDFLG,MD,XPOT)
      
      IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
       ELV=XPOT
      ELSE
       ELV=CI*W1/G*XPOT
      ENDIF

      RETURN
      END SUBROUTINE CalElevation
      
! -----------------------------------------------------------------------------------------------------
!    Compute the free-surface elevation at a point due to the separate radiation or diffraction modes
! ----------------------------------------------------------------------------------------------------- 
      SUBROUTINE WamitNondimens(VCP,PEFLG,RDFLG,MD,NVCP)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: MD
      INTEGER MEXP
      
      COMPLEX*16,INTENT(IN):: VCP
      COMPLEX*16,INTENT(OUT):: NVCP
      CHARACTER(*),INTENT(IN)::PEFLG,RDFLG
      REAL*8 NFAC

      IF (adjustl(trim(PEFLG)).EQ.'Pressure') THEN
       IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
        NFAC=RHO*AMP
       ELSE
        NFAC=RHO*G*AMP
       ENDIF
      ELSEIF (adjustl(trim(PEFLG)).EQ.'Elevation') THEN
       NFAC=AMP
      ENDIF

      IF (adjustl(trim(RDFLG)).EQ.'Diffraction') THEN
       NVCP=VCP/NFAC
      ELSEIF (adjustl(trim(RDFLG)).EQ.'Radiation') THEN
       IF (MD.LE.3) THEN
        MEXP=0
       ELSEIF (MD.GE.4) THEN
        MEXP=1
       ENDIF
       IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
        NVCP=VCP/NFAC*AMP/REFL**(MEXP+1)
       ELSE
        NVCP=VCP/NFAC*(-CI*W1)*AMP/REFL**MEXP
        !PRINT*, NFAC,REFL,MEXP
        !PAUSE
       ENDIF
      ENDIF

      IF (ABS(NVCP).LT.1.E-15) NVCP=0.D0
      
      IF (adjustl(trim(RDFLG)).EQ.'Diffraction') THEN
       NVCP=CMPLX(-IMAG(NVCP),-REAL(NVCP))
      ELSEIF (adjustl(trim(RDFLG)).EQ.'Radiation') THEN
       NVCP=CMPLX(REAL(NVCP),-IMAG(NVCP))
      ENDIF
          
      !NVCP=CONJG (NVCP)     ! Because of the Wamit format

      RETURN
      END SUBROUTINE WamitNondimens
      
! -----------------------------------------------------------
!    Output pressures and elevations into the WAMIT format
! -----------------------------------------------------------
      SUBROUTINE OutputPressureElevation_Radiation(NFILE)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: NFILE
      INTEGER IPT,MD,EMD,IHD
      REAL*8 XP(3)
      COMPLEX*16 VCP
      COMPLEX*16,ALLOCATABLE:: VCPX(:,:)
      
      ALLOCATE(VCPX(NFP,6))

      DO IPT=1,NFP
       XP=XFP(IPT,:)
       DO MD=1,6
        IF (ABS(XP(3)).GT.1.E-6) THEN
         CALL CalPressure(XP,'Radiation',MD,VCP)
         CALL WamitNondimens(VCP,'Pressure','Radiation',MD,VCPX(IPT,MD))
        ELSE
         CALL CalElevation(XP,'Radiation',MD,VCP)
         CALL WamitNondimens(VCP,'Elevation','Radiation',MD,VCPX(IPT,MD))
        ENDIF
       ENDDO
      
!   ===================================================
!    Write WAMIT-style output files
!   
      IF (ABS(XP(3)).GT.1.E-6) THEN
       WRITE(NFILE,1000) OUFR,IPT,(VCPX(IPT,MD),MD=1,6)
      ENDIF
      
      ENDDO
      
      DEALLOCATE(VCPX)
      
1000  FORMAT(ES14.6,I10,12ES14.6)
      RETURN
      END SUBROUTINE OutputPressureElevation_Radiation
      
! -----------------------------------------------------------
!    Output pressures and elevations into the WAMIT format
! -----------------------------------------------------------
      SUBROUTINE OutputPressureElevation_Diffraction(NFILE)
      IMPLICIT NONE
!
      INTEGER,INTENT(IN):: NFILE
      INTEGER IPT,MD,EMD,IHD
      REAL*8 XP(3),REL,IMG,MOD,PHS
      COMPLEX*16 VCP,NVCP

      DO IPT=1,NFP
       XP=XFP(IPT,:)
       IF (ABS(XP(3)).GT.1.E-6) THEN
        CALL CalPressure(XP,'Diffraction',7,VCP)
        CALL WamitNondimens(VCP,'Pressure','Diffraction',0,NVCP)
       ELSE
        CALL CalElevation(XP,'Diffraction',7,VCP)
        CALL WamitNondimens(VCP,'Elevation','Diffraction',0,NVCP)
       ENDIF
       
       !WRITE(NFILE,1020) OUFR,BETA*180.0D0/PI,IPT,NVCP
       REL=REAL(NVCP)
       IMG=IMAG(NVCP)
       MOD=ABS(NVCP)
       PHS=ATAN2D(IMG,REL)
       WRITE(NFILE,1020) OUFR,BETA*180.0D0/PI,IPT,MOD,PHS,REL,IMG
       
      ENDDO
      
1020  FORMAT(2ES14.6,I10,4ES14.6)
      RETURN
      END SUBROUTINE OutputPressureElevation_Diffraction
!-------------------------------------------------------------------------------
END MODULE PressureElevation
!*******************************************************************************
