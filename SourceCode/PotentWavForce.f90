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
MODULE PotentWavForce

   USE HAMS_mod
   USE Body_mod
   USE Const_mod
   USE PanelMesh_mod
   USE Inerfs_mod

   USE PatcVelct
   USE Potentials_mod
   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: EFORCE
   PUBLIC :: RFORCE
   
CONTAINS
! ------------------------------------------------------------------- 
!          Compute the wave exciting force on a 3D body
! ------------------------------------------------------------------- 
      SUBROUTINE EFORCE(WK,W1,TP,BETA,AMP,EXFC)
      IMPLICIT NONE
!
      REAL*8,INTENT(IN):: WK,W1,TP,BETA,AMP
      COMPLEX*16,INTENT(OUT):: EXFC(6)
      
      INTEGER IEL,IP,MD
      REAL*8:: XP,YP,ZP,AMFJ(6)
      REAL*8:: MOD,PHS(6),REL,IMG,NREL,NIMG
      COMPLEX*16 PHI,FORCE(6,4)

      MD=7
      FORCE=CMPLX(0.0D0, 0.0D0)
      
!!$CALL OMP_SET_NUM_THREADS(NTHREAD)
!!$OMP PARALLEL DO PRIVATE(IEL,IP,XP,YP,ZP,PHI) !$OMP REDUCTION(+:FORCE)

      DO 100  IEL=1,  NELEM
      DO 100  IP=1,  NSYS
          
       IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
        XP=RX(IP,1)*XYZ_P(IEL,1)
        YP=RX(IP,2)*XYZ_P(IEL,2)
        ZP=         XYZ_P(IEL,3)
        IF (ISOL.EQ.1) THEN
         PHI=(MXPOT(IEL,MD,IP)+VINP(XP,YP,ZP,XW(1),XW(2),BETA))*DS(IEL)
        ELSEIF (ISOL.EQ.2) THEN
         PHI=MXPOT(IEL,MD,IP)*DS(IEL)
        ENDIF
        FORCE(1,IP)=FORCE(1,IP)+PHI*RX(IP,1)*DXYZ_P(IEL,1)
        FORCE(2,IP)=FORCE(2,IP)+PHI*RX(IP,2)*DXYZ_P(IEL,2)
        FORCE(3,IP)=FORCE(3,IP)+PHI*         DXYZ_P(IEL,3)
        FORCE(4,IP)=FORCE(4,IP)+PHI*RX(IP,2)*DXYZ_P(IEL,4)
        FORCE(5,IP)=FORCE(5,IP)+PHI*RX(IP,1)*DXYZ_P(IEL,5)
        FORCE(6,IP)=FORCE(6,IP)+PHI*RX(IP,1)*RX(IP,2)*DXYZ_P(IEL,6)
       ELSE
        XP=RY(IP,1)*XYZ_P(IEL,1)
        YP=RY(IP,2)*XYZ_P(IEL,2)
        ZP=         XYZ_P(IEL,3)
        IF (ISOL.EQ.1) THEN
         PHI=(MXPOT(IEL,MD,IP)+VINP(XP,YP,ZP,XW(1),XW(2),BETA))*DS(IEL)
        ELSEIF (ISOL.EQ.2) THEN
         PHI=MXPOT(IEL,MD,IP)*DS(IEL)
        ENDIF
        FORCE(1,IP)=FORCE(1,IP)+PHI*RY(IP,1)*DXYZ_P(IEL,1)
        FORCE(2,IP)=FORCE(2,IP)+PHI*RY(IP,2)*DXYZ_P(IEL,2)
        FORCE(3,IP)=FORCE(3,IP)+PHI*          DXYZ_P(IEL,3)
        FORCE(4,IP)=FORCE(4,IP)+PHI*RY(IP,2)*DXYZ_P(IEL,4)
        FORCE(5,IP)=FORCE(5,IP)+PHI*RY(IP,1)*DXYZ_P(IEL,5)
        FORCE(6,IP)=FORCE(6,IP)+PHI*RY(IP,1)*RY(IP,2)*DXYZ_P(IEL,6)
       ENDIF

100   CONTINUE
!!$OMP END PARALLEL DO
! -------------------------------------------------
!
      DO 200 MD=1, 6
        EXFC(MD)=CMPLX(0.D0,0.D0)
      DO  200  IP=1,  NSYS
	    EXFC(MD)=EXFC(MD)+FORCE(MD,IP)
200   CONTINUE

      EXFC(:)=CI*W1*RHO*EXFC(:)
      AMFJ(:)=ABS(EXFC(:))

      DO MD=1,6
          
       REL=REAL(EXFC(MD))/(RHO*G*AMP)
       IMG=IMAG(EXFC(MD))/(RHO*G*AMP)
       MOD=ABS(EXFC(MD))/(RHO*G*AMP) !SQRT(REL**2+MDMG**2)
       NREL=-IMG
       NIMG= REL
       PHS(MD)=ATAN2D(NIMG,NREL)
       WRITE(20+MD,1010) WK,W1,REAL(EXFC(MD)),IMAG(EXFC(MD))
       
       IF (ABS(TP+1.D0).GT.1.E-6.AND.ABS(TP).GT.1.E-6) THEN
        WRITE(62,1030)  OUFR,BETA*180.0D0/PI,MD,MOD,PHS(MD),NREL,NIMG
       ENDIF
           
      ENDDO
!
!   =================================================== 
!
1010  FORMAT(F7.3,1x,F7.3,1x,6E14.6)
1030  FORMAT(2ES14.6,I6,4ES14.6)
      
      RETURN 
      END SUBROUTINE EFORCE 

! ------------------------------------------------------------------- 
!          Compute the wave radiation force on a 3D body
! ------------------------------------------------------------------- 
      SUBROUTINE RFORCE(WK,W1,TP,AMAS,BDMP)
      IMPLICIT   NONE
!
      REAL*8,INTENT(IN):: WK,W1,TP
      REAL*8,INTENT(OUT):: AMAS(6,6),BDMP(6,6)
      
      INTEGER IEL,MD,MD1,MD2,IP
      REAL*8:: AMFJ(6),NAMAS(6,6),NBDMP(6,6)
      COMPLEX*16 RPHI,IPHI

      AMAS(:,:)=0.D0
      BDMP(:,:)=0.D0
      
!!$CALL OMP_SET_NUM_THREADS(NTHREAD)
!!$OMP PARALLEL DO PRIVATE(IEL,IP,MD,RPHI,IPHI) !$OMP REDUCTION(+:AMAS,BDMP)

      DO 100  MD=1,6
        DO 100  IEL=1,  NELEM
        DO 100  IP=1,  NSYS

         IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
          RPHI=CI*RHO*MXPOT(IEL,MD,IP)*DS(IEL)
          IPHI=CI*RHO*W1*MXPOT(IEL,MD,IP)*DS(IEL)

          AMAS(1,MD)=AMAS(1,MD)+IMAG(RPHI*RX(IP,1)*DXYZ_P(IEL,1))
          AMAS(2,MD)=AMAS(2,MD)+IMAG(RPHI*RX(IP,2)*DXYZ_P(IEL,2))
          AMAS(3,MD)=AMAS(3,MD)+IMAG(RPHI*           DXYZ_P(IEL,3))
          AMAS(4,MD)=AMAS(4,MD)+IMAG(RPHI*RX(IP,2)*DXYZ_P(IEL,4))
          AMAS(5,MD)=AMAS(5,MD)+IMAG(RPHI*RX(IP,1)*DXYZ_P(IEL,5))
          AMAS(6,MD)=AMAS(6,MD)+IMAG(RPHI*RX(IP,1)*RX(IP,2)*DXYZ_P(IEL,6))
         
          BDMP(1,MD)=BDMP(1,MD)-DBLE(IPHI*RX(IP,1)*DXYZ_P(IEL,1))
          BDMP(2,MD)=BDMP(2,MD)-DBLE(IPHI*RX(IP,2)*DXYZ_P(IEL,2))
          BDMP(3,MD)=BDMP(3,MD)-DBLE(IPHI*           DXYZ_P(IEL,3))
          BDMP(4,MD)=BDMP(4,MD)-DBLE(IPHI*RX(IP,2)*DXYZ_P(IEL,4))
          BDMP(5,MD)=BDMP(5,MD)-DBLE(IPHI*RX(IP,1)*DXYZ_P(IEL,5))
          BDMP(6,MD)=BDMP(6,MD)-DBLE(IPHI*RX(IP,1)*RX(IP,2)*DXYZ_P(IEL,6))
         ELSE
          RPHI=CI*RHO*MXPOT(IEL,MD,IP)*DS(IEL)
          IPHI=CI*RHO*W1*MXPOT(IEL,MD,IP)*DS(IEL)

          AMAS(1,MD)=AMAS(1,MD)+IMAG(RPHI*RY(IP,1)*DXYZ_P(IEL,1))
          AMAS(2,MD)=AMAS(2,MD)+IMAG(RPHI*RY(IP,2)*DXYZ_P(IEL,2))
          AMAS(3,MD)=AMAS(3,MD)+IMAG(RPHI*          DXYZ_P(IEL,3))  
          AMAS(4,MD)=AMAS(4,MD)+IMAG(RPHI*RY(IP,2)*DXYZ_P(IEL,4))
          AMAS(5,MD)=AMAS(5,MD)+IMAG(RPHI*RY(IP,1)*DXYZ_P(IEL,5))
          AMAS(6,MD)=AMAS(6,MD)+IMAG(RPHI*RY(IP,1)*RY(IP,2)*DXYZ_P(IEL,6)) 
         
          BDMP(1,MD)=BDMP(1,MD)-DBLE(IPHI*RY(IP,1)*DXYZ_P(IEL,1))
          BDMP(2,MD)=BDMP(2,MD)-DBLE(IPHI*RY(IP,2)*DXYZ_P(IEL,2))
          BDMP(3,MD)=BDMP(3,MD)-DBLE(IPHI*          DXYZ_P(IEL,3))  
          BDMP(4,MD)=BDMP(4,MD)-DBLE(IPHI*RY(IP,2)*DXYZ_P(IEL,4))
          BDMP(5,MD)=BDMP(5,MD)-DBLE(IPHI*RY(IP,1)*DXYZ_P(IEL,5))
          BDMP(6,MD)=BDMP(6,MD)-DBLE(IPHI*RY(IP,1)*RY(IP,2)*DXYZ_P(IEL,6)) 
         ENDIF

100   CONTINUE
      
!!$OMP END PARALLEL DO

       DO MD1=1,6
        WRITE(30+MD1,1000) WK,W1,(AMAS(MD1,MD2),MD2=1,6)
        WRITE(40+MD1,1000) WK,W1,(BDMP(MD1,MD2),MD2=1,6)
       ENDDO  
!
!   =================================================== 
!    Write WAMIT-style output files
!
       DO MD1=1,6
          DO MD2=1,6
            IF (ABS(TP+1.D0).LT.1.E-6.OR.ABS(TP).LT.1.E-6) THEN
             WRITE(61,1020) OUFR,MD1,MD2,AMAS(MD1,MD2)/RHO
            ELSE
             WRITE(61,1020) OUFR,MD1,MD2,AMAS(MD1,MD2)/RHO,BDMP(MD1,MD2)/(RHO*W1)
            ENDIF
          ENDDO
       ENDDO
!   =================================================== 
!    Write HydroStar-style output files
!

1000    FORMAT(F7.3,1x,F7.3,1x,6E14.5)
1020    FORMAT(ES14.6,2I6,2ES14.6)

      RETURN
      END SUBROUTINE RFORCE
!-------------------------------------------------------------------------------
END MODULE PotentWavForce
!*******************************************************************************