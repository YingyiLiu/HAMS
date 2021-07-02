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
MODULE ImplementSubs

   USE HAMS_mod
   USE Const_mod 
   USE WaveDyn_mod

   IMPLICIT NONE

   ! ..... Public Subroutines ...................................................................................................
   
   PUBLIC :: CalWaveProperts

CONTAINS

! ------------------------------------------------------------------- 
!    Calculate the incident wave properties
! ------------------------------------------------------------------- 
      SUBROUTINE CalWaveProperts(KK)
      IMPLICIT  NONE
!        
      INTEGER,INTENT(IN):: KK

        IF (SYBO.EQ.1) THEN
          IF (KK.EQ.1) THEN
            V=0.D0
            WK=0.D0
            W1=0.D0
            TP=-1.D0
            WL=-1.D0
            IF (INFT.LE.3) THEN
             INFR=0.D0
            ELSE
             INFR=-1.D0
            ENDIF
            GOTO 100
            ELSEIF (KK.EQ.2) THEN
            V=-1.D0
            WK=-1.D0
            W1=-1.D0
            TP=0.D0
            WL=0.D0
            IF (INFT.LE.3) THEN
             INFR=-1.D0
            ELSE
             INFR=0.D0
            ENDIF
            GOTO 100
          ENDIF
        ENDIF
        
        INFR=WVNB(KK)
        IF (INFT .EQ. 1)  THEN
            V=WVNB(KK)
          IF (H .LE. 0.0D0) THEN
            WK=V
            W1=SQRT(G*V)
          ELSE
            W1=SQRT(G*V)
            CALL DISPERSION(WVN,NK,W1,H)
            WK=WVN(1)
          END IF
            WL=2.D0*PI/WK
            TP=2.D0*PI/W1
        ELSEIF (INFT .EQ. 2)  THEN
            WK=WVNB(KK)
          IF (H .LE. 0.0D0) THEN
            V=WK
            W1=SQRT(G*V)
          ELSE
            V=WK*TANH(WK*H)
            W1=SQRT(G*V)
            CALL DISPERSION(WVN,NK,W1,H)
          END IF
            WL=2.D0*PI/WK
            TP=2.D0*PI/W1
        ELSEIF (INFT .EQ. 3)  THEN
            W1=WVNB(KK)
          IF(H .LE. 0.0D0) THEN
            V=W1*W1/G
            WK=V
          ELSE
            V=W1*W1/G
            CALL DISPERSION(WVN,NK,W1,H)
            WK=WVN(1)
          END IF
            WL=2.D0*PI/WK
            TP=2.D0*PI/W1
        ELSEIF (INFT .EQ. 4)  THEN
            TP=WVNB(KK)
            W1=2.D0*PI/TP
          IF(H .LE. 0.0D0) THEN
            V=W1*W1/G
            WK=V
          ELSE
            V=W1*W1/G
            CALL DISPERSION(WVN,NK,W1,H)
            WK=WVN(1)
          END IF
            WL=2.D0*PI/WK
        ELSEIF (INFT .EQ. 5)  THEN
            WL=WVNB(KK)
            WK=2.D0*PI/WL
          IF(H .LE. 0.0D0) THEN
            V=WK
            W1=SQRT(G*V)
          ELSE
            V=WK*TANH(WK*H)
            W1=SQRT(G*V)
            CALL DISPERSION(WVN,NK,W1,H)
          END IF
        END IF

100     CONTINUE
        WVNB(KK)=WK
        WVFQ(KK)=W1
        
        IF (OUFT .EQ. 1)  THEN
         OUFR=V
        ELSEIF (OUFT .EQ. 2)  THEN
         OUFR=WK
        ELSEIF (OUFT .EQ. 3)  THEN
         OUFR=W1
        ELSEIF (OUFT .EQ. 4)  THEN
         OUFR=TP
        ELSEIF (OUFT .EQ. 5)  THEN
         OUFR=WL
        END IF

      END SUBROUTINE CalWaveProperts
!-------------------------------------------------------------------------------
END MODULE ImplementSubs
!*******************************************************************************