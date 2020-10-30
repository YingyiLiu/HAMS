!  ------------------------------------------------------------------------------------------------------
!                                                               
!    Program HAMS for the diffraction and radiation of waves 
!    for 3D structures by Constant Panel method.
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
      
        IF (IFWKO .EQ. 0)  THEN
            WK=WK1+(KK-1)*DWK
          IF (H .LE. 0.0D0) THEN
            W1=SQRT(G*WK)
            V=W1*W1/G
          ELSE
            W1=SQRT(G*WK*TANH(WK*H))
            V=W1*W1/G
            CALL DISPERSION(WVN,NK,DSQRT(G*V),H)
          END IF
            WL=2.D0*PI/WK
            TP=2.D0*PI/W1
        ELSEIF (IFWKO .EQ. 1)  THEN
            W1=WK1+(KK-1)*DWK
          IF(H .LE. 0.0D0) THEN
            WK=W1*W1/G
            V=WK
          ELSE
            V=W1*W1/G
            CALL DISPERSION(WVN,NK,DSQRT(G*V),H)
            WK=WVN(1)
          END IF
            WL=2.D0*PI/WK
            TP=2.D0*PI/W1
        ELSEIF (IFWKO .EQ. 2)  THEN
            TP=WK1+(KK-1)*DWK
            W1=2.D0*PI/TP
          IF(H .LE. 0.0D0) THEN
            WK=W1*W1/G
            V=WK
          ELSE
            V=W1*W1/G
            CALL DISPERSION(WVN,NK,DSQRT(G*V),H)
            WK=WVN(1)
          END IF
            WL=2.D0*PI/WK
        END IF

      END SUBROUTINE CalWaveProperts
!-------------------------------------------------------------------------------
END MODULE ImplementSubs
!*******************************************************************************