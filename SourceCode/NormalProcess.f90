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
MODULE NormalProcess

   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: CalPanelCentre
   PUBLIC :: CalPanelArea
   PUBLIC :: CalPanelArea_TransNormal
   PUBLIC :: CalPanelArea_Normal
   PUBLIC :: CalRotNormals
   PUBLIC :: CalDeltaArea
   PUBLIC :: CalDeltaArea_Improved
   
CONTAINS
!=======================================================================
!SUBROUTINE CalRotNormals


!   !   CalNormals is used to calculate normals on the immersed body surface 
!   !   and on the inner water plane


      SUBROUTINE CalRotNormals(RC,XYZ_P,DXYZ_P,NELEM)
      IMPLICIT   NONE
      
      INTEGER,INTENT(IN):: NELEM
      REAL*8,INTENT(IN):: RC(3),XYZ_P(NELEM,3)
      REAL*8,INTENT(INOUT):: DXYZ_P(NELEM,6)
      
      INTEGER IEL
      
	  DO 20 IEL=1,NELEM

          DXYZ_P(IEL,4)=(XYZ_P(IEL,2)-RC(2))*DXYZ_P(IEL,3)-(XYZ_P(IEL,3)-RC(3))*DXYZ_P(IEL,2)
          DXYZ_P(IEL,5)=(XYZ_P(IEL,3)-RC(3))*DXYZ_P(IEL,1)-(XYZ_P(IEL,1)-RC(1))*DXYZ_P(IEL,3)
          DXYZ_P(IEL,6)=(XYZ_P(IEL,1)-RC(1))*DXYZ_P(IEL,2)-(XYZ_P(IEL,2)-RC(2))*DXYZ_P(IEL,1)

20    CONTINUE 

    RETURN
    END SUBROUTINE CalRotNormals
    
!=======================================================================
!SUBROUTINE CalTransNormals


!   !   CalTransNormals is used to calculate normals on the immersed body surface 
!   !   and on the inner water plane


      SUBROUTINE CalTransNormals(XYZ,NTND,NELEM,NCN,NCON,DXYZ_P)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NTND,NELEM
      INTEGER,INTENT(IN):: NCN(NELEM),NCON(NELEM,4)
      REAL*8,INTENT(IN):: XYZ(NTND,3)
      REAL*8,INTENT(OUT):: DXYZ_P(NELEM,6)
      
      INTEGER IEL
      REAL*8::  V21(3),V23(3),V13(3),V24(3),UN(3),AUN
!
! -------------------------------------------------------------------------
!     
	  DO IEL=1,NELEM
        
        IF (NCN(IEL).EQ.3) THEN
            
         V21(:)=XYZ(NCON(IEL,1),:)-XYZ(NCON(IEL,2),:)
         V23(:)=XYZ(NCON(IEL,3),:)-XYZ(NCON(IEL,2),:)

         UN(1)=V21(2)*V23(3)-V21(3)*V23(2)
         UN(2)=V21(3)*V23(1)-V21(1)*V23(3)
         UN(3)=V21(1)*V23(2)-V21(2)*V23(1)

         AUN=SQRT(UN(1)**2+UN(2)**2+UN(3)**2)

         DXYZ_P(IEL,1)=UN(1)/AUN
         DXYZ_P(IEL,2)=UN(2)/AUN
         DXYZ_P(IEL,3)=UN(3)/AUN
                   
        ELSEIF (NCN(IEL).EQ.4) THEN
            
         V13(:)=XYZ(NCON(IEL,3),:)-XYZ(NCON(IEL,1),:)
         V24(:)=XYZ(NCON(IEL,4),:)-XYZ(NCON(IEL,2),:)
         
         UN(1)=V24(2)*V13(3)-V24(3)*V13(2)
         UN(2)=V24(3)*V13(1)-V24(1)*V13(3)
         UN(3)=V24(1)*V13(2)-V24(2)*V13(1)

         AUN=SQRT(UN(1)**2+UN(2)**2+UN(3)**2)
         
         DXYZ_P(IEL,1)=UN(1)/AUN
         DXYZ_P(IEL,2)=UN(2)/AUN
         DXYZ_P(IEL,3)=UN(3)/AUN
        
        ENDIF
        
     ENDDO
!
    RETURN
    END SUBROUTINE CalTransNormals
    
!=======================================================================
!SUBROUTINE CalPanelCentre


!   !   CalNormals is used to calculate normals on the immersed body surface 
!   !   and on the inner water plane


      SUBROUTINE CalPanelCentre(XYZ,NTND,NELEM,NCN,NCON,XYZ_P)
      IMPLICIT   NONE
	  
      INTEGER,INTENT(IN):: NTND,NELEM
      INTEGER,INTENT(IN):: NCN(NELEM),NCON(NELEM,4)
      REAL*8,INTENT(IN):: XYZ(NTND,3)
      REAL*8,INTENT(OUT):: XYZ_P(NELEM,3)
      
      INTEGER IEL,IND,J
!
! -------------------------------------------------------------------------
! 
      XYZ_P=0.D0

      DO IEL=1,NELEM
       DO J=1, NCN(IEL)
	    IND=NCON(IEL,J)
	    XYZ_P(IEL,1)=XYZ_P(IEL,1)+XYZ(IND,1)
	    XYZ_P(IEL,2)=XYZ_P(IEL,2)+XYZ(IND,2)
	    XYZ_P(IEL,3)=XYZ_P(IEL,3)+XYZ(IND,3)
	   ENDDO
	    XYZ_P(IEL,:)=XYZ_P(IEL,:)/NCN(IEL)
      ENDDO     
!
     RETURN
     END SUBROUTINE CalPanelCentre
      

    
!=======================================================================
!SUBROUTINE CalPanelArea


!   !   CalPanelArea is used to calculate panel area on the immersed body 
!   !   surface and on the inner water plane


      SUBROUTINE CalPanelArea(XYZ,NTND,NELEM,NCN,NCON,DS)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NTND,NELEM
      INTEGER,INTENT(IN):: NCN(NELEM),NCON(NELEM,4)
      REAL*8,INTENT(IN):: XYZ(NTND,3)
      REAL*8,INTENT(OUT):: DS(NELEM)
      
      INTEGER IEL
      REAL*8 ADS
!
! -------------------------------------------------------------------------
! 
	 DO IEL=1,NELEM
         
        CALL CalDeltaArea_Improved(XYZ(NCON(IEL,1),:),XYZ(NCON(IEL,2),:),XYZ(NCON(IEL,3),:),DS(IEL))
          
        IF (NCN(IEL).EQ.4) THEN
              
         CALL CalDeltaArea_Improved(XYZ(NCON(IEL,1),:),XYZ(NCON(IEL,4),:),XYZ(NCON(IEL,3),:),ADS)        
         DS(IEL)=DS(IEL)+ADS
         
        ENDIF
     ENDDO
!
     RETURN
     END SUBROUTINE CalPanelArea

    
      
!=======================================================================
!SUBROUTINE CalPanelArea_TransNormal


!   !   CalPanelArea is used to calculate panel_area + normal_components 1,2,3 
!   !   on the immersed body surface and on the inner water plane


      SUBROUTINE CalPanelArea_TransNormal(XYZ,NTND,NELEM,NCN,NCON,DS,DXYZ_P)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NTND,NELEM
      INTEGER,INTENT(IN):: NCN(NELEM),NCON(NELEM,4)
      REAL*8,INTENT(IN):: XYZ(NTND,3)
      REAL*8,INTENT(OUT):: DS(NELEM),DXYZ_P(NELEM,6)
      
      INTEGER IEL
      REAL*8::  V21(3),V23(3),V13(3),V24(3),LTE(4),UN(3),ADS,AUN
!
! -------------------------------------------------------------------------
!     
	  DO IEL=1,NELEM

        CALL CalDeltaArea_Improved(XYZ(NCON(IEL,1),:),XYZ(NCON(IEL,2),:),XYZ(NCON(IEL,3),:),DS(IEL))
        
        IF (NCN(IEL).EQ.3) THEN
            
         V21(:)=XYZ(NCON(IEL,1),:)-XYZ(NCON(IEL,2),:)
         V23(:)=XYZ(NCON(IEL,3),:)-XYZ(NCON(IEL,2),:)

         UN(1)=V21(2)*V23(3)-V21(3)*V23(2)
         UN(2)=V21(3)*V23(1)-V21(1)*V23(3)
         UN(3)=V21(1)*V23(2)-V21(2)*V23(1)

         AUN=SQRT(UN(1)**2+UN(2)**2+UN(3)**2)

         DXYZ_P(IEL,1)=UN(1)/AUN
         DXYZ_P(IEL,2)=UN(2)/AUN
         DXYZ_P(IEL,3)=UN(3)/AUN
                   
        ELSEIF (NCN(IEL).EQ.4) THEN
            
         V13(:)=XYZ(NCON(IEL,3),:)-XYZ(NCON(IEL,1),:)
         V24(:)=XYZ(NCON(IEL,4),:)-XYZ(NCON(IEL,2),:)
         
         UN(1)=V24(2)*V13(3)-V24(3)*V13(2)
         UN(2)=V24(3)*V13(1)-V24(1)*V13(3)
         UN(3)=V24(1)*V13(2)-V24(2)*V13(1)

         AUN=SQRT(UN(1)**2+UN(2)**2+UN(3)**2)
         
         DXYZ_P(IEL,1)=UN(1)/AUN
         DXYZ_P(IEL,2)=UN(2)/AUN
         DXYZ_P(IEL,3)=UN(3)/AUN
         
         CALL CalDeltaArea_Improved(XYZ(NCON(IEL,1),:),XYZ(NCON(IEL,4),:),XYZ(NCON(IEL,3),:),ADS)        
         DS(IEL)=DS(IEL)+ADS
        
        ENDIF
        
     ENDDO
!
     RETURN
     END SUBROUTINE CalPanelArea_TransNormal

     
!=======================================================================
!SUBROUTINE CalPanelArea_Normal


!   !   CalPanelArea_Normal is used to calculate panel_area + 6 normals
!   !   on the immersed body surface and on the inner water plane


      SUBROUTINE CalPanelArea_Normal(RC,XYZ,NTND,NELEM,NCN,NCON,DS,XYZ_P,DXYZ_P)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NTND,NELEM
      INTEGER,INTENT(IN):: NCN(NELEM),NCON(NELEM,4)
      REAL*8,INTENT(IN):: RC(3),XYZ(NTND,3),XYZ_P(NELEM,3)
      REAL*8,INTENT(OUT):: DS(NELEM),DXYZ_P(NELEM,6)
      
      INTEGER IEL
      REAL*8::  V21(3),V23(3),V13(3),V24(3),LTE(4),UN(3),ADS,AUN
!
! -------------------------------------------------------------------------
!     
	  DO IEL=1,NELEM

        CALL CalDeltaArea_Improved(XYZ(NCON(IEL,1),:),XYZ(NCON(IEL,2),:),XYZ(NCON(IEL,3),:),DS(IEL))
        
        IF (NCN(IEL).EQ.3) THEN
            
         V21(:)=XYZ(NCON(IEL,1),:)-XYZ(NCON(IEL,2),:)
         V23(:)=XYZ(NCON(IEL,3),:)-XYZ(NCON(IEL,2),:)

         UN(1)=V21(2)*V23(3)-V21(3)*V23(2)
         UN(2)=V21(3)*V23(1)-V21(1)*V23(3)
         UN(3)=V21(1)*V23(2)-V21(2)*V23(1)

         AUN=SQRT(UN(1)**2+UN(2)**2+UN(3)**2)

         DXYZ_P(IEL,1)=UN(1)/AUN
         DXYZ_P(IEL,2)=UN(2)/AUN
         DXYZ_P(IEL,3)=UN(3)/AUN
         DXYZ_P(IEL,4)=(XYZ_P(IEL,2)-RC(2))*DXYZ_P(IEL,3)-(XYZ_P(IEL,3)-RC(3))*DXYZ_P(IEL,2)
         DXYZ_P(IEL,5)=(XYZ_P(IEL,3)-RC(3))*DXYZ_P(IEL,1)-(XYZ_P(IEL,1)-RC(1))*DXYZ_P(IEL,3)
         DXYZ_P(IEL,6)=(XYZ_P(IEL,1)-RC(1))*DXYZ_P(IEL,2)-(XYZ_P(IEL,2)-RC(2))*DXYZ_P(IEL,1)
          
        ELSEIF (NCN(IEL).EQ.4) THEN
            
         V13(:)=XYZ(NCON(IEL,3),:)-XYZ(NCON(IEL,1),:)
         V24(:)=XYZ(NCON(IEL,4),:)-XYZ(NCON(IEL,2),:)
         
         UN(1)=V24(2)*V13(3)-V24(3)*V13(2)
         UN(2)=V24(3)*V13(1)-V24(1)*V13(3)
         UN(3)=V24(1)*V13(2)-V24(2)*V13(1)

         AUN=SQRT(UN(1)**2+UN(2)**2+UN(3)**2)
         
         DXYZ_P(IEL,1)=UN(1)/AUN
         DXYZ_P(IEL,2)=UN(2)/AUN
         DXYZ_P(IEL,3)=UN(3)/AUN
         DXYZ_P(IEL,4)=(XYZ_P(IEL,2)-RC(2))*DXYZ_P(IEL,3)-(XYZ_P(IEL,3)-RC(3))*DXYZ_P(IEL,2)
         DXYZ_P(IEL,5)=(XYZ_P(IEL,3)-RC(3))*DXYZ_P(IEL,1)-(XYZ_P(IEL,1)-RC(1))*DXYZ_P(IEL,3)
         DXYZ_P(IEL,6)=(XYZ_P(IEL,1)-RC(1))*DXYZ_P(IEL,2)-(XYZ_P(IEL,2)-RC(2))*DXYZ_P(IEL,1)
         
         CALL CalDeltaArea_Improved(XYZ(NCON(IEL,1),:),XYZ(NCON(IEL,4),:),XYZ(NCON(IEL,3),:),ADS)        
         DS(IEL)=DS(IEL)+ADS
        
        ENDIF
        
     ENDDO    
!
     RETURN
     END SUBROUTINE CalPanelArea_Normal

    
      
!=======================================================================
!SUBROUTINE CalDeltaArea


!   !   CalPanelArea is used to calculate panel_area + normal_components 1,2,3 
!   !   on the immersed body surface and on the inner water plane


      SUBROUTINE CalDeltaArea(P1,P2,P3,DS)
      IMPLICIT   NONE

      REAL*8,INTENT(IN)::  P1(3),P2(3),P3(3)
      REAL*8,INTENT(OUT)::  DS
      REAL*8::  V21(3),V13(3),V23(3),LTE(4)
   
       V21(:)=P1-P2
       V23(:)=P3-P2
       V13(:)=P3-P1

       LTE(1)=SQRT(V21(1)**2+V21(2)**2+V21(3)**2)
       LTE(2)=SQRT(V13(1)**2+V13(2)**2+V13(3)**2)
       LTE(3)=SQRT(V23(1)**2+V23(2)**2+V23(3)**2)
       LTE(4)=0.5D0*(LTE(1)+LTE(2)+LTE(3))

       DS=SQRT(LTE(4)*(LTE(4)-LTE(1))*(LTE(4)-LTE(2))*(LTE(4)-LTE(3)))
!
     RETURN
     END SUBROUTINE CalDeltaArea
    
    
!
! ========================================================
!     
      SUBROUTINE CalDeltaArea_Improved(P1,P2,P3,DS)
      IMPLICIT   NONE

      REAL*8,INTENT(IN)::  P1(3),P2(3),P3(3)
      REAL*8,INTENT(OUT)::  DS
      REAL*8::  A,B,C
      REAL*8::  V21(3),V13(3),V23(3),LTE(4)
      
       V21(:)=P1-P2
       V23(:)=P3-P2
       V13(:)=P3-P1

       A=SQRT(V21(1)**2+V21(2)**2+V21(3)**2)
       B=SQRT(V13(1)**2+V13(2)**2+V13(3)**2)
       C=SQRT(V23(1)**2+V23(2)**2+V23(3)**2)
       
       LTE(1)=A+(B+C)
       LTE(2)=C-(A-B)
       LTE(3)=C+(A-B)
       LTE(4)=A+(B-C)
       
       DS=SQRT(LTE(1)*LTE(2)*LTE(3)*LTE(4))/4.D0
!
      RETURN
      END SUBROUTINE CalDeltaArea_Improved
!
! ========================================================
!     
      SUBROUTINE CalPanelChartLength(XYZ,NTND,NELEM,NCN,NCON,PNSZ)
      IMPLICIT   NONE


      INTEGER,INTENT(IN):: NTND,NELEM
      INTEGER,INTENT(IN):: NCN(NELEM),NCON(NELEM,4)
      REAL*8,INTENT(IN):: XYZ(NTND,3)
      REAL*8,INTENT(OUT):: PNSZ(NELEM)

      INTEGER IEL
      REAL*8::  A,B,C,D
      REAL*8::  V12(3),V23(3),V31(3),V34(3),V41(3)

	  DO IEL=1,NELEM
         
      IF (NCN(IEL).EQ.3) THEN
       V12(:)=XYZ(NCON(IEL,2),:)-XYZ(NCON(IEL,1),:)
       V23(:)=XYZ(NCON(IEL,3),:)-XYZ(NCON(IEL,2),:)
       V31(:)=XYZ(NCON(IEL,1),:)-XYZ(NCON(IEL,3),:)
       A=SQRT(V12(1)**2+V12(2)**2+V12(3)**2)
       B=SQRT(V23(1)**2+V23(2)**2+V23(3)**2)
       C=SQRT(V31(1)**2+V31(2)**2+V31(3)**2)
       PNSZ(IEL)=MAX(A,B,C)
      ELSEIF (NCN(IEL).EQ.4) THEN
       V12(:)=XYZ(NCON(IEL,2),:)-XYZ(NCON(IEL,1),:)
       V23(:)=XYZ(NCON(IEL,3),:)-XYZ(NCON(IEL,2),:)
       V34(:)=XYZ(NCON(IEL,4),:)-XYZ(NCON(IEL,3),:)
       V41(:)=XYZ(NCON(IEL,1),:)-XYZ(NCON(IEL,4),:)
       A=SQRT(V12(1)**2+V12(2)**2+V12(3)**2)
       B=SQRT(V23(1)**2+V23(2)**2+V23(3)**2)
       C=SQRT(V34(1)**2+V34(2)**2+V34(3)**2)
       D=SQRT(V41(1)**2+V41(2)**2+V41(3)**2)
       PNSZ(IEL)=MAX(A,B,C,D)
      ELSE
       PRINT*,'THE NUMBER OF PANEL VERTICES IS WRONG AT PANEL NO.', IEL
      ENDIF

      ENDDO
!
      RETURN
      END SUBROUTINE CalPanelChartLength

!-------------------------------------------------------------------------------
    END MODULE NormalProcess
!******************************************************************************* 