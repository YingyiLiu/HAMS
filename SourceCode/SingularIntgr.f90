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
MODULE SingularIntgr

   USE Const_mod
   USE PanelMesh_mod
   USE Inerfs_mod

   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SGLINTBD_TRI
   PUBLIC :: SGLINTBD_QUAD
   PUBLIC :: SGLINTBD_TRI2
   PUBLIC :: SGLINTBD_QUAD2
   
CONTAINS
!-------------------------------------------------------------------------------
!          Evaluation of 1/r singular integrals in Constant Panel Method  
!      using the analytical formula derived by Newman's method (Newman,1986)
!------------------------------------------------------------------------------- 

       SUBROUTINE SGLINTBD_TRI(IS,IEL,JEL,SIJ,DIJ,IRR)     
       IMPLICIT   NONE
      
       INTEGER,INTENT(IN):: IS,IEL,JEL,IRR
       INTEGER I,J,ITNO
       
       REAL*8,INTENT(OUT):: SIJ,DIJ(3)
       REAL*8:: X,Y,Z,ALFA,DKSI,DETA,R1,R2,L,XP,YP,ZP,A1,A2,B1,B2
       REAL*8:: XV(3,3),XQ(4,3),UNITI(3),UNITJ(3),UNITK(3),KSI(4),ETA(4)
       REAL*8:: P12(3),P13(3),PO3(3),PO1(3),PO2(3),PO(3),POF(3),XF(3)
       REAL*8:: LGRN,ST,SXT,SYT,SZT,SDX,SDY,SDZ
      
       IF (IRR.EQ.1) THEN
           
        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XF(1)=SY(IS,1)*XYZ_P(IEL,1)
         XF(2)=SX(IS,1)*XYZ_P(IEL,2)
         XF(3)=          XYZ_P(IEL,3) 
        ELSE
         XF(1)=SX(IS,1)*XYZ_P(IEL,1)
         XF(2)=SY(IS,1)*XYZ_P(IEL,2)
         XF(3)=         XYZ_P(IEL,3) 
        ENDIF

        XQ(1,:)=XYZ(NCON(JEL,1),:)
        XQ(2,:)=XYZ(NCON(JEL,2),:)
        XQ(3,:)=XYZ(NCON(JEL,3),:)
        IF (NCN(JEL).EQ.4) XQ(4,:)=XYZ(NCON(JEL,4),:)
        
       ELSEIF (IRR.EQ.2) THEN
           
        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XF(1)=SY(IS,1)*XYZ_P(IEL,1)
         XF(2)=SX(IS,1)*XYZ_P(IEL,2)
         XF(3)=          XYZ_P(IEL,3) 
        ELSE
         XF(1)=SX(IS,1)*XYZ_P(IEL,1)
         XF(2)=SY(IS,1)*XYZ_P(IEL,2)
         XF(3)=         XYZ_P(IEL,3) 
        ENDIF
        
        XQ(1,:)=iXYZ(iNCON(JEL,1),:)
        XQ(2,:)=iXYZ(iNCON(JEL,2),:)
        XQ(3,:)=iXYZ(iNCON(JEL,3),:)
        IF (iNCN(JEL).EQ.4) XQ(4,:)=iXYZ(iNCON(JEL,4),:)   
        
       ELSEIF (IRR.EQ.3) THEN
           
        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XF(1)=SY(IS,1)*iXYZ_P(IEL,1)
         XF(2)=SX(IS,1)*iXYZ_P(IEL,2)
         XF(3)=          iXYZ_P(IEL,3) 
        ELSE
         XF(1)=SX(IS,1)*iXYZ_P(IEL,1)
         XF(2)=SY(IS,1)*iXYZ_P(IEL,2)
         XF(3)=         iXYZ_P(IEL,3) 
        ENDIF

        XQ(1,:)=XYZ(NCON(JEL,1),:)
        XQ(2,:)=XYZ(NCON(JEL,2),:)
        XQ(3,:)=XYZ(NCON(JEL,3),:)
        IF (NCN(JEL).EQ.4) XQ(4,:)=XYZ(NCON(JEL,4),:)
        
       ELSEIF (IRR.EQ.4) THEN
           
        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XF(1)=SY(IS,1)*iXYZ_P(IEL,1)
         XF(2)=SX(IS,1)*iXYZ_P(IEL,2)
         XF(3)=          iXYZ_P(IEL,3) 
        ELSE
         XF(1)=SX(IS,1)*iXYZ_P(IEL,1)
         XF(2)=SY(IS,1)*iXYZ_P(IEL,2)
         XF(3)=         iXYZ_P(IEL,3) 
        ENDIF
        
        XQ(1,:)=iXYZ(iNCON(JEL,1),:)
        XQ(2,:)=iXYZ(iNCON(JEL,2),:)
        XQ(3,:)=iXYZ(iNCON(JEL,3),:)
        IF (iNCN(JEL).EQ.4) XQ(4,:)=iXYZ(iNCON(JEL,4),:)    
        
       ELSE
        PRINT*, 'Please define where the field point locates.'
       ENDIF   
       
       SIJ=0.D0
       DIJ=0.D0
       SDX=0.D0
       SDY=0.D0
       SDZ=0.D0   

!      D=1,B=3,A=2  D->B->A in clockwise direction
!      This subroutine is used for the mesh with its vertices of EVERY element 
!      numbered in anti-clockwise direction

       XV(1,:)=XQ(1,:)  ! point D
       XV(2,:)=XQ(2,:)  ! point A 
       XV(3,:)=XQ(3,:)  ! point B
       
       ITNO=1
       
30     CONTINUE
       
        P12(:)=XV(2,:)-XV(1,:) !PDA
        P13(:)=XV(3,:)-XV(1,:) !PDB
        ALFA=(P12(1)*P13(1)+P12(2)*P13(2)+P12(3)*P13(3))/(P12(1)**2+P12(2)**2+P12(3)**2)

        PO(:)=XV(1,:)+ALFA*P12(:)

        PO3(:)=XV(3,:)-PO(:)   !POB

        UNITI(:)=P12(:)/SQRT(P12(1)**2+P12(2)**2+P12(3)**2)
        UNITJ(:)=PO3(:)/SQRT(PO3(1)**2+PO3(2)**2+PO3(3)**2)
        UNITK(1)=UNITI(2)*UNITJ(3)-UNITI(3)*UNITJ(2)
        UNITK(2)=UNITI(3)*UNITJ(1)-UNITI(1)*UNITJ(3)
        UNITK(3)=UNITI(1)*UNITJ(2)-UNITI(2)*UNITJ(1)

        POF(:)=XF(:)-PO(:)   ! Field point in local coordinates
        X=UNITI(1)*POF(1)+UNITI(2)*POF(2)+UNITI(3)*POF(3)
        Y=UNITJ(1)*POF(1)+UNITJ(2)*POF(2)+UNITJ(3)*POF(3)
        Z=UNITK(1)*POF(1)+UNITK(2)*POF(2)+UNITK(3)*POF(3)

        PO1(:)=XV(1,:)-PO(:)  !POD
        PO2(:)=XV(2,:)-PO(:)  !POA
        KSI(1)=-SIGN(1.0D0,ALFA)*SQRT(PO1(1)**2+PO1(2)**2+PO1(3)**2)
        ETA(1)=0.D0
        KSI(2)=0.D0
        ETA(2)=SQRT(PO3(1)**2+PO3(2)**2+PO3(3)**2)
        KSI(3)=SQRT(PO2(1)**2+PO2(2)**2+PO2(3)**2)
        ETA(3)=0.D0
        KSI(4)=KSI(1)
        ETA(4)=ETA(1)    

        ST=0.D0
        SXT=0.D0
        SYT=0.D0
        SZT=0.D0
        
        DO 100 I=1,3
          DKSI=KSI(I+1)-KSI(I)
          DETA=ETA(I+1)-ETA(I)
          R1=SQRT((KSI(I)-X)**2+(ETA(I)-Y)**2+Z**2) 
          R2=SQRT((KSI(I+1)-X)**2+(ETA(I+1)-Y)**2+Z**2)
          L=SQRT((KSI(I+1)-KSI(I))**2+(ETA(I+1)-ETA(I))**2)
		  
          B1=DETA*((KSI(I)-X)**2+Z**2)-DKSI*(KSI(I)-X)*(ETA(I)-Y)
          A1=R1*Z*DKSI
          B2=DETA*((KSI(I+1)-X)**2+Z**2)-DKSI*(KSI(I+1)-X)*(ETA(I+1)-Y)
          A2=R2*Z*DKSI

	    IF (ABS(Z).LT.1.D-6) GOTO 50
          SZT=SZT+ATAN2((B1*A2-B2*A1),(A1*A2+B1*B2))

50      IF (ABS(R1+R2-L).LT.1.D-6) GOTO 100
          LGRN=LOG((R1+R2+L)/(R1+R2-L))             
          SXT=SXT-DETA/L*LGRN 
          SYT=SYT+DKSI/L*LGRN 
          ST=ST+(DETA*(X-KSI(I))-DKSI*(Y-ETA(I)))/L*LGRN
             
100    CONTINUE

       ST=ST-Z*SZT          
       SIJ=SIJ+ST

       DIJ(1)=DIJ(1)+UNITI(1)*SXT+UNITJ(1)*SYT+UNITK(1)*SZT 
       DIJ(2)=DIJ(2)+UNITI(2)*SXT+UNITJ(2)*SYT+UNITK(2)*SZT 
       DIJ(3)=DIJ(3)+UNITI(3)*SXT+UNITJ(3)*SYT+UNITK(3)*SZT 
      
       IF (IRR.EQ.1.OR.IRR.EQ.3) THEN
           
       IF (NCN(JEL).EQ.3) THEN
       GOTO 80     
       ELSEIF (NCN(JEL).EQ.4.AND.ITNO.EQ.2) THEN
       GOTO 80            
       ELSEIF (NCN(JEL).EQ.4.AND.ITNO.EQ.1) THEN      
       ITNO=2
       XV(1,:)=XQ(1,:)  ! point D
       XV(2,:)=XQ(3,:)  ! point A
       XV(3,:)=XQ(4,:)  ! point B
       GOTO 30            
       ENDIF  
       
       ELSE
           
       IF (iNCN(JEL).EQ.3) THEN
       GOTO 80     
       ELSEIF (iNCN(JEL).EQ.4.AND.ITNO.EQ.2) THEN
       GOTO 80            
       ELSEIF (iNCN(JEL).EQ.4.AND.ITNO.EQ.1) THEN      
       ITNO=2
       XV(1,:)=XQ(1,:)  ! point D
       XV(2,:)=XQ(3,:)  ! point A
       XV(3,:)=XQ(4,:)  ! point B
       GOTO 30            
       ENDIF 
       
       ENDIF        
       
80    CONTINUE

      RETURN
      END SUBROUTINE SGLINTBD_TRI
       
       
!-------------------------------------------------------------------------------
!          Evaluation of 1/r singular integrals in Constant Panel Method  
!      using the analytical formula derived by Newman's method (Newman,1986)
!------------------------------------------------------------------------------- 

       SUBROUTINE SGLINTBD_QUAD(IS,IEL,JEL,SIJ,DIJ,IRR)     
       IMPLICIT   NONE
      
       INTEGER,INTENT(IN):: IS,IEL,JEL,IRR
       INTEGER I,J,ITNO
       
       REAL*8,INTENT(OUT):: SIJ,DIJ(3)
       REAL*8:: X,Y,Z,ALFA,DKSI,DETA,R1,R2,L,XP,YP,ZP,A1,A2,B1,B2
       REAL*8:: XV(4,3),XQ(4,3),UNITI(3),UNITJ(3),UNITK(3),KSI(5),ETA(5),DTA(5)
       REAL*8:: POV(4,3),PON(4,3),UN(3),PO(3),POF(3),XF(3),PL(3,3)
       REAL*8:: LGRN,ST,SXT,SYT,SZT,SDX,SDY,SDZ,RDET

       IF (IRR.EQ.1) THEN
           
        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XF(1)=SY(IS,1)*XYZ_P(IEL,1)
         XF(2)=SX(IS,1)*XYZ_P(IEL,2)
         XF(3)=          XYZ_P(IEL,3) 
        ELSE
         XF(1)=SX(IS,1)*XYZ_P(IEL,1)
         XF(2)=SY(IS,1)*XYZ_P(IEL,2)
         XF(3)=         XYZ_P(IEL,3) 
        ENDIF
        
        XQ(1,:)=XYZ(NCON(JEL,1),:)
        XQ(2,:)=XYZ(NCON(JEL,2),:)
        XQ(3,:)=XYZ(NCON(JEL,3),:)
        IF (NCN(JEL).EQ.4) XQ(4,:)=XYZ(NCON(JEL,4),:)
        
       ELSEIF (IRR.EQ.2) THEN
           
        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XF(1)=SY(IS,1)*XYZ_P(IEL,1)
         XF(2)=SX(IS,1)*XYZ_P(IEL,2)
         XF(3)=          XYZ_P(IEL,3) 
        ELSE
         XF(1)=SX(IS,1)*XYZ_P(IEL,1)
         XF(2)=SY(IS,1)*XYZ_P(IEL,2)
         XF(3)=         XYZ_P(IEL,3) 
        ENDIF
        
        XQ(1,:)=iXYZ(iNCON(JEL,1),:)
        XQ(2,:)=iXYZ(iNCON(JEL,2),:)
        XQ(3,:)=iXYZ(iNCON(JEL,3),:)
        IF (iNCN(JEL).EQ.4) XQ(4,:)=iXYZ(iNCON(JEL,4),:)   
        
       ELSEIF (IRR.EQ.3) THEN
           
        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XF(1)=SY(IS,1)*iXYZ_P(IEL,1)
         XF(2)=SX(IS,1)*iXYZ_P(IEL,2)
         XF(3)=          iXYZ_P(IEL,3) 
        ELSE
         XF(1)=SX(IS,1)*iXYZ_P(IEL,1)
         XF(2)=SY(IS,1)*iXYZ_P(IEL,2)
         XF(3)=         iXYZ_P(IEL,3) 
        ENDIF
        
        XQ(1,:)=XYZ(NCON(JEL,1),:)
        XQ(2,:)=XYZ(NCON(JEL,2),:)
        XQ(3,:)=XYZ(NCON(JEL,3),:)
        IF (NCN(JEL).EQ.4) XQ(4,:)=XYZ(NCON(JEL,4),:)
        
       ELSEIF (IRR.EQ.4) THEN
           
        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XF(1)=SY(IS,1)*iXYZ_P(IEL,1)
         XF(2)=SX(IS,1)*iXYZ_P(IEL,2)
         XF(3)=          iXYZ_P(IEL,3) 
        ELSE
         XF(1)=SX(IS,1)*iXYZ_P(IEL,1)
         XF(2)=SY(IS,1)*iXYZ_P(IEL,2)
         XF(3)=         iXYZ_P(IEL,3) 
        ENDIF
        
        XQ(1,:)=iXYZ(iNCON(JEL,1),:)
        XQ(2,:)=iXYZ(iNCON(JEL,2),:)
        XQ(3,:)=iXYZ(iNCON(JEL,3),:)
        IF (iNCN(JEL).EQ.4) XQ(4,:)=iXYZ(iNCON(JEL,4),:)    
        
       ELSE
        PRINT*, 'Please define where the field point locates.'
       ENDIF   
       
       SIJ=0.D0
       DIJ=0.D0
       SDX=0.D0
       SDY=0.D0
       SDZ=0.D0   

!      D=1,B=3,A=2  D->B->A in clockwise direction
!      This subroutine is used for the mesh with its vertices of EVERY element 
!      numbered in anti-clockwise direction

       XV(1,:)=XQ(1,:)  ! point D
       XV(2,:)=XQ(2,:)  ! point A
       XV(3,:)=XQ(3,:)  ! point B
       XV(4,:)=XQ(4,:)  ! point B
       
!      Judge if the quadrilateral element is a flat panel
       PL(1,:)=XV(2,:)-XV(1,:)
       PL(2,:)=XV(3,:)-XV(1,:)
       PL(3,:)=XV(4,:)-XV(1,:)

!      Calculate singularity

        PO(:)=(XV(1,:)+XV(2,:)+XV(3,:)+XV(4,:))/4.D0

        POV(1,:)=XV(1,:)-PO(:)   !PO1
        POV(2,:)=XV(2,:)-PO(:)   !PO2
        
        UN(1)=POV(1,2)*POV(2,3)-POV(1,3)*POV(2,2)    !  PO1 X PO2
        UN(2)=POV(1,3)*POV(2,1)-POV(1,1)*POV(2,3)
        UN(3)=POV(1,1)*POV(2,2)-POV(1,2)*POV(2,1)
        
        UNITI(:)=POV(1,:)/SQRT(POV(1,1)**2+POV(1,2)**2+POV(1,3)**2)
        UNITK(:)=UN(:)/SQRT(UN(1)**2+UN(2)**2+UN(3)**2)
        UNITJ(1)=UNITK(2)*UNITI(3)-UNITK(3)*UNITI(2)
        UNITJ(2)=UNITK(3)*UNITI(1)-UNITK(1)*UNITI(3)
        UNITJ(3)=UNITK(1)*UNITI(2)-UNITK(2)*UNITI(1)

!        XF(:)=XV(3,:)
        POF(:)=XF(:)-PO(:)   ! Field point in local coordinates
        X=UNITI(1)*POF(1)+UNITI(2)*POF(2)+UNITI(3)*POF(3)
        Y=UNITJ(1)*POF(1)+UNITJ(2)*POF(2)+UNITJ(3)*POF(3)
        Z=UNITK(1)*POF(1)+UNITK(2)*POF(2)+UNITK(3)*POF(3)
        
        POV(3,:)=XV(3,:)-PO(:)  !POD
        POV(4,:)=XV(4,:)-PO(:)  !POA
        
        DO I=1,4
          DO J=1,3
             PON(I,J)=POV(I,J)-(UNITK(1)*POV(I,1)+UNITK(2)*POV(I,2)+UNITK(3)*POV(I,3))*UNITK(J)
          ENDDO
        ENDDO                    
        
        KSI(1)=UNITI(1)*PON(1,1)+UNITI(2)*PON(1,2)+UNITI(3)*PON(1,3)
        ETA(1)=UNITJ(1)*PON(1,1)+UNITJ(2)*PON(1,2)+UNITJ(3)*PON(1,3)
        
        KSI(2)=UNITI(1)*PON(2,1)+UNITI(2)*PON(2,2)+UNITI(3)*PON(2,3)
        ETA(2)=UNITJ(1)*PON(2,1)+UNITJ(2)*PON(2,2)+UNITJ(3)*PON(2,3)
        
        KSI(3)=UNITI(1)*PON(3,1)+UNITI(2)*PON(3,2)+UNITI(3)*PON(3,3)
        ETA(3)=UNITJ(1)*PON(3,1)+UNITJ(2)*PON(3,2)+UNITJ(3)*PON(3,3)
        
        KSI(4)=UNITI(1)*PON(4,1)+UNITI(2)*PON(4,2)+UNITI(3)*PON(4,3)
        ETA(4)=UNITJ(1)*PON(4,1)+UNITJ(2)*PON(4,2)+UNITJ(3)*PON(4,3)
        
        KSI(5)=KSI(1)
        ETA(5)=ETA(1) 

        ST=0.D0
        SXT=0.D0
        SYT=0.D0
        SZT=0.D0
        
      DO 100 I=1,4
          
          DKSI=KSI(I+1)-KSI(I)
          DETA=ETA(I+1)-ETA(I)
          R1=SQRT((KSI(I)-X)**2+(ETA(I)-Y)**2+Z**2) 
          R2=SQRT((KSI(I+1)-X)**2+(ETA(I+1)-Y)**2+Z**2)
          L=SQRT((KSI(I+1)-KSI(I))**2+(ETA(I+1)-ETA(I))**2)
		  
          IF (DABS(L).LT.1.D-8) GOTO 100     ! to avoid overlapping vertexes
          
          B1=DETA*((KSI(I)-X)**2+Z**2)-DKSI*(KSI(I)-X)*(ETA(I)-Y)
          A1=R1*Z*DKSI
          B2=DETA*((KSI(I+1)-X)**2+Z**2)-DKSI*(KSI(I+1)-X)*(ETA(I+1)-Y)
          A2=R2*Z*DKSI

	    IF (ABS(Z).LT.1.D-6) GOTO 50
          SZT=SZT+ATAN2((B1*A2-B2*A1),(A1*A2+B1*B2))

50      IF (ABS(R1+R2-L).LT.1.D-6) GOTO 100
          LGRN=LOG((R1+R2+L)/(R1+R2-L))         
          SXT=SXT-DETA/L*LGRN 
          SYT=SYT+DKSI/L*LGRN 
          ST=ST-(DETA*(X-KSI(I))-DKSI*(Y-ETA(I)))/L*LGRN
             
100    CONTINUE

       ST=ST+Z*SZT
       SIJ=SIJ+ST

       DIJ(1)=DIJ(1)+UNITI(1)*SXT+UNITJ(1)*SYT+UNITK(1)*SZT 
       DIJ(2)=DIJ(2)+UNITI(2)*SXT+UNITJ(2)*SYT+UNITK(2)*SZT 
       DIJ(3)=DIJ(3)+UNITI(3)*SXT+UNITJ(3)*SYT+UNITK(3)*SZT
       DIJ=-DIJ
       
80    CONTINUE

      RETURN
      END SUBROUTINE SGLINTBD_QUAD
       
!--------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!          Evaluation of 1/r singular integrals in Constant Panel Method  
!      using the analytical formula derived by Newman's method (Newman,1986)
!------------------------------------------------------------------------------- 

       SUBROUTINE SGLINTBD_TRI2(XF,JEL,RK,IRR)
       IMPLICIT   NONE

       INTEGER,INTENT(IN):: JEL,IRR
       REAL*8,INTENT(IN):: XF(3)
       INTEGER I,J,ITNO

       REAL*8,INTENT(OUT):: RK(4)
       REAL*8:: X,Y,Z,ALFA,DKSI,DETA,R1,R2,L,XP,YP,ZP,A1,A2,B1,B2
       REAL*8:: XV(3,3),XQ(4,3),UNITI(3),UNITJ(3),UNITK(3),KSI(4),ETA(4)
       REAL*8:: P12(3),P13(3),PO3(3),PO1(3),PO2(3),PO(3),POF(3)
       REAL*8:: LGRN,ST,SXT,SYT,SZT,SX,SY,SZ

       IF (IRR.EQ.1) THEN

        XQ(1,:)=XYZ(NCON(JEL,1),:)
        XQ(2,:)=XYZ(NCON(JEL,2),:)
        XQ(3,:)=XYZ(NCON(JEL,3),:)
        IF (NCN(JEL).EQ.4) XQ(4,:)=XYZ(NCON(JEL,4),:)

       ELSEIF (IRR.EQ.2) THEN

        XQ(1,:)=iXYZ(iNCON(JEL,1),:)
        XQ(2,:)=iXYZ(iNCON(JEL,2),:)
        XQ(3,:)=iXYZ(iNCON(JEL,3),:)
        IF (iNCN(JEL).EQ.4) XQ(4,:)=iXYZ(iNCON(JEL,4),:)

       ELSEIF (IRR.EQ.3) THEN

        XQ(1,:)=XYZ(NCON(JEL,1),:)
        XQ(2,:)=XYZ(NCON(JEL,2),:)
        XQ(3,:)=XYZ(NCON(JEL,3),:)
        IF (NCN(JEL).EQ.4) XQ(4,:)=XYZ(NCON(JEL,4),:)

       ELSEIF (IRR.EQ.4) THEN

        XQ(1,:)=iXYZ(iNCON(JEL,1),:)
        XQ(2,:)=iXYZ(iNCON(JEL,2),:)
        XQ(3,:)=iXYZ(iNCON(JEL,3),:)
        IF (iNCN(JEL).EQ.4) XQ(4,:)=iXYZ(iNCON(JEL,4),:)

       ELSE
        PRINT*, 'Please define where the field point locates.'
       ENDIF   

       RK=0.D0
       SX=0.D0
       SY=0.D0
       SZ=0.D0

!      D=1,B=3,A=2  D->B->A in clockwise direction
!      This subroutine is used for the mesh with its vertices of EVERY element 
!      numbered in anti-clockwise direction

       XV(1,:)=XQ(1,:)  ! point D
       XV(2,:)=XQ(2,:)  ! point A
       XV(3,:)=XQ(3,:)  ! point B

       ITNO=1

30     CONTINUE

        P12(:)=XV(2,:)-XV(1,:) !PDA
        P13(:)=XV(3,:)-XV(1,:) !PDB
        ALFA=(P12(1)*P13(1)+P12(2)*P13(2)+P12(3)*P13(3))/(P12(1)**2+P12(2)**2+P12(3)**2)

        PO(:)=XV(1,:)+ALFA*P12(:)

        PO3(:)=XV(3,:)-PO(:)   !POB

        UNITI(:)=P12(:)/SQRT(P12(1)**2+P12(2)**2+P12(3)**2)
        UNITJ(:)=PO3(:)/SQRT(PO3(1)**2+PO3(2)**2+PO3(3)**2)
        UNITK(1)=UNITI(2)*UNITJ(3)-UNITI(3)*UNITJ(2)
        UNITK(2)=UNITI(3)*UNITJ(1)-UNITI(1)*UNITJ(3)
        UNITK(3)=UNITI(1)*UNITJ(2)-UNITI(2)*UNITJ(1)

        POF(:)=XF(:)-PO(:)   ! Field point in local coordinates
        X=UNITI(1)*POF(1)+UNITI(2)*POF(2)+UNITI(3)*POF(3)
        Y=UNITJ(1)*POF(1)+UNITJ(2)*POF(2)+UNITJ(3)*POF(3)
        Z=UNITK(1)*POF(1)+UNITK(2)*POF(2)+UNITK(3)*POF(3)

        PO1(:)=XV(1,:)-PO(:)  !POD
        PO2(:)=XV(2,:)-PO(:)  !POA
        KSI(1)=-SIGN(1.0D0,ALFA)*SQRT(PO1(1)**2+PO1(2)**2+PO1(3)**2)
        ETA(1)=0.D0
        KSI(2)=0.D0
        ETA(2)=SQRT(PO3(1)**2+PO3(2)**2+PO3(3)**2)
        KSI(3)=SQRT(PO2(1)**2+PO2(2)**2+PO2(3)**2)
        ETA(3)=0.D0
        KSI(4)=KSI(1)
        ETA(4)=ETA(1)

        ST=0.D0
        SXT=0.D0
        SYT=0.D0
        SZT=0.D0

        DO 100 I=1,3
          DKSI=KSI(I+1)-KSI(I)
          DETA=ETA(I+1)-ETA(I)
          R1=SQRT((KSI(I)-X)**2+(ETA(I)-Y)**2+Z**2)
          R2=SQRT((KSI(I+1)-X)**2+(ETA(I+1)-Y)**2+Z**2)
          L=SQRT((KSI(I+1)-KSI(I))**2+(ETA(I+1)-ETA(I))**2)

          B1=DETA*((KSI(I)-X)**2+Z**2)-DKSI*(KSI(I)-X)*(ETA(I)-Y)
          A1=R1*Z*DKSI
          B2=DETA*((KSI(I+1)-X)**2+Z**2)-DKSI*(KSI(I+1)-X)*(ETA(I+1)-Y)
          A2=R2*Z*DKSI

	    IF (ABS(Z).LT.1.D-6) GOTO 50
          SZT=SZT+ATAN2((B1*A2-B2*A1),(A1*A2+B1*B2))

50      IF (ABS(R1+R2-L).LT.1.D-6) GOTO 100
          LGRN=LOG((R1+R2+L)/(R1+R2-L))
          SXT=SXT-DETA/L*LGRN
          SYT=SYT+DKSI/L*LGRN
          ST=ST+(DETA*(X-KSI(I))-DKSI*(Y-ETA(I)))/L*LGRN

100    CONTINUE

       ST=ST-Z*SZT
       RK(1)=RK(1)+ST

       RK(2)=RK(2)+UNITI(1)*SXT+UNITJ(1)*SYT+UNITK(1)*SZT
       RK(3)=RK(3)+UNITI(2)*SXT+UNITJ(2)*SYT+UNITK(2)*SZT
       RK(4)=RK(4)+UNITI(3)*SXT+UNITJ(3)*SYT+UNITK(3)*SZT

       IF (IRR.EQ.1.OR.IRR.EQ.3) THEN

       IF (NCN(JEL).EQ.3) THEN
       GOTO 80
       ELSEIF (NCN(JEL).EQ.4.AND.ITNO.EQ.2) THEN
       GOTO 80
       ELSEIF (NCN(JEL).EQ.4.AND.ITNO.EQ.1) THEN
       ITNO=2
       XV(1,:)=XQ(1,:)  ! point D
       XV(2,:)=XQ(3,:)  ! point A
       XV(3,:)=XQ(4,:)  ! point B
       GOTO 30
       ENDIF

       ELSE

       IF (iNCN(JEL).EQ.3) THEN
       GOTO 80
       ELSEIF (iNCN(JEL).EQ.4.AND.ITNO.EQ.2) THEN
       GOTO 80
       ELSEIF (iNCN(JEL).EQ.4.AND.ITNO.EQ.1) THEN
       ITNO=2
       XV(1,:)=XQ(1,:)  ! point D
       XV(2,:)=XQ(3,:)  ! point A
       XV(3,:)=XQ(4,:)  ! point B
       GOTO 30
       ENDIF

       ENDIF

80    CONTINUE

      RETURN
      END SUBROUTINE SGLINTBD_TRI2
       
       
!-------------------------------------------------------------------------------
!          Evaluation of 1/r singular integrals in Constant Panel Method  
!      using the analytical formula derived by Newman's method (Newman,1986)
!------------------------------------------------------------------------------- 

       SUBROUTINE SGLINTBD_QUAD2(XF,JEL,RK,IRR)
       IMPLICIT   NONE
      
       INTEGER,INTENT(IN):: JEL,IRR
       REAL*8,INTENT(IN):: XF(3)
       INTEGER I,J,ITNO
       
       REAL*8,INTENT(OUT):: RK(4)
       REAL*8:: X,Y,Z,ALFA,DKSI,DETA,R1,R2,L,XP,YP,ZP,A1,A2,B1,B2
       REAL*8:: XV(4,3),XQ(4,3),UNITI(3),UNITJ(3),UNITK(3),KSI(5),ETA(5),DTA(5)
       REAL*8:: POV(4,3),PON(4,3),UN(3),PO(3),POF(3),PL(3,3)
       REAL*8:: LGRN,ST,SXT,SYT,SZT,SX,SY,SZ,RDET

       IF (IRR.EQ.1) THEN

        XQ(1,:)=XYZ(NCON(JEL,1),:)
        XQ(2,:)=XYZ(NCON(JEL,2),:)
        XQ(3,:)=XYZ(NCON(JEL,3),:)
        IF (NCN(JEL).EQ.4) XQ(4,:)=XYZ(NCON(JEL,4),:)

       ELSEIF (IRR.EQ.2) THEN

        XQ(1,:)=iXYZ(iNCON(JEL,1),:)
        XQ(2,:)=iXYZ(iNCON(JEL,2),:)
        XQ(3,:)=iXYZ(iNCON(JEL,3),:)
        IF (iNCN(JEL).EQ.4) XQ(4,:)=iXYZ(iNCON(JEL,4),:)

       ELSEIF (IRR.EQ.3) THEN
        
        XQ(1,:)=XYZ(NCON(JEL,1),:)
        XQ(2,:)=XYZ(NCON(JEL,2),:)
        XQ(3,:)=XYZ(NCON(JEL,3),:)
        IF (NCN(JEL).EQ.4) XQ(4,:)=XYZ(NCON(JEL,4),:)

       ELSEIF (IRR.EQ.4) THEN

        XQ(1,:)=iXYZ(iNCON(JEL,1),:)
        XQ(2,:)=iXYZ(iNCON(JEL,2),:)
        XQ(3,:)=iXYZ(iNCON(JEL,3),:)
        IF (iNCN(JEL).EQ.4) XQ(4,:)=iXYZ(iNCON(JEL,4),:)

       ELSE
        PRINT*, 'Please define where the field point locates.'
       ENDIF   

       RK=0.D0
       SX=0.D0
       SY=0.D0
       SZ=0.D0

!      D=1,B=3,A=2  D->B->A in clockwise direction
!      This subroutine is used for the mesh with its vertices of EVERY element 
!      numbered in anti-clockwise direction

       XV(1,:)=XQ(1,:)  ! point D
       XV(2,:)=XQ(2,:)  ! point A
       XV(3,:)=XQ(3,:)  ! point B
       XV(4,:)=XQ(4,:)  ! point B

!      Judge if the quadrilateral element is a flat panel
       PL(1,:)=XV(2,:)-XV(1,:)
       PL(2,:)=XV(3,:)-XV(1,:)
       PL(3,:)=XV(4,:)-XV(1,:)

       CALL BSDET(PL,3,RDET)

!       IF (ABS(RDET).GT.1.E-6) THEN
!        WRITE(*,*) 'IEL,JEL:',IEL,JEL      
!        WRITE(*,*) 'Quad coordinates are:'
!        DO I=1,4
!        WRITE(*,*) (XV(I,J),J=1,3)
!        ENDDO
!        WRITE(*,*) 'Field point is:'
!        WRITE(*,*) (XF(J),J=1,3)
!        WRITE(*,*) 'RDET is not equal to zero:', RDET
!        STOP
!       ENDIF   

!      Calculate singularity

        PO(:)=(XV(1,:)+XV(2,:)+XV(3,:)+XV(4,:))/4.D0

        POV(1,:)=XV(1,:)-PO(:)   !PO1
        POV(2,:)=XV(2,:)-PO(:)   !PO2

        UN(1)=POV(1,2)*POV(2,3)-POV(1,3)*POV(2,2)    !  PO1 X PO2
        UN(2)=POV(1,3)*POV(2,1)-POV(1,1)*POV(2,3)
        UN(3)=POV(1,1)*POV(2,2)-POV(1,2)*POV(2,1)

        UNITI(:)=POV(1,:)/SQRT(POV(1,1)**2+POV(1,2)**2+POV(1,3)**2)
        UNITK(:)=UN(:)/SQRT(UN(1)**2+UN(2)**2+UN(3)**2)
        UNITJ(1)=UNITK(2)*UNITI(3)-UNITK(3)*UNITI(2)
        UNITJ(2)=UNITK(3)*UNITI(1)-UNITK(1)*UNITI(3)
        UNITJ(3)=UNITK(1)*UNITI(2)-UNITK(2)*UNITI(1)

!        XF(:)=XV(3,:)
        POF(:)=XF(:)-PO(:)   ! Field point in local coordinates
        X=UNITI(1)*POF(1)+UNITI(2)*POF(2)+UNITI(3)*POF(3)
        Y=UNITJ(1)*POF(1)+UNITJ(2)*POF(2)+UNITJ(3)*POF(3)
        Z=UNITK(1)*POF(1)+UNITK(2)*POF(2)+UNITK(3)*POF(3)

!        IF (IEL.EQ.1.AND.JEL.EQ.44) THEN
!           WRITE(10,*) 'XYZ',X,Y,Z
!           DO I=1,4
!           WRITE(10,*) 'I=',(XV(I,J),J=1,3)
!           WRITE(10,*) 'PO',(PO(J),J=1,3)
!           ENDDO
!           PAUSE
!          WRITE(*,*) 'RDET',RDET
!          PAUSE
!        ENDIF

        POV(3,:)=XV(3,:)-PO(:)  !POD
        POV(4,:)=XV(4,:)-PO(:)  !POA

        DO I=1,4
          DO J=1,3
             PON(I,J)=POV(I,J)-(UNITK(1)*POV(I,1)+UNITK(2)*POV(I,2)+UNITK(3)*POV(I,3))*UNITK(J)
          ENDDO
        ENDDO

        KSI(1)=UNITI(1)*PON(1,1)+UNITI(2)*PON(1,2)+UNITI(3)*PON(1,3)
        ETA(1)=UNITJ(1)*PON(1,1)+UNITJ(2)*PON(1,2)+UNITJ(3)*PON(1,3)

        KSI(2)=UNITI(1)*PON(2,1)+UNITI(2)*PON(2,2)+UNITI(3)*PON(2,3)
        ETA(2)=UNITJ(1)*PON(2,1)+UNITJ(2)*PON(2,2)+UNITJ(3)*PON(2,3)

        KSI(3)=UNITI(1)*PON(3,1)+UNITI(2)*PON(3,2)+UNITI(3)*PON(3,3)
        ETA(3)=UNITJ(1)*PON(3,1)+UNITJ(2)*PON(3,2)+UNITJ(3)*PON(3,3)

        KSI(4)=UNITI(1)*PON(4,1)+UNITI(2)*PON(4,2)+UNITI(3)*PON(4,3)
        ETA(4)=UNITJ(1)*PON(4,1)+UNITJ(2)*PON(4,2)+UNITJ(3)*PON(4,3)

        KSI(5)=KSI(1)
        ETA(5)=ETA(1)

!        IF (JEL.EQ.1) THEN
!          WRITE(10,*) 'DTA',(DTA(J),J=1,4)       
!          DO I=1,4
!          WRITE(10,*) 'I=',(XV(I,J),J=1,3)
!          WRITE(10,*) 'PO',(PO(J),J=1,3)
!          ENDDO
!          PAUSE
!          WRITE(*,*) 'RDET',RDET
!          PAUSE
!        ENDIF

        ST=0.D0
        SXT=0.D0
        SYT=0.D0
        SZT=0.D0

      DO 100 I=1,4

          DKSI=KSI(I+1)-KSI(I)
          DETA=ETA(I+1)-ETA(I)
          R1=SQRT((KSI(I)-X)**2+(ETA(I)-Y)**2+Z**2)
          R2=SQRT((KSI(I+1)-X)**2+(ETA(I+1)-Y)**2+Z**2)
          L=SQRT((KSI(I+1)-KSI(I))**2+(ETA(I+1)-ETA(I))**2)
		  
          IF (DABS(L).LT.1.D-8) GOTO 100     ! to avoid overlapping vertexes
          
          B1=DETA*((KSI(I)-X)**2+Z**2)-DKSI*(KSI(I)-X)*(ETA(I)-Y)
          A1=R1*Z*DKSI
          B2=DETA*((KSI(I+1)-X)**2+Z**2)-DKSI*(KSI(I+1)-X)*(ETA(I+1)-Y)
          A2=R2*Z*DKSI

	    IF (ABS(Z).LT.1.D-6) GOTO 50
          SZT=SZT+ATAN2((B1*A2-B2*A1),(A1*A2+B1*B2))

50      IF (ABS(R1+R2-L).LT.1.D-6) GOTO 100
          LGRN=LOG((R1+R2+L)/(R1+R2-L))
          SXT=SXT-DETA/L*LGRN
          SYT=SYT+DKSI/L*LGRN
          ST=ST-(DETA*(X-KSI(I))-DKSI*(Y-ETA(I)))/L*LGRN

100    CONTINUE

       ST=ST+Z*SZT
       RK(1)=RK(1)+ST

       RK(2)=-(RK(2)+UNITI(1)*SXT+UNITJ(1)*SYT+UNITK(1)*SZT)
       RK(3)=-(RK(3)+UNITI(2)*SXT+UNITJ(2)*SYT+UNITK(2)*SZT)
       RK(4)=-(RK(4)+UNITI(3)*SXT+UNITJ(3)*SYT+UNITK(3)*SZT)

80    CONTINUE

      RETURN
      END SUBROUTINE SGLINTBD_QUAD2

!-------------------------------------------------------------------------------
!         Calculate determinant of a matrix
!------------------------------------------------------------------------------- 

    SUBROUTINE BSDET(A,N,DET)
    IMPLICIT NONE

    INTEGER,INTENT(IN):: N
    REAL*8,INTENT(INOUT):: A(N,N)
    REAL*8,INTENT(OUT):: DET

    INTEGER K,I,J,IS,JS
    REAL*8 F,D,Q

	F=1.0
	DET=1.0
	DO 100 K=1,N-1
	  Q=0.0
	  DO 10 I=K,N
	  DO 10 J=K,N
	    IF (ABS(A(I,J)).GT.Q) THEN
	      Q=ABS(A(I,J))
	      IS=I
	      JS=J
	    END IF
10	  CONTINUE
	  IF (Q+1.0.EQ.1.0) THEN
	    DET=0.0
	    RETURN
	  END IF
	  IF (IS.NE.K) THEN
	    F=-F
	    DO 20 J=K,N
	      D=A(K,J)
	      A(K,J)=A(IS,J)
	      A(IS,J)=D
20	    CONTINUE
	  END IF
	  IF (JS.NE.K) THEN
	    F=-F
	    DO 30 I=K,N
	      D=A(I,JS)
	      A(I,JS)=A(I,K)
	      A(I,K)=D
30	    CONTINUE
	  END IF
	  DET=DET*A(K,K)
	  DO 50 I=K+1,N
	    D=A(I,K)/A(K,K)
	    DO 40 J=K+1,N
40	    A(I,J)=A(I,J)-D*A(K,J)
50	  CONTINUE
100	CONTINUE
	DET=F*DET*A(N,N)
	RETURN
    END SUBROUTINE BSDET

!-------------------------------------------------------------------------------
END MODULE SingularIntgr
!*******************************************************************************