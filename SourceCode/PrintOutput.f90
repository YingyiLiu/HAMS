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
MODULE PrintOutput

   USE Const_mod
   
   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: PrintHeading
   PUBLIC :: PrintBody_RealVal
   
   PUBLIC :: PrintBody_CmplxVal
   PUBLIC :: PrintBody_Exfc
   PUBLIC :: PrintEnd
   
CONTAINS
!   --------------------------------------------------------
!       Print the headings in HydroStar output files
!   --------------------------------------------------------

      SUBROUTINE PrintHeading(NFILE,NBETA,REFL,RAOType,MD1,MD2,H,XW,XC,WVHD)
      USE Const_mod
      IMPLICIT NONE
	  
      INTEGER  NFILE,NBETA,MD1,MD2,II
      REAL*8   REFL,H,XW(2),XC(3),WVHD(NBETA),NWVHD(NBETA)
      CHARACTER(*)::RAOType
      CHARACTER(LEN=100)::FMT
      
      WRITE(NFILE, '(A11)')  '# Project :'
      WRITE(NFILE, '(A11)')  '# User    :'
      IF (RAOType.EQ.'AddedMass') THEN
	   WRITE(NFILE, '(A8,A10,A1,I1,I1,A4)')  '# File :',adjustl(trim(RAOType)),'_',MD1,MD2,'.rao'
      ELSEIF (adjustl(trim(RAOType)).EQ.'WaveDamping') THEN
	   WRITE(NFILE, '(A8,A12,A1,I1,I1,A4)')  '# File :',adjustl(trim(RAOType)),'_',MD1,MD2,'.rao'
      ELSEIF (RAOType.EQ.'Cm') THEN
	   WRITE(NFILE, '(A8,A3,A1,I1,I1,A4)')  '# File :',adjustl(trim(RAOType)),'_',MD1,MD2,'.dat'
      ELSEIF (adjustl(trim(RAOType)).EQ.'Ca') THEN
	   WRITE(NFILE, '(A8,A3,A1,I1,I1,A4)')  '# File :',adjustl(trim(RAOType)),'_',MD1,MD2,'.dat'
      ELSEIF (adjustl(trim(RAOType)).EQ.'Excitation') THEN
	   WRITE(NFILE, '(A8,A11,A1,I1,A4)')  '# File :',adjustl(trim(RAOType)),'_',MD1,'.rao' 
      ELSEIF (adjustl(trim(RAOType)).EQ.'Motion') THEN
	   WRITE(NFILE, '(A8,A7,A1,I1,A4)')  '# File :',adjustl(trim(RAOType)),'_',MD1,'.rao' 
      ENDIF      
	  WRITE(NFILE, '(A1)')  '#'       
	  WRITE(NFILE, '(A34)')  '# Constants used in computations :' 
	  WRITE(NFILE, '(A28,F11.4)')  '#     Reference length     :', REFL
	  WRITE(NFILE, '(A28,F11.4)')  '#     Water density (rho)  :', RHO
	  WRITE(NFILE, '(A28,F11.4)')  '#     Gravity acceleration :', G
      IF (H.LT.0.D0) THEN
	   WRITE(NFILE, '(A35)')  '#     Waterdepth           :  Inf. '
      ELSE
	   WRITE(NFILE, '(A28,F11.4)')  '#     Waterdepth           :',H
      ENDIF
	  WRITE(NFILE, '(A28,A2,2F9.4,A1)')  '#     Ref.pt incident wave :',' (',XW(1),XW(2),')'
	  WRITE(NFILE, '(A43)')  '#            Forward speed :   0.0000  m/s '
	  WRITE(NFILE, '(A1)')  '#'
	  WRITE(NFILE, '(A28,A2,3F9.4,A1)')  '# Reference point of body 1:' ,' (',XC(1),XC(2),XC(3),')'
	  WRITE(NFILE, '(A25)')  '# MEANVALUE :   0.0000E+00'
	  WRITE(NFILE, '(A13)')  '#   AMP/PHASE'
	  WRITE(NFILE, '(A72)')  '#------------------------------------------------------------------------'
      IF (RAOType.EQ.'AddedMass') THEN
	   WRITE(NFILE, '(A15,A9)')  '#RAOTYPE    :  ',RAOType
      ELSEIF (adjustl(trim(RAOType)).EQ.'WaveDamping') THEN
	   WRITE(NFILE, '(A15,A11)')  '#RAOTYPE    :  ',RAOType
      ELSEIF (adjustl(trim(RAOType)).EQ.'Cm') THEN
	   WRITE(NFILE, '(A15,A9)')  '#RAOTYPE    :  ',RAOType
      ELSEIF (adjustl(trim(RAOType)).EQ.'Ca') THEN
	   WRITE(NFILE, '(A15,A11)')  '#RAOTYPE    :  ',RAOType
      ELSEIF (adjustl(trim(RAOType)).EQ.'Excitation') THEN
	   WRITE(NFILE, '(A15,A10)')  '#RAOTYPE    :  ',RAOType
      ELSEIF (adjustl(trim(RAOType)).EQ.'Motion') THEN
	   WRITE(NFILE, '(A15,A6)')  '#RAOTYPE    :  ',RAOType
      ENDIF
	  WRITE(NFILE, '(A15,I1,I1)')  '#COMPONENT  :  ',MD1,MD2
      IF (adjustl(trim(RAOType)).EQ.'AddedMass') THEN
	   WRITE(NFILE, '(A15,A2)')  '#UNIT       :  ','kg'
      ELSEIF (adjustl(trim(RAOType)).EQ.'WaveDamping') THEN
	   WRITE(NFILE, '(A15,A4)')  '#UNIT       :  ','kg/s'
      ELSEIF (adjustl(trim(RAOType)).EQ.'Cm') THEN
	   WRITE(NFILE, '(A15,A2)')  '#UNIT       :  ','kg'
      ELSEIF (adjustl(trim(RAOType)).EQ.'Ca') THEN
	   WRITE(NFILE, '(A15,A4)')  '#UNIT       :  ','kg/s'
      ELSEIF (adjustl(trim(RAOType)).EQ.'Excitation') THEN
	   WRITE(NFILE, '(A15,A9)')  '#UNIT       :  ','N/m, Nm/m'
      ELSEIF (adjustl(trim(RAOType)).EQ.'Motion'.AND.MD1.LE.3) THEN
	   WRITE(NFILE, '(A15,A3)')  '#UNIT       :  ','m/m'
      ELSEIF (adjustl(trim(RAOType)).EQ.'Motion'.AND.MD1.LE.6) THEN
	   WRITE(NFILE, '(A15,A5)')  '#UNIT       :  ','deg/m'
      ENDIF
      
      WRITE(NFILE, '(A10,I6)')  '#NBHEADING',NBETA
      IF (RAOType.EQ.'Motion') THEN
       DO II=1,NBETA
        IF (WVHD(II).LT.0.D0) THEN
         NWVHD(II)=WVHD(II)+360.D0
        ELSE
         NWVHD(II)=WVHD(II)
        ENDIF
       ENDDO
       WRITE(FMT,*) '(A8,',NBETA,'(7X,F7.2))'
	   WRITE(NFILE, FMT)  '#HEADING  ',(NWVHD(II),II=1,NBETA)
!	   WRITE(NFILE, '(A8,<NBETA>(7X,F7.2))')  '#HEADING  ',(NWVHD(II),II=1,NBETA)
      ELSE
       WRITE(FMT,*) '(A8,',NBETA,'(7X,F7.2))'
	   WRITE(NFILE, FMT)  '#HEADING  ',(WVHD(II),II=1,NBETA)
!      WRITE(NFILE, '(A8,<NBETA>(7X,F7.2))')  '#HEADING  ',(WVHD(II),II=1,NBETA)
      ENDIF
	  WRITE(NFILE, '(A63)')  '#---w(r/s)-----------------------------------------------------'

      RETURN
      END SUBROUTINE PrintHeading
      
!   --------------------------------------------------------
!       Print the body in HydroStar output files
!   --------------------------------------------------------

      SUBROUTINE PrintBody_RealVal(NFILE,W1,NBETA,RAOType,VAB)
      USE Const_mod
      IMPLICIT NONE 
	  
      INTEGER  NFILE,NBETA,II
      REAL*8   W1,VAB,REL(NBETA),IMG(NBETA)
      CHARACTER(*)::RAOType
      CHARACTER(LEN=100)::FMT
      
      DO II=1,NBETA
       REL(II)=VAB
       IMG(II)=0.D0
      ENDDO

      WRITE(FMT,*) '(F8.4,',NBETA,'(ES14.6),',NBETA,'(F12.4))'
      WRITE(NFILE,FMT) W1,(REL(II),II=1,NBETA),(IMG(II),II=1,NBETA)

      RETURN
      END SUBROUTINE PrintBody_RealVal
      
!   --------------------------------------------------------
!       Print the body in HydroStar output files
!   --------------------------------------------------------

      SUBROUTINE PrintBody_CmplxVal(NFILE,W1,NBETA,RAOType,CVAB)
      USE Const_mod
      IMPLICIT NONE
	  
      INTEGER  NFILE,NBETA,II
      REAL*8   W1,NREL(NBETA),NIMG(NBETA),REL(NBETA),IMG(NBETA),MDL(NBETA),PHS(NBETA)
      COMPLEX*16 CVAB(NBETA)
      CHARACTER(*)::RAOType
      CHARACTER(LEN=100)::FMT
      
      DO II=1,NBETA
        REL(II)=DREAL(CVAB(II))
        IMG(II)=DIMAG(CVAB(II))
        MDL(II)=CDABS(CVAB(II))
        NREL(II)=-IMG(II)
        NIMG(II)=-REL(II)
        PHS(II)=ATAN2D(NIMG(II),NREL(II))
        IF (PHS(II).LT.0.D0) PHS(II)=PHS(II)+360.D0
      ENDDO
       
      WRITE(FMT,*) '(F8.4,',NBETA,'(ES14.6),',NBETA,'(F12.4))'
      WRITE(NFILE,FMT) W1,(MDL(II),II=1,NBETA),(PHS(II),II=1,NBETA)
!	  WRITE(NFILE,200) W1,(NREL(II),II=1,NBETA),(NIMG(II),II=1,NBETA) 
      
!100	  FORMAT(F8.4,<NBETA>(ES14.6),<NBETA>(F12.4))
!200	  FORMAT(F8.4,<NBETA>(ES14.6),<NBETA>(ES14.6))
      RETURN
      END SUBROUTINE PrintBody_CmplxVal

!   --------------------------------------------------------
!       Print the body in HydroStar output files
!   --------------------------------------------------------

      SUBROUTINE PrintBody_Exfc(NFILE,W1,NBETA,RAOType,CVAB)
      USE Const_mod
      IMPLICIT NONE
	  
      INTEGER  NFILE,NBETA,II
      REAL*8   W1,NREL(NBETA),NIMG(NBETA),REL(NBETA),IMG(NBETA),MDL(NBETA),PHS(NBETA)
      COMPLEX*16 CVAB(NBETA)
      CHARACTER(*)::RAOType
      CHARACTER(LEN=100)::FMT
      
      DO II=1,NBETA
        REL(II)=DREAL(CVAB(II))
        IMG(II)=DIMAG(CVAB(II))
        MDL(II)=CDABS(CVAB(II))
        PHS(II)=ATAN2D(REL(II),IMG(II))
      ENDDO
       
      WRITE(FMT,*) '(F8.4,',NBETA,'(ES14.6),',NBETA,'(F12.4))'
      WRITE(NFILE,FMT) W1,(MDL(II),II=1,NBETA),(PHS(II),II=1,NBETA)

!100	  FORMAT(F8.4,<NBETA>(ES14.6),<NBETA>(F12.4))
!200	  FORMAT(F8.4,<NBETA>(ES14.6),<NBETA>(ES14.6))
      RETURN
      END SUBROUTINE PrintBody_Exfc
      
!   --------------------------------------------------------
!       Print the end in HydroStar output files
!   --------------------------------------------------------

      SUBROUTINE PrintEnd(NFILE)
      IMPLICIT NONE
      
      INTEGER NFILE
            
      WRITE(NFILE, '(A61)')  '#------------------------------------------------------------'
      WRITE(NFILE, '(A8)')  '#ENDFILE' 
       
      RETURN
      END SUBROUTINE PrintEnd
!-------------------------------------------------------------------------------
END MODULE PrintOutput
!*******************************************************************************
