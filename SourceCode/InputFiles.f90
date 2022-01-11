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
    
!  -------------------------------------------------------------------------------------------
!      Open files, and read the main input file
!---------------------------------------------------------------------------------------------
! 
      SUBROUTINE ReadOpenFiles
      USE HAMS_mod
      USE Body_mod
      USE WaveDyn_mod
      USE FieldOutput_mod
      IMPLICIT NONE

      INTEGER I,J,err,IFS,NPET
! ======================================================
      OPEN(1, FILE='Input/ControlFile.in',      STATUS='OLD')
      OPEN(2, FILE='Input/HullMesh.pnl',        STATUS='OLD')
      OPEN(4, FILE='Input/Hydrostatic.in',     STATUS='UNKNOWN')
      OPEN(9, FILE='Output/ErrorCheck.txt' ,    STATUS='UNKNOWN')
       
! ====================================================

!  H :   Water depth; H<0: For infinite water depth; H>0: For finite water depth; 
!  AMP:  Wave amplitude
!  BETA: Wave incident angle with repect to the x-direction
!  INFT=0,  input wave number; INFT=1, input wave frequency
!  XC,YC,ZC: the coordinates of body rotation center

        READ(1,*) 
        READ(1,*) 
        READ(1,'(14x,f30.15)')     H
        READ(1,*) 
        READ(1,*) 
        READ(1,'(27x,i16)')        SYBO
        READ(1,'(25x,i16)')        INFT
        READ(1,'(25x,i17)')        OUFT
        READ(1,'(26x,i16)')        NPET
        IF (SYBO.EQ.0) THEN
          IFS=0
        ELSEIF (SYBO.EQ.1) THEN
          IFS=2
        ELSE
          PRINT*, 'Warning: SYBO must be 0 or 1.'
          PRINT*
          IFS=0
        ENDIF
        IF (NPET.GE.0) THEN
         NPER=IFS+NPET
         ALLOCATE(WVNB(NPER))
         READ(1,*) (WVNB(I),I=IFS+1,NPER)
        ELSEIF (NPET.LT.0) THEN
         NPET=ABS(NPET)
         NPER=IFS+NPET
         ALLOCATE(WVNB(NPER))
         READ(1,'(27x,f30.15)')     WK1
         READ(1,'(19x,f30.15)')     DWK
         DO I=IFS+1,NPER
          WVNB(I)=WK1+(I-IFS-1)*DWK
         ENDDO
        ENDIF
        READ(1,*) 
        READ(1,*) 
        READ(1,*) 
        READ(1,'(23x,i16)')        NBETA
        IF (NBETA.GT.0) THEN
         ALLOCATE(WVHD(NBETA))
         READ(1,*) (WVHD(I),I=1,NBETA)
        ELSEIF (NBETA.LT.0) THEN
         NBETA=ABS(NBETA)
         ALLOCATE(WVHD(NBETA))
         READ(1,'(20x,f30.15)')     BETA1
         READ(1,'(17x,f30.15)')     DBETA
         DO I=1,NBETA
          WVHD(I)=BETA1+(I-1)*DBETA
         ENDDO
        ENDIF
        
        READ(1,*) 
        READ(1,*) 
        READ(1,'(26x,3f12.3)')      (XR(I), I=1,3)
        WRITE(9,*) "The rotation center is input as (please confirm if it is correct):"
        WRITE(9,'(3f12.3)') (XR(I), I=1,3)
        READ(1,'(26x,f30.15)')     REFL
        READ(1,'(26x,i16)')        ISOL
        READ(1,'(23x,i16)')        IRSP
        READ(1,'(23x,i16)')        NTHREAD
        IF (IRSP.NE.0) THEN
           OPEN(5, FILE='Input/WaterplaneMesh.pnl', STATUS='OLD', &
                IOSTAT=err)
          if (err/=0) then
          PRINT*, 'Error: The waterplane mesh file does not exist.'
          PRINT*
          stop
          endif
        ENDIF

        READ(1,*) 
        READ(1,*) 
        READ(1,'(27x,i16)')        NFP
        ALLOCATE(XFP(NFP,3))
        DO I=1,NFP
!          READ(1,'(26x,3(1x,f10.4))')     (XFP(I,J), J=1,3)
           READ(1,*)     (XFP(I,J), J=1,3)
        ENDDO
             
! ====================================================
!
! Exciting forces
!
        OPEN(21,FILE='Output/Hams_format/OEXFOR1.txt',&
             STATUS='UNKNOWN')
        OPEN(22,FILE='Output/Hams_format/OEXFOR2.txt',&
             STATUS='UNKNOWN')
        OPEN(23,FILE='Output/Hams_format/OEXFOR3.txt',&
             STATUS='UNKNOWN')
        OPEN(24,FILE='Output/Hams_format/OEXFOR4.txt',&
             STATUS='UNKNOWN')
        OPEN(25,FILE='Output/Hams_format/OEXFOR5.txt',&
             STATUS='UNKNOWN')
        OPEN(26,FILE='Output/Hams_format/OEXFOR6.txt',&
             STATUS='UNKNOWN')
!                                           
! Radiation forces and generated waves
!
        OPEN(31,FILE='Output/Hams_format/OAMASS1.txt',&
             STATUS='UNKNOWN')
        OPEN(32,FILE='Output/Hams_format/OAMASS2.txt',&
             STATUS='UNKNOWN')
        OPEN(33,FILE='Output/Hams_format/OAMASS3.txt',&
             STATUS='UNKNOWN')
        OPEN(34,FILE='Output/Hams_format/OAMASS4.txt',&
             STATUS='UNKNOWN')
        OPEN(35,FILE='Output/Hams_format/OAMASS5.txt',&
             STATUS='UNKNOWN')
        OPEN(36,FILE='Output/Hams_format/OAMASS6.txt',&
             STATUS='UNKNOWN')
!        
        OPEN(41,FILE='Output/Hams_format/ODAMPING1.txt',&
             STATUS='UNKNOWN')
        OPEN(42,FILE='Output/Hams_format/ODAMPING2.txt',&
             STATUS='UNKNOWN')
        OPEN(43,FILE='Output/Hams_format/ODAMPING3.txt',&
             STATUS='UNKNOWN')
        OPEN(44,FILE='Output/Hams_format/ODAMPING4.txt',&
             STATUS='UNKNOWN')
        OPEN(45,FILE='Output/Hams_format/ODAMPING5.txt',&
             STATUS='UNKNOWN')
        OPEN(46,FILE='Output/Hams_format/ODAMPING6.txt',&
             STATUS='UNKNOWN')

! Outputs in WAMIT style
!        
        !OPEN(61,FILE='Output/Wamit_format/AmssDamp.1',&
        !     STATUS='UNKNOWN')
        !OPEN(62,FILE='Output/Wamit_format/ExcForce.3',&
        !     STATUS='UNKNOWN')
        !OPEN(63,FILE='Output/Wamit_format/Motion.4',&
        !     STATUS='UNKNOWN')
        !OPEN(64,FILE='Output/Wamit_format/PressureElevation.6p',&
        !     STATUS='UNKNOWN')
        !OPEN(65,FILE='Output/Wamit_format/Hydrostat.hst',&
        !     STATUS='UNKNOWN')
        
        OPEN(61,FILE='Output/Wamit_format/Buoy.1',&
             STATUS='UNKNOWN')
        OPEN(62,FILE='Output/Wamit_format/Buoy.3',&
             STATUS='UNKNOWN')
        OPEN(63,FILE='Output/Wamit_format/Buoy.4',&
             STATUS='UNKNOWN')
        OPEN(64,FILE='Output/Wamit_format/Buoy.6p',&
             STATUS='UNKNOWN')
        OPEN(65,FILE='Output/Wamit_format/Buoy.hst',&
             STATUS='UNKNOWN')
! Outputs in HydroStar style
!        
        OPEN(71,FILE='Output/Hydrostar_format/AddedMass_11.rao',&
             STATUS='UNKNOWN')
        OPEN(72,FILE='Output/Hydrostar_format/AddedMass_12.rao',&
             STATUS='UNKNOWN')
        OPEN(73,FILE='Output/Hydrostar_format/AddedMass_13.rao',&
             STATUS='UNKNOWN')
        OPEN(74,FILE='Output/Hydrostar_format/AddedMass_14.rao',&
             STATUS='UNKNOWN')
        OPEN(75,FILE='Output/Hydrostar_format/AddedMass_15.rao',&
             STATUS='UNKNOWN')
        OPEN(76,FILE='Output/Hydrostar_format/AddedMass_16.rao',&
             STATUS='UNKNOWN')

        OPEN(81,FILE='Output/Hydrostar_format/AddedMass_21.rao',&
             STATUS='UNKNOWN')
        OPEN(82,FILE='Output/Hydrostar_format/AddedMass_22.rao',&
             STATUS='UNKNOWN')
        OPEN(83,FILE='Output/Hydrostar_format/AddedMass_23.rao',&
             STATUS='UNKNOWN')
        OPEN(84,FILE='Output/Hydrostar_format/AddedMass_24.rao',&
             STATUS='UNKNOWN')
        OPEN(85,FILE='Output/Hydrostar_format/AddedMass_25.rao',&
             STATUS='UNKNOWN')
        OPEN(86,FILE='Output/Hydrostar_format/AddedMass_26.rao',&
             STATUS='UNKNOWN')
        
        OPEN(91,FILE='Output/Hydrostar_format/AddedMass_31.rao',&
             STATUS='UNKNOWN')
        OPEN(92,FILE='Output/Hydrostar_format/AddedMass_32.rao',&
             STATUS='UNKNOWN')
        OPEN(93,FILE='Output/Hydrostar_format/AddedMass_33.rao',&
             STATUS='UNKNOWN')
        OPEN(94,FILE='Output/Hydrostar_format/AddedMass_34.rao',&
             STATUS='UNKNOWN')
        OPEN(95,FILE='Output/Hydrostar_format/AddedMass_35.rao',&
             STATUS='UNKNOWN')
        OPEN(96,FILE='Output/Hydrostar_format/AddedMass_36.rao',&
             STATUS='UNKNOWN')
        
        OPEN(101,FILE='Output/Hydrostar_format/AddedMass_41.rao',&
             STATUS='UNKNOWN')
        OPEN(102,FILE='Output/Hydrostar_format/AddedMass_42.rao',&
             STATUS='UNKNOWN')
        OPEN(103,FILE='Output/Hydrostar_format/AddedMass_43.rao',&
             STATUS='UNKNOWN')
        OPEN(104,FILE='Output/Hydrostar_format/AddedMass_44.rao',&
             STATUS='UNKNOWN')
        OPEN(105,FILE='Output/Hydrostar_format/AddedMass_45.rao',&
             STATUS='UNKNOWN')
        OPEN(106,FILE='Output/Hydrostar_format/AddedMass_46.rao',&
             STATUS='UNKNOWN')
        
        OPEN(111,FILE='Output/Hydrostar_format/AddedMass_51.rao',&
             STATUS='UNKNOWN')
        OPEN(112,FILE='Output/Hydrostar_format/AddedMass_52.rao',&
             STATUS='UNKNOWN')
        OPEN(113,FILE='Output/Hydrostar_format/AddedMass_53.rao',&
             STATUS='UNKNOWN')
        OPEN(114,FILE='Output/Hydrostar_format/AddedMass_54.rao',&
             STATUS='UNKNOWN')
        OPEN(115,FILE='Output/Hydrostar_format/AddedMass_55.rao',&
             STATUS='UNKNOWN')
        OPEN(116,FILE='Output/Hydrostar_format/AddedMass_56.rao',&
             STATUS='UNKNOWN')
        
        OPEN(121,FILE='Output/Hydrostar_format/AddedMass_61.rao',&
             STATUS='UNKNOWN')
        OPEN(122,FILE='Output/Hydrostar_format/AddedMass_62.rao',&
             STATUS='UNKNOWN')
        OPEN(123,FILE='Output/Hydrostar_format/AddedMass_63.rao',&
             STATUS='UNKNOWN')
        OPEN(124,FILE='Output/Hydrostar_format/AddedMass_64.rao',&
             STATUS='UNKNOWN')
        OPEN(125,FILE='Output/Hydrostar_format/AddedMass_65.rao',&
             STATUS='UNKNOWN')
        OPEN(126,FILE='Output/Hydrostar_format/AddedMass_66.rao',&
             STATUS='UNKNOWN')

        OPEN(131,FILE='Output/Hydrostar_format/WaveDamping_11.rao',&
             STATUS='UNKNOWN')
        OPEN(132,FILE='Output/Hydrostar_format/WaveDamping_12.rao',&
             STATUS='UNKNOWN')
        OPEN(133,FILE='Output/Hydrostar_format/WaveDamping_13.rao',&
             STATUS='UNKNOWN')
        OPEN(134,FILE='Output/Hydrostar_format/WaveDamping_14.rao',&
             STATUS='UNKNOWN')
        OPEN(135,FILE='Output/Hydrostar_format/WaveDamping_15.rao',&
             STATUS='UNKNOWN')
        OPEN(136,FILE='Output/Hydrostar_format/WaveDamping_16.rao',&
             STATUS='UNKNOWN')

        OPEN(141,FILE='Output/Hydrostar_format/WaveDamping_21.rao',&
             STATUS='UNKNOWN')
        OPEN(142,FILE='Output/Hydrostar_format/WaveDamping_22.rao',&
             STATUS='UNKNOWN')
        OPEN(143,FILE='Output/Hydrostar_format/WaveDamping_23.rao',&
             STATUS='UNKNOWN')
        OPEN(144,FILE='Output/Hydrostar_format/WaveDamping_24.rao',&
             STATUS='UNKNOWN')
        OPEN(145,FILE='Output/Hydrostar_format/WaveDamping_25.rao',&
             STATUS='UNKNOWN')
        OPEN(146,FILE='Output/Hydrostar_format/WaveDamping_26.rao',&
             STATUS='UNKNOWN')

        OPEN(151,FILE='Output/Hydrostar_format/WaveDamping_31.rao',&
             STATUS='UNKNOWN')
        OPEN(152,FILE='Output/Hydrostar_format/WaveDamping_32.rao',&
             STATUS='UNKNOWN')
        OPEN(153,FILE='Output/Hydrostar_format/WaveDamping_33.rao',&
             STATUS='UNKNOWN')
        OPEN(154,FILE='Output/Hydrostar_format/WaveDamping_34.rao',&
             STATUS='UNKNOWN')
        OPEN(155,FILE='Output/Hydrostar_format/WaveDamping_35.rao',&
             STATUS='UNKNOWN')
        OPEN(156,FILE='Output/Hydrostar_format/WaveDamping_36.rao',&
             STATUS='UNKNOWN')

        OPEN(161,FILE='Output/Hydrostar_format/WaveDamping_41.rao',&
             STATUS='UNKNOWN')
        OPEN(162,FILE='Output/Hydrostar_format/WaveDamping_42.rao',&
             STATUS='UNKNOWN')
        OPEN(163,FILE='Output/Hydrostar_format/WaveDamping_43.rao',&
             STATUS='UNKNOWN')
        OPEN(164,FILE='Output/Hydrostar_format/WaveDamping_44.rao',&
             STATUS='UNKNOWN')
        OPEN(165,FILE='Output/Hydrostar_format/WaveDamping_45.rao',&
             STATUS='UNKNOWN')
        OPEN(166,FILE='Output/Hydrostar_format/WaveDamping_46.rao',&
             STATUS='UNKNOWN')

        OPEN(171,FILE='Output/Hydrostar_format/WaveDamping_51.rao',&
             STATUS='UNKNOWN')
        OPEN(172,FILE='Output/Hydrostar_format/WaveDamping_52.rao',&
             STATUS='UNKNOWN')
        OPEN(173,FILE='Output/Hydrostar_format/WaveDamping_53.rao',&
             STATUS='UNKNOWN')
        OPEN(174,FILE='Output/Hydrostar_format/WaveDamping_54.rao',&
             STATUS='UNKNOWN')
        OPEN(175,FILE='Output/Hydrostar_format/WaveDamping_55.rao',&
             STATUS='UNKNOWN')
        OPEN(176,FILE='Output/Hydrostar_format/WaveDamping_56.rao',&
             STATUS='UNKNOWN')

        OPEN(181,FILE='Output/Hydrostar_format/WaveDamping_61.rao',&
             STATUS='UNKNOWN')
        OPEN(182,FILE='Output/Hydrostar_format/WaveDamping_62.rao',&
             STATUS='UNKNOWN')
        OPEN(183,FILE='Output/Hydrostar_format/WaveDamping_63.rao',&
             STATUS='UNKNOWN')
        OPEN(184,FILE='Output/Hydrostar_format/WaveDamping_64.rao',&
             STATUS='UNKNOWN')
        OPEN(185,FILE='Output/Hydrostar_format/WaveDamping_65.rao',&
             STATUS='UNKNOWN')
        OPEN(186,FILE='Output/Hydrostar_format/WaveDamping_66.rao',&
             STATUS='UNKNOWN')
        
        OPEN(191,FILE='Output/Hydrostar_format/Excitation_1.rao',&
             STATUS='UNKNOWN')
        OPEN(192,FILE='Output/Hydrostar_format/Excitation_2.rao',&
             STATUS='UNKNOWN')
        OPEN(193,FILE='Output/Hydrostar_format/Excitation_3.rao',&
             STATUS='UNKNOWN')
        OPEN(194,FILE='Output/Hydrostar_format/Excitation_4.rao',&
             STATUS='UNKNOWN')
        OPEN(195,FILE='Output/Hydrostar_format/Excitation_5.rao',&
             STATUS='UNKNOWN')
        OPEN(196,FILE='Output/Hydrostar_format/Excitation_6.rao',&
             STATUS='UNKNOWN')
        
        OPEN(201,FILE='Output/Hydrostar_format/Motion_1.rao',&
             STATUS='UNKNOWN')
        OPEN(202,FILE='Output/Hydrostar_format/Motion_2.rao',&
             STATUS='UNKNOWN')
        OPEN(203,FILE='Output/Hydrostar_format/Motion_3.rao',&
             STATUS='UNKNOWN')
        OPEN(204,FILE='Output/Hydrostar_format/Motion_4.rao',&
             STATUS='UNKNOWN')
        OPEN(205,FILE='Output/Hydrostar_format/Motion_5.rao',&
             STATUS='UNKNOWN')
        OPEN(206,FILE='Output/Hydrostar_format/Motion_6.rao',&
             STATUS='UNKNOWN')
!          
      END Subroutine ReadOpenFiles
 

