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
!---------------------------------------------------------------------------------------------
!        Data module for declaring input variables in HAMS
!---------------------------------------------------------------------------------------------					
!
        MODULE HAMS_mod
!
        INTEGER,PUBLIC:: IRSP,INFT,OUFT,SYBO,ISOL
        INTEGER,PUBLIC:: NTHREAD

        REAL*8,PUBLIC::  WK1,DWK,BETA1,DBETA

        END MODULE HAMS_mod

        
!---------------------------------------------------------------------------------------------
!        Data module for declaring wave relevant variables in Morison_stick numerical model
!---------------------------------------------------------------------------------------------	
!
        MODULE WaveDyn_mod
!
        INTEGER,PUBLIC:: NPER,NBETA,NK
        PARAMETER(NK=200)
        REAL*8,PUBLIC:: A,H,AMP,V,BETA,BLNR(6,6),BQDR(6,6)
        REAL*8,PUBLIC:: WK,W1,WL,TP,WVN(NK),INFR,OUFR
!
        REAL*8,ALLOCATABLE,PUBLIC:: WVNB(:),WVFQ(:),WVHD(:)
        REAL*8,ALLOCATABLE,PUBLIC:: AMAS(:,:,:),BDMP(:,:,:)

        COMPLEX*16,ALLOCATABLE,PUBLIC:: EXFC(:,:,:),DSPL(:,:,:)
!
        END MODULE WaveDyn_mod

        
!---------------------------------------------------------------------------------------------
!        Data module for declaring body relevant variables in Morison_stick numerical model
!---------------------------------------------------------------------------------------------	
!
        MODULE Body_mod
        
        REAL*8,PUBLIC:: XR(3),XG(3),XB(3),XW(2)
        REAL*8,PUBLIC:: VOL,MASS,REFL
!
        REAL*8,PUBLIC:: IB(3,3),MATX(6,6),CRS(6,6),KSTF(6,6),RAO(6)
!
        END MODULE Body_mod
        
        
!---------------------------------------------------------------------------------------------
!        Data module for declaring variables in panel method
!---------------------------------------------------------------------------------------------
!
        MODULE PanelMesh_mod
!
        INTEGER NELEM,NTND,NTNDD
        INTEGER ISYS,NSYS,ISX,ISY
!
        INTEGER,ALLOCATABLE,PUBLIC:: NCN(:),NCON(:,:),NCOND(:,:)
        INTEGER,ALLOCATABLE,PUBLIC:: IPIV(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: XYZ(:,:),DXYZ_P(:,:),XYZ_P(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: DS(:),PNSZ(:)
!
        END MODULE PanelMesh_mod
        
!---------------------------------------------------------------------------------------------
!        Data module for declaring variables of general computation
!---------------------------------------------------------------------------------------------						
!
        MODULE Const_mod
!
        REAL*8,PUBLIC:: G,RHO,PI
        REAL*8,PUBLIC:: RXY(2,2),RX(2,2),RY(2,2),SY(2,2),SX(2,2)
!
        COMPLEX*16,PUBLIC:: CI
!
        DATA G,PI,RHO/9.80665D0,3.141592653589793D0, 1025.D0/ 
        DATA CI/(0.0D0, 1.0D0)/
        DATA RXY  /  1.D0,  1.D0,  1.D0, -1.D0  /
        DATA RX   /  1.D0, -1.D0,  1.D0,  1.D0  /
        DATA RY   /  1.D0,  1.D0,  1.D0, -1.D0  /
        DATA SX   /  1.D0,  1.D0,  1.D0,  1.D0  /
        DATA SY  /   1.D0, -1.D0, -1.D0,  1.D0  /

        END MODULE Const_mod
        

!---------------------------------------------------------------------------------------------
!        Data module for declaring variables for removing irregular frequencies
!---------------------------------------------------------------------------------------------			
!
        MODULE Inerfs_mod
!
        INTEGER,PUBLIC:: iNELEM,iNTND,tNELEM,tNTND
!     
        INTEGER,ALLOCATABLE,PUBLIC:: iNCN(:),iNCON(:,:),iNCOND(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: iXYZ(:,:),iDXYZ_P(:,:),iXYZ_P(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: iDS(:),iPNSZ(:)
!
        END MODULE Inerfs_mod
    
!---------------------------------------------------------------------------------------------
!        Data module for declaring variables in linear algeraic system
!---------------------------------------------------------------------------------------------						
!
        MODULE LinearMatrix_mod
!
        COMPLEX*16,ALLOCATABLE,PUBLIC:: AMAT(:,:,:),BRMAT(:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: BDMAT(:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: CMAT(:,:,:),DRMAT(:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: DDMAT(:,:)
!
        END MODULE LinearMatrix_mod

!---------------------------------------------------------------------------------------------
!        Data module for declaring variables for potentials
!---------------------------------------------------------------------------------------------						
!
        MODULE Potentials_mod
!
        COMPLEX*16,ALLOCATABLE,PUBLIC:: MXPOT(:,:,:),DPOT(:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: CGRN(:,:,:,:),RKBN(:,:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: DGRN(:,:,:,:),PKBN(:,:,:,:)
!
        END MODULE Potentials_mod


!---------------------------------------------------------------------------------------------
!        Data module for declaring variables for field points
!---------------------------------------------------------------------------------------------
!
        MODULE FieldOutput_mod
        !USE Precision
!
        INTEGER NFP
!
        REAL*8,ALLOCATABLE,PUBLIC:: XFP(:,:)
!
        END MODULE FieldOutput_mod
        
