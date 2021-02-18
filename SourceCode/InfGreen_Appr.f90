! 
!==================================================================================
! 
!     Purpose: This program evaluates the three-dimensional free-surface
!              Green function and derivatives in deep water
! 
!              Code Original Author: Hui Liang       created on  2016.12.26    
! 
!  License:
!
!
!    This routine is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option) 
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3  
!    or later), along with this routine. If not, see <http://www.gnu.org/licenses/>.
!
!  Modified on:
! 
!    January 02, 2018
! 
!  Reference:
! 
!
!	[1]	H. Wu, C. Zhang, Y. Zhu, W. Li, D. Wan, F. Noblesse, 
!		A global approximation to the Green function for 
!		diffraction radiation of water waves, 
!		Eur. J. Mech. B Fluids 65 (2017) 54-64.
!
!	[2]	H. Liang, H. Wu, F. Noblesse,
!		Validation of a global approximation for 
!		wave diffraction-radiation in deep water,
!		Appl. Ocean Res. 74 (2018) 80-86.
! 
!   Remarks
! 
!	The local-flow component is approximated by mean of the global 
!	approximations [1]. The computations reported in [2] provides
!	strong evidence that the global approximations are sufficiently 
!	accurate to compute linear and second-order wave loads in practice.
!	
!
!	It should be noted that the Rankine source term -2/d appeared in L_z 
!	given by (8a) in Ref. [2] is not evaluated here (there is a typo in (8a)).  
!	These Rankine source components are required to evaluate in another routine  
!	due to strong singular behavior.
!
!	For any questions, please contact: lianghuistar@gmail.com
!                   
!   We define the flow-field point (x,y,z) and source point (xi,eta,zeta).               
!   Please note that all variables are non-dimensionalized with respect to 
!	the wavenumber k0
!                   
!   The Green function is defiend as                
!                   
!   G = -1/r-1/d+GF,
!
!   where   
!                                          
!   r = sqrt((x-xi)^2+(y-eta)^2+(z-zeta)^2);                  
!   d = sqrt((x-xi)^2+(y-eta)^2+(z+zeta)^2).
!
!   Parameters:
!
!       Input:      dx  --- x-xi
!                   dy  --- y-eta
!                   vv  --- z+zeta
! 
!       Output:     GF(0)  --- free surface itself GF
!                   GF(1)  --- x-derivative GF_x
!                   GF(2)  --- y-derivative GF_y
!                   GF(3)  --- part of z-derivative GF_z = GF(3) - 2/d
!
! ==================================================================================

module INFG3D_Open

      USE Const_mod
      implicit none

      PUBLIC :: INFGREEN3D

!      REAL*8            ::	pi	=	4.0D0*datan(1.0D0)
      REAL*8			::	gama=	0.5772156649D0
      COMPLEX*16		::  Im	=	dcmplx(0.0D0,1.0D0)
	
      REAL*8			::	GF_alpha
      REAL*8			::	GF_beta
      REAL*8			::	GF_sigma
      REAL*8			::	GF_rho
      REAL*8			::	GF_dd

      REAL*8			::	GF_hh
      REAL*8			::	GF_vv

contains

!***********************************************************************
!        FREE-SURFACE GREEN FUNCTION AND ITS NORMAL DERIVATIVES
!***********************************************************************

SUBROUTINE INFGREEN3D(XF,XP,YF,YP,ZF,ZP,V,GRN,FLAG)

        IMPLICIT NONE
        
        INTEGER,INTENT(IN):: FLAG
        REAL*8,INTENT(IN):: XF,XP,YF,YP,ZF,ZP,V
        COMPLEX*16,INTENT(OUT):: GRN(4)
        
        REAL*8:: R,R1,RR,XPQ,YPQ,ZPQ,PI
        COMPLEX*16:: GF,GFH

        DATA  PI /3.141592653589793D0/

        XPQ=XF-XP
        YPQ=YF-YP
        ZPQ=ZP+ZF
        
        RR  =SQRT(XPQ**2+YPQ**2)
        R =SQRT(RR**2+(ZP-ZF)**2)
        R1 =SQRT(RR**2+ZPQ*ZPQ)
        
        GRN(:)=CMPLX(0.D0,0.D0)
        
        CALL HavelockGF(V*RR,V*ZPQ,GF,GFH)


        IF (ABS(V).LT.1.E-8) THEN

         IF (FLAG.EQ.1) THEN
          GRN(1)=GRN(1)+1.D0/R1
          GRN(2)=GRN(2)-XPQ/R1**3
          GRN(3)=GRN(3)-YPQ/R1**3
          GRN(4)=GRN(4)-ZPQ/R1**3
         ELSE
          GRN(1)=GRN(1)+1.D0/R1+1.D0/R
          GRN(2)=GRN(2)-XPQ/R1**3-XPQ/R**3
          GRN(3)=GRN(3)-YPQ/R1**3-YPQ/R**3
          GRN(4)=GRN(4)-ZPQ/R1**3-(ZF-ZP)/R**3
         ENDIF

        ELSEIF (ABS(V+1.D0).LT.1.E-8) THEN
            
         IF (FLAG.EQ.1) THEN
          GRN(1)=GRN(1)-1.D0/R1
          GRN(2)=GRN(2)+XPQ/R1**3
          GRN(3)=GRN(3)+YPQ/R1**3
          GRN(4)=GRN(4)+ZPQ/R1**3
         ELSE
          GRN(1)=GRN(1)-1.D0/R1+1.D0/R
          GRN(2)=GRN(2)+XPQ/R1**3-XPQ/R**3
          GRN(3)=GRN(3)+YPQ/R1**3-YPQ/R**3
          GRN(4)=GRN(4)+ZPQ/R1**3-(ZF-ZP)/R**3
         ENDIF

        ELSE

         IF(RR.GT.1.E-6) THEN
          GRN(2)=-V**2.D0*GFH*XPQ/RR
          GRN(3)=-V**2.D0*GFH*YPQ/RR
         ELSE
          GRN(2)=CMPLX(0.D0,0.D0)
          GRN(3)=CMPLX(0.D0,0.D0)
         ENDIF
          
          GRN(1)=-V*GF
          GRN(4)=-V**2.D0*GF+2.D0*V/R1
           
         IF (FLAG.EQ.1) THEN
          GRN(1)=GRN(1)+1.D0/R1
          GRN(2)=GRN(2)-XPQ/R1**3
          GRN(3)=GRN(3)-YPQ/R1**3
          GRN(4)=GRN(4)-ZPQ/R1**3
         ELSE
          GRN(1)=GRN(1)+1.D0/R1+1.D0/R
          GRN(2)=GRN(2)-XPQ/R1**3-XPQ/R**3
          GRN(3)=GRN(3)-YPQ/R1**3-YPQ/R**3
          GRN(4)=GRN(4)-ZPQ/R1**3-(ZF-ZP)/R**3
         ENDIF

        ENDIF

      RETURN
END SUBROUTINE INFGREEN3D    
    
!===============================================================
subroutine HavelockGF(hh,vv,GF,GFh)
!
    implicit none
! --- Variables -------------------------------------------
    REAL*8,intent(in)        ::    hh,vv
    COMPLEX*16,intent(out)    ::    GF, GFh

! --- Local variables -------------------------------------
!    REAL*8                    ::    hh

!    hh    =    dsqrt(dx*dx+dy*dy)
    call    GF_DivParameters(hh,vv)
 
    GF    =    GF_Func_L0(hh,vv)+GF_Func_W(hh,vv)
    GFh =    GF_Func_Ls(hh,vv)+GF_Func_Wh(hh,vv)

end subroutine HavelockGF
    
!=============================================================
! Calculate parameters
subroutine GF_DivParameters(hh,vv)
    implicit none

! --- Variables -------------------------------------------
    REAL*8,intent(in)        ::    hh,vv

! --- Local variables -------------------------------------

    GF_dd		=	dsqrt(hh*hh+vv*vv)     ! r1=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)
	GF_alpha	=	-vv/GF_dd
	GF_beta		=	hh/GF_dd
	GF_sigma	=	hh/(GF_dd-vv)
	GF_rho		=	GF_dd/(1.0D0+GF_dd)
    
	return
end subroutine GF_DivParameters

!=============================================================
function GF_Func_L0(hh,vv)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	hh,vv

! --- Local variables -------------------------------------
	integer					::	ii
	REAL*8					::	PP,Lp,Cal_L0
	REAL*8					::	GF_Func_L0

	PP	=	dlog(0.5D0*(GF_dd-vv))+gama-2.0D0*GF_dd*GF_dd
	PP	=	dexp(vv)*PP
	PP	=	PP+GF_dd*GF_dd-vv      ! Eq. (22b) in [1] 
    
	Lp	=	GF_Func_Lp(hh,vv)   ! L', i.e. Eq.(22c) in [1]
    
	GF_Func_L0	=	2.0D0*PP/(1.0D0+GF_dd**3)+2.0D0*Lp  ! Eq. (22a) in [1] excluding the term -1/d

	return
end function GF_Func_L0

!=============================================================
function GF_Func_Lp(hh,vv)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	hh,vv

! --- Local variables -------------------------------------
	REAL*8					::	A,B,C,D,RR,Lp
	REAL*8					::	GF_Func_Lp

	A	=	GF_FuncA(GF_rho)    ! below: Eq.(33d) ~ (33g) in [1]
	B	=	GF_FuncB(GF_rho)
	C	=	GF_FuncC(GF_rho)
	D	=	GF_FuncD(GF_rho)

	RR	=	(1.0D0-GF_beta)*A    ! below: Eq (26b) in [1]
	RR	=	RR-GF_beta*B
	RR	=	RR-GF_alpha*C/(1.0D0+6.0D0*GF_alpha*GF_rho*(1.0D0-GF_rho))
	RR	=	RR+GF_beta*(1.0D0-GF_beta)*D

	GF_Func_Lp	=	GF_rho*(1.0D0-GF_rho)**3*RR     ! Eq.(26a) in [1]

	return
end function GF_Func_Lp

!=============================================================
function GF_Func_W(hh,vv)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	hh,vv

! --- Local variables -------------------------------------
	REAL*8					::	H0,J0
	COMPLEX*16				::	GF_Func_W

	H0	=	StruveH0(hh)
	J0	=	BesselJ0(hh)

	GF_Func_W	=	2.0D0*pi*(H0-Im*J0)*dexp(vv)

	return
end function GF_Func_W

!=============================================================
function GF_FuncA(tt)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	tt

! --- Local variables -------------------------------------
	REAL*8					::	A0,A1,A2,A3,A4
	REAL*8					::	A5,A6,A7,A8,A9
	REAL*8					::	GF_FuncA

	data	A0,A1,A2,A3,A4,A5,A6,A7,A8,A9			&
		/	+1.21D0,	-13.328D0,	+215.896D0,		&
			-1763.96D0,	+8418.94D0,	-24314.21D0,	&
			+42002.57D0,-41592.9D0,	21859.0D0,		&
			-4838.6D0	/

	GF_FuncA	=	(A8+A9*tt)*tt
	GF_FuncA	=	(A7+GF_FuncA)*tt
	GF_FuncA	=	(A6+GF_FuncA)*tt
	GF_FuncA	=	(A5+GF_FuncA)*tt
	GF_FuncA	=	(A4+GF_FuncA)*tt
	GF_FuncA	=	(A3+GF_FuncA)*tt
	GF_FuncA	=	(A2+GF_FuncA)*tt
	GF_FuncA	=	(A1+GF_FuncA)*tt
	GF_FuncA	=	A0+GF_FuncA

	return
end function GF_FuncA

!=============================================================
function GF_FuncB(tt)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	tt

! --- Local variables -------------------------------------
	REAL*8					::	B0,B1,B2,B3,B4
	REAL*8					::	B5,B6,B7,B8,B9
	REAL*8					::	GF_FuncB

	data	B0,B1,B2,B3,B4,B5,B6,B7,B8,B9			&
		/	+0.938D0,	+5.373D0,	-67.92D0,		&
			+796.534D0,	-4780.77D0,	+17137.74D0,	&
			-36618.81D0,+44894.06D0,-29030.24D0,	&
			+7671.22D0	/

	GF_FuncB	=	(B8+B9*tt)*tt
	GF_FuncB	=	(B7+GF_FuncB)*tt
	GF_FuncB	=	(B6+GF_FuncB)*tt
	GF_FuncB	=	(B5+GF_FuncB)*tt
	GF_FuncB	=	(B4+GF_FuncB)*tt
	GF_FuncB	=	(B3+GF_FuncB)*tt
	GF_FuncB	=	(B2+GF_FuncB)*tt
	GF_FuncB	=	(B1+GF_FuncB)*tt
	GF_FuncB	=	B0+GF_FuncB

	return
end function GF_FuncB

!=============================================================
function GF_FuncC(tt)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	tt

! --- Local variables -------------------------------------
	REAL*8					::	C0,C1,C2,C3,C4
	REAL*8					::	C5,C6,C7
	REAL*8					::	GF_FuncC

	data	C0,C1,C2,C3,C4,C5,C6,C7					&
		/	+1.268D0,	-9.747D0,	+209.653D0,		&
			-1397.89D0,	+5155.67D0,	-9844.35D0,		&
			+9136.4D0,	-3272.62D0		/

	GF_FuncC	=	(C6+C7*tt)*tt
	GF_FuncC	=	(C5+GF_FuncC)*tt
	GF_FuncC	=	(C4+GF_FuncC)*tt
	GF_FuncC	=	(C3+GF_FuncC)*tt
	GF_FuncC	=	(C2+GF_FuncC)*tt
	GF_FuncC	=	(C1+GF_FuncC)*tt
	GF_FuncC	=	C0+GF_FuncC

	return
end function GF_FuncC

!=============================================================
function GF_FuncD(tt)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	tt

! --- Local variables -------------------------------------
	REAL*8					::	D0,D1,D2,D3,D4
	REAL*8					::	D5,D6,D7,D8,D9
	REAL*8					::	GF_FuncD
			
	data	D0,D1,D2,D3,D4,D5,D6,D7,D8,D9			&
		/	+0.632D0,	-40.97D0,	+667.16D0,		&
			-6072.07D0,	+31127.39D0,-96293.05D0,	&
			+181856.75D0,			-205690.43D0,	&
			+128170.2D0,-33744.6D0	/
    
	GF_FuncD	=	(D8+D9*tt)*tt
	GF_FuncD	=	(D7+GF_FuncD)*tt
	GF_FuncD	=	(D6+GF_FuncD)*tt
	GF_FuncD	=	(D5+GF_FuncD)*tt
	GF_FuncD	=	(D4+GF_FuncD)*tt
	GF_FuncD	=	(D3+GF_FuncD)*tt
	GF_FuncD	=	(D2+GF_FuncD)*tt
	GF_FuncD	=	(D1+GF_FuncD)*tt
	GF_FuncD	=	D0+GF_FuncD

	return
end function GF_FuncD
    
!=============================================================
function GF_Func_Ls(hh,vv)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	hh,vv

! --- Local variables -------------------------------------
	REAL*8					::	PS,QS,Lsp
	REAL*8					::	GF_Func_Ls

	PS	=	(GF_beta+hh)/(GF_dd-vv)
	PS	=	PS-2.0D0*GF_beta+2.0D0*dexp(vv)*GF_dd-hh   ! the first equation in Eq.(27b) in [1]
 
	QS	=	dexp(-GF_dd)*(1.0D0-GF_beta)
	QS	=	QS*(1.0D0+GF_dd/(1.0D0+GF_dd**3))    ! the second equation in Eq.(27b) in [1]
 
	Lsp	=	GF_Func_Lsp(hh,vv)
 
	GF_Func_Ls	=	2.0D0*PS/(1.0D0+GF_dd**3)-4.0D0*QS+2.0D0*Lsp    ! Eq.(27a) in [1]

	return
end function GF_Func_Ls

!=============================================================
function GF_Func_Lsp(hh,vv)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	hh,vv

! --- Local variables -------------------------------------
	REAL*8					::	A,B,C,RR
	REAL*8					::	GF_Func_Lsp

	A	=	GF_dFuncA(GF_rho)   ! below: Eq.(34e) ~ (34g) in [1]
	B	=	GF_dFuncB(GF_rho)
	C	=	GF_dFuncC(GF_rho)

	RR	=	GF_beta*A           ! below: Eq (31b) in [1]
	RR	=	RR-(1.0D0-GF_alpha)*B
	RR	=	RR+GF_beta*(1.0D0-GF_beta)*GF_rho*(1.0D0-2.0D0*GF_rho)*C

	GF_Func_Lsp	=	GF_rho*(1.0D0-GF_rho)**3*RR      ! Eq.(31a) in [1]

	return
end function GF_Func_Lsp

!=============================================================
function GF_Func_Wh(hh,vv)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	hh,vv

! --- Local variables -------------------------------------
	REAL*8					::	H1,J1
	COMPLEX*16				::	GF_Func_Wh

	H1	=	StruveH1(hh)
	J1	=	BesselJ1(hh)
    
	GF_Func_Wh	=	2.0D0*pi*(2.0D0/pi-H1+Im*J1)*dexp(vv)
!	GF_Func_Wh	=	2.0D0*pi*(-H1+Im*J1)*dexp(vv)
    
	return
end function GF_Func_Wh

!=============================================================
function GF_dFuncA(tt)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	tt

! --- Local variables -------------------------------------
	REAL*8					::	A0,A1,A2,A3,A4
	REAL*8					::	A5,A6,A7,A8,A9
	REAL*8					::	GF_dFuncA
			
	data	A0,A1,A2,A3,A4,A5,A6,A7,A8,A9			&
		/	+2.948D0,	-24.53D0,	+249.69D0,		&
			-754.85D0,	-1187.71D0,	+16370.75D0,	&
			-48811.41D0,+68220.87D0,-46688.0D0,		&
			+12622.25D0	/

	GF_dFuncA	=	(A8+A9*tt)*tt
	GF_dFuncA	=	(A7+GF_dFuncA)*tt
	GF_dFuncA	=	(A6+GF_dFuncA)*tt
	GF_dFuncA	=	(A5+GF_dFuncA)*tt
	GF_dFuncA	=	(A4+GF_dFuncA)*tt
	GF_dFuncA	=	(A3+GF_dFuncA)*tt
	GF_dFuncA	=	(A2+GF_dFuncA)*tt
	GF_dFuncA	=	(A1+GF_dFuncA)*tt
	GF_dFuncA	=	A0+GF_dFuncA

	return
end function GF_dFuncA

!=============================================================
function GF_dFuncB(tt)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	tt

! --- Local variables -------------------------------------
	REAL*8					::	B0,B1,B2,B3,B4
	REAL*8					::	B5,B6,B7,B8,B9
	REAL*8					::	GF_dFuncB
			
	data	B0,B1,B2,B3,B4,B5,B6,B7,B8,B9			&
		/	+1.11D0,	+2.894D0,	-76.765D0,		&
			+1565.35D0,	-11336.19D0,+44270.15D0,	&
			-97014.11D0,+118879.26D0,-76209.82D0,	&
			+19923.28D0	/
    
	GF_dFuncB	=	(B8+B9*tt)*tt
	GF_dFuncB	=	(B7+GF_dFuncB)*tt
	GF_dFuncB	=	(B6+GF_dFuncB)*tt
	GF_dFuncB	=	(B5+GF_dFuncB)*tt
	GF_dFuncB	=	(B4+GF_dFuncB)*tt
	GF_dFuncB	=	(B3+GF_dFuncB)*tt
	GF_dFuncB	=	(B2+GF_dFuncB)*tt
	GF_dFuncB	=	(B1+GF_dFuncB)*tt
	GF_dFuncB	=	B0+GF_dFuncB

	return
end function GF_dFuncB

!=============================================================
function GF_dFuncC(tt)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	tt

! --- Local variables -------------------------------------
	REAL*8					::	C0,C1,C2,C3,C4,C5
	REAL*8					::	GF_dFuncC
			
	data	C0,C1,C2,C3,C4,C5						&
		/	+14.19D0,	-148.24D0,	+847.8D0,		&
			-2318.58D0,	+3168.35D0,	-1590.27D0		/
    
	GF_dFuncC	=	(C4+C5*tt)*tt
	GF_dFuncC	=	(C3+GF_dFuncC)*tt
	GF_dFuncC	=	(C2+GF_dFuncC)*tt
	GF_dFuncC	=	(C1+GF_dFuncC)*tt
	GF_dFuncC	=	C0+GF_dFuncC

	return
end function GF_dFuncC

!**************************************
function BesselJ0(xx)
!
	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	xx

! --- Local variables -------------------------------------
	REAL*8					::	yy,y2,f0,theta0
	REAL*8					::	BesselJ0				
	REAL*8					::	P0,P1,P2,P3,P4,P5,P6
	REAL*8					::	R0,R1,R2,R3,R4,R5
	REAL*8					::  S1,S2,S3,S4,S5
    
	data	P0,P1,P2,P3,P4,P5,P6		&
		/	+0.999999999D0,	-2.249999879D0,	+1.265623060D0,	&
			-0.316394552D0,	+0.044460948D0,	-0.003954479D0,	&
			+0.000212950D0	/
    
	data	R0,R1,R2,R3,R4,R5			&
		/	+0.79788454D0,	-0.00553897D0,	+0.00099336D0,	&
			-0.00044346D0,	+0.00020445D0,	-0.00004959D0	/
    
	data	S1,S2,S3,S4,S5				&
		/	-0.04166592D0,	+0.00239399D0,	-0.00073984D0,	&
			+0.00031099D0,	-0.00007605D0	/


	if(xx <= 3.0D0) then
		yy	=	(xx/3.0D0)**2
		BesselJ0	=	 P0+(P1+(P2+(P3+(P4+(P5+P6*yy)*yy)*yy)*yy)*yy)*yy
	else
		yy	=	3.0D0/xx
		y2	=	yy**2
        
		f0			=	R0+(R1+(R2+(R3+(R4+R5*y2)*y2)*y2)*y2)*y2

		theta0		=	 xx	-	0.25D0*pi		&
						+(S1+(S2+(S3+(S4+S5*y2)*y2)*y2)*y2)*yy

		BesselJ0	=	f0*dcos(theta0)/dsqrt(xx)
	end if

	return
end function BesselJ0

!===============================================================
function BesselJ1(xx)
!
	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	xx

! --- Local variables -------------------------------------
	REAL*8					::	yy,y2,f1,theta1
	REAL*8					::	BesselJ1
    
	REAL*8					::	P0,P1,P2,P3,P4,P5,P6
	REAL*8					::	R0,R1,R2,R3,R4,R5
	REAL*8					::  S1,S2,S3,S4,S5
    
	data	P0,P1,P2,P3,P4,P5,P6		&
		/	+0.500000000D0,	-0.562499992D0,	+0.210937377D0,	&
			-0.039550040D0,	+0.004447331D0,	-0.000330547D0,	&
			+0.000015525D0	/
    
	data	R0,R1,R2,R3,R4,R5			&
		/	+0.79788459D0,	+0.01662008D0,	-0.00187002D0,	&
			+0.00068519D0,	-0.00029440D0,	+0.00006952D0	/
    
	data	S1,S2,S3,S4,S5				&
		/	+0.12499895D0,	-0.00605240D0,	+0.00135825D0,	&
			-0.00049616D0,	+0.00011531D0	/

	if(xx <= 3.0D0) then
		yy	=	xx/3.0D0
		y2	=	yy*yy
		BesselJ1	=	P0+(P1+(P2+(P3+(P4+(P5+P6*y2)*y2)*y2)*y2)*y2)*y2

		BesselJ1	=	BesselJ1*xx
	else
		yy	=	3.0D0/xx
		y2	=	yy*yy
		f1			=	R0+(R1+(R2+(R3+(R4+R5*y2)*y2)*y2)*y2)*y2

		theta1		=	 xx	-	0.75D0*pi		&
						+(S1+(S2+(S3+(S4+S5*y2)*y2)*y2)*y2)*yy

		BesselJ1	=	f1*dcos(theta1)/dsqrt(xx)
	end if

	return
end function BesselJ1

!===============================================================
function BesselY0(xx)
!
	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	xx

! --- Local variables -------------------------------------
	REAL*8					::	f0,theta0
	REAL*8					::	BesselY0				

	if(xx <= 3.0D0) then
		BesselY0	=	(2.0D0/pi)*dlog(xx/2.0D0)*BesselJ0(xx)	&	
						+0.367466907D0							&
						+0.605593797D0*(xx/3.0D0)**2			&														
						-0.743505078D0*(xx/3.0D0)**4			&														
						+0.253005481D0*(xx/3.0D0)**6			&														
						-0.042619616D0*(xx/3.0D0)**8			&														
						+0.004285691D0*(xx/3.0D0)**10			&														
						-0.000250716D0*(xx/3.0D0)**12														
	else
		f0			=	 0.79788454D0					&
						-0.00553897D0*(3.0D0/xx)**2		&														
						+0.00099336D0*(3.0D0/xx)**4		&														
						-0.00044346D0*(3.0D0/xx)**6		&														
						+0.00020445D0*(3.0D0/xx)**8		&														
						-0.00004959D0*(3.0D0/xx)**10														

		theta0		=	 xx		-	0.25D0*pi			&
						-0.04166592D0*(3.0D0/xx)**1		&														
						+0.00239399D0*(3.0D0/xx)**3		&														
						-0.00073984D0*(3.0D0/xx)**5		&														
						+0.00031099D0*(3.0D0/xx)**7		&														
						-0.00007605D0*(3.0D0/xx)**9														

		BesselY0	=	f0*dsin(theta0)/dsqrt(xx)																		

	end if

	return
end function BesselY0

!===============================================================
	function BesselY1(xx)
!
	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	xx

! --- Local variables -------------------------------------
	REAL*8					::	f1,theta1
	REAL*8					::	BesselY1				

	if(xx <= 3.0D0) then
		BesselY1	=	(2.0D0/pi)*(dlog(xx/2.0D0)*BesselJ1(xx)-1.0D0/xx)	&	
						+0.07373571D0*(xx/3.0D0)**1				&														
						+0.72276433D0*(xx/3.0D0)**3				&														
						-0.43885620D0*(xx/3.0D0)**5				&														
						+0.10418264D0*(xx/3.0D0)**7				&														
						-0.01340825D0*(xx/3.0D0)**9				&														
						+0.00094249D0*(xx/3.0D0)**11														
	else
		f1			=	 0.79788459D0						&
						+0.01662008D0*(3.0D0/xx)**2		&														
						-0.00187002D0*(3.0D0/xx)**4		&														
						+0.00068519D0*(3.0D0/xx)**6		&														
						-0.00029440D0*(3.0D0/xx)**8		&														
						+0.00006952D0*(3.0D0/xx)**10														

		theta1		=	 xx		-	3.0D0*pi/4.0D0		&
						+0.12499895D0*(3.0D0/xx)**1		&														
						-0.00605240D0*(3.0D0/xx)**3		&														
						+0.00135825D0*(3.0D0/xx)**5		&														
						-0.00049616D0*(3.0D0/xx)**7		&														
						+0.00011531D0*(3.0D0/xx)**9														

		BesselY1	=	f1*dsin(theta1)/dsqrt(xx)																		

	end if

	return
end function BesselY1


!===============================================================
function StruveH0(xx)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	xx

! --- Local variables -------------------------------------

	integer					::	ii
	REAL*8					::	P0,P1,P2
	REAL*8					::	P3,P4,P5
	REAL*8					::	a0,a1,a2,a3
	REAL*8					::	b1,b2,b3
	REAL*8					::	c1,c2
	REAL*8					::	yy,StruveH0

	if(xx <= 3.0D0) then
        
		yy	=	(xx/3.0D0)**2
        
		P0	=	+1.909859164D0
		P1	=	-1.909855001D0
		P2	=	+0.687514637D0
		P3	=	-0.126164557D0
 		P4	=	+0.013828813D0
		P5	=	-0.000876918D0
        
		StruveH0	=	P0+(P1+(P2+(P3+(P4+P5*yy)*yy)*yy)*yy)*yy
        
		StruveH0	=	StruveH0*(xx/3.0D0)
        
    else
        
		yy	=	(3.0D0/xx)**2

		a0	=	0.99999906D0
		a1	=	4.77228920D0
		a2	=	3.85542044D0
		a3	=	0.32303607D0

		b1	=	4.88331068D0
		b2	=	4.28957333D0
		b3	=	0.52120508D0

		c1	=	2.0D0*(a0	+	(a1+(a2+a3*yy)*yy)*yy) 
		c2	=	pi*xx*(1.0D0+	(b1+(b2+b3*yy)*yy)*yy) 
																
		StruveH0	=	c1/c2+	BesselY0(xx)																
																				
	end if

	return
end function StruveH0

!===============================================================
function StruveH1(xx)

	implicit none

! --- Variables -------------------------------------------
	REAL*8,intent(in)		::	xx

! --- Local variables -------------------------------------

	integer					::	ii
	REAL*8					::	P1,P2,P3
	REAL*8					::	P4,P5,P6
	REAL*8					::	a0,a1,a2,a3
	REAL*8					::	b1,b2,b3
	REAL*8					::	c1,c2,yy
	REAL*8					::	StruveH1				

	if(xx <= 3.0D0) then
        
		yy	=	(xx/3.0D0)**2
        
		P1	=	+1.909859286D0
		P2	=	-1.145914713D0
		P3	=	+0.294656958D0
		P4	=	-0.042070508D0
 		P5	=	+0.003785727D0
		P6	=	-0.000207183D0
      
		StruveH1	=	(P1+(P2+(P3+(P4+(P5+P6*yy)*yy)*yy)*yy)*yy)*yy

    else
        
		yy	=	(3.0D0/xx)**2

		a0	=	1.00000004D0
		a1	=	3.92205313D0
		a2	=	2.64893033D0
		a3	=	0.27450895D0

		b1	=	3.81095112D0
		b2	=	2.26216956D0
		b3	=	0.10885141D0

		c1	=	2.0D0*(a0	+	(a1+(a2+a3*yy)*yy)*yy) 
		c2	=	pi*(1.0D0	+	(b1+(b2+b3*yy)*yy)*yy) 

		StruveH1	=	c1/c2	+	BesselY1(xx)

	end if

	return
end function StruveH1

end module INFG3D_Open
