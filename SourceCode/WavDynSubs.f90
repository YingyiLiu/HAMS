! ==================================================================================
!   Purpose: This subroutine computes roots of the water-wave dispersion equation 
!             in finite water depth, by using a higher-order iterative method
!
!  License:
! 
!    This routine is part of FinGreen3D.
!
!    FinGreen3D is a free software package: you can redistribute it and/or modify it
!    under the terms of the GNU Lesser General Public License as published by the
!    Free Software Foundation, either version 3 of the License, or (at your option) 
!    any later version.
!
!    You should have received a copy of the GNU General Public License (either V3  
!    or later), along with FinGreen3D. If not, see <http://www.gnu.org/licenses/>.
!
!  Author:
! 
!    Yingyi Liu on Mar.23, 2017
! 
!  Reference:
! 
!    Yingyi Liu, Shigeo Yoshida, Liang Sun, Junliang Gao
!    FinGreen3D: An efficient open-source package for computing free-surface 
!    Green's function in finite water depth using new region-decomposition 
!    strategy with applications to wave-structure interactions
!    Computer Physics Communications, 2016
! 
!    J.N. Newman
!    Numerical solutions of the water-wave dispersion relation
!    Applied Ocean Research 12 (1990) 14-18
!
!  Parameters:
!      Input:   NRT --- Integer, the number of roots required
!                W   --- Real, wave angular frequency
!                H   --- Real, water depth (h>0)
!      Output:  WVN --- Real, an array storing roots of the dispersion equation
! ==================================================================================
  
      SUBROUTINE DISPERSION(WVN,NRT,W,H)
!
!   Evaluation of the roots of the following equations 
!   by higher-order iterative method
!   first root stored in WVN is from Eq. (i)
!   the rest roots are from Eq. (ii)
!   i) w*w/g = k tanh ( kh )
!   ii) -w*w/g = Um tan ( Umh )
!

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NRT
      REAL*8,INTENT(IN):: W,H
      REAL*8,INTENT(OUT):: WVN(1:NRT)
      INTEGER I, M
      REAL*8 T,X,U,Y,DNM,G,PI
      REAL*8 FUN,DFUN,D2FUN,TRIAL,EXX

      DATA G,PI/9.807d0,3.141592653589793d0/

!------------------------------------------------------------------
! I. calculation of wave number (root of Eq. (i))
!------------------------------------------------------------------
!
!   initialize iteration by an accurate Chebyshev approximation
!   if y=x, use the approximation directly insteady of iteration
!   to avoid the singularity in the denomenator of the transcendental
!   function; otherwise, do the iterative procedure. 
!
      X=W*W*H/G
      IF (X.GT.0.D0.AND.X.LE.2.D0) THEN
       Y=DSQRT(X)*(0.9994D0+0.1701D0*X+0.0305*X*X)
      ELSE
       T=X*DEXP(-2.D0*X)
       Y=X+2.D0*T-6.D0*T*T
      ENDIF

      IF (DABS(Y-X).LT.1.E-10) THEN
       WVN(1)=X/H
      ELSE
       M=0
       EXX=1.D0
       DO WHILE (EXX.GT.1.0D-10)
        TRIAL=Y
        DNM=TRIAL*TRIAL-X*X
        FUN=DLOG((TRIAL+X)/(TRIAL-X))/2.D0-TRIAL
        DFUN=-X/DNM-1.D0
        D2FUN=2.D0*X*TRIAL/(DNM*DNM)
        Y=TRIAL-FUN/DFUN*(1.D0+(FUN/DFUN)*(D2FUN/DFUN)/2.D0)
        EXX=DABS(Y-TRIAL)
        M=M+1
       ENDDO
       WVN(1)=Y/H
      ENDIF

!------------------------------------------------------------------
! II. calcultion of roots of Eq. (ii), which characterizes
!     the evanescene modes in eigenfucntion
!------------------------------------------------------------------
!
!   initialize iteration by a suitable starting approximation
!
      U=3.D0*X/(7.D0+3.D0*X)
      T=0.0159D0+0.1032D0*U+4.3152D0*U*U-2.8768D0*U*U*U
!
!   perform iterative procedure to find exact solution of Um (m=1,..NRT-1)
!   of the transcendental equation Eq. (ii)
!
      DO I=2,NRT
       M=0
       EXX=1.D0
       DO WHILE (EXX.GT.1.0D-10)
        TRIAL=T
        Y=(I-1)*PI-TRIAL
        DNM=Y*Y+X*X
        FUN=ATAN2(X,Y)-TRIAL
        DFUN=X/DNM-1.D0
        D2FUN=2.D0*X*TRIAL/(DNM*DNM)
        T=TRIAL-FUN/DFUN*(1.D0+(FUN/DFUN)*(D2FUN/DFUN)/2.D0)
        EXX=DABS(T-TRIAL)
        M=M+1
       ENDDO
       Y=(I-1)*PI-T
       WVN(I)=Y/H
       T=T-PI*X/(X*X+PI*I*(PI*(I-1)-T))
      ENDDO

      END SUBROUTINE DISPERSION